/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2024/25
* 
* Progetto dell'algoritmo Predizione struttura terziaria proteine 221 231 a
* in linguaggio assembly x86-32 + SSE
* 
* F. Angiulli F. Fassetti S. Nisticò, novembre 2024
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf32 pst32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o pst32.o pst32c.c -o pst32c -lm && ./pst32c $pars
* 
* oppure
* 
* ./runpst32
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

#define random() (((type) rand())/RAND_MAX)

// Definizione delle nostre funzioni
type distance(VECTOR, VECTOR);
type energy(char*, VECTOR, VECTOR, int);
type packing_energy(char*, MATRIX, int, type*);
type electrostatic_energy(char*, MATRIX, int, type*);
type hydrophobic_energy(char*, MATRIX, int, type*);
type rama_energy(VECTOR, VECTOR, int);
void normalize(VECTOR, int);
type cosine(type);
type sine(type);
MATRIX rotation(VECTOR, type);
void apply_rotation(VECTOR, MATRIX);
MATRIX backbone(char*, VECTOR, VECTOR, int);

int amino_index(char);

type hydrophobicity[] = {1.8, -1, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -1, -3.9, 3.8, 1.9, -3.5, -1, -1.6, -3.5, -4.5, -0.8, -0.7, -1, 4.2, -0.9, -1, -1.3, -1};		// hydrophobicity
type volume[] = {88.6, -1, 108.5, 111.1, 138.4, 189.9, 60.1, 153.2, 166.7, -1, 168.6, 166.7, 162.9, 114.1, -1, 112.7, 143.8, 173.4, 89.0, 116.1, -1, 140.0, 227.8, -1, 193.6, -1};		// volume
type charge[] = {0, -1, 0, -1, -1, 0, 0, 0.5, 0, -1, 1, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0, -1, 0, -1};		// charge

typedef struct {
	char* seq;		// sequenza di amminoacidi
	int N;			// lunghezza sequenza
	unsigned int sd; 	// seed per la generazione casuale
	type to;		// temperatura INIZIALE
	type alpha;		// tasso di raffredamento
	type k;		// costante
	VECTOR hydrophobicity; // hydrophobicity
	VECTOR volume;		// volume
	VECTOR charge;		// charge
	VECTOR phi;		// vettore angoli phi
	VECTOR psi;		// vettore angoli psi
	type e;		// energy
	int display;
	int silent;

} params;


/*
* 
*	Le funzioni sono state scritte assumendo che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,16); 
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}

char* alloc_char_matrix(int rows, int cols) {
	return (char*) get_block(sizeof(char),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri floating-point a precisione singola
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 
* 	load_seq
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*1 byte: matrix data in row-major order --> charatteri che compongono la stringa
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
char* load_seq(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	
	char* data = alloc_char_matrix(rows,cols);
	status = fread(data, sizeof(char), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o floating-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

/*
* 	save_out
* 	=========
* 
*	Salva su file un array lineare composto da k elementi.
* 
* 	Codifica del file:
* 	primi 4 byte: contenenti l'intero 1 		--> numero intero a 32 bit
* 	successivi 4 byte: numero di elementi k     --> numero intero a 32 bit
* 	successivi byte: elementi del vettore 		--> k numero floating-point a precisione singola
*/
void save_out(char* filename, MATRIX X, int k) {
	FILE* fp;
	int i;
	int n = 1;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(X, sizeof(type), k, fp);
	}
	fclose(fp);
}

/*
* 	gen_rnd_mat
* 	=========
* 
*	Genera in maniera casuale numeri reali tra -pi e pi
*	per riempire una struttura dati di dimensione Nx1
* 
*/
void gen_rnd_mat(VECTOR v, int N){
	int i;

	for(i=0; i<N; i++){
		// Campionamento del valore + scalatura
		v[i] = (random()*2 * M_PI) - M_PI;
	}
}

// PROCEDURE ASSEMBLY
// extern void prova(params* input);

type distance(VECTOR v, VECTOR w) {
    type dist = 0.0;
    for (int i = 0; i < 3; i++) {
        type diff = v[i] - w[i];
		dist += diff * diff;
    }
	return (type)sqrtf(dist);
}

type energy(char* seq, VECTOR phi, VECTOR psi, int N) {
	
	MATRIX coords = alloc_matrix(3 * N, 3);

	coords = backbone(seq, phi, psi, N);
			

	//MATRIX *coords = backbone(input);

	type rama_e = rama_energy(phi, psi, N);
	printf("Rama_e: %.3f\n", rama_e);
	type hydro_e = hydrophobic_energy(seq, coords, N, hydrophobicity);
	printf("Hydro %.3f\n", hydro_e);
	type elec_e = electrostatic_energy(seq, coords, N, charge);
	printf("Elec: %.3f\n", hydro_e);
	type pack_e = packing_energy(seq, coords, N, volume);
	printf("Pack_e: %.3f\n", pack_e);

	type w_rama = 1.0f;
	type w_hydro = 0.5f;
	type w_elec = 0.2f;
	type w_pack = 0.3f;

	type tot_e = (w_rama * rama_e) + (w_hydro * hydro_e) + (w_elec * elec_e) + (w_pack * pack_e);
	printf("TOT_E: %.3f\n", tot_e);
	return tot_e;
}

/*
type packing_energy(char* seq, MATRIX coords, int N, type* volume) {
	type E = 0;
	// printf("\nPACKING ENERGY:\n");
	VECTOR v = alloc_matrix(1, 3);
	VECTOR w = alloc_matrix(1, 3);
	for(int i = 0; i < N; i++) {
		type density = 0;
		
		for (int k = 0; k < 3; k++) {
        	v[k] = coords[(k + i + 1)*3];
    	}
		
		//v[i] = coords[(i +1)*3];
		
		
		for(int j = 0; j < N; j++) {
			if(i != j) {
				//w[j] = coords[((j+1) * 3)];				// !
				for (int k = 0; k < 3; k++) {
        			w[k] = coords[(k + j + 1)*3];
    			}
				//printf("V: %.3f, coord[%f]+ %d\n", v[i], coords[(i +1)*3], i);
				//printf("W: %.3f coord[%f]+ %d\n", w[j], coords[(j +1)*3], j);
				type dist = distance(v, w);
				//printf("\nDistance: %.3f\n", dist);		
				if ((dist > 1e-6) && (dist < 10.0f)) {
					int aminoacido = amino_index(seq[j]);
						if(aminoacido >=0){
							density += (volume[aminoacido] / (dist * dist * dist));
							// printf("\n");
							//printf("Volume %d: %.3f\n", j, volume[seq[j]]);
							printf("%c",  aminoacido);
							//printf("\n");
							//printf("Density: %.3f\n", density);
					}
				}
			}
		}
		// printf("\nDistance: %.3f\n", distance);
		type diff = volume[(unsigned char)seq[i]] - density;
		// printf("\n");
		// printf("Volume %d: %.3f\n", i, volume[seq[i]]);
		E += diff * diff;
		// E = E + pow((volume[seq[i]] - density), 2);
		// printf("\n");
		// printf("Energy: %.3f\n", E);
		// ! volume[seq[j]] e volume[seq[i]] stampano le stesse cose 
	}
	return E;
}
*/
type packing_energy(char* seq, MATRIX coords, int N, type* volume) {
    type E = 0;
    VECTOR v = alloc_matrix(1, 3); // Alloca un vettore 3D per il punto v
    VECTOR w = alloc_matrix(1, 3); // Alloca un vettore 3D per il punto w

    for (int i = 0; i < N; i++) {
        type density = 0;

        // Copia le coordinate del punto `i` in `v`
        for (int k = 0; k < 3; k++) {
            v[k] = coords[(k + i + 1)*3];
        }

        for (int j = 0; j < N; j++) {
            if (i != j) {
                // Copia le coordinate del punto `j` in `w`
                for (int k = 0; k < 3; k++) {
                    w[k] = coords[(k + j + 1)*3];
                }

                type dist = distance(v, w); // Calcola la distanza tra `v` e `w`

                if ((dist > 1e-6) && (dist < 10.0f)) {
                    int aminoacido = amino_index(seq[j]);
                    if (aminoacido >= 0) {
                        density += (volume[aminoacido] / (dist * dist * dist));
                        // Debug: stampa informazioni utili
                        printf("Amminoacido: %c, Indice: %d, Volume: %.3f, Densità: %.3f\n",
                               seq[j], aminoacido, volume[aminoacido], density);
                    }
                }
            }
        }

        // Calcola la differenza di densità
        int aminoacido_i = amino_index(seq[i]);
        if (aminoacido_i >= 0) {
            type diff = volume[aminoacido_i] - density;
            E += diff * diff; // Aggiungi il contributo all'energia totale
        }
    }

    free(v); // Libera la memoria allocata per `v`
    free(w); // Libera la memoria allocata per `w`

    return E;
}

int amino_index(char amino) {
    // Array di amminoacidi standard in ordine
    const char amino_acids[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const int num_amino_acids = 20; // Numero di amminoacidi standard

    for (int i = 0; i < num_amino_acids; i++) {
        if (amino_acids[i] == amino  && volume[i] !=-1) {
            return i; // Ritorna l'indice corrispondente
        }
    }

    // Se non trovato, stampa un messaggio di errore e ritorna -1
    return -1;
}

type electrostatic_energy(char* seq, MATRIX coords, int N, type* charge) {
	type E = 0;

	for(int i = 0; i < N; i++) {
		VECTOR v = &coords[i * 3];
		type density = 0;
		for(int j = i + 1; j < N; j++) {
			VECTOR w = &coords[j * 3];
			// Abbiamo tolto l'if che controlla i != j
			type dist = distance(v, w);
			if ((dist < 10.0f) && (charge[seq[i]] != 0.0f) && (charge[seq[j]] != 0.0f)) {
				E = E + ((charge[seq[i]] * charge[seq[j]]) / (dist * 4.0f));
			}
		}
	}

	return E;
}

type hydrophobic_energy(char* seq, MATRIX coords, int N, type* hydrophobicity) {
	type E = 0;

	for(int i = 0; i < N; i++) {
		VECTOR v = &coords[i * 3];			// Primo atomo C_alpha
		for(int j = i + 1; j < N; j++) {
			VECTOR w = &coords[j * 3];		// Secondo atomo C_alpha
			//type dist = sqrt(pow((coords[i] - coords[j]), 2));
			type dist = distance(v, w);

			if (dist < 10.0f) {
				E = E + ((hydrophobicity[seq[i]] * hydrophobicity[seq[j]]) / (dist));
			}
		}
	}

	return E;
}

type rama_energy(VECTOR phi, VECTOR psi, int N) {
	type alpha_psi = -47.0f;
	type alpha_phi = -57.8f;
	type beta_psi = 113.0f;
	type beta_phi = -119.0f;

	type E = 0;

	for(int i = 0; i < N; i++) {
		type alpha_dist = sqrtf((pow((phi[i] - alpha_phi), 2) + (pow((psi[i] - alpha_psi), 2))));
		type beta_dist = sqrtf((pow((phi[i] - beta_phi), 2) + (pow((psi[i] - beta_psi), 2))));

		type min = alpha_dist;
		if (alpha_dist > beta_dist)
			min = beta_dist;

		E = E + (0.5f * min);
	}

	return E;
}

// Funzione per normalizzare un vettore
void normalize(VECTOR v, int n) {
    type norm = 0;
    for (int i = 0; i < n; i++) {
        norm += v[i] * v[i];
    }
    norm = sqrtf(norm);

	if (norm != 0) {
    	for (int i = 0; i < n; i++) {
        	v[i] /= norm;
    	}
	}
}

type cosine(type theta) {
	type ret = 0.0;
	ret = 1 - (pow(theta, 2) / 2.0f) + (pow(theta, 4) / 24.0f) - (pow(theta, 6) / 720.0f);
	return ret;
}

type sine(type theta) {
	type ret = 0.0;
	ret = theta - (pow(theta, 3) / 6.0f) + (pow(theta, 5) / 120.0f) - (pow(theta, 7) / 5040.0f);
	return ret; 
}

MATRIX rotation(VECTOR axis, type theta) {
	MATRIX rot = alloc_matrix(3, 3);

	for(int i = 0; i < 9; i++)
		rot[i] = 0;

	type scalar = sqrtf(pow(axis[0], 2) + pow(axis[1], 2) + pow(axis[2], 2));
	for(int i = 0; i < 3; i++) {
		axis[i] = axis[i] / scalar;
	}

	type a = cosine(theta / 2.0f);

	VECTOR bcd = alloc_matrix(1, 3);

	for(int i = 0; i < 3; i++) {
		bcd[i] = (-1.0f) * (axis[i]) * (sine(theta / 2.0f));
	}

	rot[0] = (pow(a, 2)) + (pow(bcd[0], 2)) - (pow(bcd[1], 2)) - (pow(bcd[2], 2));
	rot[1] = (2.0f) * ((bcd[0]) * (bcd[1]) + (a * bcd[2]));
	rot[2] = (2.0f) * ((bcd[0]) * (bcd[2]) - (a * bcd[1]));
	rot[3] = (2.0f) * ((bcd[0]) * (bcd[1]) - (a * bcd[2]));
	rot[4] = (pow(a, 2)) + (pow(bcd[1], 2)) - (pow(bcd[0], 2)) - (pow(bcd[2], 2));
	rot[5] = (2.0f) * ((bcd[1]) * (bcd[2]) + (a * bcd[0]));
	rot[6] = (2.0f) * ((bcd[0]) * (bcd[2]) + (a * bcd[1]));
	rot[7] = (2.0f) * ((bcd[1]) * (bcd[2]) - (a * bcd[0]));
	rot[8] = (pow(a, 2)) + (pow(bcd[2], 2)) - (pow(bcd[0], 2)) - (pow(bcd[1], 2));

	dealloc_matrix(bcd);

	return rot;
}

// Funzione per applicare la matrice di rotazione ad un vettore
void apply_rotation(VECTOR vec, MATRIX rot) {
    vec[0] = (rot[0] * vec[0]) + (rot[1] * vec[1]) + (rot[2] * vec[2]);
    vec[1] = (rot[3] * vec[0]) + (rot[4] * vec[1]) + (rot[5] * vec[2]);
    vec[2] = (rot[6] * vec[0]) + (rot[7] * vec[1]) + (rot[8] * vec[2]);
}

// MATRIX backbone(char* s, VECTOR phi, VECTOR psi, int N) {
MATRIX backbone(char* seq, VECTOR phi, VECTOR psi, int N) {
	type r_ca_n = 1.46f;
	type r_ca_c = 1.52f;
	type r_c_n = 1.33f;

	type theta_ca_c_n = 2.028f;
	type theta_c_n_ca = 2.124f;
	type theta_n_ca_c = 1.940f;

	// Da rivedere dimensione
	// MATRIX coords[3 * input->N * 3];
	MATRIX coords = alloc_matrix(3 * N, 3);

	// Da rivedere l'assegnamento
	coords[0] = 0;
	coords[1] = 0;
	coords[2] = 0;

	coords[3] = r_ca_n;
	coords[4] = 0;
	coords[5] = 0;
	

	// coords[6] = r_ca_n + r_ca_c;
	// coords[7] = 0;
	// coords[8] = 0;

	// coords[1] = (r_ca_n, 0, 0);

	// Support vectors
	VECTOR v1 = alloc_matrix(1, 3);
	VECTOR v2 = alloc_matrix(1, 3);
	VECTOR v3 = alloc_matrix(1, 3);

	MATRIX rot = alloc_matrix(3, 3);

	VECTOR newv = alloc_matrix(1, 3);

	newv[0] = 0;
	newv[1] = r_c_n;		// !
	newv[2] = 0;

	for(int i = 0; i < N; i++) {
		int idx = i * 3;

		if (i > 0) {
			// Parte di aggiornamento per l'atomo N
			for (int j = 0; j < 3; j++) {
				v1[j] = coords[idx - 1 + j] - coords[idx - 2 + j];
				// printf("For iteration [%d, %d] v1[%d] is: %.3f\n", i, j, j, v1[j]);
				// printf("idx: %d\n", idx);
			}
			normalize(v1, 3);
			//for (int j = 0; j < 3; j++)
				// printf("After normalization v1[%d, %d] is: %.3f\n", i, j, v1[j]);
			rot = rotation(v1, theta_c_n_ca);
			apply_rotation(newv, rot);
			// for (int j = 0; j < 3; j++)
				// printf("After applying rotation newv[%d, %d] is: %.3f\n", i, j, newv[j]);
			for(int k = 0; k < 3; k++) {
				coords[(idx * 3) + k] = coords[(idx - 1) * 3 + k] + (newv[k]);
				// printf("New coords at position [%d, %d] is: %.3f\n", i, k, coords[(idx * 3) + k]);
			}
			
			// Parte di aggiornamento per l'atomo C_a
			for (int j = 0; j < 3; j++) 
				v2[j] = coords[idx + j] - coords[(idx - 1) + j];
			normalize(v2, 3);
			rot = rotation(v2, phi[i]);
			newv[0] = 0;
			newv[1] = r_ca_n;		// !
			newv[2] = 0;
			apply_rotation(newv, rot);
			for(int k = 0; k < 3; k++) {
				coords[(idx + 1) * 3 + k] = coords[(idx * 3) + k] + newv[k];
				// printf("New coords at position [%d, %d] is: %.3f\n", i, k, coords[k]);
				// printf("New coords at position [%d, %d] is: %.3f\n", i, k, coords[(idx + 1) * 3 + k]);
			}
		}

		// Parte di aggiornamento per l'atomo C
		for (int j = 0; j < 3; j++) 
			v3[j] = coords[idx + 1 + j] - coords[idx + j];
		normalize(v3, 3);
		rot = rotation(v3, psi[i]);
		newv[0] = 0;
		newv[1] = r_ca_c;
		newv[2] = 0;
		apply_rotation(newv, rot);
		for(int k = 0; k < 3; k++) {
			coords[((idx + 2) * 3) + k] = coords[((idx + 1) * 3) + k] + newv[k];
		
		}
	

	

	}

	return coords;
}

void pst(params* input){
	/*
	char* seq = input->seq;
	int N = input->N;
	int to = input->to;
	int T = to;
	int k = input->k;
	type E = input->e;
	int alpha = input->alpha;
	VECTOR phi = input->phi;
	VECTOR psi = input->psi;
	*/

	VECTOR phi = input->phi;
	VECTOR psi = input->psi;

	int to = input->to;
	int T = to;

	gen_rnd_mat(input->phi, input->N);
	gen_rnd_mat(input->psi, input->N);

	type E = energy(input->seq, phi, psi, input->N);

	int t = 0;

	while(T <= 0) {
		int i = rand() % (input->N + 1);						// Da controllare se è tra 0 e 1

		type delta_phi = (random()*2 * M_PI) - M_PI;
		type delta_psi = (random()*2 * M_PI) - M_PI;
		
		phi[i] = phi[i] + delta_phi;
		psi[i] = psi[i] + delta_psi;

		type delta_energy = energy(input->seq, phi, psi, input->N) - E;

		if (delta_energy <= 0) {
			E = energy(input->seq, phi, psi, input->N);
		} else {
			type P = (type)(exp((-delta_energy) / (input->k * T)));		// Attenzione alla funzione divisione
			type r = random();								// Da controllare se è tra 0 e 1

			if (r <= P) {
				E = energy(input->seq, phi, psi, input->N);
			} else {
				phi[i] = phi[i] - delta_phi;
				psi[i] = psi[i] - delta_psi;
			}
		}

		t += 1;
		T = to - sqrtf(input->alpha * t);
	}
}

int main(int argc, char** argv) {
	char fname_phi[256];
	char fname_psi[256];
	char* seqfilename = NULL;
	clock_t t;
	float time;
	int d;
	
	//
	// Imposta i valori di default dei parametri
	//
	params* input = malloc(sizeof(params));
	input->seq = NULL;	
	input->N = -1;			
	input->to = -1;
	input->alpha = -1;
	input->k = -1;		
	input->sd = -1;		
	input->phi = NULL;		
	input->psi = NULL;
	input->silent = 0;
	input->display = 0;
	input->e = -1;
	input->hydrophobicity = hydrophobicity;
	input->volume = volume;
	input->charge = charge;


	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//
	if(argc <= 1){
		printf("%s -seq <SEQ> -to <to> -alpha <alpha> -k <k> -sd <sd> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tSEQ: il nome del file ds2 contenente la sequenza amminoacidica\n");
		printf("\tto: parametro di temperatura\n");
		printf("\talpha: tasso di raffredamento\n");
		printf("\tk: costante\n");
		printf("\tsd: seed per la generazione casuale\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-seq") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing dataset file name!\n");
				exit(1);
			}
			seqfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-to") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing to value!\n");
				exit(1);
			}
			input->to = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-alpha") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing alpha value!\n");
				exit(1);
			}
			input->alpha = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sd") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing seed value!\n");
				exit(1);
			}
			input->sd = atoi(argv[par]);
			par++;
		}else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//
	if(seqfilename == NULL || strlen(seqfilename) == 0){
		printf("Missing ds file name!\n");
		exit(1);
	}

	input->seq = load_seq(seqfilename, &input->N, &d);

	
	if(d != 1){
		printf("Invalid size of sequence file, should be %ix1!\n", input->N);
		exit(1);
	} 

	if(input->to <= 0){
		printf("Invalid value of to parameter!\n");
		exit(1);
	}

	if(input->k <= 0){
		printf("Invalid value of k parameter!\n");
		exit(1);
	}

	if(input->alpha <= 0){
		printf("Invalid value of alpha parameter!\n");
		exit(1);
	}

	input->phi = alloc_matrix(input->N, 1);
	input->psi = alloc_matrix(input->N, 1);
	// Impostazione seed 
	srand(input->sd);
	// Inizializzazione dei valori
	gen_rnd_mat(input->phi, input->N);
	gen_rnd_mat(input->psi, input->N);

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Dataset file name: '%s'\n", seqfilename);
		printf("Sequence lenght: %d\n", input->N);
	}

	// COMMENTARE QUESTA RIGA!
	// prova(input);
	//

	//
	// Predizione struttura terziaria
	//
	t = clock();
	pst(input);
	t = clock() - t;
	time = ((float)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("PST time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato
	//
	sprintf(fname_phi, "out32_%d_%d_phi.ds2", input->N, input->sd);
	save_out(fname_phi, input->phi, input->N);
	sprintf(fname_psi, "out32_%d_%d_psi.ds2", input->N, input->sd);
	save_out(fname_psi, input->psi, input->N);
	if(input->display){
		if(input->phi == NULL || input->psi)
			printf("out: NULL\n");
		else{
			int i,j;
			printf("energy: %f, phi: [", input->e);
			for(i=0; i<input->N; i++){
				printf("%f,", input->phi[i]);
			}
			printf("]\n");
			printf("psi: [");
			for(i=0; i<input->N; i++){
				printf("%f,", input->psi[i]);
			}
			printf("]\n");
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	dealloc_matrix(input->phi);
	dealloc_matrix(input->psi);
	free(input);

	return 0;
}
