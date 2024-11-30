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
type energy(char*, VECTOR, VECTOR, MATRIX, int);
type packing_energy(char*, MATRIX, int, type*);
type electrostatic_energy(char*, MATRIX, int, type*);
type hydrophobic_energy(char*, MATRIX, int, type*);
type rama_energy(VECTOR, VECTOR, int);
void normalize(VECTOR, int);
type cosine(type);
type sine(type);
MATRIX rotation(VECTOR, type);
VECTOR apply_rotation(VECTOR, MATRIX);
void backbone(char*, MATRIX, VECTOR, VECTOR, int);

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


///! Funzioni di test:
/*
void test_amino_index() {
    char test_seq[] = "ACDEFGHIKLMNPQRSTVWY";
    char invalid_seq[] = "XYZ123";
    
    printf("Testing `amino_index` with valid amino acids:\n");
    for (int i = 0; i < strlen(test_seq); i++) {
        int index = amino_index(test_seq[i]);
        printf("Amino: %c, Index: %d\n", test_seq[i], index);
        if (index == -1) {
            printf("Error: Valid amino acid '%c' returned -1\n", test_seq[i]);
        }
    }
    
    printf("\nTesting `amino_index` with invalid characters:\n");
    for (int i = 0; i < strlen(invalid_seq); i++) {
        int index = amino_index(invalid_seq[i]);
        printf("Invalid Char: %c, Index: %d\n", invalid_seq[i], index);
        if (index != -1) {
            printf("Error: Invalid character '%c' returned valid index %d\n", invalid_seq[i], index);
        }
    }
}

void test_distance() {
    VECTOR v1 = alloc_matrix(1, 3);
    VECTOR v2 = alloc_matrix(1, 3);

    v1[0] = 0; v1[1] = 0; v1[2] = 0;
    v2[0] = 3; v2[1] = 4; v2[2] = 0;

    type dist = distance(v1, v2);
    printf("Distance between [0,0,0] and [3,4,0]: %.3f\n", dist);
    if (fabs(dist - 5.0) > 1e-6) {
        printf("Error: Distance calculation is incorrect.\n");
    }

    dealloc_matrix(v1);
    dealloc_matrix(v2);
}


void test_normalize() {
    VECTOR v = alloc_matrix(1, 3);
    v[0] = 3; v[1] = 4; v[2] = 0;

    normalize(v, 3);
    printf("Normalized vector: [%.3f, %.3f, %.3f]\n", v[0], v[1], v[2]);

    type norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    printf("Norm after normalization: %.3f\n", norm);

    if (fabs(norm - 1.0) > 1e-6) {
        printf("Error: Normalization is incorrect.\n");
    }

    dealloc_matrix(v);
}


void test_rotation() {
    VECTOR axis = alloc_matrix(1, 3);
    axis[0] = 0; axis[1] = 0; axis[2] = 1; // Rotazione attorno a Z
    type theta = M_PI / 2; // 90 gradi

    MATRIX rot = rotation(axis, theta);
    printf("Rotation Matrix for 90 degrees around Z-axis:\n");
    for (int i = 0; i < 3; i++) {
        printf("[%.3f, %.3f, %.3f]\n", rot[i * 3 + 0], rot[i * 3 + 1], rot[i * 3 + 2]);
    }

    // Test per verificare il determinante
    type det = rot[0] * (rot[4] * rot[8] - rot[5] * rot[7])
             - rot[1] * (rot[3] * rot[8] - rot[5] * rot[6])
             + rot[2] * (rot[3] * rot[7] - rot[4] * rot[6]);
    printf("Determinant of the rotation matrix: %.3f\n", det);

    if (fabs(det - 1.0) > 1e-6) {
        printf("Error: Rotation matrix determinant is incorrect.\n");
    }

    dealloc_matrix(axis);
    dealloc_matrix(rot);
}


void test_backbone() {
    int N = 5; // Numero di residui nella catena
    VECTOR phi = (VECTOR) malloc(sizeof(float) * N);
    VECTOR psi = (VECTOR) malloc(sizeof(float) * N);

    // Definiamo una sequenza casuale e angoli noti
    char seq[] = "ACDEFG";
    for (int i = 0; i < N; i++) {
        phi[i] = M_PI / 2; // 90° in radianti
        psi[i] = M_PI / 2; // 90° in radianti
    }

    // Esegui la funzione backbone
    MATRIX coords = backbone(seq, phi, psi, N);

    // Stampa i risultati
    printf("Backbone coordinates:\n");
    for (int i = 0; i < 3 * N; i++) {
        printf("Atom %d: [%f, %f, %f]\n", i,
               coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]);
    }

    // Libera memoria
    free(phi);
    free(psi);
    free(coords);
}
*/

void pst(params* input){
	VECTOR phi = input->phi;
	VECTOR psi = input->psi;
  	MATRIX coords = alloc_matrix(3 * input->N, 3);

	type to = input->to;
	type T = to;

	type E = energy(input->seq, phi, psi, coords, input->N);
	//printf("Energia iniziale: %.3f\n", E);

	int t = 0;
  	//int cnt = 1;

	type new_E;
	type delta_energy;
	type P, r;

	while(T > 0.0f) {

		int i = rand() % (input->N);

		type delta_phi = (random()*2 * M_PI) - M_PI;
		type delta_psi = (random()*2 * M_PI) - M_PI;

    	//printf("\n");
    	//if(T<20) printf("Inizio %d iterazione da qui:\n", cnt);
    	//printf("phi_i[%d]: %.3f, psi_i[%d]: %.3f\n", i, phi[i], i, psi[i]);
		//printf("delta_phi: %.3f , delta_psi: %.3f, T: %.3f\n", delta_phi, delta_psi, T);
		
		phi[i] = phi[i] + delta_phi;
		psi[i] = psi[i] + delta_psi;

		//printf("phi[%d] %.3f, delta_phi %.3f, T: %.3f\n", i, phi[i], delta_phi, T);
		//printf("psi[%d] %.3f, delta_psi %.3f, T: %.3f\n", i, psi[i], delta_psi, T);

		new_E = energy(input->seq, phi, psi, coords, input->N);
		//printf("new_E: %.3f, T: %.3f\n", new_E, T);
		// type delta_energy = energy(input->seq, phi, psi, input->N) - E;

		delta_energy = new_E - E;
		//printf("delta_energy %.3f, T: %.3f\n",  delta_energy, T);

		if (delta_energy <= 0) {
			E = new_E;
      		//printf("Configurazione accettata per decremento energia!\n");
		} else {
			//printf("delta_energy_else: %.3f, T: %.3f\n",  delta_energy, T);
			P = (exp((-delta_energy) / (input->k * T)));  	//! Attenzione alla funzione divisione
			r = random();        							//! Da controllare se è tra 0 e 1
			//printf("p %.3f, r: %.3f  \n", P, r);
			if (r <= P) {
				E = new_E;
				//printf("Configurazione accettata con prob: %.3f, E: %.3f, T: %.3f\n", r, E, T);
			}else {
        		//printf("Configurazione rifiutata! Reset dei vettori:\n");
				phi[i] = phi[i] - delta_phi;
				psi[i] = psi[i] - delta_psi;
				//printf("phi reset: %.3f  , delta_phi: %.3f, i: %d \n",  phi[i], delta_phi, i);
				//printf("psi reset: %.3f  , delta_psi: %.3f, i: %d \n",  psi[i], delta_psi, i);
			}
		}

		t += 1;
		T = to - sqrtf(input->alpha * t);
    	// cnt++;
	}

	input->e = E;
	input->phi = phi;
	input->psi = psi;

	dealloc_matrix(coords);
}

type energy(char* seq, VECTOR phi, VECTOR psi, MATRIX coords, int N) {
	// La matrice coords viene passata nei parametri di energy che la
	// distribuisce alle varie energy

	type w_rama = 1.0f;
	type w_hydro = 0.5f;
	type w_elec = 0.2f;
	type w_pack = 0.3f;

	backbone(seq, coords, phi, psi, N);

	type rama_e = rama_energy(phi, psi, N);
	//printf("Rama_e: %.3f\n", rama_e);
	type hydro_e = hydrophobic_energy(seq, coords, N, hydrophobicity);
	//printf("Hydro %.3f\n", hydro_e);
	type elec_e = electrostatic_energy(seq, coords, N, charge);
	//printf("Elec: %.3f\n", elec_e);
	type pack_e = packing_energy(seq, coords, N, volume);
	//printf("Pack_e: %.3f\n", pack_e);

	type tot_e = (w_rama * rama_e) + (w_hydro * hydro_e) + (w_elec * elec_e) + (w_pack * pack_e);
	//printf("TOT_E: %.3f\n", tot_e);

	return tot_e;
}

type packing_energy(char* seq, MATRIX coords, int N, type* volume) {
    type E = 0.0f;
    VECTOR v = alloc_matrix(1, 3); 				// Alloca un vettore 3D per il punto v
    VECTOR w = alloc_matrix(1, 3); 				// Alloca un vettore 3D per il punto w

    for (int i = 0; i < N; i++) {
        type density = 0;

        // Caricamento dei C_alpha nel vettore v
        for (int k = 0; k < 3; k++) {			
            v[k] = coords[(i * 9) + 3 + k];
			//printf("V: %.3f, coords[%d]: %.3f\n", v[k], i, coords[(i * 9) + 3 + k]);
        }
		
        for (int j = 0; j < N; j++) {
            if (i != j) {
				// Caricamento dei C_alpha nel vettore w
                for (int k = 0; k < 3; k++) {
                    w[k] = coords[(j * 9) + 3 + k];
					//printf("W: %.3f, coords[%d]: %.3f\n", w[k], j, coords[(j * 9) + 3 + k]);
                }
				
                type dist = distance(v, w); 			// Calcola la distanza tra v e w
                //printf("Distance: %.3f, i: %d\n", dist, i);
				//if ((dist > 1e-6) && (dist < 10.0f)) {
                if ((dist < 10.0f)) {
                    int amminoacido_j = amino_index(seq[j]);
                    if (volume[amminoacido_j] > 0) {
                        density += ((volume[amminoacido_j]) / (dist * dist * dist));
                        //printf("Amminoacido: %c, Indice: %d, Volume: %.3f, Densità: %.3f, Distance: %.3f\n",
                        //       seq[j], aminoacido_j, volume[aminoacido_j], density, dist);
                    }
                }
            }
        }
        // Calcolo della differenza di densità
        int amminoacido_i = amino_index(seq[i]);
        if (volume[amminoacido_i] > 0) {
            type diff = volume[amminoacido_i] - density;
            E += diff * diff;
        }
    }

    dealloc_matrix(v);
    dealloc_matrix(w);

    return E;
}

type electrostatic_energy(char* seq, MATRIX coords, int N, type* charge) {
	type E = 0.0f;
	VECTOR v = alloc_matrix(1, 3);
    VECTOR w = alloc_matrix(1, 3);

	for(int i = 0; i < N; i++) {
		type density = 0;

		// Caricamento dei C_alpha nel vettore v
		for (int k = 0; k < 3; k++) {
            v[k] = coords[(i * 9) + 3 + k];
        }
		
		for(int j = i + 1; j < N; j++) {
			// Caricamento dei C_alpha nel vettore w
			for (int k = 0; k < 3; k++) {
                    w[k] = coords[(j * 9) + 3 + k];
                }
			//! Abbiamo tolto l'if che controlla i != j
			type dist = distance(v, w);

			int amminoacido_i = amino_index(seq[i]);
			int amminoacido_j = amino_index(seq[j]);
			
			//vedere se inserire le chiamate direttamente nell'if
			if ((dist < 10.0f) && ((charge[amminoacido_i] != 0.0f) && (charge[amminoacido_j] != 0.0f))
				&& ((volume[amminoacido_i] != -1) && (volume[amminoacido_j] != -1))) {
				//la carica può essere anche -1, quindi verifico che l'amminoacido esista facendo riferimento a volume
				E += ((charge[amminoacido_i] * charge[amminoacido_j]) / (dist * 4.0f));
			}
		}
	}
	return E;

	dealloc_matrix(v);
    dealloc_matrix(w);
}

type hydrophobic_energy(char* seq, MATRIX coords, int N, type* hydrophobicity) {
	type E = 0.0f;

	VECTOR v = alloc_matrix(1, 3);
    VECTOR w = alloc_matrix(1, 3);


	for(int i = 0; i < N; i++) {
		// Caricamento dei C_alpha nel vettore v
		for (int k = 0; k < 3; k++) {
            v[k] = coords[(i * 9) + 3 + k];
        }			

		for(int j = i + 1; j < N; j++) {
			// Caricamento dei C_alpha nel vettore v
			for (int k = 0; k < 3; k++) {
                    w[k] = coords[(j * 9) + 3 + k];
            }		

			type dist = distance(v, w);

			int amminoacido_i = amino_index(seq[i]);
			int amminoacido_j = amino_index(seq[j]);

			if ((dist < 10.0f) &&
				((hydrophobicity[amminoacido_i] != -1) && (volume[amminoacido_j] != -1))) {
				E += ((hydrophobicity[amminoacido_i] * hydrophobicity[amminoacido_j]) / (dist));
			}
		}
	}

	return E;

	dealloc_matrix(v);
    dealloc_matrix(w);
}

type rama_energy(VECTOR phi, VECTOR psi, int N) {
    type alpha_psi = -47.0f;
    type alpha_phi = -57.8f;
    type beta_psi = 113.0f;
    type beta_phi = -119.0f;

    type E = 0.0f;

    for(int i = 0; i < N; i++) {
        type alpha_dist = sqrtf(( (phi[i] - alpha_phi) * (phi[i] - alpha_phi) ) + 
                                ( (psi[i] - alpha_psi) * (psi[i] - alpha_psi) ));
        
        type beta_dist = sqrtf(( (phi[i] - beta_phi) * (phi[i] - beta_phi) ) + 
                                ( (psi[i] - beta_psi) * (psi[i] - beta_psi) ));

        type min = alpha_dist;
        if (alpha_dist > beta_dist)
            min = beta_dist;

        E += (0.5f * min);
    }

    return E;
}

int amino_index(char amino) {
    // Array di amminoacidi standard in ordine
    const char aminoacids[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    for (int i = 0; i < 26; i++) {
        if (aminoacids[i] == amino) {
            return i; // Ritorna l'indice corrispondente
        }
    }
    return -1;
}

void backbone(char* seq, MATRIX coords, VECTOR phi, VECTOR psi, int N) {
	type r_ca_n = 1.46f;
	type r_ca_c = 1.52f;
	type r_c_n = 1.33f;

	type theta_ca_c_n = 2.028f;
	type theta_c_n_ca = 2.124f;
	type theta_n_ca_c = 1.940f;

	coords[0] = 0.0f;
	coords[1] = 0.0f;
	coords[2] = 0.0f;

	coords[3] = r_ca_n;
	coords[4] = 0.0f;
	coords[5] = 0.0f;

	// Vettori di direzione
	VECTOR v1 = alloc_matrix(1, 3);
	VECTOR v2 = alloc_matrix(1, 3);
	VECTOR v3 = alloc_matrix(1, 3);

	MATRIX rot = alloc_matrix(3, 3);

	VECTOR newv = alloc_matrix(1, 3);

	for(int i = 0; i < N; i++) {
		int idx = i * 3;

		if (i > 0) {
			// Parte di aggiornamento per l'atomo N
			for (int j = 0; j < 3; j++) {
				//v1[j] = coords[((idx - 1) * 3) + j] - coords[((idx - 2) * 3) + j];
				v1[j] = coords[((idx - 2) * 3) + j] - coords[((idx - 1) * 3) + j];
				//printf("For iteration [%d, %d] v1[%d] is: %.3f\n", i, j, j, v1[j]);
				// printf("idx: %d\n", idx);
			}
			normalize(v1, 3);
			
			rot = rotation(v1, theta_c_n_ca);

			//! Prima newv era esterno al for quindi apply_rotation nnon veniva fatto correttamente

			newv[0] = 0.0f;
			newv[1] = r_c_n;
			newv[2] = 0.0f;

			newv = apply_rotation(newv, rot);
			
			for(int k = 0; k < 3; k++)
				coords[(idx * 3) + k] = coords[((idx - 1) * 3) + k] + (newv[k]);
			
			// Parte di aggiornamento per l'atomo C_a
			for (int j = 0; j < 3; j++)
				v2[j] = coords[(idx * 3) + j] - coords[((idx - 1) * 3) + j];
				
			normalize(v2, 3);
			
			rot = rotation(v2, phi[i]);
	
			newv[0] = 0.0f;
			newv[1] = r_ca_n;		// !
			newv[2] = 0.0f;
			newv = apply_rotation(newv, rot);
			
			for(int k = 0; k < 3; k++) {
				coords[((idx + 1) * 3) + k] = coords[((idx * 3)) + k] + newv[k];
				//printf("New coords at position [%d, %d] is: %.3f\n", i, k, coords[k]);
				//printf("New coords at position [%d, %d] is: %.3f\n", i, k, coords[((idx + 1) * 3) + k]);
			}
		}

		// Parte di aggiornamento per l'atomo C
		for (int j = 0; j < 3; j++) {
			v3[j] = coords[((idx + 1) * 3) + j] - coords[(idx * 3) + j];
			//printf("v3: %.3f, coords[idx + 1 + j]: %.3f - coords[idx + j]: %.3f\n", v3[j], coords[(idx + 1)*3 + j], coords[idx*3 + j] );
		}
		normalize(v3, 3);
		
		rot = rotation(v3, psi[i]);

		newv[0] = 0.0f;
		newv[1] = r_ca_c;
		newv[2] = 0.0f;

		newv = apply_rotation(newv, rot);

		for(int k = 0; k < 3; k++) {
			coords[((idx + 2) * 3) + k] = coords[((idx + 1) * 3) + k] + newv[k];
			//printf(" coords[((idx + 2) * 3) + k]: %.3f - coords[((idx + 1) * 3) + k]: %.3f  i: %d\n", coords[((idx + 2) * 3) + k], coords[((idx + 1) * 3) + k] , i);
		
		}
	
	}

	dealloc_matrix(v1);
	dealloc_matrix(v2);
	dealloc_matrix(v3);
	dealloc_matrix(rot);
	dealloc_matrix(newv);
}

type distance(VECTOR v, VECTOR w) {
    type dist = 0.0;
    for (int i = 0; i < 3; i++) {
        type diff = v[i] - w[i];
		    dist += diff * diff;
    }
	return sqrtf(dist);
}

void normalize(VECTOR v, int n) {
    type norm = 0;
    for (int i = 0; i < n; i++) {
        norm += v[i] * v[i];
    }
    norm = sqrtf(norm);

	if (norm != 0) {
    	for (int i = 0; i < n; i++)
        	v[i] = v[i] / norm;
	}
}

type cosine(type theta) {
	type theta2 = theta * theta;
    return 1.0f - (theta2 / 2.0f) + (theta2 * theta2 / 24.0f) - (theta2 * theta2 * theta2 / 720.0f);
}

type sine(type theta) {
	type theta2 = theta * theta;
    return theta - (theta * theta2 / 6.0f) + (theta * theta2 * theta2 / 120.0f) - (theta * theta2 * theta2 * theta2 / 5040.0f);
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

VECTOR apply_rotation(VECTOR vec, MATRIX rot) {

	VECTOR ris= alloc_matrix(1,3);

    ris[0] = (rot[0] * vec[0]) + (rot[3] * vec[1]) + (rot[6] * vec[2]);
    ris[1] = (rot[1] * vec[0]) + (rot[4] * vec[1]) + (rot[7] * vec[2]);
    ris[2] = (rot[2] * vec[0]) + (rot[5] * vec[1]) + (rot[8] * vec[2]);

	return ris;
	
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

	/*
	printf("Testing `amino_index`...\n");
    test_amino_index();
    printf("\nTesting `distance`...\n");
    test_distance();
    printf("\nTesting `normalize`...\n");
    test_normalize();
    printf("\nTesting `rotation`...\n");
    test_rotation();
    printf("\nTesting `backbone`...\n");
    test_backbone();
	*/

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
		if(input->phi == NULL || input->psi==NULL)				//!!! abbiamo inserito input->psi==NULL, prima era input->psi
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
