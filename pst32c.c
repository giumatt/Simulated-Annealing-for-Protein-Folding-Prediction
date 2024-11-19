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

type energy(char* seq, VECTOR phi, VECTOR psi) {
	MATRIX *coords = backbone(seq, phi, psi);

	type rama_e = rama_energy(phi, psi, );
	type hydro_e = hydrophobic_energy(seq, coords);
	type elec_e = electrostatic_energy(seq, coords);
	type pack_e = packing_energy(seq, coords);

	type w_rama = 1.0;
	type w_hydro = 0.5;
	type w_elec = 0.2;
	type w_pack = 0.3;

	type tot_e = (w_rama * rama_e) + (w_hydro * hydro_e) + (w_elec * elec_e) + (w_pack * pack_e);

	return tot_e;
}

type packing_energy(char* seq, VECTOR coords, int N, VECTOR volume) {
	type E = 0;

	for(int i = 0; i < N; i++) {
		type density = 0;
		for(int j = 0; j < N; j++) {
			if(i != j) {
				type dist = sqrt(pow((coords[i] - coords[j]), 2));
				if (dist < 10.0) {
					density = density + (volume[seq[j]] / pow(dist, 3));
				}
			}
		}
		E = E + pow((volume[seq[i]] - density), 2);
	}

	return E;
}

type electrostatic_energy(char* seq, VECTOR coords, int N, VECTOR charge) {
	type E = 0;

	for(int i = 0; i < N; i++) {
		type density = 0;
		for(int j = i + 1; j < N; j++) {
			// Abbiamo tolto l'if che controlla i != j
			type dist = sqrt(pow((coords[i] - coords[j]), 2));
			if ((dist < 10.0) && (charge[seq[i]] != 0.0) && (charge[seq[j]] != 0.0)) {
				E = E + ((charge[seq[i]] * charge[seq[j]]) / (dist * 4.0));
			}
		}
	}

	return E;
}

type hydrophobic_energy(char* seq, VECTOR coords, int N, VECTOR hydrophobicity) {
	type E = 0;

	for(int i = 0; i < N; i++) {
		for(int j = i + 1; j < N; j++) {
			type dist = sqrt(pow((coords[i] - coords[j]), 2));
			if (dist < 10.0) {
				E = E + ((hydrophobicity[seq[i]] * hydrophobicity[seq[j]]) / (dist));
			}
		}
	}

	return E;
}

type rama_energy(VECTOR phi, VECTOR psi, int N) {
	type alpha_psi = -47.0;
	type alpha_phi = -57.8;
	type beta_psi = 113.0;
	type beta_phi = -119.0;

	type E = 0;

	for(int i = 0; i < N; i++) {
		type alpha_dist = sqrt((pow((phi[i] - alpha_phi), 2) + (pow((psi[i] - alpha_psi), 2))));
		type beta_dist = sqrt((pow((phi[i] - beta_phi), 2) + (pow((psi[i] - beta_psi), 2))));

		type min = alpha_dist;
		if (alpha_dist > beta_dist)
			min = beta_dist;

		E = E + (0.5 * min);
	}

	return E;
}

// Funzione per normalizzare un vettore
void normalize(VECTOR v, int n) {
    type norm = 0.0;
    for (int i = 0; i < n; i++) {
        norm += v[i] * v[i];
    }
    norm = sqrt(norm);
    for (int i = 0; i < n; i++) {
        v[i] /= norm;
    }
}

// Funzione per applicare la matrice di rotazione ad un vettore
void apply_rotation(type* vec, type* rot) {
    vec[0] = rot[0] * vec[0] + rot[1] * vec[1] + rot[2] * vec[2];
    vec[1] = rot[3] * vec[0] + rot[4] * vec[1] + rot[5] * vec[2];
    vec[2] = rot[6] * vec[0] + rot[7] * vec[1] + rot[8] * vec[2];
}

MATRIX backbone(char* s, VECTOR phi, VECTOR psi, int N) {
	type r_ca_n = 1.46;
	type r_ca_c = 1.52;
	type r_c_n = 1.33;

	type theta_ca_c_n = 2.028;
	type theta_c_n_ca = 2.124;
	type theta_n_ca_c = 1.940;

	// Da rivedere dimensione
	type coords[3 * N * 3];

	// Da rivedere l'assegnamento
	coords[0] = 0;
	coords[1] = 0;
	coords[2] = 0;

	coords[3] = r_ca_n;
	coords[4] = 0;
	coords[5] = 0;

	coords[6] = r_ca_n + r_ca_c;
	coords[7] = 0;
	coords[8] = 0;

	// coords[1] = (r_ca_n, 0, 0);

	// Support vectors
	VECTOR v1;
	VECTOR v2;
	VECTOR v3;
	MATRIX rot;

	type newv[3];

	newv[0] = 0.0;
	newv[1] = r_c_n;
	newv[2] = 0.0;

	for(int i = 0; i < N; i++) {
		int idx = i * 3;

		if (i > 0) {
			// Parte di aggiornamento per l'atomo N
			for (int j = 0; j < 3; j++) 
				v1[j] = coords[idx - 1 + j] - coords[idx - 2 + j];
			normalize(v1, 3);
			rot = rotation();
			apply_rotation(newv, rot);
			for(int k = 0; k < 3; k++)
				coords[idx * 3 + k] = coords[(idx - 1) * 3 + k] + newv[k];
			
			// Parte di aggiornamento per l'atomo C_a
			for (int j = 0; j < 3; j++) 
				v2[j] = coords[idx + j] - coords[idx - 1 + j];
			normalize(v2, 3);
			rot = rotation();
			newv[0] = 0;
			newv[1] = r_ca_n;
			newv[2] = 0;
			apply_rotation(newv, rot);
			for(int k = 0; k < 3; k++)
				coords[(idx + 1) * 3 + k] = coords[idx * 3 + k] + newv[k];
			
			// Parte di aggiornamento per l'atomo C
			for (int j = 0; j < 3; j++) 
				v3[j] = coords[idx + 1 + j] - coords[idx + j];
			normalize(v3, 3);
			rot = rotation();
			newv[0] = 0;
			newv[1] = r_ca_c;
			newv[2] = 0;
			apply_rotation(newv, rot);
			for(int k = 0; k < 3; k++)
				coords[(idx + 1) * 3 + k] = coords[idx * 3 + k] + newv[k];
		}
		return coords;
	}
}

void pst(params* input){
	char* seq = input->seq;
	int N = input->N;
	int to = input->to;
	int T = to;
	int k = input->k;
	type E = input->e;
	int alpha = input->alpha;
	VECTOR phi = input->phi;
	VECTOR psi = input->psi;

	gen_rnd_mat(phi, N);
	gen_rnd_mat(psi, N);

	type E = energy(seq, phi, psi);

	int t = 0;

	while(T <= 0) {
		int i = rand() % (N + 1);						// Da controllare se è tra 0 e 1

		type delta_phi = (random()*2 * M_PI) - M_PI;
		type delta_psi = (random()*2 * M_PI) - M_PI;
		
		phi[i] = phi[i] + delta_phi;
		psi[i] = psi[i] + delta_psi;

		type delta_energy = energy(seq, phi, psi) - E;

		if (delta_energy <= 0) {
			E = energy(seq, phi, psi);
		} else {
			type P = exp((-delta_energy) / (k * T));		// Attenzione alla funzione divisione
			type r = random();								// Da controllare se è tra 0 e 1

			if (r <= P) {
				E = energy(seq, phi, psi);
			} else {
				phi[i] = phi[i] - delta_phi;
				psi[i] = psi[i] - delta_psi;
			}
		}

		t += 1;
		T = to - sqrt(alpha * t);
	}
	return phi, psi;
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
