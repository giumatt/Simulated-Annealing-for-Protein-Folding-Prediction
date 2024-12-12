#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv) {
    // Controllo degli argomenti
    if (argc != 3) {
        printf("Usage: %s <file1> <file2>\n", argv[0]);
        return 1;
    }

    char* file1 = argv[1];
    char* file2 = argv[2];
    int count = 0;       // Contatore per i valori letti
    double sum = 0.0;    // Somma delle differenze al quadrato

    // Apertura dei file
    FILE *fp1 = fopen(file1, "r");
    if (fp1 == NULL) {
        perror("Error opening file1");
        return 1;
    }

    FILE *fp2 = fopen(file2, "r");
    if (fp2 == NULL) {
        perror("Error opening file2");
        fclose(fp1);
        return 1;
    }

    double val1, val2;
    // Leggi i file valore per valore
    while (fscanf(fp1, "%lf", &val1) == 1 && fscanf(fp2, "%lf", &val2) == 1) {
        double diff = val1 - val2;
        sum += diff * diff;
        count++;
    }  

    fclose(fp1);
    fclose(fp2);

    // Controlla se sono stati letti valori
    if (count == 0) {
        printf("Error: No valid data found in the files.\n");
        return 1;
    }

    // Calcola lo scarto quadratico medio
    double rmsd = sqrt(sum / count);
    printf("RMSD: %.8lf\n", rmsd);

    return 0;
}
