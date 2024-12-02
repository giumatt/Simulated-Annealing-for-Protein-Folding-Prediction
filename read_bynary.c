#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define x86 float
#define x64 double

int main(int argc, char** argv) {

    // Controllo degli argomenti
    if (argc < 3) {
        printf("Usage: ./read_binary <-32|-64|-char> <filename>\n");
        return 1;
    }

    int count = 0; // Contatore per i valori letti

    char* filename = argv[2];
    //printf("Reading %s\n", filename);

    // Apertura del file
    FILE *fp = fopen(filename, "rb");
    if (fp == NULL) {
        printf("Error opening file: %s\n", filename);
        return 1;
    }
    
    

    // Se non è l'opzione -char, leggi i due valori interi iniziali
   
    int firstInt;
    if (fread(&firstInt, sizeof(int), 1, fp) != 1) {
        printf("Error reading the first int value.\n");
        fclose(fp);
        return 1;
    }
    //printf("First int value: %d\n", firstInt);

    int secondInt;
    if (fread(&secondInt, sizeof(int), 1, fp) != 1) {
        printf("Error reading the second int value.\n");
        fclose(fp);
        return 1;
    }
    //printf("Second int value: %d\n", secondInt);
    

    // Determina il tipo di dati da leggere
    if (strcmp(argv[1], "-32") == 0) {
        x86 value;
        //printf("Values (float, 32-bit):\n");
        while (fread(&value, sizeof(x86), 1, fp)) {
            printf("%12.6f ", value); // Stampa con 6 cifre decimali
            count++;
            if (count % 10 == 0) {
                printf("\n");
            }
        }
    } else if (strcmp(argv[1], "-64") == 0) {
        x64 value;
        //printf("Values (double, 64-bit):\n");
        while (fread(&value, sizeof(x64), 1, fp)) {
            printf("%15.8lf ", value); // Stampa con 8 cifre decimali
            count++;
            if (count % 10 == 0) {
                printf("\n");
            }
        }
    } else if (strcmp(argv[1], "-char") == 0) {
        char value;
        //printf("Values (char):\n");
        while (fread(&value, sizeof(char), 1, fp)) {
            printf("%c", value); // Stampa il carattere
            count++;
            if (count % 80 == 0) { // Aggiunge un ritorno a capo ogni 80 caratteri
                printf("\n");
            }
        }
    } else {
        printf("Invalid option: %s. Use -32 for float, -64 for double, or -char for characters.\n", argv[1]);
        fclose(fp);
        return 1;
    }

    // Aggiunge un ritorno a capo finale se necessario
    if (count % 10 != 0) {
        printf("\n");
    }

    //printf("Total values read: %d\n", count); // Stampa il conteggio totale dei valori

    // Chiusura del file
    fclose(fp);
    return 0;
}
