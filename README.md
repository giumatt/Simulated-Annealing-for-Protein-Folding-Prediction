
# **Protein Folding Prediction**

## **Project Description**
This project implements an algorithm for predicting the tertiary structure of proteins using the simulated annealing method.
The solution was developed in C, with subsequent optimizations in Assembly and OpenMP to improve performance.

The prediction of the three-dimensional structure of a protein is based on the optimization of the dihedral angles φ and ψ that define the folding
of the amino acid chain, minimizing the overall energy of the system.

---

## **Main Features**
- Calculation of the total energy of a protein structure through various contributions (Ramachandran, hydrophobic, electrostatic, packing).
- Optimization of the configuration using the simulated annealing algorithm.
- Modularity for subsequent optimizations in Assembly and OpenMP.

---

## **System Requirements**
### Required Software:
- **GCC** (C compiler)
- **NASM** (assembler for Assembly language)
- **System Libraries**:
  - `libc6-dev-i386`
  - `lib32gcc-<version>-dev`

### Installation (on Ubuntu):
```bash
sudo apt-get install gcc nasm libc6-dev-i386 lib32gcc-10-dev
```

---

## **Compilation Instructions**
### To compile and execute:
```bash
nasm -f elf32 pst32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o pst32.o pst32c.c -o pst32c -lm
./pst32c -seq <sequence_file> -to <initial_temperature> -alpha <cooling_rate> -k <constant> -sd <seed>
```

### Example Execution:
```bash
./pst32c -seq input.ds2 -to 100 -alpha 0.9 -k 1 -sd 42
```

---

## **Repository Structure**
- **`pst32c.c`**: C source code.
- **`pst32.nasm`**: Assembly source code for optimizations.
- **`README.md`**: Project documentation.
- **`seq_256.ds2`**: Example file containing an amino acid sequence.
- **`output/`**: Folder for execution results (predicted φ and ψ).

---

# **Predizione della Struttura delle Proteine**

## **Descrizione del Progetto**
Questo progetto implementa un algoritmo per la predizione della struttura terziaria delle proteine utilizzando il metodo del *simulated annealing*.
La soluzione è stata sviluppata in linguaggio C, con ottimizzazioni successive in Assembly e OpenMP per migliorare le prestazioni.

La predizione della struttura tridimensionale di una proteina si basa sull'ottimizzazione degli angoli diedri φ e ψ che definiscono il ripiegamento della catena
aminoacidica, minimizzando l'energia complessiva del sistema.

---

## **Funzionalità Principali**
- Calcolo dell'energia complessiva di una struttura proteica tramite diversi contributi (Ramachandran, idrofobico, elettrostatico, impaccamento).
- Ottimizzazione della configurazione tramite l'algoritmo di *simulated annealing*.
- Modularità per successive ottimizzazioni in linguaggio Assembly e OpenMP.

---

## **Requisiti di Sistema**
### Software richiesto:
- **GCC** (compilatore C)
- **NASM** (assembler per linguaggio Assembly)
- **Librerie di sistema**:
  - `libc6-dev-i386`
  - `lib32gcc-<versione>-dev`

### Installazione (su Ubuntu):
```bash
sudo apt-get install gcc nasm libc6-dev-i386 lib32gcc-10-dev
```

---

## **Istruzioni di Compilazione**
### Per compilare ed eseguire:
```bash
nasm -f elf32 pst32.nasm && gcc -m32 -msse -O0 -no-pie sseutils32.o pst32.o pst32c.c -o pst32c -lm
./pst32c -seq <file_sequenza> -to <temperatura_iniziale> -alpha <tasso_raffreddamento> -k <costante> -sd <seed>
```

### Esempio di esecuzione:
```bash
./pst32c -seq input.ds2 -to 100 -alpha 0.9 -k 1 -sd 42
```

---

## **Struttura della Repository**
- **`pst32c.c`**: Codice sorgente in linguaggio C.
- **`pst32.nasm`**: Codice sorgente in linguaggio Assembly per ottimizzazioni.
- **`README.md`**: Documentazione del progetto.
- **`seq_256.ds2`**: File di esempio contenente una sequenza aminoacidica.
- **`output/`**: Cartella per i risultati dell'esecuzione (`φ` e `ψ` predetti).
