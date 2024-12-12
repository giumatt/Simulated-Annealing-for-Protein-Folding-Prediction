#!/bin/bash

read_binary -64 ../Outputs/phi_float.ds2 > out-loro-phi.txt
read_binary -32 ../out32_256_3_phi.ds2 > out-nostro-phi.txt

file1="out-loro-phi.txt"
file2="out-nostro-phi.txt"

# Variabile per contare i numeri uguali
count=0

# Confronta i due file riga per riga e cerca i numeri uguali con gli indici (arrotondati alla seconda cifra decimale)
paste -d' ' "$file1" "$file2" | awk '{ 
    for(i=1; i<=NF/2; i++) {
        # Arrotonda i numeri alla seconda cifra decimale
        num1 = sprintf("%f", $i)
        num2 = sprintf("%f", $(i+NF/2))
        
        # Confronta i numeri arrotondati
        if (num1 == num2) {
            # Stampa riga, colonna e numero
            #se volessi stampare tutte le posizioni decommento la prossima print
            print "Riga: " NR-1 ", Colonna: " i-1 ", Numero: " num1
            # Incrementa il contatore per i numeri uguali
            count++
        }
    }
}

# Stampa il totale dei numeri uguali
END { 
    print "Totale numeri uguali: " count
}
'
