#!/bin/bash

echo "Choose the architecture to launch:"
echo "1) 32-bit (runpst32)"
echo "2) 64-bit (runpst64)"
echo "3) 32-bit (runpst32_omp)"
echo "4) 64-bit (runpst64_omp)"
read -p "Type the number of your choice (1 or 2): " choice

if [[ "$choice" == "1" ]]; then
    command="./runpst32"
elif [[ "$choice" == "2" ]]; then
    command="./runpst64"
elif [[ "$choice" == "3" ]]; then
    command="./compile32omp"
elif [[ "$choice" == "4" ]]; then
    command="./compile64_omp"
else
    echo "Invalid choice. Exiting."
    exit 1
fi

# Parameters
SEQ="seq_256.ds2"
TO=20
K=1
ALPHA=1
SD=3

if [[ "$command" == "./runpst64" ]]; then
    "$command" -seq "$SEQ" -to "$TO" -k "$K" -alpha "$ALPHA" -sd "$SD" -d
elif [[ "$choice" == "./compile32omp" ]]; then
    "$command" -seq "$SEQ" -to "$TO" -k "$K" -alpha "$ALPHA" -sd "$SD" -d
elif [[ "$choice" == "./compile64_omp" ]]; then
    "$command" -seq "$SEQ" -to "$TO" -k "$K" -alpha "$ALPHA" -sd "$SD" -d
else
    "$command" -seq "$SEQ" -to "$TO" -k "$K" -alpha "$ALPHA" -sd "$SD" -d
fi
