#!/bin/bash

echo "Choose the architecture to run:"
echo "1) 32-bit"
echo "2) 64-bit"
echo "3) 32-bit + SSE"
echo "4) 64-bit + AVX"
echo "5) 32-bit + OpenMP"
echo "6) 64-bit + OpenMP"
read -p "Type the number of your choice (1 or 2): " choice

commands=(
    "./C/runpst32"
    "./C/runpst64"
    "./C + SSE_AVX/runpst32"
    "./C + SSE_AVX/runpst64"
    "./OpenMP/runpst32_omp"
    "./OpenMP/runpst64_omp"
)

# Parameters
SEQ="seq_256.ds2"
TO=20
K=1
ALPHA=1
SD=3

if [[ "$choice" -ge 1 && "$choice" -le 6 ]]; then
    command="${commands[choice-1]}"
    echo "Executing: $command -seq $SEQ -to $TO -k $K -alpha $ALPHA -sd $SD -d"
    eval "$command -seq $SEQ -to $TO -k $K -alpha $ALPHA -sd $SD -d"
else
    echo "Invalid choice. Exiting."
    exit 1
fi