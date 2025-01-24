%include "sseutils32.nasm"

extern get_block
extern free_block

%macro getmem 2
    mov eax, %1
    push eax
    mov eax, %2
    push eax
    call get_block
    add esp, 8
%endmacro

%macro fremem 1
    push %1
    call free_block
    add esp, 4
%endmacro

section .data

section .bss
    alignb 16
    i resd 1
    j resd 1
    N resd 1

section .text
    global normalize_sse, apply_rotation_sse

; -------------------------------------
; Funzione normalize_sse
; Normalizza un vettore di 4 float
; -------------------------------------

normalize_sse:
    push ebp
    mov ebp, esp
    push ebx
    push esi
    push edi

    ; Carica i parametri
    mov esi, [ebp + 8]              

    ; Calcolo della norma
    movaps xmm0, [esi]            ; Carica 4 float da vector1 in xmm0
    mulps xmm0, xmm0              ; Calcola il quadrato di ciascun elemento
    haddps xmm0, xmm0             ; Somma i primi due elementi
    haddps xmm0, xmm0             ; Somma i risultati rimanenti
    sqrtps xmm0, xmm0             ; Calcola la radice quadrata

    ; Normalizzazione del vettore
    movaps xmm1, [esi]            ; Ricarica il vettore originale
    divps xmm1, xmm0              ; Divide ciascun elemento per la norma
    movaps [esi], xmm1            ; Salva il risultato normalizzato in vector1

    pop edi
    pop esi
    pop ebx
    mov esp, ebp
    pop ebp
    ret

; -------------------------------------
; Funzione apply_rotation_sse
; Rotazione un vettore di 9 float
; -------------------------------------

apply_rotation_sse:
    push ebp
    mov ebp, esp
    push esi
    push edi

    ; Carica i parametri
    mov esi, [ebp + 8]    ; Indirizzo del vettore 
    mov edi, [ebp + 12]   ; Indirizzo della matrice di rotazione

    ; Carica il vettore in xmm0
    movaps xmm0, [esi]  
    movaps xmm1, [edi]             ; Carica la prima riga della matrice 
   
    mulps xmm1, xmm0               ; Moltiplica elemento per elemento
    haddps xmm1, xmm1              ; Somma orizzontale
    haddps xmm1, xmm1              ; Somma finale
    movss [esi], xmm1              ; Salva il risultato in ris[0]

    ; Calcola il secondo elemento del risultato (ris[1])
    movaps xmm2, [edi + 16]          ; Carica la seconda riga della matrice
    mulps xmm2, xmm0               ; Moltiplica elemento per elemento
    haddps xmm2, xmm2              ; Somma orizzontale
    haddps xmm2, xmm2              ; Somma finale
    movss [esi + 4], xmm2            ; Salva il risultato in ris[1]

    ; Calcola il terzo elemento del risultato (ris[2])
    movaps xmm3, [edi + 32]          ; Carica la terza riga della matrice
    mulps xmm3, xmm0               ; Moltiplica elemento per elemento
    haddps xmm3, xmm3              ; Somma orizzontale
    haddps xmm3, xmm3              ; Somma finale
    movss [esi + 8], xmm3            ; Salva il risultato in ris[2]

    pop edi
    pop esi
    mov esp, ebp
    pop ebp
    ret
