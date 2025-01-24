%include "sseutils32.nasm"

extern get_block
extern free_block
extern get_distance_index
extern charge
extern hydrophobicity

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
    alignb 16
    neg_uno dd -1.0, -1.0, -1.0, -1.0
    ten_float dd 10.0, 10.0, 10.0, 10.0
    four_float dd 4.0, 4.0, 4.0, 4.0
    alignb 16
    zero_float dd 0.0
    w_elec dd 0.2
    w_hydro dd 0.5
    alignb 16
    alpha_psi dd -47.0, -47.0, -47.0, -47.0             
    alpha_phi dd -57.8, -57.8, -57.8, -57.8             
    beta_psi dd 113.0, 113.0, 113.0, 113.0                  
    beta_phi dd -119.0, -119.0, -119.0, -119.0

section .bss
    alignb 16
    i resd 1
    j resd 1
    N resd 1
    

section .text
    global combined_energy_sse, normalize_sse, apply_rotation_sse

; -------------------------------------
; Funzione combined_energy_sse
; Calcola il valore combinato delle energie elettrostatica e idrofobica
; -------------------------------------

combined_energy_sse:
    push ebp
    mov ebp, esp
    push esi
    push edi
    push ebx

    ; Carica i parametri
    mov edi, [ebp + 8]          ; seq
    mov esi, [ebp + 12]         ; distances
    mov ecx, [ebp + 16]         ; N
    mov eax, [ebp + 20]         ; &comb
    mov dword [N], ecx

    xorps xmm0, xmm0          ; E_elec = 0.0
    xorps xmm1, xmm1          ; E_hydro = 0.0

    mov dword [i], 0          ; i = 0

.for_i:
    cmp dword [i], ecx
    jge .end                  ; Se i >= N, termina

    mov ebx, dword [i]
    mov dword [j], ebx         ; j = i
    inc dword [j]              ; j = i + 1

.inner_loop:
    cmp dword [j], ecx
    jge .next_i         ; Se j >= N, passo al prossimo i

    mov edx, dword [j]

    push ecx
    push edx
    push ebx
    call get_distance_index
    add esp, 12

    mov ecx, dword [N]

    shl eax, 2
    movss xmm3, [esi + eax]   ; dist = distances[get_distance_index(i, j)]

    movss xmm4, [ten_float]
    comiss xmm3, xmm4
    jae .next_j               ; Salta se dist >= 10.0

    ; Calcolo energia idrofobica
    mov eax, dword [i]
    movzx eax, byte [edi + eax]
    sub eax, 'A'
    shl eax, 2
    movss xmm4, [hydrophobicity + eax]
    comiss xmm4, [neg_uno]
    je .next_i               ; Salta se hydrophobicity[i] = -1

    mov eax, dword [j]
    movzx eax, byte [edi + eax]
    sub eax, 'A'
    shl eax, 2
    movss xmm5, [hydrophobicity + eax]
    comiss xmm5, [neg_uno]
    je .next_j               ; Salta se hydrophobicity[j] = -1

    mulss xmm4, xmm5          ; hydrophobicity[i] * hydrophobicity[j]
    divss xmm4, xmm3          ; / dist
    addss xmm1, xmm4           

    ; Calcolo energia elettrostatica
    mov eax, dword [i]
    movzx eax, byte [edi + eax]
    sub eax, 'A'
    shl eax, 2
    movss xmm4, [charge + eax]

    mov eax, dword [j]
    movzx eax, byte [edi + eax]
    sub eax, 'A'
    shl eax, 2
    movss xmm5, [charge + eax]

    mulss xmm4, xmm5          ; charge[seq[i]-65] * charge[seq[j]-65]
    comiss xmm4, [zero_float]
    je .next_j                ; Salta se risultato == 0

    movss xmm6, xmm3          ; Salva dist in xmm6 prima di moltiplicare
    mulss xmm6, [four_float]  ; dist * 4.0f
    divss xmm4, xmm6          ; (charge[i]*charge[j]) / (dist * 4.0)
    addss xmm0, xmm4         

.next_j:
    inc dword [j]             ; j++
    jmp .inner_loop

.next_i:
    inc dword [i]             ; i++
    jmp .for_i

.end:
    ; Salva i risultati
    mov eax, [ebp + 20]
    mulss xmm0, [w_elec]
    mulss xmm1, [w_hydro]
    addss xmm0, xmm1
    movss [eax], xmm0         ; comb

    pop ebx
    pop edi
    pop esi
    mov esp, ebp
    pop ebp
    ret

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
    movss [esi+4], xmm2            ; Salva il risultato in ris[1]

    ; Calcola il terzo elemento del risultato (ris[2])
    movaps xmm3, [edi + 32]          ; Carica la terza riga della matrice
    mulps xmm3, xmm0               ; Moltiplica elemento per elemento
    haddps xmm3, xmm3              ; Somma orizzontale
    haddps xmm3, xmm3              ; Somma finale
    movss [esi+8], xmm3            ; Salva il risultato in ris[2]

    pop edi
    pop esi
    mov esp, ebp
    pop ebp
    ret
