%include "sseutils64.nasm"

extern get_block
extern free_block

%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

section .data

section .bss
	alignb 32
    i resq 1
    j resq 1
    N resq 1
    alignb 32
    ris resq 1
    seq resq 1

section .text
global apply_rotation_avx, normalize_avx

; -------------------------------------
; Funzione apply_rotation_avx
; -------------------------------------

apply_rotation_avx:
    push rbp
    mov rbp, rsp

    vzeroall
    vzeroupper  
    
    ; Carica il vettore da rdi in ymm0
    vmovapd ymm0, [rdi]           ; ymm0 contiene il vettore da ruotare

    ; Calcolo di ris[0]
    vmovapd ymm1, [rsi]           ; Carica la prima riga della matrice in ymm1
    vmulpd ymm1, ymm1, ymm0       ; Moltiplicazione elemento per elemento
    vperm2f128 ymm2, ymm1, ymm1, 0x01 ; Scambia i 128 bit alti e bassi
    vaddpd ymm1, ymm1, ymm2       ; Somma le due metà
    vhaddpd ymm1, ymm1, ymm1      ; Somma i due elementi rimanenti
    vmovsd [rdi], xmm1            ; Salva ris[0] nel vettore risultato

    ; Calcolo di ris[1]
    vmovapd ymm1, [rsi + 32]        ; Carica la seconda riga della matrice in ymm1
    vmulpd ymm1, ymm1, ymm0       ; Moltiplicazione elemento per elemento
    vperm2f128 ymm2, ymm1, ymm1, 0x01 ; Scambia i 128 bit alti e bassi
    vaddpd ymm1, ymm1, ymm2       ; Somma le due metà
    vhaddpd ymm1, ymm1, ymm1      ; Somma i due elementi rimanenti
    vmovsd [rdi+8], xmm1          ; Salva ris[1] nel vettore risultato

    ; Calcolo di ris[2]
    vmovapd ymm1, [rsi + 64]        ; Carica la terza riga della matrice in ymm1
    vmulpd ymm1, ymm1, ymm0       ; Moltiplicazione elemento per elemento
    vperm2f128 ymm2, ymm1, ymm1, 0x01 ; Scambia i 128 bit alti e bassi
    vaddpd ymm1, ymm1, ymm2       ; Somma le due metà
    vhaddpd ymm1, ymm1, ymm1      ; Somma i due elementi rimanenti
    vmovsd [rdi+16], xmm1         ; Salva ris[2] nel vettore risultato

    vzeroupper                   
    pop rbp
    ret

; -------------------------------------
; Funzione normalize_avx
; -------------------------------------

normalize_avx:
    push rbp
    mov rbp, rsp
    sub rsp, 32             

    ; Carica i parametri
    vmovapd ymm0, [rdi]          ; Carica il vettore in ymm0

    ; Calcolo della norma
    vmulpd ymm1, ymm0, ymm0      ; Quadrato di ciascun elemento

    ; Somma gli elementi del vettore
    vperm2f128 ymm2, ymm1, ymm1, 0x01  ; Scambia i 128 bit alti e bassi
    vaddpd ymm1, ymm1, ymm2           ; Somma i registri
    vhaddpd ymm1, ymm1, ymm1          ; Somma orizzontalmente i due rimanenti

    ; Calcola la radice quadrata della norma
    vsqrtsd xmm1, xmm1, xmm1          ; Radice quadrata della norma

    ; Normalizzazione del vettore
    vbroadcastsd ymm2, xmm1           ; Duplica la norma in tutti gli elementi di ymm2
    vdivpd ymm0, ymm0, ymm2           ; Divide ciascun elemento del vettore per la norma
    vmovapd [rdi], ymm0               ; Salva il vettore normalizzato in memoria

    add rsp, 32                       
    pop rbp
    vzeroupper                       
    ret