%include "sseutils64.nasm"

extern get_block
extern free_block
extern volume
extern get_distance_element
extern charge
extern hydrophobicity
extern distances

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
	alignb 32
    neg_uno dq -1.0, -1.0, -1.0, -1.0
    ten_double dq 10.0, 10.0, 10.0, 10.0
    four_float dq 4.0, 4.0, 4.0, 4.0
    alignb 32
    zero_float dq 0.0
    alignb 32
    w_hydro dq 0.5, 0.5, 0.5, 0.5
    w_elec dq 0.2, 0.2, 0.2, 0.2
    eps dq 1e-6, 1e-6, 1e-6, 1e-6
    half dq 0.5, 0.5, 0.5, 0.5
    alpha_psi dq -47.0, -47.0, -47.0, -47.0
    alpha_phi dq -57.8, -57.8, -57.8, -57.8
    beta_psi dq 113.0, 113.0, 113.0, 113.0
    beta_phi dq -119.0, -119.0, -119.0, -119.0

section .bss
	alignb 32
    i resq 1
    j resq 1
    N resq 1
    alignb 32
    ris resq 1
    seq resq 1

section .text
global combined_energy_avx, apply_rotation_avx, normalize_avx

; -------------------------------------
; Funzione combined_energy_avx
; -------------------------------------

combined_energy_avx:
    push rbp
    mov rbp, rsp
    sub rsp, 32

    mov [N], rdx               ; Salva N nello stack
    mov [ris], rcx
    mov [seq], rdi

    vxorps ymm2, ymm2, ymm2    ; E_elec = 0.0
    vxorps ymm1, ymm1, ymm1    ; E_hydro = 0.0
    mov qword [i], 0           ; i = 0
    mov rax, [i]

.for_i:
    cmp rax, [N]               ; Confronta con N
    jge .end                   ; Se i >= N, termina

    mov rbx, [i]               ; j = i
    add rbx, 1
    mov [j], rbx

.inner_loop:
    mov rax, [j]
    cmp rax, [N]               ; Confronta con N
    jge .next_i                ; Se j >= N, passa al prossimo i

    ; Passa i parametri per get_distance_element
    mov rdi, rbx               ; Primo parametro: i (in rdi)
    mov rsi, rax               ; Secondo parametro: j (in rsi)
    mov rdx, [N]               ; Terzo parametro: N (in rdx)

    ; Chiamata della funzione
    call get_distance_element 

    mov rcx, [ris]
    mov rdi, [seq]

    ; se dist >= 10.0, salta
    vcomisd xmm0, [ten_double]
    jae .next_j

    ; Calcolo energia idrofobica
    movzx rax, byte [rdi + rbx]
    sub rax, 'A'
    vmovsd xmm4, [hydrophobicity + rax * 8]

    vcomisd xmm4, [neg_uno]
    je .next_j
    
    vmovsd xmm5, [hydrophobicity + rax * 8]

    vcomisd xmm5, [neg_uno]
    je .next_j

    vmulsd xmm4, xmm4, xmm5
    vdivsd xmm4, xmm4, xmm0
    vaddsd xmm1, xmm1, xmm4

    ; Calcolo energia elettrostatica
    
    vmovsd xmm4, [charge + rax * 8]

    vmovsd xmm5, [charge + rax * 8]

    vmulsd xmm4, xmm4, xmm5
    vcomisd xmm4, [zero_float]
    je .next_j

    vmovsd xmm6, xmm0          ; Salva dist
    vmulsd xmm6, xmm6, [four_float]
    vdivsd xmm4, xmm4, xmm6
    vaddsd xmm2, xmm2, xmm4

.next_j:
    add qword [j], 1
    jmp .inner_loop

.next_i:
    add qword [i], 1
    jmp .for_i

.end:
    ; Salva i risultati
    vmulsd xmm2, xmm2, [w_elec]
    vmulsd xmm1, xmm1, [w_hydro]
    vaddsd xmm2, xmm2, xmm1
    vmovsd [rcx], xmm2          ; Salva comb

    add rsp, 32
    pop rbp
    vzeroupper
    ret

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