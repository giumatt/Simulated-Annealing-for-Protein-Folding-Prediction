; --------------------------------------------------
; Macro per il linking di file oggetto
; NASM a 64-bit con la C standard library e 
; alcune utility per l'uso di SSE e AVX
;
; Uso:
;	%include "sseutils64.nasm"
;	(da specificare nel file sorgente NASM
;	prima della section .data)
; --------------------------------------------------
; Fabrizio Angiulli, 2/3/2013
;

extern	printf

section	.bss

align 32
dbuf:	resq	1

section	.data

;dmask:	db		'%18.16f ',0
dmask:	db		'%f ',0
cr:	db		10,0
br1:	db		'( ',0
br2:	db		')',10,0

%macro	vpush	1
	sub	rsp, 32
	vmovups	[rsp], %1
%endmacro

%macro	vpop	1
	vmovups	%1, [rsp]
	add	rsp, 32
%endmacro

%macro	pushaq	0
	push	rax
	push	rbx
	push	rcx
	push	rdx
;	push	rbp
	push	rsi
	push	rdi
;	push	rsp
	push	r8
	push	r9
	push	r10
	push	r11
	push	r12
	push	r13
	push	r14
	push	r15
	
	vpush	ymm0
	vpush	ymm1
	vpush	ymm2
	vpush	ymm3
	vpush	ymm4
	vpush	ymm5
	vpush	ymm6
	vpush	ymm7
	vpush	ymm8
	vpush	ymm9
	vpush	ymm10
	vpush	ymm11
	vpush	ymm12
	vpush	ymm13
	vpush	ymm14
	vpush	ymm15
%endmacro

%macro	popaq	0
	vpop	ymm15
	vpop	ymm14
	vpop	ymm13
	vpop	ymm12
	vpop	ymm11
	vpop	ymm10
	vpop	ymm9
	vpop	ymm8
	vpop	ymm7
	vpop	ymm6
	vpop	ymm5
	vpop	ymm4
	vpop	ymm3
	vpop	ymm2
	vpop	ymm1
	vpop	ymm0
	
	pop	r15
	pop	r14
	pop	r13
	pop	r12
	pop	r11
	pop	r10
	pop	r9
	pop	r8
;	pop	rsp
	pop	rdi
	pop	rsi
;	pop	rbp
	pop	rdx
	pop	rcx
	pop	rbx
	pop	rax
%endmacro

%macro	start	0
		push	rbp
		mov	rbp, rsp
%endmacro

%macro	stop	0
		mov	rsp, rbp
		pop	rbp
		ret
%endmacro

%macro	prints	1
		pushaq
		xor	rax, rax
		mov	rdi, %1
		call	printf
		popaq
%endmacro

%macro	dprint	1
		pushaq
		mov	rax, 1
		mov	rdi, dmask
		movsd	xmm0, [%1]
		call	printf
		popaq
%endmacro

%macro	sprint	1
		pushaq
		mov		rax, 1
		mov		rdi, dmask
		cvtss2sd	xmm0, [%1]
		call		printf
		popaq
%endmacro

%macro	printss	1
		sprint	%1
		prints	cr
%endmacro

%macro	printps	2
		prints	br1
		push	rdx
		push	rcx
		mov	rdx, %1
		mov	rcx, %2
%%loopps:
		sprint	rdx
		sprint	rdx+4
		sprint	rdx+8
		sprint	rdx+12
		add	rdx, 16
		dec	rcx
		jnz	%%loopps
		pop	rcx
		pop	rdx
		prints	br2
%endmacro

%macro	vprintps	2
		prints	br1
		push	rdx
		push	rcx
		mov	rdx, %1
		mov	rcx, %2
%%vloopps:
		sprint	rdx
		sprint	rdx+4
		sprint	rdx+8
		sprint	rdx+12
		sprint	rdx+16
		sprint	rdx+20
		sprint	rdx+24
		sprint	rdx+28
		add	rdx, 32
		dec	rcx
		jnz	%%vloopps
		pop	rcx
		pop	rdx
		prints	br2
%endmacro

%macro	printsd	1
		dprint	%1
		prints	cr
%endmacro

%macro	printpd	2
		prints	br1
		push	rdx
		push	rcx
		mov	rdx, %1
		mov	rcx, %2
%%looppd:
		dprint	rdx
		dprint	rdx+8
		add	rdx, 16
		dec	rcx
		jnz	%%looppd
		pop	rcx
		pop	rdx
		prints	br2
%endmacro

%macro	vprintpd	2
		prints	br1
		push	rdx
		push	rcx
		mov	rdx, %1
		mov	rcx, %2
%%vlooppd:
		dprint	rdx
		dprint	rdx+8
		dprint	rdx+16
		dprint	rdx+24
		add	rdx, 32
		dec	rcx
		jnz	%%vlooppd
		pop	rcx
		pop	rdx
		prints	br2
%endmacro
