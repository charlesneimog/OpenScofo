	.file	"mfcc.cpp"
	.text
	.p2align 4
	.globl	_Z3fooPKdS0_Pdi
	.type	_Z3fooPKdS0_Pdi, @function
_Z3fooPKdS0_Pdi:
.LFB0:
	.cfi_startproc
	testl	%ecx, %ecx
	movq	%rdi, %r8
	movq	%rsi, %r9
	movq	%rdx, %r10
	jle	.L7
	leal	-1(%rcx), %eax
	vpxor	%xmm0, %xmm0, %xmm0
	cmpl	$14, %eax
	jbe	.L8
	movl	%ecx, %esi
	movq	%rdi, %rax
	vmovapd	%ymm0, %ymm3
	movq	%r9, %rdx
	shrl	$4, %esi
	vmovapd	%ymm0, %ymm1
	vmovapd	%ymm0, %ymm2
	movl	%esi, %edi
	salq	$7, %rdi
	addq	%r8, %rdi
	.p2align 6
	.p2align 4,,10
	.p2align 3
.L4:
	vmovupd	(%rax), %ymm4
	vmovupd	32(%rax), %ymm5
	subq	$-128, %rax
	subq	$-128, %rdx
	vmovupd	-64(%rax), %ymm6
	vmovupd	-32(%rax), %ymm7
	vfmadd231pd	-128(%rdx), %ymm4, %ymm0
	vfmadd231pd	-96(%rdx), %ymm5, %ymm2
	vfmadd231pd	-64(%rdx), %ymm6, %ymm1
	vfmadd231pd	-32(%rdx), %ymm7, %ymm3
	cmpq	%rdi, %rax
	jne	.L4
	vaddpd	%ymm2, %ymm0, %ymm0
	vaddpd	%ymm3, %ymm1, %ymm1
	movl	%esi, %eax
	sall	$4, %eax
	cmpl	%eax, %ecx
	vaddpd	%ymm0, %ymm1, %ymm1
	vextractf128	$0x1, %ymm1, %xmm0
	vaddpd	%xmm1, %xmm0, %xmm0
	vunpckhpd	%xmm0, %xmm0, %xmm1
	vaddpd	%xmm0, %xmm1, %xmm1
	je	.L31
	movl	%eax, %edi
	vzeroupper
.L3:
	subl	%eax, %ecx
	cmpl	$1, %ecx
	je	.L5
	movl	%ecx, %edx
	salq	$3, %rax
	shrl	%edx
	leaq	(%r8,%rax), %rsi
	addq	%r9, %rax
	cmpl	$1, %edx
	vmovupd	(%rsi), %xmm2
	vfmadd231pd	(%rax), %xmm2, %xmm0
	je	.L6
	cmpl	$2, %edx
	vmovupd	16(%rax), %xmm3
	vfmadd231pd	16(%rsi), %xmm3, %xmm0
	je	.L6
	cmpl	$3, %edx
	vmovupd	32(%rax), %xmm2
	vfmadd231pd	32(%rsi), %xmm2, %xmm0
	je	.L6
	cmpl	$4, %edx
	vmovupd	48(%rsi), %xmm3
	vfmadd231pd	48(%rax), %xmm3, %xmm0
	je	.L6
	cmpl	$5, %edx
	vmovupd	64(%rsi), %xmm2
	vfmadd231pd	64(%rax), %xmm2, %xmm0
	je	.L6
	cmpl	$6, %edx
	vmovupd	80(%rsi), %xmm3
	vfmadd231pd	80(%rax), %xmm3, %xmm0
	je	.L6
	vmovupd	96(%rsi), %xmm2
	vfmadd231pd	96(%rax), %xmm2, %xmm0
.L6:
	addl	%edx, %edx
	vunpckhpd	%xmm0, %xmm0, %xmm1
	cmpl	%edx, %ecx
	vaddpd	%xmm0, %xmm1, %xmm1
	je	.L2
	addl	%edx, %edi
.L5:
	vmovsd	(%r8,%rdi,8), %xmm3
	vfmadd231sd	(%r9,%rdi,8), %xmm3, %xmm1
.L2:
	vmovsd	%xmm1, (%r10)
	ret
	.p2align 4,,10
	.p2align 3
.L7:
	vxorpd	%xmm1, %xmm1, %xmm1
	vmovsd	%xmm1, (%r10)
	ret
.L8:
	xorl	%eax, %eax
	xorl	%edi, %edi
	vxorpd	%xmm1, %xmm1, %xmm1
	jmp	.L3
.L31:
	vzeroupper
	jmp	.L2
	.cfi_endproc
.LFE0:
	.size	_Z3fooPKdS0_Pdi, .-_Z3fooPKdS0_Pdi
	.ident	"GCC: (GNU) 16.1.1 20260430"
	.section	.note.GNU-stack,"",@progbits
