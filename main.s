	.file	"main.cpp"
	.text
	.section	.text._ZNKSt5ctypeIcE8do_widenEc,"axG",@progbits,_ZNKSt5ctypeIcE8do_widenEc,comdat
	.align 2
	.p2align 4
	.weak	_ZNKSt5ctypeIcE8do_widenEc
	.type	_ZNKSt5ctypeIcE8do_widenEc, @function
_ZNKSt5ctypeIcE8do_widenEc:
.LFB1535:
	.cfi_startproc
	movl	%esi, %eax
	ret
	.cfi_endproc
.LFE1535:
	.size	_ZNKSt5ctypeIcE8do_widenEc, .-_ZNKSt5ctypeIcE8do_widenEc
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC0:
	.ascii	"Eigen::CwiseNull"
	.string	"aryOp<NullaryOp, MatrixType>::CwiseNullaryOp(Eigen::Index, Eigen::Index, const NullaryOp&) [with NullaryOp = Eigen::internal::scalar_constant_op<std::complex<double> >; PlainObjectType = Eigen::Matrix<std::complex<double>, -1, -1>; Eigen::Index = long int]"
	.align 8
.LC1:
	.string	"/usr/local/include/Eigen/src/Core/CwiseNullaryOp.h"
	.align 8
.LC2:
	.string	"rows >= 0 && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == rows) && cols >= 0 && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == cols)"
	.text
	.align 2
	.p2align 4
	.type	_ZN5Eigen9DenseBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8ConstantEllRKS3_.part.0, @function
_ZN5Eigen9DenseBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8ConstantEllRKS3_.part.0:
.LFB14825:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	leaq	.LC0(%rip), %rcx
	movl	$71, %edx
	leaq	.LC1(%rip), %rsi
	leaq	.LC2(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE14825:
	.size	_ZN5Eigen9DenseBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8ConstantEllRKS3_.part.0, .-_ZN5Eigen9DenseBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8ConstantEllRKS3_.part.0
	.section	.rodata.str1.8
	.align 8
.LC3:
	.ascii	"Eigen::DenseCoeffsBase<Derived, 0>::CoeffReturnType Eigen::D"
	.ascii	"enseCoeffsBase<Derived, 0>::operator()(Eigen::Index, Eigen::"
	.ascii	"Index) const [with Derived = Eigen::CwiseBinaryOp<Eigen::int"
	.ascii	"ernal::scalar_product_op<double, double>, const Eigen::Cwise"
	.ascii	"NullaryOp<Eigen::internal::scalar_constant_op<double>, const"
	.ascii	" Eigen::Matri"
	.string	"x<double, -1, -1> >, const Eigen::CwiseBinaryOp<Eigen::internal::scalar_difference_op<double, double>, const Eigen::Matrix<double, -1, -1>, const Eigen::Transpose<const Eigen::Matrix<double, -1, -1> > > >; CoeffReturnType = double; Eigen::Index = long int]"
	.align 8
.LC4:
	.string	"/usr/local/include/Eigen/src/Core/DenseCoeffsBase.h"
	.align 8
.LC5:
	.string	"row >= 0 && row < rows() && col >= 0 && col < cols()"
	.text
	.align 2
	.p2align 4
	.type	_ZNK5Eigen15DenseCoeffsBaseINS_13CwiseBinaryOpINS_8internal17scalar_product_opIddEEKNS_14CwiseNullaryOpINS2_18scalar_constant_opIdEEKNS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEEKNS1_INS2_20scalar_difference_opIddEESA_KNS_9TransposeISA_EEEEEELi0EEclEll, @function
_ZNK5Eigen15DenseCoeffsBaseINS_13CwiseBinaryOpINS_8internal17scalar_product_opIddEEKNS_14CwiseNullaryOpINS2_18scalar_constant_opIdEEKNS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEEKNS1_INS2_20scalar_difference_opIddEESA_KNS_9TransposeISA_EEEEEELi0EEclEll:
.LFB10606:
	.cfi_startproc
	movq	40(%rdi), %rcx
	cmpq	16(%rcx), %rsi
	jge	.L6
	movq	8(%rcx), %r8
	cmpq	%r8, %rdx
	jge	.L6
	movq	32(%rdi), %r9
	movq	%rdi, %rax
	movq	8(%r9), %rdi
	movq	(%r9), %r9
	imulq	%rdx, %rdi
	addq	%rsi, %rdi
	imulq	%r8, %rsi
	movsd	(%r9,%rdi,8), %xmm0
	addq	%rdx, %rsi
	movq	(%rcx), %rdx
	subsd	(%rdx,%rsi,8), %xmm0
	mulsd	24(%rax), %xmm0
	ret
.L6:
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	leaq	.LC3(%rip), %rcx
	movl	$118, %edx
	leaq	.LC4(%rip), %rsi
	leaq	.LC5(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE10606:
	.size	_ZNK5Eigen15DenseCoeffsBaseINS_13CwiseBinaryOpINS_8internal17scalar_product_opIddEEKNS_14CwiseNullaryOpINS2_18scalar_constant_opIdEEKNS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEEKNS1_INS2_20scalar_difference_opIddEESA_KNS_9TransposeISA_EEEEEELi0EEclEll, .-_ZNK5Eigen15DenseCoeffsBaseINS_13CwiseBinaryOpINS_8internal17scalar_product_opIddEEKNS_14CwiseNullaryOpINS2_18scalar_constant_opIdEEKNS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEEKNS1_INS2_20scalar_difference_opIddEESA_KNS_9TransposeISA_EEEEEELi0EEclEll
	.align 2
	.p2align 4
	.type	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb1ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_.constprop.0, @function
_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb1ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_.constprop.0:
.LFB14843:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	leaq	-1(%rdi), %rax
	movapd	%xmm0, %xmm15
	movapd	%xmm1, %xmm14
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%rsi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$264, %rsp
	.cfi_def_cfa_offset 320
	movq	(%rdx), %rsi
	movq	8(%rdx), %rdx
	movq	%rcx, 56(%rsp)
	leaq	-3(%rdi), %rcx
	movq	%rdx, %r10
	movq	%rax, 104(%rsp)
	xorl	%eax, %eax
	salq	$4, %r10
	movq	%rdi, 80(%rsp)
	movq	%r8, 96(%rsp)
	movq	%rsi, 88(%rsp)
	cmpq	$32000, %r10
	ja	.L11
	subq	$7, %rdi
	movq	%rdi, 64(%rsp)
	testq	%rdi, %rdi
	jle	.L11
	leaq	(%rdx,%rdx,2), %r14
	leaq	(%rdx,%rdx,4), %r15
	movq	%rdx, %r11
	movq	%rdx, %r9
	salq	$4, %r15
	movq	%r8, %rbx
	movq	%r10, 120(%rsp)
	movq	%r14, %r8
	salq	$5, %r9
	salq	$4, %r8
	movq	%rdx, %rdi
	subq	%r10, %rsi
	salq	$5, %r14
	movq	%rdx, 112(%rsp)
	negq	%r11
	salq	$6, %rdi
	leaq	0(,%rdx,8), %r13
	salq	$7, %r11
	movq	%r10, %r12
	movq	%rcx, 128(%rsp)
	subq	%rdx, %r13
	movq	%rbp, %rcx
	movq	%r15, %rdx
	movq	%rsi, 72(%rsp)
	salq	$4, %r13
	movq	%r14, %rbp
	xorl	%esi, %esi
	movq	%r9, %r14
	movapd	.LC7(%rip), %xmm12
	movq	%r13, %r15
	movq	%rcx, %r9
	movq	%r8, %r13
	movq	.LC8(%rip), %xmm13
	movq	%rdx, %r8
	.p2align 4,,10
	.p2align 3
.L12:
	testq	%r9, %r9
	jle	.L56
	movq	56(%rsp), %rcx
	movq	%rax, (%rsp)
	pxor	%xmm8, %xmm8
	xorl	%edx, %edx
	movapd	%xmm8, %xmm7
	movapd	%xmm8, %xmm6
	movapd	%xmm8, %xmm0
	movq	(%rcx), %r10
	movq	72(%rsp), %rcx
	movapd	%xmm8, %xmm10
	movapd	%xmm8, %xmm4
	movapd	%xmm8, %xmm5
	movapd	%xmm8, %xmm2
	addq	%r12, %rcx
	.p2align 4,,10
	.p2align 3
.L13:
	movq	%rdx, %rax
	movupd	(%rcx), %xmm9
	addq	$1, %rdx
	salq	$4, %rax
	movupd	(%r10,%rax), %xmm1
	xorpd	%xmm12, %xmm9
	leaq	(%rsi,%rcx), %rax
	addq	$16, %rcx
	pshufd	$238, %xmm9, %xmm11
	pshufd	$68, %xmm9, %xmm9
	mulpd	%xmm1, %xmm9
	pshufd	$78, %xmm1, %xmm3
	mulpd	%xmm3, %xmm11
	xorpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm9
	addpd	%xmm9, %xmm2
	movupd	(%rax,%r12), %xmm9
	xorpd	%xmm12, %xmm9
	pshufd	$238, %xmm9, %xmm11
	pshufd	$68, %xmm9, %xmm9
	mulpd	%xmm3, %xmm11
	mulpd	%xmm1, %xmm9
	xorpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm9
	addpd	%xmm9, %xmm5
	movupd	(%rax,%r14), %xmm9
	xorpd	%xmm12, %xmm9
	pshufd	$238, %xmm9, %xmm11
	pshufd	$68, %xmm9, %xmm9
	mulpd	%xmm3, %xmm11
	mulpd	%xmm1, %xmm9
	xorpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm9
	addpd	%xmm9, %xmm4
	movupd	(%rax,%r13), %xmm9
	xorpd	%xmm12, %xmm9
	pshufd	$238, %xmm9, %xmm11
	pshufd	$68, %xmm9, %xmm9
	mulpd	%xmm3, %xmm11
	mulpd	%xmm1, %xmm9
	xorpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm9
	addpd	%xmm9, %xmm10
	movupd	(%rax,%rdi), %xmm9
	xorpd	%xmm12, %xmm9
	pshufd	$238, %xmm9, %xmm11
	pshufd	$68, %xmm9, %xmm9
	mulpd	%xmm3, %xmm11
	mulpd	%xmm1, %xmm9
	xorpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm9
	addpd	%xmm9, %xmm0
	movupd	(%rax,%r8), %xmm9
	xorpd	%xmm12, %xmm9
	pshufd	$238, %xmm9, %xmm11
	pshufd	$68, %xmm9, %xmm9
	mulpd	%xmm3, %xmm11
	mulpd	%xmm1, %xmm9
	xorpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm9
	addpd	%xmm9, %xmm6
	movupd	(%rax,%rbp), %xmm9
	xorpd	%xmm12, %xmm9
	pshufd	$238, %xmm9, %xmm11
	pshufd	$68, %xmm9, %xmm9
	mulpd	%xmm3, %xmm11
	mulpd	%xmm1, %xmm9
	xorpd	%xmm13, %xmm11
	addpd	%xmm11, %xmm9
	addpd	%xmm9, %xmm7
	movupd	(%rax,%r15), %xmm9
	xorpd	%xmm12, %xmm9
	pshufd	$238, %xmm9, %xmm11
	pshufd	$68, %xmm9, %xmm9
	mulpd	%xmm11, %xmm3
	mulpd	%xmm9, %xmm1
	xorpd	%xmm13, %xmm3
	addpd	%xmm3, %xmm1
	addpd	%xmm1, %xmm8
	cmpq	%r9, %rdx
	jne	.L13
	movapd	%xmm2, %xmm3
	movq	(%rsp), %rax
	movapd	%xmm4, %xmm9
	unpckhpd	%xmm4, %xmm4
	unpckhpd	%xmm3, %xmm3
	movhpd	%xmm10, (%rsp)
	movapd	%xmm0, %xmm11
	movhpd	%xmm5, 8(%rsp)
	movhpd	%xmm0, 16(%rsp)
	movhpd	%xmm6, 24(%rsp)
	movhpd	%xmm7, 32(%rsp)
	movhpd	%xmm8, 40(%rsp)
.L23:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	mulsd	%xmm3, %xmm1
	mulsd	%xmm2, %xmm0
	subsd	%xmm1, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm2, %xmm1
	movsd	%xmm1, 48(%rsp)
	movapd	%xmm15, %xmm1
	mulsd	%xmm3, %xmm1
	addsd	48(%rsp), %xmm1
	ucomisd	%xmm0, %xmm1
	jp	.L57
.L14:
	addsd	8(%rbx), %xmm1
	addsd	(%rbx), %xmm0
	movsd	8(%rsp), %xmm2
	movsd	%xmm0, (%rbx)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 8(%rbx)
	mulsd	%xmm5, %xmm0
	movapd	%xmm2, %xmm1
	mulsd	%xmm14, %xmm1
	mulsd	%xmm15, %xmm2
	subsd	%xmm1, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm5, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L58
.L15:
	addsd	24(%rbx), %xmm1
	addsd	16(%rbx), %xmm0
	movapd	%xmm14, %xmm2
	mulsd	%xmm9, %xmm2
	movsd	%xmm0, 16(%rbx)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 24(%rbx)
	mulsd	%xmm9, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm4, %xmm1
	subsd	%xmm1, %xmm0
	movapd	%xmm15, %xmm1
	mulsd	%xmm4, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L59
.L16:
	addsd	40(%rbx), %xmm1
	addsd	32(%rbx), %xmm0
	movapd	%xmm14, %xmm2
	movsd	(%rsp), %xmm5
	mulsd	%xmm10, %xmm2
	movsd	%xmm0, 32(%rbx)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 40(%rbx)
	mulsd	%xmm10, %xmm0
	movapd	%xmm5, %xmm1
	mulsd	%xmm14, %xmm1
	mulsd	%xmm15, %xmm5
	subsd	%xmm1, %xmm0
	movapd	%xmm5, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L60
.L17:
	addsd	56(%rbx), %xmm1
	addsd	48(%rbx), %xmm0
	movapd	%xmm14, %xmm2
	movsd	16(%rsp), %xmm4
	mulsd	%xmm11, %xmm2
	movsd	%xmm0, 48(%rbx)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 56(%rbx)
	mulsd	%xmm11, %xmm0
	movapd	%xmm4, %xmm1
	mulsd	%xmm14, %xmm1
	mulsd	%xmm15, %xmm4
	subsd	%xmm1, %xmm0
	movapd	%xmm4, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L61
.L18:
	addsd	72(%rbx), %xmm1
	addsd	64(%rbx), %xmm0
	movapd	%xmm14, %xmm2
	movsd	24(%rsp), %xmm5
	mulsd	%xmm6, %xmm2
	movsd	%xmm0, 64(%rbx)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 72(%rbx)
	mulsd	%xmm6, %xmm0
	movapd	%xmm5, %xmm1
	mulsd	%xmm14, %xmm1
	mulsd	%xmm15, %xmm5
	subsd	%xmm1, %xmm0
	movapd	%xmm5, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L62
.L19:
	addsd	88(%rbx), %xmm1
	addsd	80(%rbx), %xmm0
	movapd	%xmm14, %xmm2
	movsd	32(%rsp), %xmm6
	mulsd	%xmm7, %xmm2
	movsd	%xmm0, 80(%rbx)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 88(%rbx)
	mulsd	%xmm7, %xmm0
	movapd	%xmm6, %xmm1
	mulsd	%xmm14, %xmm1
	mulsd	%xmm15, %xmm6
	subsd	%xmm1, %xmm0
	movapd	%xmm6, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L63
.L20:
	addsd	104(%rbx), %xmm1
	addsd	96(%rbx), %xmm0
	movapd	%xmm14, %xmm2
	movsd	40(%rsp), %xmm7
	mulsd	%xmm8, %xmm2
	movsd	%xmm0, 96(%rbx)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 104(%rbx)
	mulsd	%xmm8, %xmm0
	movapd	%xmm7, %xmm1
	mulsd	%xmm14, %xmm1
	mulsd	%xmm15, %xmm7
	subsd	%xmm1, %xmm0
	movapd	%xmm7, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L64
.L21:
	addsd	120(%rbx), %xmm1
	addsd	112(%rbx), %xmm0
	addq	%r11, %rsi
	subq	%r11, %r12
	addq	$8, %rax
	subq	$-128, %rbx
	subq	%r11, %r14
	subq	%r11, %r13
	subq	%r11, %rdi
	subq	%r11, %r8
	subq	%r11, %rbp
	subq	%r11, %r15
	movsd	%xmm0, -16(%rbx)
	movsd	%xmm1, -8(%rbx)
	cmpq	%rax, 64(%rsp)
	jg	.L12
	movq	80(%rsp), %rax
	movq	112(%rsp), %rdx
	movq	%r9, %rbp
	movq	120(%rsp), %r10
	movq	128(%rsp), %rcx
	subq	$8, %rax
	shrq	$3, %rax
	leaq	8(,%rax,8), %rax
.L11:
	cmpq	%rcx, %rax
	jge	.L24
	movq	%r10, %r15
	leaq	1(%rax), %r14
	leaq	2(%rax), %r13
	movq	96(%rsp), %rdi
	leaq	3(%rax), %r12
	movq	%rax, %rbx
	movq	88(%rsp), %rsi
	movq	%rdx, (%rsp)
	imulq	%rax, %r15
	salq	$4, %rbx
	movq	56(%rsp), %r11
	movq	%rax, %r8
	imulq	%r10, %r14
	addq	%rdi, %rbx
	movq	%rdx, %rdi
	movapd	.LC7(%rip), %xmm5
	imulq	%r10, %r13
	movq	.LC8(%rip), %xmm4
	salq	$6, %rdi
	imulq	%r10, %r12
	addq	%rsi, %r15
	addq	%rsi, %r14
	addq	%rsi, %r13
	addq	%rsi, %r12
	.p2align 4,,10
	.p2align 3
.L25:
	testq	%rbp, %rbp
	jle	.L65
	pxor	%xmm9, %xmm9
	movq	(%r11), %r9
	xorl	%edx, %edx
	xorl	%esi, %esi
	movapd	%xmm9, %xmm8
	movapd	%xmm9, %xmm7
	movapd	%xmm9, %xmm6
	.p2align 4,,10
	.p2align 3
.L26:
	movupd	(%r15,%rdx), %xmm0
	movupd	(%r9,%rdx), %xmm1
	addq	$1, %rsi
	xorpd	%xmm5, %xmm0
	pshufd	$78, %xmm1, %xmm2
	pshufd	$238, %xmm0, %xmm3
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm2, %xmm3
	mulpd	%xmm1, %xmm0
	xorpd	%xmm4, %xmm3
	addpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm6
	movupd	(%r14,%rdx), %xmm0
	xorpd	%xmm5, %xmm0
	pshufd	$238, %xmm0, %xmm3
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm2, %xmm3
	mulpd	%xmm1, %xmm0
	xorpd	%xmm4, %xmm3
	addpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm7
	movupd	0(%r13,%rdx), %xmm0
	xorpd	%xmm5, %xmm0
	pshufd	$238, %xmm0, %xmm3
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm2, %xmm3
	mulpd	%xmm1, %xmm0
	xorpd	%xmm4, %xmm3
	addpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm8
	movupd	(%r12,%rdx), %xmm0
	addq	$16, %rdx
	xorpd	%xmm5, %xmm0
	pshufd	$238, %xmm0, %xmm3
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm3, %xmm2
	mulpd	%xmm0, %xmm1
	xorpd	%xmm4, %xmm2
	addpd	%xmm2, %xmm1
	addpd	%xmm1, %xmm9
	cmpq	%rsi, %rbp
	jne	.L26
	movapd	%xmm6, %xmm2
	unpckhpd	%xmm2, %xmm2
	movapd	%xmm2, %xmm3
	movapd	%xmm7, %xmm2
	unpckhpd	%xmm2, %xmm2
	movapd	%xmm2, %xmm12
	movapd	%xmm8, %xmm2
	unpckhpd	%xmm2, %xmm2
	movapd	%xmm2, %xmm11
	movapd	%xmm9, %xmm2
	unpckhpd	%xmm2, %xmm2
	movapd	%xmm2, %xmm10
.L32:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm14, %xmm2
	mulsd	%xmm3, %xmm1
	mulsd	%xmm6, %xmm0
	mulsd	%xmm6, %xmm2
	subsd	%xmm1, %xmm0
	movapd	%xmm15, %xmm1
	mulsd	%xmm3, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L66
.L27:
	addsd	8(%rbx), %xmm1
	addsd	(%rbx), %xmm0
	movapd	%xmm14, %xmm2
	mulsd	%xmm7, %xmm2
	movsd	%xmm0, (%rbx)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 8(%rbx)
	mulsd	%xmm7, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm12, %xmm1
	subsd	%xmm1, %xmm0
	movapd	%xmm15, %xmm1
	mulsd	%xmm12, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L67
.L28:
	addsd	24(%rbx), %xmm1
	addsd	16(%rbx), %xmm0
	movapd	%xmm14, %xmm2
	mulsd	%xmm8, %xmm2
	movsd	%xmm0, 16(%rbx)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 24(%rbx)
	mulsd	%xmm8, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm11, %xmm1
	subsd	%xmm1, %xmm0
	movapd	%xmm15, %xmm1
	mulsd	%xmm11, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L68
.L29:
	addsd	40(%rbx), %xmm1
	addsd	32(%rbx), %xmm0
	movapd	%xmm14, %xmm2
	mulsd	%xmm9, %xmm2
	movsd	%xmm0, 32(%rbx)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 40(%rbx)
	mulsd	%xmm9, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm10, %xmm1
	subsd	%xmm1, %xmm0
	movapd	%xmm15, %xmm1
	mulsd	%xmm10, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L69
.L30:
	addsd	56(%rbx), %xmm1
	addsd	48(%rbx), %xmm0
	addq	%rdi, %r15
	addq	%rdi, %r14
	addq	$4, %r8
	addq	$64, %rbx
	addq	%rdi, %r13
	addq	%rdi, %r12
	movsd	%xmm0, -16(%rbx)
	movsd	%xmm1, -8(%rbx)
	cmpq	%rcx, %r8
	jl	.L25
	movq	80(%rsp), %rdi
	movq	(%rsp), %rdx
	leaq	-4(%rdi), %rcx
	subq	%rax, %rcx
	andq	$-4, %rcx
	leaq	4(%rax,%rcx), %rax
.L24:
	movq	104(%rsp), %rdi
	cmpq	%rdi, %rax
	jge	.L33
	movq	96(%rsp), %rsi
	movq	%rax, %r13
	movq	%r10, %r12
	leaq	1(%rax), %rbx
	imulq	%rax, %r12
	salq	$4, %r13
	salq	$5, %rdx
	movq	%rax, %r15
	imulq	%r10, %rbx
	addq	%rsi, %r13
	movq	88(%rsp), %rsi
	movapd	.LC7(%rip), %xmm5
	movq	.LC8(%rip), %xmm4
	movq	%rdx, %r14
	addq	%rsi, %r12
	addq	%rsi, %rbx
.L34:
	testq	%rbp, %rbp
	jle	.L70
	movq	56(%rsp), %rsi
	pxor	%xmm7, %xmm7
	xorl	%edx, %edx
	xorl	%ecx, %ecx
	movapd	%xmm7, %xmm6
	movq	(%rsi), %rsi
	.p2align 4,,10
	.p2align 3
.L35:
	movupd	(%r12,%rdx), %xmm0
	movupd	(%rsi,%rdx), %xmm1
	addq	$1, %rcx
	xorpd	%xmm5, %xmm0
	pshufd	$78, %xmm1, %xmm2
	pshufd	$238, %xmm0, %xmm3
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm2, %xmm3
	mulpd	%xmm1, %xmm0
	xorpd	%xmm4, %xmm3
	addpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm6
	movupd	(%rbx,%rdx), %xmm0
	addq	$16, %rdx
	xorpd	%xmm5, %xmm0
	pshufd	$238, %xmm0, %xmm3
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm3, %xmm2
	mulpd	%xmm0, %xmm1
	xorpd	%xmm4, %xmm2
	addpd	%xmm2, %xmm1
	addpd	%xmm1, %xmm7
	cmpq	%rcx, %rbp
	jne	.L35
	movapd	%xmm6, %xmm2
	unpckhpd	%xmm2, %xmm2
	movapd	%xmm2, %xmm3
	movapd	%xmm6, %xmm2
	movapd	%xmm7, %xmm6
	unpckhpd	%xmm6, %xmm6
.L39:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm14, %xmm8
	mulsd	%xmm3, %xmm1
	mulsd	%xmm2, %xmm0
	mulsd	%xmm2, %xmm8
	subsd	%xmm1, %xmm0
	movapd	%xmm15, %xmm1
	mulsd	%xmm3, %xmm1
	addsd	%xmm8, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L71
.L36:
	addsd	8(%r13), %xmm1
	addsd	0(%r13), %xmm0
	movapd	%xmm14, %xmm2
	mulsd	%xmm7, %xmm2
	movsd	%xmm0, 0(%r13)
	movapd	%xmm15, %xmm0
	movsd	%xmm1, 8(%r13)
	mulsd	%xmm7, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm6, %xmm1
	subsd	%xmm1, %xmm0
	movapd	%xmm15, %xmm1
	mulsd	%xmm6, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L72
.L37:
	addsd	24(%r13), %xmm1
	addq	$2, %r15
	addq	%r14, %r12
	addq	%r14, %rbx
	addsd	16(%r13), %xmm0
	addq	$32, %r13
	movsd	%xmm1, -8(%r13)
	movsd	%xmm0, -16(%r13)
	cmpq	%rdi, %r15
	jl	.L34
	movq	80(%rsp), %rdi
	leaq	-2(%rdi), %rdx
	subq	%rax, %rdx
	andq	$-2, %rdx
	leaq	2(%rax,%rdx), %rax
.L33:
	movq	80(%rsp), %r13
	cmpq	%rax, %r13
	jle	.L10
	movq	96(%rsp), %rdi
	movq	%rax, %r12
	movq	%r10, %rbx
	movq	56(%rsp), %r14
	salq	$4, %r12
	imulq	%rax, %rbx
	movq	%r10, %r15
	movapd	.LC7(%rip), %xmm5
	addq	%r12, %rdi
	movq	.LC8(%rip), %xmm4
	movq	%rdi, %r12
	movq	88(%rsp), %rdi
	addq	%rbx, %rdi
	movq	%rdi, %rbx
.L41:
	pxor	%xmm2, %xmm2
	movapd	%xmm2, %xmm3
	testq	%rbp, %rbp
	jle	.L45
	movq	(%r14), %rsi
	xorl	%edx, %edx
	pxor	%xmm2, %xmm2
	xorl	%ecx, %ecx
	.p2align 4,,10
	.p2align 3
.L42:
	movupd	(%rbx,%rdx), %xmm0
	movdqu	(%rsi,%rdx), %xmm7
	addq	$1, %rcx
	movupd	(%rsi,%rdx), %xmm6
	addq	$16, %rdx
	xorpd	%xmm5, %xmm0
	pshufd	$78, %xmm7, %xmm1
	pshufd	$238, %xmm0, %xmm3
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm3, %xmm1
	mulpd	%xmm6, %xmm0
	xorpd	%xmm4, %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm0, %xmm2
	cmpq	%rcx, %rbp
	jne	.L42
	movapd	%xmm2, %xmm7
	unpckhpd	%xmm7, %xmm7
	movapd	%xmm7, %xmm3
.L45:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm14, %xmm6
	mulsd	%xmm3, %xmm1
	mulsd	%xmm2, %xmm0
	mulsd	%xmm2, %xmm6
	subsd	%xmm1, %xmm0
	movapd	%xmm15, %xmm1
	mulsd	%xmm3, %xmm1
	addsd	%xmm6, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L73
.L43:
	addsd	8(%r12), %xmm1
	addq	$1, %rax
	addq	$16, %r12
	addq	%r15, %rbx
	addsd	-16(%r12), %xmm0
	movsd	%xmm1, -8(%r12)
	movsd	%xmm0, -16(%r12)
	cmpq	%rax, %r13
	jne	.L41
.L10:
	addq	$264, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L65:
	.cfi_restore_state
	pxor	%xmm9, %xmm9
	movapd	%xmm9, %xmm10
	movapd	%xmm9, %xmm8
	movapd	%xmm9, %xmm11
	movapd	%xmm9, %xmm7
	movapd	%xmm9, %xmm12
	movapd	%xmm9, %xmm6
	movapd	%xmm9, %xmm3
	jmp	.L32
	.p2align 4,,10
	.p2align 3
.L56:
	pxor	%xmm8, %xmm8
	movapd	%xmm8, %xmm7
	movapd	%xmm8, %xmm6
	movsd	%xmm8, 40(%rsp)
	movapd	%xmm8, %xmm11
	movsd	%xmm8, 32(%rsp)
	movapd	%xmm8, %xmm5
	movapd	%xmm8, %xmm2
	movapd	%xmm8, %xmm10
	movsd	%xmm8, (%rsp)
	movapd	%xmm8, %xmm3
	movapd	%xmm8, %xmm4
	movapd	%xmm8, %xmm9
	movsd	%xmm8, 24(%rsp)
	movsd	%xmm8, 16(%rsp)
	movsd	%xmm8, 8(%rsp)
	jmp	.L23
.L70:
	pxor	%xmm7, %xmm7
	movapd	%xmm7, %xmm6
	movapd	%xmm7, %xmm2
	movapd	%xmm7, %xmm3
	jmp	.L39
.L64:
	movsd	40(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm8, %xmm2
	movq	%r11, 144(%rsp)
	movq	%r9, 136(%rsp)
	movq	%rsi, 48(%rsp)
	movq	%rdi, 32(%rsp)
	movq	%r8, 24(%rsp)
	movq	%rax, 16(%rsp)
	movsd	%xmm14, 8(%rsp)
	movsd	%xmm15, (%rsp)
	call	__muldc3@PLT
	movq	48(%rsp), %rsi
	movq	32(%rsp), %rdi
	movq	144(%rsp), %r11
	movq	136(%rsp), %r9
	movq	24(%rsp), %r8
	movq	16(%rsp), %rax
	movq	.LC8(%rip), %xmm13
	movsd	8(%rsp), %xmm14
	movapd	.LC7(%rip), %xmm12
	movsd	(%rsp), %xmm15
	jmp	.L21
.L69:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm10, %xmm3
	movq	%r11, 72(%rsp)
	movapd	%xmm9, %xmm2
	movq	%rcx, 64(%rsp)
	movq	%rdi, 48(%rsp)
	movq	%r8, 40(%rsp)
	movq	%rax, 32(%rsp)
	movq	%r10, 24(%rsp)
	movsd	%xmm14, 16(%rsp)
	movsd	%xmm15, 8(%rsp)
	call	__muldc3@PLT
	movq	72(%rsp), %r11
	movq	64(%rsp), %rcx
	movq	48(%rsp), %rdi
	movq	40(%rsp), %r8
	movq	32(%rsp), %rax
	movq	24(%rsp), %r10
	movq	.LC8(%rip), %xmm4
	movapd	.LC7(%rip), %xmm5
	movsd	16(%rsp), %xmm14
	movsd	8(%rsp), %xmm15
	jmp	.L30
.L68:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm11, %xmm3
	movq	%r11, 120(%rsp)
	movapd	%xmm8, %xmm2
	movq	%rcx, 112(%rsp)
	movq	%rdi, 72(%rsp)
	movq	%r8, 64(%rsp)
	movq	%rax, 48(%rsp)
	movq	%r10, 40(%rsp)
	movsd	%xmm9, 32(%rsp)
	movsd	%xmm10, 24(%rsp)
	movsd	%xmm14, 16(%rsp)
	movsd	%xmm15, 8(%rsp)
	call	__muldc3@PLT
	movq	120(%rsp), %r11
	movq	112(%rsp), %rcx
	movq	72(%rsp), %rdi
	movq	64(%rsp), %r8
	movq	48(%rsp), %rax
	movq	40(%rsp), %r10
	movq	.LC8(%rip), %xmm4
	movapd	.LC7(%rip), %xmm5
	movsd	32(%rsp), %xmm9
	movsd	24(%rsp), %xmm10
	movsd	16(%rsp), %xmm14
	movsd	8(%rsp), %xmm15
	jmp	.L29
.L67:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm7, %xmm2
	movq	%r11, 136(%rsp)
	movapd	%xmm12, %xmm3
	movq	%rcx, 128(%rsp)
	movq	%rdi, 120(%rsp)
	movq	%r8, 112(%rsp)
	movq	%rax, 72(%rsp)
	movq	%r10, 64(%rsp)
	movsd	%xmm9, 48(%rsp)
	movsd	%xmm10, 40(%rsp)
	movsd	%xmm8, 32(%rsp)
	movsd	%xmm11, 24(%rsp)
	movsd	%xmm14, 16(%rsp)
	movsd	%xmm15, 8(%rsp)
	call	__muldc3@PLT
	movq	120(%rsp), %rdi
	movq	112(%rsp), %r8
	movq	136(%rsp), %r11
	movq	128(%rsp), %rcx
	movq	72(%rsp), %rax
	movq	64(%rsp), %r10
	movq	.LC8(%rip), %xmm4
	movapd	.LC7(%rip), %xmm5
	movsd	48(%rsp), %xmm9
	movsd	40(%rsp), %xmm10
	movsd	32(%rsp), %xmm8
	movsd	24(%rsp), %xmm11
	movsd	16(%rsp), %xmm14
	movsd	8(%rsp), %xmm15
	jmp	.L28
.L66:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm6, %xmm2
	movq	%r11, 152(%rsp)
	movq	%rcx, 144(%rsp)
	movq	%rdi, 136(%rsp)
	movq	%r8, 128(%rsp)
	movq	%rax, 120(%rsp)
	movq	%r10, 112(%rsp)
	movsd	%xmm9, 72(%rsp)
	movsd	%xmm10, 64(%rsp)
	movsd	%xmm8, 48(%rsp)
	movsd	%xmm11, 40(%rsp)
	movsd	%xmm7, 32(%rsp)
	movsd	%xmm12, 24(%rsp)
	movsd	%xmm14, 16(%rsp)
	movsd	%xmm15, 8(%rsp)
	call	__muldc3@PLT
	movq	152(%rsp), %r11
	movq	144(%rsp), %rcx
	movq	136(%rsp), %rdi
	movq	128(%rsp), %r8
	movq	120(%rsp), %rax
	movq	112(%rsp), %r10
	movq	.LC8(%rip), %xmm4
	movapd	.LC7(%rip), %xmm5
	movsd	72(%rsp), %xmm9
	movsd	64(%rsp), %xmm10
	movsd	48(%rsp), %xmm8
	movsd	40(%rsp), %xmm11
	movsd	32(%rsp), %xmm7
	movsd	24(%rsp), %xmm12
	movsd	16(%rsp), %xmm14
	movsd	8(%rsp), %xmm15
	jmp	.L27
.L57:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movq	%r11, 248(%rsp)
	movq	%r9, 240(%rsp)
	movq	%rsi, 232(%rsp)
	movq	%rdi, 224(%rsp)
	movq	%r8, 216(%rsp)
	movq	%rax, 208(%rsp)
	movsd	%xmm8, 200(%rsp)
	movsd	%xmm7, 192(%rsp)
	movsd	%xmm6, 184(%rsp)
	movsd	%xmm11, 176(%rsp)
	movsd	%xmm5, 168(%rsp)
	movsd	%xmm10, 160(%rsp)
	movsd	%xmm4, 152(%rsp)
	movsd	%xmm9, 144(%rsp)
	movsd	%xmm14, 136(%rsp)
	movsd	%xmm15, 48(%rsp)
	call	__muldc3@PLT
	movq	248(%rsp), %r11
	movq	240(%rsp), %r9
	movq	232(%rsp), %rsi
	movq	224(%rsp), %rdi
	movq	216(%rsp), %r8
	movsd	48(%rsp), %xmm15
	movq	208(%rsp), %rax
	movq	.LC8(%rip), %xmm13
	movapd	.LC7(%rip), %xmm12
	movsd	200(%rsp), %xmm8
	movsd	192(%rsp), %xmm7
	movsd	184(%rsp), %xmm6
	movsd	176(%rsp), %xmm11
	movsd	168(%rsp), %xmm5
	movsd	160(%rsp), %xmm10
	movsd	152(%rsp), %xmm4
	movsd	144(%rsp), %xmm9
	movsd	136(%rsp), %xmm14
	jmp	.L14
.L58:
	movsd	8(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm5, %xmm2
	movapd	%xmm15, %xmm0
	movq	%r11, 232(%rsp)
	movq	%r9, 224(%rsp)
	movq	%rsi, 216(%rsp)
	movq	%rdi, 208(%rsp)
	movq	%r8, 200(%rsp)
	movq	%rax, 192(%rsp)
	movsd	%xmm8, 184(%rsp)
	movsd	%xmm7, 176(%rsp)
	movsd	%xmm6, 168(%rsp)
	movsd	%xmm11, 160(%rsp)
	movsd	%xmm10, 152(%rsp)
	movsd	%xmm4, 144(%rsp)
	movsd	%xmm9, 136(%rsp)
	movsd	%xmm14, 48(%rsp)
	movsd	%xmm15, 8(%rsp)
	call	__muldc3@PLT
	movq	232(%rsp), %r11
	movq	224(%rsp), %r9
	movq	216(%rsp), %rsi
	movq	208(%rsp), %rdi
	movq	200(%rsp), %r8
	movsd	48(%rsp), %xmm14
	movq	192(%rsp), %rax
	movq	.LC8(%rip), %xmm13
	movapd	.LC7(%rip), %xmm12
	movsd	8(%rsp), %xmm15
	movsd	184(%rsp), %xmm8
	movsd	176(%rsp), %xmm7
	movsd	168(%rsp), %xmm6
	movsd	160(%rsp), %xmm11
	movsd	152(%rsp), %xmm10
	movsd	144(%rsp), %xmm4
	movsd	136(%rsp), %xmm9
	jmp	.L15
.L59:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm4, %xmm3
	movq	%r11, 216(%rsp)
	movapd	%xmm9, %xmm2
	movq	%r9, 208(%rsp)
	movq	%rsi, 200(%rsp)
	movq	%rdi, 192(%rsp)
	movq	%r8, 184(%rsp)
	movq	%rax, 176(%rsp)
	movsd	%xmm8, 168(%rsp)
	movsd	%xmm7, 160(%rsp)
	movsd	%xmm6, 152(%rsp)
	movsd	%xmm11, 144(%rsp)
	movsd	%xmm10, 136(%rsp)
	movsd	%xmm14, 48(%rsp)
	movsd	%xmm15, 8(%rsp)
	call	__muldc3@PLT
	movq	216(%rsp), %r11
	movq	208(%rsp), %r9
	movq	200(%rsp), %rsi
	movq	192(%rsp), %rdi
	movq	184(%rsp), %r8
	movsd	48(%rsp), %xmm14
	movq	176(%rsp), %rax
	movq	.LC8(%rip), %xmm13
	movapd	.LC7(%rip), %xmm12
	movsd	8(%rsp), %xmm15
	movsd	168(%rsp), %xmm8
	movsd	160(%rsp), %xmm7
	movsd	152(%rsp), %xmm6
	movsd	144(%rsp), %xmm11
	movsd	136(%rsp), %xmm10
	jmp	.L16
.L60:
	movsd	(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm10, %xmm2
	movq	%r11, 200(%rsp)
	movq	%r9, 192(%rsp)
	movq	%rsi, 184(%rsp)
	movq	%rdi, 176(%rsp)
	movq	%r8, 168(%rsp)
	movq	%rax, 160(%rsp)
	movsd	%xmm8, 152(%rsp)
	movsd	%xmm7, 144(%rsp)
	movsd	%xmm6, 136(%rsp)
	movsd	%xmm11, 48(%rsp)
	movsd	%xmm14, 8(%rsp)
	movsd	%xmm15, (%rsp)
	call	__muldc3@PLT
	movq	200(%rsp), %r11
	movq	192(%rsp), %r9
	movq	184(%rsp), %rsi
	movq	176(%rsp), %rdi
	movq	168(%rsp), %r8
	movsd	48(%rsp), %xmm11
	movq	160(%rsp), %rax
	movsd	(%rsp), %xmm15
	movq	.LC8(%rip), %xmm13
	movsd	8(%rsp), %xmm14
	movapd	.LC7(%rip), %xmm12
	movsd	152(%rsp), %xmm8
	movsd	144(%rsp), %xmm7
	movsd	136(%rsp), %xmm6
	jmp	.L17
.L61:
	movsd	16(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm11, %xmm2
	movq	%r11, 192(%rsp)
	movq	%r9, 184(%rsp)
	movq	%rsi, 176(%rsp)
	movq	%rdi, 168(%rsp)
	movq	%r8, 160(%rsp)
	movq	%rax, 152(%rsp)
	movsd	%xmm8, 144(%rsp)
	movsd	%xmm7, 136(%rsp)
	movsd	%xmm6, 48(%rsp)
	movsd	%xmm14, 8(%rsp)
	movsd	%xmm15, (%rsp)
	call	__muldc3@PLT
	movq	192(%rsp), %r11
	movq	184(%rsp), %r9
	movq	176(%rsp), %rsi
	movq	168(%rsp), %rdi
	movq	160(%rsp), %r8
	movsd	48(%rsp), %xmm6
	movq	152(%rsp), %rax
	movsd	(%rsp), %xmm15
	movq	.LC8(%rip), %xmm13
	movsd	8(%rsp), %xmm14
	movapd	.LC7(%rip), %xmm12
	movsd	144(%rsp), %xmm8
	movsd	136(%rsp), %xmm7
	jmp	.L18
.L62:
	movsd	24(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm6, %xmm2
	movapd	%xmm15, %xmm0
	movq	%r11, 176(%rsp)
	movq	%r9, 168(%rsp)
	movq	%rsi, 160(%rsp)
	movq	%rdi, 152(%rsp)
	movq	%r8, 144(%rsp)
	movq	%rax, 136(%rsp)
	movsd	%xmm8, 48(%rsp)
	movsd	%xmm7, 16(%rsp)
	movsd	%xmm14, 8(%rsp)
	movsd	%xmm15, (%rsp)
	call	__muldc3@PLT
	movq	176(%rsp), %r11
	movq	168(%rsp), %r9
	movq	160(%rsp), %rsi
	movq	152(%rsp), %rdi
	movq	144(%rsp), %r8
	movsd	48(%rsp), %xmm8
	movq	136(%rsp), %rax
	movsd	16(%rsp), %xmm7
	movq	.LC8(%rip), %xmm13
	movsd	8(%rsp), %xmm14
	movapd	.LC7(%rip), %xmm12
	movsd	(%rsp), %xmm15
	jmp	.L19
.L63:
	movsd	32(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm7, %xmm2
	movapd	%xmm15, %xmm0
	movq	%r11, 160(%rsp)
	movq	%r9, 152(%rsp)
	movq	%rsi, 144(%rsp)
	movq	%rdi, 136(%rsp)
	movq	%r8, 48(%rsp)
	movq	%rax, 24(%rsp)
	movsd	%xmm8, 16(%rsp)
	movsd	%xmm14, 8(%rsp)
	movsd	%xmm15, (%rsp)
	call	__muldc3@PLT
	movq	160(%rsp), %r11
	movq	152(%rsp), %r9
	movq	144(%rsp), %rsi
	movq	136(%rsp), %rdi
	movq	48(%rsp), %r8
	movq	24(%rsp), %rax
	movq	.LC8(%rip), %xmm13
	movsd	16(%rsp), %xmm8
	movapd	.LC7(%rip), %xmm12
	movsd	8(%rsp), %xmm14
	movsd	(%rsp), %xmm15
	jmp	.L20
.L72:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movapd	%xmm6, %xmm3
	movq	%rdi, 32(%rsp)
	movapd	%xmm7, %xmm2
	movq	%rax, 24(%rsp)
	movq	%r10, 16(%rsp)
	movsd	%xmm14, 8(%rsp)
	movsd	%xmm15, (%rsp)
	call	__muldc3@PLT
	movq	32(%rsp), %rdi
	movq	24(%rsp), %rax
	movq	16(%rsp), %r10
	movq	.LC8(%rip), %xmm4
	movapd	.LC7(%rip), %xmm5
	movsd	8(%rsp), %xmm14
	movsd	(%rsp), %xmm15
	jmp	.L37
.L71:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movq	%rdi, 48(%rsp)
	movq	%rax, 40(%rsp)
	movq	%r10, 32(%rsp)
	movsd	%xmm7, 24(%rsp)
	movsd	%xmm6, 16(%rsp)
	movsd	%xmm14, 8(%rsp)
	movsd	%xmm15, (%rsp)
	call	__muldc3@PLT
	movq	48(%rsp), %rdi
	movq	40(%rsp), %rax
	movq	32(%rsp), %r10
	movq	.LC8(%rip), %xmm4
	movapd	.LC7(%rip), %xmm5
	movsd	24(%rsp), %xmm7
	movsd	16(%rsp), %xmm6
	movsd	8(%rsp), %xmm14
	movsd	(%rsp), %xmm15
	jmp	.L36
.L73:
	movapd	%xmm14, %xmm1
	movapd	%xmm15, %xmm0
	movq	%rax, 16(%rsp)
	movsd	%xmm14, 8(%rsp)
	movsd	%xmm15, (%rsp)
	call	__muldc3@PLT
	movq	16(%rsp), %rax
	movq	.LC8(%rip), %xmm4
	movapd	.LC7(%rip), %xmm5
	movsd	8(%rsp), %xmm14
	movsd	(%rsp), %xmm15
	jmp	.L43
	.cfi_endproc
.LFE14843:
	.size	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb1ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_.constprop.0, .-_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb1ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_.constprop.0
	.align 2
	.p2align 4
	.type	_ZN5Eigen8internal13gemm_pack_lhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi1EEELi1ELi1ENS0_9Packet1cdELi1ELb0ELb0EEclEPS3_RKS5_llll.constprop.0, @function
_ZN5Eigen8internal13gemm_pack_lhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi1EEELi1ELi1ENS0_9Packet1cdELi1ELb0ELb0EEclEPS3_RKS5_llll.constprop.0:
.LFB14844:
	.cfi_startproc
	movq	%rdi, %r10
	movq	%rsi, %r8
#APP
# 2263 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#EIGEN PRODUCT PACK LHS
# 0 "" 2
#NO_APP
	testq	%rcx, %rcx
	jle	.L74
	xorl	%r11d, %r11d
	xorl	%r9d, %r9d
	testq	%rdx, %rdx
	jle	.L81
	.p2align 4,,10
	.p2align 3
.L76:
	movq	%r11, %rdi
	xorl	%esi, %esi
	salq	$4, %rdi
	addq	%r10, %rdi
	.p2align 4,,10
	.p2align 3
.L77:
	movq	8(%r8), %rax
	addq	$16, %rdi
	imulq	%r9, %rax
	addq	%rsi, %rax
	addq	$1, %rsi
	salq	$4, %rax
	addq	(%r8), %rax
	movupd	(%rax), %xmm0
	movaps	%xmm0, -16(%rdi)
	cmpq	%rdx, %rsi
	jne	.L77
	addq	$1, %r9
	addq	%rdx, %r11
	cmpq	%r9, %rcx
	jne	.L76
.L74:
	ret
.L81:
	ret
	.cfi_endproc
.LFE14844:
	.size	_ZN5Eigen8internal13gemm_pack_lhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi1EEELi1ELi1ENS0_9Packet1cdELi1ELb0ELb0EEclEPS3_RKS5_llll.constprop.0, .-_ZN5Eigen8internal13gemm_pack_lhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi1EEELi1ELi1ENS0_9Packet1cdELi1ELb0ELb0EEclEPS3_RKS5_llll.constprop.0
	.align 2
	.p2align 4
	.type	_ZN5Eigen8internal13gemm_pack_lhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi1ELi1ENS0_9Packet1cdELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0, @function
_ZN5Eigen8internal13gemm_pack_lhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi1ELi1ENS0_9Packet1cdELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0:
.LFB14845:
	.cfi_startproc
	movq	%rdi, %r10
	movq	%rsi, %r8
#APP
# 2107 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#EIGEN PRODUCT PACK LHS
# 0 "" 2
#NO_APP
	testq	%rcx, %rcx
	jle	.L82
	xorl	%r11d, %r11d
	xorl	%r9d, %r9d
	testq	%rdx, %rdx
	jle	.L89
	.p2align 4,,10
	.p2align 3
.L84:
	movq	%r11, %rdi
	xorl	%esi, %esi
	salq	$4, %rdi
	addq	%r10, %rdi
	.p2align 4,,10
	.p2align 3
.L85:
	movq	8(%r8), %rax
	addq	$16, %rdi
	imulq	%rsi, %rax
	addq	$1, %rsi
	addq	%r9, %rax
	salq	$4, %rax
	addq	(%r8), %rax
	movupd	(%rax), %xmm0
	movaps	%xmm0, -16(%rdi)
	cmpq	%rdx, %rsi
	jne	.L85
	addq	$1, %r9
	addq	%rdx, %r11
	cmpq	%rcx, %r9
	jne	.L84
.L82:
	ret
.L89:
	ret
	.cfi_endproc
.LFE14845:
	.size	_ZN5Eigen8internal13gemm_pack_lhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi1ELi1ENS0_9Packet1cdELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0, .-_ZN5Eigen8internal13gemm_pack_lhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi1ELi1ENS0_9Packet1cdELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0
	.align 2
	.p2align 4
	.type	_ZN5Eigen8internal13gemm_pack_rhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi4ELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0, @function
_ZN5Eigen8internal13gemm_pack_rhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi4ELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0:
.LFB14846:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movq	%rdx, %r10
	movq	%rcx, %rdx
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	%rdi, -32(%rsp)
	movq	%rsi, %rdi
#APP
# 2389 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#EIGEN PRODUCT PACK RHS COLMAJOR
# 0 "" 2
#NO_APP
	sarq	$63, %rcx
	movq	%rdx, %rbx
	shrq	$62, %rcx
	leaq	(%rdx,%rcx), %rax
	andl	$3, %eax
	subq	%rcx, %rax
	subq	%rax, %rbx
	testq	%rbx, %rbx
	jle	.L98
	movq	8(%rsi), %r14
	movq	(%rdi), %r15
	movq	%r10, %r11
	movq	%rdi, -16(%rsp)
	leaq	-4(,%r10,4), %rax
	xorl	%r9d, %r9d
	xorl	%r13d, %r13d
	xorl	%ebp, %ebp
	movq	%rax, -24(%rsp)
	movq	%r14, %r12
	movq	%r14, %rsi
	salq	$4, %r11
	movq	%rdx, -8(%rsp)
	salq	$6, %r12
	salq	$5, %rsi
	addq	%r15, %r11
	salq	$4, %r14
	.p2align 4,,10
	.p2align 3
.L94:
	leaq	(%r14,%r9), %rdi
	leaq	(%r14,%rsi), %r8
	testq	%r10, %r10
	jle	.L92
	movq	-32(%rsp), %rdx
	movq	%r13, %rax
	leaq	(%r15,%r9), %rcx
	salq	$4, %rax
	addq	%rdx, %rax
	.p2align 4,,10
	.p2align 3
.L93:
	movsd	(%rcx), %xmm0
	movq	%rcx, %rdx
	addq	$16, %rcx
	addq	$64, %rax
	subq	%r9, %rdx
	movsd	%xmm0, -64(%rax)
	movsd	-8(%rcx), %xmm0
	movsd	%xmm0, -56(%rax)
	movsd	(%rdx,%rdi), %xmm0
	movsd	%xmm0, -48(%rax)
	movsd	8(%rdx,%rdi), %xmm0
	movsd	%xmm0, -40(%rax)
	movsd	(%rdx,%rsi), %xmm0
	movsd	%xmm0, -32(%rax)
	movsd	8(%rdx,%rsi), %xmm0
	movsd	%xmm0, -24(%rax)
	movsd	(%rdx,%r8), %xmm0
	movsd	%xmm0, -16(%rax)
	movsd	8(%rdx,%r8), %xmm0
	movsd	%xmm0, -8(%rax)
	cmpq	%rcx, %r11
	jne	.L93
	movq	-24(%rsp), %rax
	leaq	4(%r13,%rax), %r13
.L92:
	addq	$4, %rbp
	addq	%r12, %r9
	addq	%r12, %rsi
	addq	%r12, %r11
	cmpq	%rbp, %rbx
	jg	.L94
	movq	-16(%rsp), %rdi
	movq	-8(%rsp), %rdx
.L91:
	cmpq	%rbx, %rdx
	jle	.L90
	movq	(%rdi), %r9
	movq	8(%rdi), %r8
	testq	%r10, %r10
	jle	.L90
	movq	-32(%rsp), %rcx
	salq	$4, %r13
	movq	%rbx, %rdi
	salq	$4, %r10
	imulq	%r8, %rdi
	addq	%r13, %rcx
	.p2align 4,,10
	.p2align 3
.L97:
	movq	%rdi, %rsi
	xorl	%eax, %eax
	salq	$4, %rsi
	addq	%r9, %rsi
	.p2align 4,,10
	.p2align 3
.L96:
	movsd	(%rsi,%rax), %xmm0
	movsd	%xmm0, (%rcx,%rax)
	movsd	8(%rsi,%rax), %xmm0
	movsd	%xmm0, 8(%rcx,%rax)
	addq	$16, %rax
	cmpq	%rax, %r10
	jne	.L96
	addq	$1, %rbx
	addq	%r8, %rdi
	addq	%r10, %rcx
	cmpq	%rbx, %rdx
	jne	.L97
.L90:
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
.L98:
	.cfi_restore_state
	xorl	%r13d, %r13d
	jmp	.L91
	.cfi_endproc
.LFE14846:
	.size	_ZN5Eigen8internal13gemm_pack_rhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi4ELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0, .-_ZN5Eigen8internal13gemm_pack_rhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi4ELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0
	.align 2
	.p2align 4
	.type	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi0EEELi0ELb0ES3_NS4_IS3_lLi1EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_.isra.0, @function
_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi0EEELi0ELb0ES3_NS4_IS3_lLi1EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_.isra.0:
.LFB14847:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	movq	%rcx, %r13
	leaq	-7(%rdi), %rcx
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$264, %rsp
	.cfi_def_cfa_offset 320
	movq	(%rdx), %rbx
	movq	8(%rdx), %rax
	movsd	%xmm0, 240(%rsp)
	movq	%rbx, 16(%rsp)
	leaq	-3(%rdi), %rbx
	movq	%rbx, 72(%rsp)
	leaq	-2(%rdi), %rbx
	movsd	%xmm1, 248(%rsp)
	movapd	240(%rsp), %xmm6
	movq	%rbx, 80(%rsp)
	leaq	-1(%rdi), %rbx
	movq	%rdi, 40(%rsp)
	movq	%rsi, 24(%rsp)
	movq	%r8, 32(%rsp)
	movq	%rax, 112(%rsp)
	movq	%rbx, 88(%rsp)
	movaps	%xmm6, (%rsp)
	cmpq	$127, %rsi
	jle	.L105
	salq	$4, %rax
	cmpq	$32000, %rax
	sbbq	%rax, %rax
	andl	$12, %eax
	addq	$4, %rax
	movq	%rax, 64(%rsp)
.L106:
	movq	112(%rsp), %rax
	movq	64(%rsp), %rsi
	movq	%rcx, 96(%rsp)
	xorl	%edi, %edi
	movq	.LC8(%rip), %xmm4
	pshufd	$78, (%rsp), %xmm12
	movdqa	%xmm12, %xmm8
	imulq	%rax, %rsi
	salq	$4, %rax
	movq	%rax, %rbx
	movq	%rsi, 56(%rsp)
	movq	40(%rsp), %rsi
	leaq	-8(%rsi), %rax
	salq	$4, %rsi
	shrq	$3, %rax
	movq	%rsi, %r9
	movq	%rdi, %rsi
	leaq	8(,%rax,8), %rax
	movq	%rax, 104(%rsp)
	movq	32(%rsp), %rax
	addq	%rax, %r9
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L149:
	movq	64(%rsp), %rcx
	movq	%rsi, 48(%rsp)
	movq	%rsi, %r14
	addq	%rcx, %rsi
	movq	24(%rsp), %rcx
	cmpq	%rcx, %rsi
	movq	%rcx, %rbp
	movq	96(%rsp), %rcx
	cmovle	%rsi, %rbp
	testq	%rcx, %rcx
	jle	.L153
	movq	16(%rsp), %rdi
	movq	%rax, %r11
	movq	32(%rsp), %r8
	xorl	%r12d, %r12d
	salq	$4, %r11
	addq	%rdi, %r11
	.p2align 4,,10
	.p2align 3
.L109:
	cmpq	%rbp, %r14
	jge	.L177
	movq	8(%r13), %rdi
	movq	0(%r13), %r15
	cmpq	$1, %rdi
	jne	.L178
	pxor	%xmm11, %xmm11
	movq	%r11, %rdx
	movq	%r14, %rdi
	movapd	%xmm11, %xmm10
	movapd	%xmm11, %xmm9
	movapd	%xmm11, %xmm7
	movapd	%xmm11, %xmm6
	movapd	%xmm11, %xmm5
	movapd	%xmm11, %xmm3
	movapd	%xmm11, %xmm2
	.p2align 4,,10
	.p2align 3
.L110:
	movq	%rdi, %r10
	movdqu	(%rdx), %xmm15
	addq	$1, %rdi
	salq	$4, %r10
	movupd	(%r15,%r10), %xmm0
	pshufd	$238, %xmm15, %xmm14
	pshufd	$68, %xmm15, %xmm13
	movdqu	16(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	pshufd	$78, %xmm0, %xmm1
	mulpd	%xmm1, %xmm14
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm2
	pshufd	$68, %xmm15, %xmm13
	movdqu	32(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm3
	pshufd	$68, %xmm15, %xmm13
	movdqu	48(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm5
	pshufd	$68, %xmm15, %xmm13
	movdqu	64(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm6
	pshufd	$68, %xmm15, %xmm13
	movdqu	80(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm7
	pshufd	$68, %xmm15, %xmm13
	movdqu	96(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm9
	pshufd	$68, %xmm15, %xmm13
	movdqu	112(%rdx), %xmm15
	addq	%rbx, %rdx
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	addpd	%xmm13, %xmm10
	pshufd	$238, %xmm15, %xmm13
	mulpd	%xmm13, %xmm1
	pshufd	$68, %xmm15, %xmm13
	mulpd	%xmm13, %xmm0
	xorpd	%xmm4, %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm0, %xmm11
	cmpq	%rdi, %rbp
	jne	.L110
.L114:
	pshufd	$238, %xmm2, %xmm0
	addq	$8, %r12
	subq	$-128, %r8
	subq	$-128, %r11
	mulpd	%xmm12, %xmm0
	pshufd	$68, %xmm2, %xmm2
	movupd	-128(%r8), %xmm1
	mulpd	(%rsp), %xmm2
	xorpd	%xmm4, %xmm0
	addpd	%xmm0, %xmm2
	pshufd	$238, %xmm3, %xmm0
	pshufd	$68, %xmm3, %xmm3
	mulpd	%xmm12, %xmm0
	mulpd	(%rsp), %xmm3
	addpd	%xmm1, %xmm2
	pshufd	$238, %xmm6, %xmm1
	mulpd	%xmm12, %xmm1
	xorpd	%xmm4, %xmm0
	movups	%xmm2, -128(%r8)
	movupd	-112(%r8), %xmm2
	addpd	%xmm0, %xmm3
	pshufd	$238, %xmm5, %xmm0
	pshufd	$68, %xmm5, %xmm5
	mulpd	%xmm12, %xmm0
	mulpd	(%rsp), %xmm5
	xorpd	%xmm4, %xmm1
	addpd	%xmm2, %xmm3
	pshufd	$68, %xmm9, %xmm2
	mulpd	(%rsp), %xmm2
	xorpd	%xmm4, %xmm0
	movups	%xmm3, -112(%r8)
	movupd	-96(%r8), %xmm3
	addpd	%xmm0, %xmm5
	pshufd	$68, %xmm6, %xmm0
	movupd	-80(%r8), %xmm6
	mulpd	(%rsp), %xmm0
	addpd	%xmm3, %xmm5
	pshufd	$68, %xmm7, %xmm3
	mulpd	(%rsp), %xmm3
	addpd	%xmm1, %xmm0
	pshufd	$68, %xmm10, %xmm1
	movups	%xmm5, -96(%r8)
	movupd	-48(%r8), %xmm5
	mulpd	(%rsp), %xmm1
	addpd	%xmm6, %xmm0
	movupd	-64(%r8), %xmm6
	movups	%xmm0, -80(%r8)
	pshufd	$238, %xmm7, %xmm0
	mulpd	%xmm12, %xmm0
	xorpd	%xmm4, %xmm0
	addpd	%xmm0, %xmm3
	pshufd	$238, %xmm9, %xmm0
	mulpd	%xmm12, %xmm0
	addpd	%xmm6, %xmm3
	xorpd	%xmm4, %xmm0
	movups	%xmm3, -64(%r8)
	addpd	%xmm0, %xmm2
	pshufd	$238, %xmm10, %xmm0
	mulpd	%xmm12, %xmm0
	addpd	%xmm5, %xmm2
	movupd	-32(%r8), %xmm5
	xorpd	%xmm4, %xmm0
	movups	%xmm2, -48(%r8)
	addpd	%xmm0, %xmm1
	pshufd	$68, %xmm11, %xmm0
	mulpd	(%rsp), %xmm0
	addpd	%xmm5, %xmm1
	movupd	-16(%r8), %xmm5
	movups	%xmm1, -32(%r8)
	pshufd	$238, %xmm11, %xmm1
	mulpd	%xmm12, %xmm1
	xorpd	%xmm4, %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm5, %xmm0
	movups	%xmm0, -16(%r8)
	cmpq	%r12, %rcx
	jg	.L109
	movq	104(%rsp), %rdi
.L108:
	cmpq	%rdi, 72(%rsp)
	jg	.L179
.L115:
	cmpq	%rdi, 80(%rsp)
	jg	.L180
.L121:
	cmpq	%rdi, 88(%rsp)
	jg	.L181
.L127:
	cmpq	%rdi, 40(%rsp)
	jg	.L182
.L133:
	movq	56(%rsp), %rcx
	addq	%rcx, %rax
	movq	24(%rsp), %rcx
	cmpq	%rcx, %rsi
	jl	.L149
.L104:
	addq	$264, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L178:
	.cfi_restore_state
	movq	%rdi, %r10
	imulq	%r14, %rdi
	pxor	%xmm11, %xmm11
	movq	%r11, %rdx
	salq	$4, %r10
	movapd	%xmm11, %xmm10
	movapd	%xmm11, %xmm9
	movapd	%xmm11, %xmm7
	movapd	%xmm11, %xmm6
	movapd	%xmm11, %xmm5
	salq	$4, %rdi
	movapd	%xmm11, %xmm3
	movapd	%xmm11, %xmm2
	addq	%r15, %rdi
	movq	%r14, %r15
	.p2align 4,,10
	.p2align 3
.L111:
	movupd	(%rdi), %xmm0
	movdqu	(%rdx), %xmm15
	addq	$1, %r15
	addq	%r10, %rdi
	pshufd	$78, %xmm0, %xmm1
	pshufd	$238, %xmm15, %xmm14
	pshufd	$68, %xmm15, %xmm13
	movdqu	16(%rdx), %xmm15
	mulpd	%xmm1, %xmm14
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm2
	pshufd	$68, %xmm15, %xmm13
	movdqu	32(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm3
	pshufd	$68, %xmm15, %xmm13
	movdqu	48(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm5
	pshufd	$68, %xmm15, %xmm13
	movdqu	64(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm6
	pshufd	$68, %xmm15, %xmm13
	movdqu	80(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm7
	pshufd	$68, %xmm15, %xmm13
	movdqu	96(%rdx), %xmm15
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	pshufd	$238, %xmm15, %xmm14
	mulpd	%xmm1, %xmm14
	addpd	%xmm13, %xmm9
	pshufd	$68, %xmm15, %xmm13
	movdqu	112(%rdx), %xmm15
	addq	%rbx, %rdx
	mulpd	%xmm0, %xmm13
	xorpd	%xmm4, %xmm14
	addpd	%xmm14, %xmm13
	addpd	%xmm13, %xmm10
	pshufd	$238, %xmm15, %xmm13
	mulpd	%xmm13, %xmm1
	pshufd	$68, %xmm15, %xmm13
	mulpd	%xmm13, %xmm0
	xorpd	%xmm4, %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm0, %xmm11
	cmpq	%rbp, %r15
	jne	.L111
	jmp	.L114
	.p2align 4,,10
	.p2align 3
.L177:
	pxor	%xmm11, %xmm11
	movdqa	%xmm11, %xmm10
	movdqa	%xmm11, %xmm9
	movdqa	%xmm11, %xmm7
	movdqa	%xmm11, %xmm6
	movdqa	%xmm11, %xmm5
	movdqa	%xmm11, %xmm3
	movdqa	%xmm11, %xmm2
	jmp	.L114
.L182:
	cmpq	%rbp, %r14
	jge	.L159
	movq	8(%r13), %rdx
	movq	0(%r13), %r11
	cmpq	$1, %rdx
	jne	.L160
	cmpq	$1, 112(%rsp)
	jne	.L160
	movq	16(%rsp), %rcx
	movq	%rdi, %r10
	movq	%r14, %rdx
	pxor	%xmm0, %xmm0
	salq	$4, %r10
	movapd	%xmm4, %xmm11
	addq	%rcx, %r10
	.p2align 4,,10
	.p2align 3
.L139:
	movq	%rdx, %r8
	addq	$1, %rdx
	salq	$4, %r8
	movupd	(%r11,%r8), %xmm2
	movupd	(%r10,%r8), %xmm5
	pshufd	$78, %xmm2, %xmm6
	pshufd	$238, %xmm5, %xmm3
	pshufd	$68, %xmm5, %xmm5
	mulpd	%xmm6, %xmm3
	mulpd	%xmm5, %xmm2
	xorpd	%xmm4, %xmm3
	addpd	%xmm3, %xmm2
	addpd	%xmm2, %xmm0
	cmpq	%rdx, %rbp
	jne	.L139
.L134:
	pshufd	$238, %xmm0, %xmm2
	pshufd	$68, %xmm0, %xmm0
	movq	32(%rsp), %rcx
	movq	%rdi, %rdx
	mulpd	%xmm8, %xmm2
	movapd	%xmm11, %xmm1
	salq	$4, %rdx
	mulpd	(%rsp), %xmm0
	leaq	(%rcx,%rdx), %r8
	movupd	(%r8), %xmm6
	xorpd	%xmm2, %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm6, %xmm0
	movups	%xmm0, (%r8)
	leaq	1(%rdi), %r8
	cmpq	%r8, 40(%rsp)
	jle	.L133
	leaq	1(%rdi,%rax), %r15
	leaq	16(%rcx,%rdx), %r12
	movq	16(%rsp), %rcx
	movq	%rbp, %rdi
	salq	$4, %r15
	movsd	8(%rsp), %xmm3
	movsd	(%rsp), %xmm2
	addq	%rcx, %r15
	movq	48(%rsp), %rcx
	movq	%rsi, 48(%rsp)
	movapd	%xmm3, %xmm6
	movapd	%xmm2, %xmm5
	subq	%rcx, %rdi
	salq	$4, %rcx
	unpcklpd	%xmm6, %xmm6
	unpcklpd	%xmm5, %xmm5
	movq	%rdi, %r11
	movq	%rcx, %r10
	salq	$4, %r11
	.p2align 4,,10
	.p2align 3
.L140:
	pxor	%xmm9, %xmm9
	movapd	%xmm9, %xmm1
	cmpq	%rbp, %r14
	jge	.L148
	movq	8(%r13), %rcx
	movq	0(%r13), %rdx
	cmpq	$1, %rcx
	jne	.L183
	addq	%r10, %rdx
	movq	%r15, %rcx
	pxor	%xmm9, %xmm9
	leaq	(%r11,%rdx), %rsi
	.p2align 4,,10
	.p2align 3
.L144:
	movupd	(%rcx), %xmm7
	movupd	(%rdx), %xmm1
	addq	$16, %rdx
	addq	%rbx, %rcx
	movapd	%xmm7, %xmm0
	unpcklpd	%xmm1, %xmm1
	shufpd	$1, %xmm7, %xmm0
	mulpd	%xmm0, %xmm1
	movupd	-16(%rdx), %xmm0
	unpckhpd	%xmm0, %xmm0
	mulpd	%xmm7, %xmm0
	movapd	%xmm0, %xmm7
	addpd	%xmm1, %xmm7
	subpd	%xmm0, %xmm1
	movsd	%xmm7, %xmm1
	addpd	%xmm1, %xmm9
	cmpq	%rsi, %rdx
	jne	.L144
.L176:
	movapd	%xmm9, %xmm1
	shufpd	$1, %xmm9, %xmm1
.L148:
	mulpd	%xmm6, %xmm9
	movapd	%xmm1, %xmm0
	mulpd	%xmm5, %xmm0
	movapd	%xmm0, %xmm7
	subpd	%xmm9, %xmm0
	addpd	%xmm9, %xmm7
	movapd	%xmm7, %xmm9
	unpckhpd	%xmm7, %xmm7
	ucomisd	%xmm7, %xmm0
	movsd	%xmm0, %xmm9
	jp	.L184
	movupd	(%r12), %xmm0
	addq	$16, %r12
	addq	$16, %r15
	addpd	%xmm9, %xmm0
	movups	%xmm0, -16(%r12)
	cmpq	%r12, %r9
	jne	.L140
.L174:
	movq	48(%rsp), %rsi
	jmp	.L133
	.p2align 4,,10
	.p2align 3
.L183:
	movq	%rcx, %r8
	imulq	%r14, %rcx
	pxor	%xmm9, %xmm9
	movq	%r15, %rsi
	salq	$4, %r8
	salq	$4, %rcx
	addq	%rcx, %rdx
	xorl	%ecx, %ecx
	.p2align 4,,10
	.p2align 3
.L142:
	movupd	(%rsi), %xmm7
	movupd	(%rdx), %xmm1
	addq	$1, %rcx
	addq	%rbx, %rsi
	movapd	%xmm7, %xmm0
	unpcklpd	%xmm1, %xmm1
	shufpd	$1, %xmm7, %xmm0
	mulpd	%xmm0, %xmm1
	movupd	(%rdx), %xmm0
	addq	%r8, %rdx
	unpckhpd	%xmm0, %xmm0
	mulpd	%xmm7, %xmm0
	movapd	%xmm0, %xmm7
	addpd	%xmm1, %xmm7
	subpd	%xmm0, %xmm1
	movsd	%xmm7, %xmm1
	addpd	%xmm1, %xmm9
	cmpq	%rcx, %rdi
	jne	.L142
	jmp	.L176
.L181:
	cmpq	%rbp, %r14
	jge	.L158
	movq	8(%r13), %rdx
	movq	0(%r13), %r11
	cmpq	$1, %rdx
	jne	.L185
	movq	16(%rsp), %rcx
	leaq	(%rax,%rdi), %rdx
	pxor	%xmm6, %xmm6
	movq	%r14, %r8
	salq	$4, %rdx
	movapd	%xmm6, %xmm7
	movapd	%xmm4, %xmm11
	addq	%rcx, %rdx
	.p2align 4,,10
	.p2align 3
.L132:
	movq	%r8, %r10
	movdqu	(%rdx), %xmm5
	movdqu	(%rdx), %xmm1
	addq	$1, %r8
	salq	$4, %r10
	movupd	(%r11,%r10), %xmm2
	pshufd	$238, %xmm5, %xmm5
	pshufd	$68, %xmm1, %xmm0
	mulpd	%xmm2, %xmm0
	pshufd	$78, %xmm2, %xmm3
	mulpd	%xmm3, %xmm5
	xorpd	%xmm4, %xmm5
	addpd	%xmm5, %xmm0
	movdqu	16(%rdx), %xmm5
	addq	%rbx, %rdx
	addpd	%xmm0, %xmm7
	pshufd	$238, %xmm5, %xmm0
	mulpd	%xmm0, %xmm3
	pshufd	$68, %xmm5, %xmm0
	mulpd	%xmm0, %xmm2
	xorpd	%xmm4, %xmm3
	addpd	%xmm3, %xmm2
	addpd	%xmm2, %xmm6
	cmpq	%r8, %rbp
	jne	.L132
.L128:
	pshufd	$238, %xmm7, %xmm2
	pshufd	$68, %xmm7, %xmm0
	movq	32(%rsp), %rcx
	movq	%rdi, %rdx
	mulpd	%xmm8, %xmm2
	salq	$4, %rdx
	movapd	%xmm11, %xmm1
	addq	$2, %rdi
	mulpd	(%rsp), %xmm0
	leaq	(%rcx,%rdx), %r8
	leaq	16(%rcx,%rdx), %rdx
	movupd	(%r8), %xmm5
	xorpd	%xmm11, %xmm2
	addpd	%xmm2, %xmm0
	addpd	%xmm5, %xmm0
	movups	%xmm0, (%r8)
	pshufd	$238, %xmm6, %xmm0
	mulpd	%xmm8, %xmm0
	xorpd	%xmm0, %xmm1
	pshufd	$68, %xmm6, %xmm0
	movupd	(%rdx), %xmm6
	mulpd	(%rsp), %xmm0
	addpd	%xmm1, %xmm0
	addpd	%xmm6, %xmm0
	movups	%xmm0, (%rdx)
	jmp	.L127
.L180:
	cmpq	%rbp, %r14
	jge	.L157
	movq	8(%r13), %rdx
	movq	0(%r13), %r11
	cmpq	$1, %rdx
	jne	.L186
	movq	16(%rsp), %rcx
	leaq	(%rax,%rdi), %rdx
	pxor	%xmm3, %xmm3
	movq	%r14, %r8
	salq	$4, %rdx
	movapd	%xmm3, %xmm5
	movapd	%xmm3, %xmm6
	addq	%rcx, %rdx
	movapd	%xmm4, %xmm11
	.p2align 4,,10
	.p2align 3
.L126:
	movq	%r8, %r10
	movdqu	(%rdx), %xmm2
	movdqu	16(%rdx), %xmm1
	addq	$1, %r8
	salq	$4, %r10
	movupd	(%r11,%r10), %xmm0
	pshufd	$238, %xmm2, %xmm9
	pshufd	$68, %xmm2, %xmm2
	pshufd	$68, %xmm1, %xmm1
	mulpd	%xmm0, %xmm2
	pshufd	$78, %xmm0, %xmm7
	mulpd	%xmm7, %xmm9
	mulpd	%xmm0, %xmm1
	xorpd	%xmm4, %xmm9
	addpd	%xmm9, %xmm2
	addpd	%xmm2, %xmm6
	movdqu	16(%rdx), %xmm2
	pshufd	$238, %xmm2, %xmm2
	mulpd	%xmm7, %xmm2
	xorpd	%xmm4, %xmm2
	addpd	%xmm2, %xmm1
	movdqu	32(%rdx), %xmm2
	pshufd	$238, %xmm2, %xmm2
	mulpd	%xmm2, %xmm7
	movdqu	32(%rdx), %xmm2
	addpd	%xmm1, %xmm5
	addq	%rbx, %rdx
	pshufd	$68, %xmm2, %xmm2
	mulpd	%xmm2, %xmm0
	xorpd	%xmm4, %xmm7
	addpd	%xmm7, %xmm0
	addpd	%xmm0, %xmm3
	cmpq	%r8, %rbp
	jne	.L126
.L122:
	pshufd	$238, %xmm6, %xmm2
	pshufd	$68, %xmm6, %xmm0
	movq	32(%rsp), %rcx
	movq	%rdi, %rdx
	mulpd	%xmm8, %xmm2
	salq	$4, %rdx
	movapd	%xmm11, %xmm1
	addq	$3, %rdi
	mulpd	(%rsp), %xmm0
	leaq	(%rcx,%rdx), %r8
	movupd	(%r8), %xmm6
	xorpd	%xmm11, %xmm2
	addpd	%xmm2, %xmm0
	pshufd	$238, %xmm5, %xmm2
	mulpd	%xmm8, %xmm2
	addpd	%xmm6, %xmm0
	movups	%xmm0, (%r8)
	pshufd	$68, %xmm5, %xmm0
	xorpd	%xmm11, %xmm2
	leaq	16(%rcx,%rdx), %r8
	mulpd	(%rsp), %xmm0
	movupd	(%r8), %xmm6
	leaq	32(%rcx,%rdx), %rdx
	addpd	%xmm2, %xmm0
	addpd	%xmm6, %xmm0
	movups	%xmm0, (%r8)
	pshufd	$238, %xmm3, %xmm0
	movupd	(%rdx), %xmm6
	mulpd	%xmm8, %xmm0
	xorpd	%xmm0, %xmm1
	pshufd	$68, %xmm3, %xmm0
	mulpd	(%rsp), %xmm0
	addpd	%xmm1, %xmm0
	addpd	%xmm6, %xmm0
	movups	%xmm0, (%rdx)
	jmp	.L121
.L179:
	cmpq	%rbp, %r14
	jge	.L156
	movq	8(%r13), %rdx
	movq	0(%r13), %r11
	cmpq	$1, %rdx
	jne	.L187
	movq	16(%rsp), %rcx
	leaq	(%rax,%rdi), %rdx
	pxor	%xmm5, %xmm5
	movq	%r14, %r8
	salq	$4, %rdx
	movapd	%xmm5, %xmm6
	movapd	%xmm5, %xmm7
	addq	%rcx, %rdx
	movapd	%xmm5, %xmm9
	movapd	%xmm4, %xmm11
	.p2align 4,,10
	.p2align 3
.L120:
	movq	%r8, %r10
	movdqu	(%rdx), %xmm3
	addq	$1, %r8
	salq	$4, %r10
	movupd	(%r11,%r10), %xmm0
	pshufd	$238, %xmm3, %xmm10
	pshufd	$68, %xmm3, %xmm3
	mulpd	%xmm0, %xmm3
	pshufd	$78, %xmm0, %xmm2
	mulpd	%xmm2, %xmm10
	xorpd	%xmm4, %xmm10
	addpd	%xmm10, %xmm3
	addpd	%xmm3, %xmm9
	movdqu	16(%rdx), %xmm3
	pshufd	$238, %xmm3, %xmm10
	pshufd	$68, %xmm3, %xmm3
	mulpd	%xmm2, %xmm10
	mulpd	%xmm0, %xmm3
	xorpd	%xmm4, %xmm10
	addpd	%xmm10, %xmm3
	addpd	%xmm3, %xmm7
	movdqu	32(%rdx), %xmm3
	pshufd	$238, %xmm3, %xmm10
	pshufd	$68, %xmm3, %xmm3
	mulpd	%xmm2, %xmm10
	mulpd	%xmm0, %xmm3
	xorpd	%xmm4, %xmm10
	addpd	%xmm10, %xmm3
	addpd	%xmm3, %xmm6
	movdqu	48(%rdx), %xmm3
	pshufd	$238, %xmm3, %xmm3
	mulpd	%xmm3, %xmm2
	movdqu	48(%rdx), %xmm3
	addq	%rbx, %rdx
	pshufd	$68, %xmm3, %xmm3
	mulpd	%xmm3, %xmm0
	xorpd	%xmm4, %xmm2
	addpd	%xmm2, %xmm0
	addpd	%xmm0, %xmm5
	cmpq	%r8, %rbp
	jne	.L120
.L116:
	pshufd	$238, %xmm9, %xmm2
	movq	32(%rsp), %rcx
	movq	%rdi, %rdx
	movapd	%xmm11, %xmm1
	mulpd	%xmm8, %xmm2
	pshufd	$68, %xmm9, %xmm0
	salq	$4, %rdx
	addq	$4, %rdi
	mulpd	(%rsp), %xmm0
	leaq	(%rcx,%rdx), %r8
	movupd	(%r8), %xmm3
	xorpd	%xmm11, %xmm2
	addpd	%xmm2, %xmm0
	pshufd	$238, %xmm7, %xmm2
	mulpd	%xmm8, %xmm2
	addpd	%xmm3, %xmm0
	movups	%xmm0, (%r8)
	pshufd	$68, %xmm7, %xmm0
	xorpd	%xmm11, %xmm2
	leaq	16(%rcx,%rdx), %r8
	mulpd	(%rsp), %xmm0
	movupd	(%r8), %xmm3
	addpd	%xmm2, %xmm0
	pshufd	$238, %xmm6, %xmm2
	mulpd	%xmm8, %xmm2
	addpd	%xmm3, %xmm0
	movups	%xmm0, (%r8)
	pshufd	$68, %xmm6, %xmm0
	xorpd	%xmm11, %xmm2
	leaq	32(%rcx,%rdx), %r8
	mulpd	(%rsp), %xmm0
	movupd	(%r8), %xmm6
	leaq	48(%rcx,%rdx), %rdx
	addpd	%xmm2, %xmm0
	addpd	%xmm6, %xmm0
	movups	%xmm0, (%r8)
	pshufd	$238, %xmm5, %xmm0
	movupd	(%rdx), %xmm6
	mulpd	%xmm8, %xmm0
	xorpd	%xmm0, %xmm1
	pshufd	$68, %xmm5, %xmm0
	mulpd	(%rsp), %xmm0
	addpd	%xmm1, %xmm0
	addpd	%xmm6, %xmm0
	movups	%xmm0, (%rdx)
	jmp	.L115
.L153:
	xorl	%edi, %edi
	jmp	.L108
.L105:
	movq	%rsi, 64(%rsp)
	testq	%rsi, %rsi
	jg	.L106
	jmp	.L104
.L187:
	movq	%rdx, %r12
	imulq	%r14, %rdx
	movq	16(%rsp), %rcx
	pxor	%xmm5, %xmm5
	salq	$4, %r12
	movapd	%xmm5, %xmm6
	movapd	%xmm5, %xmm7
	movq	%r14, %r10
	movq	.LC8(%rip), %xmm11
	movapd	%xmm5, %xmm9
	movq	%rdx, %r8
	leaq	(%rax,%rdi), %rdx
	salq	$4, %r8
	salq	$4, %rdx
	addq	%r11, %r8
	addq	%rcx, %rdx
.L118:
	movupd	(%r8), %xmm0
	movdqu	(%rdx), %xmm3
	addq	$1, %r10
	addq	%r12, %r8
	pshufd	$78, %xmm0, %xmm2
	pshufd	$238, %xmm3, %xmm10
	pshufd	$68, %xmm3, %xmm3
	mulpd	%xmm2, %xmm10
	mulpd	%xmm0, %xmm3
	xorpd	%xmm11, %xmm10
	addpd	%xmm10, %xmm3
	addpd	%xmm3, %xmm9
	movdqu	16(%rdx), %xmm3
	pshufd	$238, %xmm3, %xmm10
	pshufd	$68, %xmm3, %xmm3
	mulpd	%xmm2, %xmm10
	mulpd	%xmm0, %xmm3
	xorpd	%xmm11, %xmm10
	addpd	%xmm10, %xmm3
	addpd	%xmm3, %xmm7
	movdqu	32(%rdx), %xmm3
	pshufd	$238, %xmm3, %xmm10
	pshufd	$68, %xmm3, %xmm3
	mulpd	%xmm2, %xmm10
	mulpd	%xmm0, %xmm3
	xorpd	%xmm11, %xmm10
	addpd	%xmm10, %xmm3
	addpd	%xmm3, %xmm6
	movdqu	48(%rdx), %xmm3
	pshufd	$238, %xmm3, %xmm3
	mulpd	%xmm3, %xmm2
	movdqu	48(%rdx), %xmm3
	addq	%rbx, %rdx
	pshufd	$68, %xmm3, %xmm3
	mulpd	%xmm3, %xmm0
	xorpd	%xmm11, %xmm2
	addpd	%xmm2, %xmm0
	addpd	%xmm0, %xmm5
	cmpq	%rbp, %r10
	jne	.L118
	jmp	.L116
.L186:
	movq	%rdx, %r12
	imulq	%r14, %rdx
	movq	16(%rsp), %rcx
	pxor	%xmm3, %xmm3
	salq	$4, %r12
	movapd	%xmm3, %xmm5
	movapd	%xmm3, %xmm6
	movq	%r14, %r10
	movq	.LC8(%rip), %xmm11
	movq	%rdx, %r8
	leaq	(%rax,%rdi), %rdx
	salq	$4, %r8
	salq	$4, %rdx
	addq	%r11, %r8
	addq	%rcx, %rdx
.L124:
	movupd	(%r8), %xmm0
	movdqu	(%rdx), %xmm2
	addq	$1, %r10
	addq	%r12, %r8
	pshufd	$78, %xmm0, %xmm7
	pshufd	$238, %xmm2, %xmm9
	pshufd	$68, %xmm2, %xmm2
	mulpd	%xmm7, %xmm9
	mulpd	%xmm0, %xmm2
	xorpd	%xmm11, %xmm9
	addpd	%xmm9, %xmm2
	addpd	%xmm2, %xmm6
	movdqu	16(%rdx), %xmm2
	pshufd	$238, %xmm2, %xmm9
	pshufd	$68, %xmm2, %xmm2
	mulpd	%xmm7, %xmm9
	mulpd	%xmm0, %xmm2
	xorpd	%xmm11, %xmm9
	addpd	%xmm9, %xmm2
	addpd	%xmm2, %xmm5
	movdqu	32(%rdx), %xmm2
	pshufd	$238, %xmm2, %xmm2
	mulpd	%xmm2, %xmm7
	movdqu	32(%rdx), %xmm2
	addq	%rbx, %rdx
	pshufd	$68, %xmm2, %xmm2
	mulpd	%xmm2, %xmm0
	xorpd	%xmm11, %xmm7
	addpd	%xmm7, %xmm0
	addpd	%xmm0, %xmm3
	cmpq	%rbp, %r10
	jne	.L124
	jmp	.L122
.L185:
	movq	%rdx, %r12
	imulq	%r14, %rdx
	movq	16(%rsp), %rcx
	leaq	(%rax,%rdi), %r8
	salq	$4, %r8
	pxor	%xmm6, %xmm6
	salq	$4, %r12
	movq	%r14, %r10
	movq	.LC8(%rip), %xmm11
	addq	%rcx, %r8
	movapd	%xmm6, %xmm7
	salq	$4, %rdx
	addq	%r11, %rdx
.L130:
	movupd	(%rdx), %xmm2
	movdqu	(%r8), %xmm5
	addq	$1, %r10
	addq	%r12, %rdx
	movdqu	(%r8), %xmm1
	pshufd	$78, %xmm2, %xmm3
	pshufd	$238, %xmm5, %xmm5
	mulpd	%xmm3, %xmm5
	pshufd	$68, %xmm1, %xmm0
	mulpd	%xmm2, %xmm0
	xorpd	%xmm11, %xmm5
	addpd	%xmm5, %xmm0
	movdqu	16(%r8), %xmm5
	addq	%rbx, %r8
	addpd	%xmm0, %xmm7
	pshufd	$238, %xmm5, %xmm0
	mulpd	%xmm0, %xmm3
	pshufd	$68, %xmm5, %xmm0
	mulpd	%xmm0, %xmm2
	xorpd	%xmm11, %xmm3
	addpd	%xmm3, %xmm2
	addpd	%xmm2, %xmm6
	cmpq	%rbp, %r10
	jne	.L130
	jmp	.L128
.L160:
	movq	%rdx, %r12
	imulq	%r14, %rdx
	movq	16(%rsp), %rcx
	leaq	(%rax,%rdi), %r8
	salq	$4, %r8
	salq	$4, %r12
	movq	%r14, %r10
	pxor	%xmm0, %xmm0
	movq	.LC8(%rip), %xmm11
	addq	%rcx, %r8
	salq	$4, %rdx
	addq	%r11, %rdx
.L137:
	movdqu	(%rdx), %xmm6
	addq	$1, %r10
	pshufd	$78, %xmm6, %xmm3
	movdqu	(%r8), %xmm6
	addq	%rbx, %r8
	pshufd	$238, %xmm6, %xmm2
	mulpd	%xmm2, %xmm3
	pshufd	$68, %xmm6, %xmm2
	movupd	(%rdx), %xmm6
	addq	%r12, %rdx
	mulpd	%xmm6, %xmm2
	xorpd	%xmm11, %xmm3
	addpd	%xmm3, %xmm2
	addpd	%xmm2, %xmm0
	cmpq	%rbp, %r10
	jne	.L137
	jmp	.L134
.L159:
	movq	.LC8(%rip), %xmm11
	pxor	%xmm0, %xmm0
	jmp	.L134
.L156:
	movq	.LC8(%rip), %xmm11
	pxor	%xmm5, %xmm5
	movdqa	%xmm5, %xmm6
	movdqa	%xmm5, %xmm7
	movdqa	%xmm5, %xmm9
	jmp	.L116
.L157:
	movq	.LC8(%rip), %xmm11
	pxor	%xmm3, %xmm3
	movdqa	%xmm3, %xmm5
	movdqa	%xmm3, %xmm6
	jmp	.L122
.L158:
	movq	.LC8(%rip), %xmm11
	pxor	%xmm6, %xmm6
	movdqa	%xmm6, %xmm7
	jmp	.L128
.L184:
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm1
	movq	%r10, 232(%rsp)
	addq	$16, %r12
	movq	%rax, 224(%rsp)
	addq	$16, %r15
	movq	%r11, 216(%rsp)
	movq	%rdi, 208(%rsp)
	movq	%r9, 200(%rsp)
	movaps	%xmm8, 176(%rsp)
	movaps	%xmm12, 160(%rsp)
	movaps	%xmm6, 144(%rsp)
	movaps	%xmm5, 128(%rsp)
	movsd	%xmm3, 192(%rsp)
	movsd	%xmm2, 120(%rsp)
	call	__muldc3@PLT
	movupd	-16(%r12), %xmm6
	movq	200(%rsp), %r9
	unpcklpd	%xmm1, %xmm0
	movsd	120(%rsp), %xmm2
	movsd	192(%rsp), %xmm3
	addpd	%xmm6, %xmm0
	movapd	128(%rsp), %xmm5
	movapd	144(%rsp), %xmm6
	movq	.LC8(%rip), %xmm4
	movq	208(%rsp), %rdi
	movq	216(%rsp), %r11
	movq	224(%rsp), %rax
	movups	%xmm0, -16(%r12)
	cmpq	%r12, %r9
	movdqa	160(%rsp), %xmm12
	movdqa	176(%rsp), %xmm8
	movq	232(%rsp), %r10
	jne	.L140
	jmp	.L174
	.cfi_endproc
.LFE14847:
	.size	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi0EEELi0ELb0ES3_NS4_IS3_lLi1EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_.isra.0, .-_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi0EEELi0ELb0ES3_NS4_IS3_lLi1EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_.isra.0
	.section	.rodata.str1.8
	.align 8
.LC9:
	.string	"XprType& Eigen::CommaInitializer<MatrixType>::finished() [with XprType = Eigen::Matrix<std::complex<double>, -1, -1>]"
	.align 8
.LC10:
	.string	"/usr/local/include/Eigen/src/Core/CommaInitializer.h"
	.align 8
.LC11:
	.string	"((m_row+m_currentBlockRows) == m_xpr.rows() || m_xpr.cols() == 0) && m_col == m_xpr.cols() && \"Too few coefficients passed to comma initializer (operator<<)\""
	.text
	.align 2
	.p2align 4
	.type	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8finishedEv.isra.0, @function
_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8finishedEv.isra.0:
.LFB14848:
	.cfi_startproc
	movq	(%rdi), %rdx
	movq	24(%rdi), %rax
	addq	8(%rdi), %rax
	movq	16(%rdx), %rcx
	cmpq	8(%rdx), %rax
	je	.L189
	testq	%rcx, %rcx
	jne	.L190
.L189:
	cmpq	%rcx, 16(%rdi)
	jne	.L190
	ret
.L190:
	pushq	%rax
	.cfi_def_cfa_offset 16
	leaq	.LC9(%rip), %rcx
	movl	$122, %edx
	leaq	.LC10(%rip), %rsi
	leaq	.LC11(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE14848:
	.size	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8finishedEv.isra.0, .-_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8finishedEv.isra.0
	.align 2
	.p2align 4
	.type	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0, @function
_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0:
.LFB14849:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	subq	%rsi, %rdx
	movq	%rsi, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	movq	%rdi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movq	%rdx, %rbx
	subq	$16, %rsp
	.cfi_def_cfa_offset 48
	movq	%fs:40, %rax
	movq	%rax, 8(%rsp)
	xorl	%eax, %eax
	movq	%rdx, (%rsp)
	cmpq	$15, %rdx
	ja	.L204
	movq	(%rdi), %rdi
	cmpq	$1, %rdx
	jne	.L197
	movzbl	(%rsi), %eax
	movb	%al, (%rdi)
	movq	(%rsp), %rbx
	movq	0(%rbp), %rdi
.L198:
	movq	%rbx, 8(%rbp)
	movb	$0, (%rdi,%rbx)
	movq	8(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L205
	addq	$16, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L197:
	.cfi_restore_state
	testq	%rdx, %rdx
	je	.L198
	jmp	.L196
	.p2align 4,,10
	.p2align 3
.L204:
	movq	%rsp, %rsi
	xorl	%edx, %edx
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_createERmm@PLT
	movq	%rax, 0(%rbp)
	movq	%rax, %rdi
	movq	(%rsp), %rax
	movq	%rax, 16(%rbp)
.L196:
	movq	%rbx, %rdx
	movq	%r12, %rsi
	call	memcpy@PLT
	movq	(%rsp), %rbx
	movq	0(%rbp), %rdi
	jmp	.L198
.L205:
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE14849:
	.size	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0, .-_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
	.section	.rodata.str1.8
	.align 8
.LC12:
	.string	"Eigen::CommaInitializer<MatrixType>& Eigen::CommaInitializer<MatrixType>::operator,(const Eigen::DenseBase<OtherDerived>&) [with OtherDerived = Eigen::Matrix<std::complex<double>, -1, -1>; XprType = Eigen::Matrix<std::complex<double>, -1, -1>]"
	.align 8
.LC13:
	.string	"m_row+m_currentBlockRows<=m_xpr.rows() && \"Too many rows passed to comma initializer (operator<<)\""
	.align 8
.LC14:
	.string	"(m_col + other.cols() <= m_xpr.cols()) && \"Too many coefficients passed to comma initializer (operator<<)\""
	.align 8
.LC15:
	.string	"m_currentBlockRows==other.rows()"
	.align 8
.LC16:
	.string	"Eigen::MapBase<Derived, 0>::MapBase(PointerType, Eigen::Index, Eigen::Index) [with Derived = Eigen::Block<Eigen::Matrix<std::complex<double>, -1, -1>, -1, -1, false>; PointerType = std::complex<double>*; Eigen::Index = long int]"
	.align 8
.LC17:
	.string	"/usr/local/include/Eigen/src/Core/MapBase.h"
	.align 8
.LC18:
	.string	"(dataPtr == 0) || ( rows >= 0 && (RowsAtCompileTime == Dynamic || RowsAtCompileTime == rows) && cols >= 0 && (ColsAtCompileTime == Dynamic || ColsAtCompileTime == cols))"
	.align 8
.LC19:
	.ascii	"Eigen::Block<XprType, B"
	.string	"lockRows, BlockCols, InnerPanel>::Block(XprType&, Eigen::Index, Eigen::Index, Eigen::Index, Eigen::Index) [with XprType = Eigen::Matrix<std::complex<double>, -1, -1>; int BlockRows = -1; int BlockCols = -1; bool InnerPanel = false; Eigen::Index = long int]"
	.align 8
.LC20:
	.string	"/usr/local/include/Eigen/src/Core/Block.h"
	.align 8
.LC21:
	.string	"startRow >= 0 && blockRows >= 0 && startRow <= xpr.rows() - blockRows && startCol >= 0 && blockCols >= 0 && startCol <= xpr.cols() - blockCols"
	.text
	.align 2
	.p2align 4
	.type	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmIS4_EERS5_RKNS_9DenseBaseIT_EE.isra.0, @function
_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmIS4_EERS5_RKNS_9DenseBaseIT_EE.isra.0:
.LFB14851:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rdi, %r10
	movq	%rsi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$8, %rsp
	.cfi_def_cfa_offset 32
	movq	(%r10), %rcx
	movq	16(%rdi), %rdi
	movq	16(%rsi), %r9
	movq	16(%rcx), %rax
	cmpq	%rax, %rdi
	je	.L240
	leaq	(%rdi,%r9), %rdx
	cmpq	%rdx, %rax
	jl	.L227
	movq	8(%rsi), %rsi
	cmpq	%rsi, 24(%r10)
	jne	.L241
.L209:
	movq	8(%rcx), %r11
	movq	%rdi, %r8
	movq	8(%r10), %rdx
	imulq	%r11, %r8
.L226:
	movq	%rsi, %rbx
	addq	%rdx, %r8
	notq	%rbx
	salq	$4, %r8
	shrq	$63, %rbx
	addq	(%rcx), %r8
	je	.L212
	movq	%rsi, %rcx
	orq	%r9, %rcx
	js	.L242
.L212:
	movq	%rdx, %rcx
	orq	%rdi, %rcx
	orq	%r9, %rcx
	js	.L213
	testb	%bl, %bl
	je	.L213
	movq	%r11, %rcx
	subq	%rsi, %rcx
	cmpq	%rdx, %rcx
	jl	.L213
	subq	%r9, %rax
	cmpq	%rdi, %rax
	jl	.L213
	movq	0(%rbp), %rdx
	testb	$15, %r8b
	jne	.L215
	testq	%r9, %r9
	jle	.L222
	testq	%rsi, %rsi
	jle	.L225
	salq	$4, %rsi
	xorl	%ebx, %ebx
	xorl	%edi, %edi
	.p2align 4,,10
	.p2align 3
.L224:
	movq	%rbx, %rcx
	xorl	%eax, %eax
	salq	$4, %rcx
	addq	%r8, %rcx
	.p2align 4,,10
	.p2align 3
.L223:
	movupd	(%rdx,%rax), %xmm1
	movups	%xmm1, (%rcx,%rax)
	addq	$16, %rax
	cmpq	%rax, %rsi
	jne	.L223
	addq	$1, %rdi
	addq	%rsi, %rdx
	addq	%r11, %rbx
	cmpq	%r9, %rdi
	jne	.L224
	movq	16(%rbp), %r9
.L225:
	movq	16(%r10), %rdi
.L222:
	addq	%r9, %rdi
	movq	%rdi, 16(%r10)
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
.L240:
	.cfi_restore_state
	movq	8(%rsi), %rsi
	movq	24(%r10), %rdx
	testq	%r9, %r9
	je	.L243
.L208:
	addq	8(%r10), %rdx
	movq	8(%rcx), %r11
	movq	$0, 16(%r10)
	leaq	(%rdx,%rsi), %rdi
	movq	%rdx, 8(%r10)
	movq	%rsi, 24(%r10)
	cmpq	%r11, %rdi
	jg	.L244
	cmpq	%r9, %rax
	jl	.L227
	xorl	%r8d, %r8d
	xorl	%edi, %edi
	jmp	.L226
.L215:
	testq	%r9, %r9
	jle	.L222
	testq	%rsi, %rsi
	jle	.L222
	salq	$4, %rsi
	xorl	%ebp, %ebp
	xorl	%ebx, %ebx
	.p2align 4,,10
	.p2align 3
.L220:
	movq	%rbp, %rcx
	xorl	%eax, %eax
	salq	$4, %rcx
	addq	%r8, %rcx
	.p2align 4,,10
	.p2align 3
.L221:
	movsd	(%rdx,%rax), %xmm0
	movsd	%xmm0, (%rcx,%rax)
	movsd	8(%rdx,%rax), %xmm0
	movsd	%xmm0, 8(%rcx,%rax)
	addq	$16, %rax
	cmpq	%rsi, %rax
	jne	.L221
	addq	$1, %rbx
	addq	%r11, %rbp
	addq	%rsi, %rdx
	cmpq	%r9, %rbx
	jne	.L220
	addq	%r9, %rdi
	movq	%rdi, 16(%r10)
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
.L243:
	.cfi_restore_state
	cmpq	%rdx, %rsi
	jne	.L208
	jmp	.L209
.L244:
	leaq	.LC12(%rip), %rcx
	movl	$92, %edx
	leaq	.LC10(%rip), %rsi
	leaq	.LC13(%rip), %rdi
	call	__assert_fail@PLT
	.p2align 4,,10
	.p2align 3
.L213:
	leaq	.LC19(%rip), %rcx
	movl	$146, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC21(%rip), %rdi
	call	__assert_fail@PLT
.L242:
	leaq	.LC16(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L241:
	leaq	.LC12(%rip), %rcx
	movl	$97, %edx
	leaq	.LC10(%rip), %rsi
	leaq	.LC15(%rip), %rdi
	call	__assert_fail@PLT
.L227:
	leaq	.LC12(%rip), %rcx
	movl	$95, %edx
	leaq	.LC10(%rip), %rsi
	leaq	.LC14(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE14851:
	.size	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmIS4_EERS5_RKNS_9DenseBaseIT_EE.isra.0, .-_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmIS4_EERS5_RKNS_9DenseBaseIT_EE.isra.0
	.p2align 4
	.type	_ZSt5flushIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0, @function
_ZSt5flushIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0:
.LFB14853:
	.cfi_startproc
	jmp	_ZNSo5flushEv@PLT
	.cfi_endproc
.LFE14853:
	.size	_ZSt5flushIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0, .-_ZSt5flushIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0
	.p2align 4
	.type	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0, @function
_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0:
.LFB14854:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rdi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$8, %rsp
	.cfi_def_cfa_offset 32
	testq	%rsi, %rsi
	je	.L249
	movq	%rsi, %rdi
	movq	%rsi, %rbx
	call	strlen@PLT
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	movq	%rbx, %rsi
	movq	%rbp, %rdi
	popq	%rbx
	.cfi_def_cfa_offset 16
	movq	%rax, %rdx
	popq	%rbp
	.cfi_def_cfa_offset 8
	jmp	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	.p2align 4,,10
	.p2align 3
.L249:
	.cfi_restore_state
	movq	(%rdi), %rax
	movq	-24(%rax), %rdi
	addq	%rbp, %rdi
	movl	32(%rdi), %esi
	addq	$8, %rsp
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	orl	$1, %esi
	jmp	_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate@PLT
	.cfi_endproc
.LFE14854:
	.size	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0, .-_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
	.section	.rodata.str1.8
	.align 8
.LC22:
	.string	"basic_string: construction from null is not valid"
	.text
	.align 2
	.p2align 4
	.type	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC2IS3_EEPKcRKS3_.constprop.3, @function
_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC2IS3_EEPKcRKS3_.constprop.3:
.LFB14855:
	.cfi_startproc
	pushq	%r13
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	leaq	16(%rdi), %r13
	pushq	%r12
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	pushq	%rbp
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	pushq	%rbx
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
	subq	$24, %rsp
	.cfi_def_cfa_offset 64
	movq	%fs:40, %rax
	movq	%rax, 8(%rsp)
	xorl	%eax, %eax
	movq	%r13, (%rdi)
	testq	%rsi, %rsi
	je	.L261
	movq	%rdi, %rbx
	movq	%rsi, %rdi
	movq	%rsi, %r12
	call	strlen@PLT
	movq	%rax, (%rsp)
	movq	%rax, %rbp
	cmpq	$15, %rax
	ja	.L262
	cmpq	$1, %rax
	jne	.L254
	movzbl	(%r12), %eax
	movb	%al, 16(%rbx)
.L255:
	movq	(%rsp), %rax
	movq	(%rbx), %rdx
	movq	%rax, 8(%rbx)
	movb	$0, (%rdx,%rax)
	movq	8(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L263
	addq	$24, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	popq	%rbx
	.cfi_def_cfa_offset 32
	popq	%rbp
	.cfi_def_cfa_offset 24
	popq	%r12
	.cfi_def_cfa_offset 16
	popq	%r13
	.cfi_def_cfa_offset 8
	ret
.L254:
	.cfi_restore_state
	testq	%rax, %rax
	je	.L255
	jmp	.L253
.L262:
	movq	%rsp, %rsi
	xorl	%edx, %edx
	movq	%rbx, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_createERmm@PLT
	movq	%rax, (%rbx)
	movq	%rax, %r13
	movq	(%rsp), %rax
	movq	%rax, 16(%rbx)
.L253:
	movq	%rbp, %rdx
	movq	%r12, %rsi
	movq	%r13, %rdi
	call	memcpy@PLT
	jmp	.L255
.L263:
	call	__stack_chk_fail@PLT
.L261:
	leaq	.LC22(%rip), %rdi
	call	_ZSt19__throw_logic_errorPKc@PLT
	.cfi_endproc
.LFE14855:
	.size	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC2IS3_EEPKcRKS3_.constprop.3, .-_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC2IS3_EEPKcRKS3_.constprop.3
	.align 2
	.p2align 4
	.type	_ZN5Eigen8internal22lhs_process_one_packetILi4ELl1ELl1ESt7complexIdES3_S3_NS0_12DoublePacketIDv2_dEES5_S6_NS0_9Packet1cdENS0_11gebp_traitsIS3_S3_Lb1ELb0ELi1ELi0EEENS0_16BlasLinearMapperIS3_lLi0ELi1EEENS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEEEclERKSD_PKS3_SI_S3_llllllilllll.constprop.0.isra.0, @function
_ZN5Eigen8internal22lhs_process_one_packetILi4ELl1ELl1ESt7complexIdES3_S3_NS0_12DoublePacketIDv2_dEES5_S6_NS0_9Packet1cdENS0_11gebp_traitsIS3_S3_Lb1ELb0ELi1ELi0EEENS0_16BlasLinearMapperIS3_lLi0ELi1EEENS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEEEclERKSD_PKS3_SI_S3_llllllilllll.constprop.0.isra.0:
.LFB14859:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movq	%rdi, %r15
	movq	%rsi, %rdi
	movq	%rcx, %rsi
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$592, %rsp
	.cfi_def_cfa_offset 648
	movq	656(%rsp), %rcx
	movq	664(%rsp), %r13
	movsd	%xmm0, 568(%rsp)
	movsd	%xmm1, 576(%rsp)
	movq	680(%rsp), %r11
	movapd	568(%rsp), %xmm9
	testq	%rsi, %rsi
	jle	.L264
	movq	%rdx, %rax
	movq	%r9, %rdx
	salq	$4, %r8
	movq	648(%rsp), %r9
	movq	%rdx, %rbp
	movq	%r8, 552(%rsp)
	pshufd	$78, %xmm9, %xmm10
	salq	$4, %r9
	salq	$4, %rbp
	movq	%rsi, 560(%rsp)
	leaq	(%rdi,%r9), %r14
	movq	%rdx, %rdi
	xorl	%r9d, %r9d
	imulq	688(%rsp), %rdx
	salq	$6, %rdi
	movq	%r9, %r10
	movq	%rdi, 504(%rsp)
	movq	%rcx, %rdi
	salq	$6, %rdi
	addq	%rax, %rdi
	addq	%rcx, %rdx
	movq	%rdi, 544(%rsp)
	leaq	-1(%r13), %rdi
	salq	$4, %rdx
	shrq	$3, %rdi
	addq	%rdx, %rax
	addq	$1, %rdi
	movq	%rax, 536(%rsp)
	movq	%rdi, %rbx
	salq	$7, %rdi
	salq	$9, %rbx
	movq	%rdi, 496(%rsp)
	movq	%rbx, 520(%rsp)
	movq	%r13, %rbx
	movq	%rbp, %r13
	movq	688(%rsp), %rbp
	.p2align 4,,10
	.p2align 3
.L266:
	movq	496(%rsp), %rax
	movq	544(%rsp), %r8
	xorl	%edi, %edi
	addq	%r14, %rax
	movq	%rax, 512(%rsp)
	testq	%rbp, %rbp
	jle	.L280
	movq	%r13, 528(%rsp)
	movq	%rbx, %r13
	.p2align 4,,10
	.p2align 3
.L271:
	movq	8(%r15), %rbx
	leaq	1(%rdi), %rdx
	leaq	2(%rdi), %rcx
	movq	(%r15), %r9
	leaq	3(%rdi), %rsi
	prefetcht0	(%r14)
	prefetcht0	(%r8)
	movq	%rbx, %rax
	imulq	%rbx, %rdx
	imulq	%rdi, %rax
	imulq	%rbx, %rcx
	imulq	%rbx, %rsi
	addq	%r10, %rdx
	addq	%r10, %rax
	salq	$4, %rdx
	addq	%r10, %rcx
	salq	$4, %rax
	addq	%r9, %rdx
	addq	%r10, %rsi
	salq	$4, %rcx
	addq	%r9, %rax
	salq	$4, %rsi
	addq	%r9, %rcx
	addq	%r9, %rsi
	movq	%r8, %r9
	testq	%r13, %r13
	jle	.L281
	xorl	%r9d, %r9d
	movq	%r11, -120(%rsp)
	movq	%r14, %rbx
	movq	%r8, %r11
	leaq	768(%r8), %r12
	movq	%rdi, %r8
	movq	%r9, %rdi
	movq	%r10, %r9
	movaps	%xmm9, 456(%rsp)
	movq	%r14, %r10
	movq	%rbp, %r14
	pxor	%xmm3, %xmm3
	movq	%r15, %rbp
	movq	%r14, %r15
	movq	-120(%rsp), %r14
	movaps	%xmm3, 56(%rsp)
	movaps	%xmm3, 40(%rsp)
	movapd	%xmm3, %xmm12
	movapd	%xmm3, %xmm5
	movapd	%xmm3, %xmm9
	movaps	%xmm3, 24(%rsp)
	movaps	%xmm3, 8(%rsp)
	movaps	%xmm3, -104(%rsp)
	movaps	%xmm3, -88(%rsp)
	movaps	%xmm3, -8(%rsp)
	movaps	%xmm3, -24(%rsp)
	movaps	%xmm3, -56(%rsp)
	movaps	%xmm3, -40(%rsp)
	movaps	%xmm3, -72(%rsp)
	movaps	%xmm3, 88(%rsp)
	movaps	%xmm3, 72(%rsp)
	movaps	%xmm10, 472(%rsp)
	.p2align 4,,10
	.p2align 3
.L268:
#APP
# 1264 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin gebp micro kernel 1/half/quarterX4
# 0 "" 2
#NO_APP
	prefetcht0	(%r12)
#APP
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	-760(%r12), %xmm3
	movapd	(%rbx), %xmm15
	movsd	-728(%r12), %xmm6
	movsd	-720(%r12), %xmm4
	movsd	-712(%r12), %xmm7
	movaps	%xmm15, -120(%rsp)
	movsd	-768(%r12), %xmm2
	movsd	%xmm3, 424(%rsp)
	movsd	-752(%r12), %xmm3
	movsd	-744(%r12), %xmm11
	movsd	-736(%r12), %xmm8
	movsd	%xmm6, 440(%rsp)
	movsd	%xmm3, 432(%rsp)
	movsd	%xmm4, 448(%rsp)
	movsd	%xmm7, 488(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	-704(%r12), %xmm1
	movapd	16(%rbx), %xmm6
	movsd	-696(%r12), %xmm0
	movsd	-688(%r12), %xmm10
	movsd	-680(%r12), %xmm13
	movsd	-672(%r12), %xmm14
	movaps	%xmm6, 344(%rsp)
	movsd	-664(%r12), %xmm3
	movsd	-656(%r12), %xmm4
	movsd	%xmm1, 360(%rsp)
	movsd	-648(%r12), %xmm7
	movsd	%xmm0, 368(%rsp)
	movsd	%xmm10, 376(%rsp)
	movsd	%xmm13, 384(%rsp)
	movsd	%xmm14, 392(%rsp)
	movsd	%xmm3, 400(%rsp)
	movsd	%xmm4, 408(%rsp)
	movsd	%xmm7, 416(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	-632(%r12), %xmm10
	movapd	32(%rbx), %xmm0
	movsd	-640(%r12), %xmm1
	movsd	-600(%r12), %xmm15
	movsd	-584(%r12), %xmm14
	movsd	%xmm10, 320(%rsp)
	movsd	-624(%r12), %xmm13
	movsd	-616(%r12), %xmm10
	movsd	-608(%r12), %xmm7
	movsd	%xmm1, 312(%rsp)
	movsd	-592(%r12), %xmm3
	movsd	%xmm15, 328(%rsp)
	movsd	%xmm14, 336(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	-576(%r12), %xmm1
	movapd	48(%rbx), %xmm4
	movsd	-568(%r12), %xmm15
	movsd	-560(%r12), %xmm14
	movsd	%xmm1, 248(%rsp)
	movsd	-552(%r12), %xmm1
	movsd	-528(%r12), %xmm6
	movsd	%xmm15, 256(%rsp)
	movsd	-544(%r12), %xmm15
	movsd	%xmm14, 264(%rsp)
	movsd	-536(%r12), %xmm14
	movsd	%xmm1, 272(%rsp)
	movsd	-520(%r12), %xmm1
	movaps	%xmm4, 232(%rsp)
	movsd	%xmm15, 280(%rsp)
	movsd	%xmm14, 288(%rsp)
	movsd	%xmm6, 296(%rsp)
	movsd	%xmm1, 304(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
#NO_APP
	prefetcht0	256(%r12)
#APP
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	-496(%r12), %xmm14
	movapd	64(%rbx), %xmm1
	movsd	-512(%r12), %xmm15
	movsd	-488(%r12), %xmm6
	movsd	%xmm14, 200(%rsp)
	movsd	-464(%r12), %xmm14
	movsd	-472(%r12), %xmm4
	movsd	%xmm15, 192(%rsp)
	movsd	-504(%r12), %xmm15
	movsd	%xmm14, 216(%rsp)
	movsd	-456(%r12), %xmm14
	movsd	%xmm6, 208(%rsp)
	movsd	-480(%r12), %xmm6
	movsd	%xmm14, 224(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movapd	80(%rbx), %xmm14
	movaps	%xmm14, 104(%rsp)
	movsd	-448(%r12), %xmm14
	movsd	%xmm14, 128(%rsp)
	movsd	-440(%r12), %xmm14
	movsd	%xmm14, 136(%rsp)
	movsd	-432(%r12), %xmm14
	movsd	%xmm14, 144(%rsp)
	movsd	-424(%r12), %xmm14
	movsd	%xmm14, 152(%rsp)
	movsd	-416(%r12), %xmm14
	movsd	%xmm14, 160(%rsp)
	movsd	-408(%r12), %xmm14
	movsd	%xmm14, 168(%rsp)
	movsd	-400(%r12), %xmm14
	movsd	%xmm14, 176(%rsp)
	movsd	-392(%r12), %xmm14
	movsd	%xmm14, 184(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movapd	%xmm2, %xmm14
	unpcklpd	%xmm15, %xmm15
	unpcklpd	%xmm7, %xmm7
	movsd	312(%rsp), %xmm2
	unpcklpd	%xmm14, %xmm14
	unpcklpd	%xmm8, %xmm8
	unpcklpd	%xmm6, %xmm6
	unpcklpd	%xmm2, %xmm2
	unpcklpd	%xmm11, %xmm11
	unpcklpd	%xmm13, %xmm13
	mulpd	%xmm0, %xmm2
	unpcklpd	%xmm10, %xmm10
	unpcklpd	%xmm4, %xmm4
	mulpd	-120(%rsp), %xmm14
	unpcklpd	%xmm3, %xmm3
	addpd	72(%rsp), %xmm14
	mulpd	%xmm1, %xmm15
	mulpd	%xmm0, %xmm7
	mulpd	%xmm1, %xmm6
	addpd	%xmm2, %xmm14
	mulpd	%xmm0, %xmm13
	movsd	192(%rsp), %xmm2
	mulpd	%xmm0, %xmm10
	mulpd	%xmm1, %xmm4
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm1, %xmm2
	mulpd	%xmm0, %xmm3
	addpd	%xmm2, %xmm14
	movsd	-384(%r12), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	96(%rbx), %xmm2
	addpd	%xmm2, %xmm14
	movsd	424(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	-120(%rsp), %xmm2
	movaps	%xmm14, 72(%rsp)
	addpd	88(%rsp), %xmm2
	movapd	%xmm2, %xmm14
	movsd	320(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm0, %xmm2
	movaps	%xmm2, 88(%rsp)
	movapd	88(%rsp), %xmm2
	addpd	%xmm14, %xmm2
	addpd	%xmm2, %xmm15
	movsd	-376(%r12), %xmm2
	unpcklpd	%xmm2, %xmm2
	movapd	%xmm15, %xmm14
	movapd	96(%rbx), %xmm15
	mulpd	%xmm2, %xmm15
	movsd	-352(%r12), %xmm2
	unpcklpd	%xmm2, %xmm2
	addpd	%xmm14, %xmm15
	movsd	432(%rsp), %xmm14
	unpcklpd	%xmm14, %xmm14
	movaps	%xmm15, 88(%rsp)
	movapd	-120(%rsp), %xmm15
	mulpd	%xmm15, %xmm8
	addpd	-72(%rsp), %xmm8
	mulpd	%xmm15, %xmm14
	mulpd	%xmm15, %xmm11
	addpd	%xmm8, %xmm7
	addpd	%xmm12, %xmm14
	movsd	200(%rsp), %xmm12
	addpd	%xmm6, %xmm7
	movapd	96(%rbx), %xmm6
	addpd	%xmm9, %xmm11
	movsd	208(%rsp), %xmm9
	unpcklpd	%xmm12, %xmm12
	mulpd	%xmm2, %xmm6
	unpcklpd	%xmm9, %xmm9
	addpd	%xmm14, %xmm13
	mulpd	%xmm1, %xmm12
	movapd	%xmm7, %xmm2
	addpd	%xmm11, %xmm10
	mulpd	%xmm1, %xmm9
	addpd	%xmm6, %xmm2
	addpd	%xmm12, %xmm13
	movsd	-368(%r12), %xmm12
	addpd	%xmm9, %xmm10
	movsd	-360(%r12), %xmm9
	movaps	%xmm2, -72(%rsp)
	unpcklpd	%xmm12, %xmm12
	movsd	440(%rsp), %xmm2
	mulpd	96(%rbx), %xmm12
	unpcklpd	%xmm9, %xmm9
	mulpd	96(%rbx), %xmm9
	unpcklpd	%xmm2, %xmm2
	movapd	%xmm2, %xmm6
	movsd	-344(%r12), %xmm2
	mulpd	%xmm15, %xmm6
	unpcklpd	%xmm2, %xmm2
	addpd	%xmm13, %xmm12
	addpd	%xmm10, %xmm9
	addpd	%xmm5, %xmm6
	movsd	328(%rsp), %xmm5
	unpcklpd	%xmm5, %xmm5
	mulpd	%xmm0, %xmm5
	addpd	%xmm6, %xmm5
	addpd	%xmm4, %xmm5
	movapd	96(%rbx), %xmm4
	mulpd	%xmm2, %xmm4
	movsd	448(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	addpd	%xmm4, %xmm5
	movapd	%xmm2, %xmm4
	movsd	216(%rsp), %xmm2
	mulpd	%xmm15, %xmm4
	addpd	-40(%rsp), %xmm4
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm1, %xmm2
	addpd	%xmm4, %xmm3
	addpd	%xmm2, %xmm3
	movsd	-336(%r12), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	96(%rbx), %xmm2
	addpd	%xmm2, %xmm3
	movsd	488(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm15, %xmm2
	movaps	%xmm3, -40(%rsp)
	movapd	-56(%rsp), %xmm3
	addpd	%xmm2, %xmm3
	movsd	336(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm0, %xmm2
	movsd	224(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm0, %xmm1
	addpd	%xmm3, %xmm2
	addpd	%xmm1, %xmm2
	movsd	-328(%r12), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	96(%rbx), %xmm1
	addpd	%xmm1, %xmm2
	movaps	%xmm2, -56(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	360(%rsp), %xmm1
	movapd	-24(%rsp), %xmm2
	movapd	344(%rsp), %xmm6
	movapd	104(%rsp), %xmm7
	unpcklpd	%xmm1, %xmm1
	movapd	112(%rbx), %xmm3
	movapd	232(%rsp), %xmm4
	movsd	128(%rsp), %xmm0
	mulpd	%xmm6, %xmm1
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm7, %xmm0
	addpd	%xmm1, %xmm2
	movsd	248(%rsp), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm4, %xmm1
	addpd	%xmm2, %xmm1
	addpd	%xmm0, %xmm1
	movsd	-320(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	368(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm6, %xmm0
	movaps	%xmm1, -24(%rsp)
	movapd	-8(%rsp), %xmm1
	addpd	%xmm0, %xmm1
	movsd	256(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm4, %xmm0
	addpd	%xmm0, %xmm1
	movsd	136(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm7, %xmm0
	addpd	%xmm0, %xmm1
	movsd	-312(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	376(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm6, %xmm0
	movaps	%xmm1, -8(%rsp)
	movapd	-88(%rsp), %xmm1
	addpd	%xmm0, %xmm1
	movsd	264(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm4, %xmm0
	addpd	%xmm0, %xmm1
	movsd	144(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm7, %xmm0
	addpd	%xmm0, %xmm1
	movsd	-304(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	384(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	movaps	%xmm1, -88(%rsp)
	mulpd	%xmm6, %xmm0
	movsd	272(%rsp), %xmm1
	addpd	-104(%rsp), %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm4, %xmm1
	addpd	%xmm1, %xmm0
	movsd	152(%rsp), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm7, %xmm1
	addpd	%xmm0, %xmm1
	movsd	-296(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	392(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	movaps	%xmm1, -104(%rsp)
	mulpd	%xmm6, %xmm0
	movsd	280(%rsp), %xmm1
	addpd	8(%rsp), %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm4, %xmm1
	addpd	%xmm1, %xmm0
	movsd	160(%rsp), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm7, %xmm1
	addpd	%xmm0, %xmm1
	movsd	-288(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	400(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	movaps	%xmm1, 8(%rsp)
	mulpd	%xmm6, %xmm0
	movsd	288(%rsp), %xmm1
	addpd	24(%rsp), %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm4, %xmm1
	addpd	%xmm1, %xmm0
	movsd	168(%rsp), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm7, %xmm1
	addpd	%xmm0, %xmm1
	movsd	-280(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	408(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	movaps	%xmm1, 24(%rsp)
	mulpd	%xmm6, %xmm0
	movsd	296(%rsp), %xmm1
	addpd	40(%rsp), %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm4, %xmm1
	addpd	%xmm1, %xmm0
	movsd	176(%rsp), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm7, %xmm1
	addpd	%xmm0, %xmm1
	movsd	-272(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	416(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	movaps	%xmm1, 40(%rsp)
	mulpd	%xmm6, %xmm0
	movsd	304(%rsp), %xmm1
	addpd	56(%rsp), %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm4, %xmm1
	addpd	%xmm0, %xmm1
	movsd	184(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm7, %xmm0
	addpd	%xmm1, %xmm0
	movsd	-264(%r12), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm3, %xmm1
	addpd	%xmm1, %xmm0
	movaps	%xmm0, 56(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
#NO_APP
	subq	$-128, %rbx
#APP
# 1282 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end gebp micro kernel 1/half/quarterX4
# 0 "" 2
#NO_APP
	addq	$8, %rdi
	addq	$512, %r12
	cmpq	%r13, %rdi
	jl	.L268
	movq	%r15, %rbx
	movapd	%xmm12, %xmm15
	movq	%rbp, %r15
	movq	%r8, %rdi
	movapd	-88(%rsp), %xmm3
	movapd	-104(%rsp), %xmm4
	movq	%rbx, %rbp
	movq	%r11, %r8
	movq	520(%rsp), %rbx
	movapd	%xmm9, %xmm12
	movapd	%xmm5, %xmm7
	movapd	-24(%rsp), %xmm0
	movapd	-8(%rsp), %xmm1
	movapd	8(%rsp), %xmm6
	movq	%r14, %r11
	addpd	%xmm15, %xmm3
	movapd	40(%rsp), %xmm8
	movapd	56(%rsp), %xmm11
	addpd	%xmm12, %xmm4
	movq	%r10, %r14
	addpd	72(%rsp), %xmm0
	addpd	88(%rsp), %xmm1
	movq	%r9, %r10
	addpd	-72(%rsp), %xmm6
	addpd	24(%rsp), %xmm7
	addpd	-40(%rsp), %xmm8
	leaq	(%r8,%rbx), %r9
	movapd	456(%rsp), %xmm9
	addpd	-56(%rsp), %xmm11
	movq	512(%rsp), %rbx
	movdqa	472(%rsp), %xmm10
.L267:
	cmpq	%r11, %r13
	jge	.L269
	movq	%r13, %r12
	.p2align 4,,10
	.p2align 3
.L270:
#APP
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	(%r9), %xmm5
	movapd	(%rbx), %xmm2
	unpcklpd	%xmm5, %xmm5
	mulpd	%xmm2, %xmm5
	addpd	%xmm5, %xmm0
	movsd	8(%r9), %xmm5
	unpcklpd	%xmm5, %xmm5
	mulpd	%xmm2, %xmm5
	addpd	%xmm5, %xmm1
	movsd	16(%r9), %xmm5
	unpcklpd	%xmm5, %xmm5
	mulpd	%xmm2, %xmm5
	addpd	%xmm5, %xmm3
	movsd	24(%r9), %xmm5
	unpcklpd	%xmm5, %xmm5
	mulpd	%xmm2, %xmm5
	addpd	%xmm5, %xmm4
	movsd	32(%r9), %xmm5
	unpcklpd	%xmm5, %xmm5
	mulpd	%xmm2, %xmm5
	addpd	%xmm5, %xmm6
	movsd	40(%r9), %xmm5
	unpcklpd	%xmm5, %xmm5
	mulpd	%xmm2, %xmm5
	addpd	%xmm5, %xmm7
	movsd	48(%r9), %xmm5
	unpcklpd	%xmm5, %xmm5
	mulpd	%xmm2, %xmm5
	addpd	%xmm5, %xmm8
	movsd	56(%r9), %xmm5
	unpcklpd	%xmm5, %xmm5
	mulpd	%xmm5, %xmm2
	addpd	%xmm2, %xmm11
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
#NO_APP
	addq	$1, %r12
	addq	$64, %r9
	addq	$16, %rbx
	cmpq	%r11, %r12
	jne	.L270
.L269:
	shufpd	$1, %xmm4, %xmm4
	shufpd	$1, %xmm1, %xmm1
	xorpd	.LC7(%rip), %xmm3
	xorpd	.LC7(%rip), %xmm0
	shufpd	$1, %xmm11, %xmm11
	shufpd	$1, %xmm7, %xmm7
	addq	$4, %rdi
	xorpd	.LC7(%rip), %xmm8
	addpd	%xmm3, %xmm4
	addpd	%xmm0, %xmm1
	movupd	(%rdx), %xmm3
	xorpd	.LC7(%rip), %xmm6
	addpd	%xmm8, %xmm11
	addpd	%xmm6, %xmm7
	pshufd	$238, %xmm4, %xmm2
	pshufd	$68, %xmm4, %xmm4
	pshufd	$238, %xmm1, %xmm0
	mulpd	%xmm9, %xmm4
	pshufd	$68, %xmm1, %xmm1
	mulpd	%xmm10, %xmm2
	mulpd	%xmm9, %xmm1
	xorpd	.LC8(%rip), %xmm2
	mulpd	%xmm10, %xmm0
	xorpd	.LC8(%rip), %xmm0
	addpd	%xmm2, %xmm4
	addpd	%xmm0, %xmm1
	pshufd	$238, %xmm7, %xmm0
	pshufd	$68, %xmm7, %xmm7
	addpd	%xmm3, %xmm4
	movupd	(%rax), %xmm3
	mulpd	%xmm9, %xmm7
	mulpd	%xmm10, %xmm0
	addpd	%xmm3, %xmm1
	xorpd	.LC8(%rip), %xmm0
	movups	%xmm1, (%rax)
	pshufd	$238, %xmm11, %xmm1
	addpd	%xmm0, %xmm7
	pshufd	$68, %xmm11, %xmm11
	mulpd	%xmm9, %xmm11
	movups	%xmm4, (%rdx)
	movupd	(%rsi), %xmm3
	movq	504(%rsp), %rax
	mulpd	%xmm10, %xmm1
	xorpd	.LC8(%rip), %xmm1
	addq	%rax, %r8
	addpd	%xmm1, %xmm11
	addpd	%xmm3, %xmm11
	movupd	(%rcx), %xmm3
	addpd	%xmm3, %xmm7
	movups	%xmm7, (%rcx)
	movups	%xmm11, (%rsi)
	cmpq	%rbp, %rdi
	jl	.L271
	movq	%r13, %rbx
	movq	528(%rsp), %r13
.L280:
	movq	496(%rsp), %rax
	movq	536(%rsp), %r8
	movq	%rbp, %r9
	leaq	(%r14,%rax), %r12
	cmpq	672(%rsp), %rbp
	jge	.L278
	movq	%rbp, 688(%rsp)
	movq	672(%rsp), %rbp
	.p2align 4,,10
	.p2align 3
.L277:
	movq	8(%r15), %rax
	movq	%r14, %rdx
	prefetcht0	(%r14)
	imulq	%r9, %rax
	addq	%r10, %rax
	salq	$4, %rax
	addq	(%r15), %rax
	movq	%rax, %rsi
	movq	%r8, %rax
	testq	%rbx, %rbx
	jle	.L282
	pxor	%xmm0, %xmm0
	movq	%r12, %rcx
	xorl	%edi, %edi
	movq	%r10, %r12
	movapd	%xmm0, %xmm1
	movq	%r9, %r10
	movaps	%xmm9, -8(%rsp)
	movq	%rsi, %r9
	movaps	%xmm10, 8(%rsp)
	movq	%r8, %rsi
	movq	%r13, %r8
	.p2align 4,,10
	.p2align 3
.L274:
#APP
# 1336 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1350 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1350 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	8(%rax), %xmm3
	movupd	(%rdx), %xmm11
	movq	(%rax), %r13
	movsd	%xmm3, -120(%rsp)
#APP
# 1350 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1351 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1351 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	16(%rax), %xmm3
	movsd	24(%rax), %xmm4
	movupd	16(%rdx), %xmm8
	movsd	%xmm3, -104(%rsp)
	movsd	%xmm4, -88(%rsp)
#APP
# 1351 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1352 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1352 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	32(%rax), %xmm6
	movsd	40(%rax), %xmm5
	movupd	32(%rdx), %xmm7
	movsd	%xmm6, -72(%rsp)
	movsd	%xmm5, -56(%rsp)
#APP
# 1352 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1353 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1353 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	48(%rax), %xmm2
	movupd	48(%rdx), %xmm6
	movsd	56(%rax), %xmm15
	movsd	%xmm2, -40(%rsp)
#APP
# 1353 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1354 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1354 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	64(%rax), %xmm10
	movupd	64(%rdx), %xmm5
	movsd	72(%rax), %xmm14
	movsd	%xmm10, -24(%rsp)
#APP
# 1354 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1355 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1355 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movupd	80(%rdx), %xmm4
	movsd	80(%rax), %xmm10
	movsd	88(%rax), %xmm13
#APP
# 1355 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1356 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1356 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movupd	96(%rdx), %xmm3
	movsd	96(%rax), %xmm9
	movsd	104(%rax), %xmm12
#APP
# 1356 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1357 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1357 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movq	%r13, %xmm2
	unpcklpd	%xmm10, %xmm10
	unpcklpd	%xmm9, %xmm9
	unpcklpd	%xmm2, %xmm2
	unpcklpd	%xmm15, %xmm15
	unpcklpd	%xmm14, %xmm14
	mulpd	%xmm11, %xmm2
	unpcklpd	%xmm13, %xmm13
	unpcklpd	%xmm12, %xmm12
	mulpd	%xmm4, %xmm10
	mulpd	%xmm3, %xmm9
	mulpd	%xmm6, %xmm15
	mulpd	%xmm5, %xmm14
	addpd	%xmm2, %xmm1
	movsd	-104(%rsp), %xmm2
	mulpd	%xmm4, %xmm13
	mulpd	%xmm3, %xmm12
	movupd	112(%rdx), %xmm3
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm8, %xmm2
	addpd	%xmm2, %xmm1
	movsd	-72(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm7, %xmm2
	addpd	%xmm2, %xmm1
	movsd	-40(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm6, %xmm2
	addpd	%xmm2, %xmm1
	movsd	-24(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm5, %xmm2
	addpd	%xmm2, %xmm1
	movupd	112(%rdx), %xmm2
	addpd	%xmm10, %xmm1
	addpd	%xmm9, %xmm1
	movsd	112(%rax), %xmm9
	unpcklpd	%xmm9, %xmm9
	mulpd	%xmm2, %xmm9
	movsd	-120(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm11, %xmm2
	addpd	%xmm9, %xmm1
	addpd	%xmm2, %xmm0
	movsd	-88(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm8, %xmm2
	addpd	%xmm2, %xmm0
	movsd	-56(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm7, %xmm2
	addpd	%xmm2, %xmm0
	movsd	120(%rax), %xmm2
	unpcklpd	%xmm2, %xmm2
	addpd	%xmm15, %xmm0
	mulpd	%xmm3, %xmm2
	addpd	%xmm14, %xmm0
	addpd	%xmm13, %xmm0
	addpd	%xmm12, %xmm0
	addpd	%xmm2, %xmm0
#APP
# 1357 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
#NO_APP
	subq	$-128, %rax
	subq	$-128, %rdx
#APP
# 1362 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end gebp micro kernel 1/half/quarterX1
# 0 "" 2
#NO_APP
	addq	$8, %rdi
	cmpq	%rdi, %rbx
	jg	.L274
	movq	496(%rsp), %rax
	movq	%r8, %r13
	movq	%rsi, %r8
	movq	%rcx, %rdx
	movapd	-8(%rsp), %xmm9
	movdqa	8(%rsp), %xmm10
	movq	%r9, %rsi
	movq	%r10, %r9
	addq	%r8, %rax
	movq	%r12, %r10
	movq	%rcx, %r12
.L273:
	cmpq	%r11, %rbx
	jge	.L275
	movq	%rbx, %rcx
	.p2align 4,,10
	.p2align 3
.L276:
#APP
# 1369 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1369 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	(%rax), %xmm2
	movupd	(%rdx), %xmm3
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm3, %xmm2
	addpd	%xmm2, %xmm1
	movsd	8(%rax), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm3, %xmm2
	addpd	%xmm2, %xmm0
#APP
# 1369 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
#NO_APP
	addq	$1, %rcx
	addq	$16, %rax
	addq	$16, %rdx
	cmpq	%rcx, %r11
	jne	.L276
.L275:
	shufpd	$1, %xmm0, %xmm0
	movupd	(%rsi), %xmm3
	addq	$1, %r9
	addq	%r13, %r8
	xorpd	.LC7(%rip), %xmm1
	addpd	%xmm1, %xmm0
	pshufd	$238, %xmm0, %xmm1
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm9, %xmm0
	mulpd	%xmm10, %xmm1
	xorpd	.LC8(%rip), %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm3, %xmm0
	movups	%xmm0, (%rsi)
	cmpq	%rbp, %r9
	jne	.L277
	movq	688(%rsp), %rbp
.L278:
	movq	552(%rsp), %rax
	addq	$1, %r10
	addq	%rax, %r14
	movq	560(%rsp), %rax
	cmpq	%rax, %r10
	jne	.L266
.L264:
	addq	$592, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L281:
	.cfi_restore_state
	pxor	%xmm11, %xmm11
	movq	%r14, %rbx
	movapd	%xmm11, %xmm8
	movapd	%xmm11, %xmm7
	movapd	%xmm11, %xmm6
	movapd	%xmm11, %xmm4
	movapd	%xmm11, %xmm3
	movapd	%xmm11, %xmm1
	movapd	%xmm11, %xmm0
	jmp	.L267
	.p2align 4,,10
	.p2align 3
.L282:
	pxor	%xmm0, %xmm0
	movapd	%xmm0, %xmm1
	jmp	.L273
	.cfi_endproc
.LFE14859:
	.size	_ZN5Eigen8internal22lhs_process_one_packetILi4ELl1ELl1ESt7complexIdES3_S3_NS0_12DoublePacketIDv2_dEES5_S6_NS0_9Packet1cdENS0_11gebp_traitsIS3_S3_Lb1ELb0ELi1ELi0EEENS0_16BlasLinearMapperIS3_lLi0ELi1EEENS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEEEclERKSD_PKS3_SI_S3_llllllilllll.constprop.0.isra.0, .-_ZN5Eigen8internal22lhs_process_one_packetILi4ELl1ELl1ESt7complexIdES3_S3_NS0_12DoublePacketIDv2_dEES5_S6_NS0_9Packet1cdENS0_11gebp_traitsIS3_S3_Lb1ELb0ELi1ELi0EEENS0_16BlasLinearMapperIS3_lLi0ELi1EEENS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEEEclERKSD_PKS3_SI_S3_llllllilllll.constprop.0.isra.0
	.align 2
	.p2align 4
	.type	_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb1ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0, @function
_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb1ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0:
.LFB14868:
	.cfi_startproc
	movq	%rdx, %r10
	movq	%r9, %rdx
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	sarq	$63, %rdx
	shrq	$62, %rdx
	leaq	(%r9,%rdx), %rax
	andl	$3, %eax
	subq	%rdx, %rax
	movq	%r9, %rdx
	subq	%rax, %rdx
	movq	%r8, %rax
	andq	$-8, %rax
	pushq	%rdx
	.cfi_def_cfa_offset 24
	movq	%r10, %rdx
	pushq	%r8
	.cfi_def_cfa_offset 32
	pushq	%r9
	.cfi_def_cfa_offset 40
	movq	%r8, %r9
	pushq	%rax
	.cfi_def_cfa_offset 48
	pushq	$0
	.cfi_def_cfa_offset 56
	pushq	$0
	.cfi_def_cfa_offset 64
	call	_ZN5Eigen8internal22lhs_process_one_packetILi4ELl1ELl1ESt7complexIdES3_S3_NS0_12DoublePacketIDv2_dEES5_S6_NS0_9Packet1cdENS0_11gebp_traitsIS3_S3_Lb1ELb0ELi1ELi0EEENS0_16BlasLinearMapperIS3_lLi0ELi1EEENS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEEEclERKSD_PKS3_SI_S3_llllllilllll.constprop.0.isra.0
	addq	$56, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE14868:
	.size	_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb1ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0, .-_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb1ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0
	.align 2
	.p2align 4
	.type	_ZN5Eigen8internal22lhs_process_one_packetILi4ELl1ELl1ESt7complexIdES3_S3_NS0_12DoublePacketIDv2_dEES5_S6_NS0_9Packet1cdENS0_11gebp_traitsIS3_S3_Lb0ELb0ELi1ELi0EEENS0_16BlasLinearMapperIS3_lLi0ELi1EEENS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEEEclERKSD_PKS3_SI_S3_llllllilllll.constprop.0.isra.0, @function
_ZN5Eigen8internal22lhs_process_one_packetILi4ELl1ELl1ESt7complexIdES3_S3_NS0_12DoublePacketIDv2_dEES5_S6_NS0_9Packet1cdENS0_11gebp_traitsIS3_S3_Lb0ELb0ELi1ELi0EEENS0_16BlasLinearMapperIS3_lLi0ELi1EEENS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEEEclERKSD_PKS3_SI_S3_llllllilllll.constprop.0.isra.0:
.LFB14870:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movq	%rdi, %r15
	movq	%rsi, %rdi
	movq	%rcx, %rsi
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$592, %rsp
	.cfi_def_cfa_offset 648
	movq	656(%rsp), %rcx
	movq	664(%rsp), %r13
	movsd	%xmm0, 568(%rsp)
	movsd	%xmm1, 576(%rsp)
	movq	680(%rsp), %r11
	movapd	568(%rsp), %xmm9
	testq	%rsi, %rsi
	jle	.L295
	movq	%rdx, %rax
	movq	%r9, %rdx
	salq	$4, %r8
	movq	648(%rsp), %r9
	movq	%rdx, %rbp
	movq	%r8, 552(%rsp)
	pshufd	$78, %xmm9, %xmm10
	salq	$4, %r9
	salq	$4, %rbp
	movq	%rsi, 560(%rsp)
	leaq	(%rdi,%r9), %r14
	movq	%rdx, %rdi
	xorl	%r9d, %r9d
	imulq	688(%rsp), %rdx
	salq	$6, %rdi
	movq	%r9, %r10
	movq	%rdi, 504(%rsp)
	movq	%rcx, %rdi
	salq	$6, %rdi
	addq	%rax, %rdi
	addq	%rcx, %rdx
	movq	%rdi, 544(%rsp)
	leaq	-1(%r13), %rdi
	salq	$4, %rdx
	shrq	$3, %rdi
	addq	%rdx, %rax
	addq	$1, %rdi
	movq	%rax, 536(%rsp)
	movq	%rdi, %rbx
	salq	$7, %rdi
	salq	$9, %rbx
	movq	%rdi, 496(%rsp)
	movq	%rbx, 520(%rsp)
	movq	%r13, %rbx
	movq	%rbp, %r13
	movq	688(%rsp), %rbp
	.p2align 4,,10
	.p2align 3
.L297:
	movq	496(%rsp), %rax
	movq	544(%rsp), %r8
	xorl	%edi, %edi
	addq	%r14, %rax
	movq	%rax, 512(%rsp)
	testq	%rbp, %rbp
	jle	.L311
	movq	%r13, 528(%rsp)
	movq	%rbx, %r13
	.p2align 4,,10
	.p2align 3
.L302:
	movq	8(%r15), %rbx
	leaq	1(%rdi), %rdx
	leaq	2(%rdi), %rcx
	movq	(%r15), %r9
	leaq	3(%rdi), %rsi
	prefetcht0	(%r14)
	prefetcht0	(%r8)
	movq	%rbx, %rax
	imulq	%rbx, %rdx
	imulq	%rdi, %rax
	imulq	%rbx, %rcx
	imulq	%rbx, %rsi
	addq	%r10, %rdx
	addq	%r10, %rax
	salq	$4, %rdx
	addq	%r10, %rcx
	salq	$4, %rax
	addq	%r9, %rdx
	addq	%r10, %rsi
	salq	$4, %rcx
	addq	%r9, %rax
	salq	$4, %rsi
	addq	%r9, %rcx
	addq	%r9, %rsi
	movq	%r8, %r9
	testq	%r13, %r13
	jle	.L312
	xorl	%r9d, %r9d
	movq	%r11, -120(%rsp)
	movq	%r14, %rbx
	movq	%r8, %r11
	leaq	768(%r8), %r12
	movq	%rdi, %r8
	movq	%r9, %rdi
	movq	%r10, %r9
	movaps	%xmm9, 456(%rsp)
	movq	%r14, %r10
	movq	%rbp, %r14
	pxor	%xmm5, %xmm5
	movq	%r15, %rbp
	movq	%r14, %r15
	movq	-120(%rsp), %r14
	movaps	%xmm5, 56(%rsp)
	movaps	%xmm5, 40(%rsp)
	movapd	%xmm5, %xmm12
	movapd	%xmm5, %xmm9
	movaps	%xmm5, 24(%rsp)
	movaps	%xmm5, 8(%rsp)
	movaps	%xmm5, -104(%rsp)
	movaps	%xmm5, -88(%rsp)
	movaps	%xmm5, -8(%rsp)
	movaps	%xmm5, -24(%rsp)
	movaps	%xmm5, -56(%rsp)
	movaps	%xmm5, -40(%rsp)
	movaps	%xmm5, -72(%rsp)
	movaps	%xmm5, 88(%rsp)
	movaps	%xmm5, 72(%rsp)
	movaps	%xmm10, 472(%rsp)
	.p2align 4,,10
	.p2align 3
.L299:
#APP
# 1264 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin gebp micro kernel 1/half/quarterX4
# 0 "" 2
#NO_APP
	prefetcht0	(%r12)
#APP
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	-760(%r12), %xmm3
	movapd	(%rbx), %xmm15
	movsd	-728(%r12), %xmm4
	movsd	-720(%r12), %xmm6
	movsd	-712(%r12), %xmm7
	movaps	%xmm15, -120(%rsp)
	movsd	-768(%r12), %xmm2
	movsd	%xmm3, 424(%rsp)
	movsd	-752(%r12), %xmm3
	movsd	-744(%r12), %xmm11
	movsd	-736(%r12), %xmm8
	movsd	%xmm4, 440(%rsp)
	movsd	%xmm3, 432(%rsp)
	movsd	%xmm6, 448(%rsp)
	movsd	%xmm7, 488(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	-704(%r12), %xmm1
	movapd	16(%rbx), %xmm4
	movsd	-696(%r12), %xmm0
	movsd	-688(%r12), %xmm10
	movsd	-680(%r12), %xmm13
	movsd	-672(%r12), %xmm14
	movaps	%xmm4, 344(%rsp)
	movsd	-664(%r12), %xmm3
	movsd	-656(%r12), %xmm6
	movsd	%xmm1, 360(%rsp)
	movsd	-648(%r12), %xmm7
	movsd	%xmm0, 368(%rsp)
	movsd	%xmm10, 376(%rsp)
	movsd	%xmm13, 384(%rsp)
	movsd	%xmm14, 392(%rsp)
	movsd	%xmm3, 400(%rsp)
	movsd	%xmm6, 408(%rsp)
	movsd	%xmm7, 416(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	-632(%r12), %xmm10
	movapd	32(%rbx), %xmm0
	movsd	-640(%r12), %xmm1
	movsd	-600(%r12), %xmm15
	movsd	-584(%r12), %xmm14
	movsd	%xmm10, 320(%rsp)
	movsd	-624(%r12), %xmm13
	movsd	-616(%r12), %xmm10
	movsd	-608(%r12), %xmm7
	movsd	%xmm1, 312(%rsp)
	movsd	-592(%r12), %xmm3
	movsd	%xmm15, 328(%rsp)
	movsd	%xmm14, 336(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	-576(%r12), %xmm1
	movapd	48(%rbx), %xmm6
	movsd	-568(%r12), %xmm15
	movsd	-560(%r12), %xmm14
	movsd	%xmm1, 248(%rsp)
	movsd	-552(%r12), %xmm1
	movsd	-528(%r12), %xmm4
	movsd	%xmm15, 256(%rsp)
	movsd	-544(%r12), %xmm15
	movsd	%xmm14, 264(%rsp)
	movsd	-536(%r12), %xmm14
	movsd	%xmm1, 272(%rsp)
	movsd	-520(%r12), %xmm1
	movaps	%xmm6, 232(%rsp)
	movsd	%xmm15, 280(%rsp)
	movsd	%xmm14, 288(%rsp)
	movsd	%xmm4, 296(%rsp)
	movsd	%xmm1, 304(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
#NO_APP
	prefetcht0	256(%r12)
#APP
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	-496(%r12), %xmm14
	movapd	64(%rbx), %xmm1
	movsd	-512(%r12), %xmm15
	movsd	-488(%r12), %xmm4
	movsd	%xmm14, 200(%rsp)
	movsd	-464(%r12), %xmm14
	movsd	-480(%r12), %xmm6
	movsd	%xmm15, 192(%rsp)
	movsd	-504(%r12), %xmm15
	movsd	%xmm14, 216(%rsp)
	movsd	-456(%r12), %xmm14
	movsd	%xmm4, 208(%rsp)
	movsd	-472(%r12), %xmm4
	movsd	%xmm14, 224(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movapd	80(%rbx), %xmm14
	movaps	%xmm14, 104(%rsp)
	movsd	-448(%r12), %xmm14
	movsd	%xmm14, 128(%rsp)
	movsd	-440(%r12), %xmm14
	movsd	%xmm14, 136(%rsp)
	movsd	-432(%r12), %xmm14
	movsd	%xmm14, 144(%rsp)
	movsd	-424(%r12), %xmm14
	movsd	%xmm14, 152(%rsp)
	movsd	-416(%r12), %xmm14
	movsd	%xmm14, 160(%rsp)
	movsd	-408(%r12), %xmm14
	movsd	%xmm14, 168(%rsp)
	movsd	-400(%r12), %xmm14
	movsd	%xmm14, 176(%rsp)
	movsd	-392(%r12), %xmm14
	movsd	%xmm14, 184(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movapd	%xmm2, %xmm14
	unpcklpd	%xmm15, %xmm15
	unpcklpd	%xmm7, %xmm7
	movsd	312(%rsp), %xmm2
	unpcklpd	%xmm14, %xmm14
	unpcklpd	%xmm8, %xmm8
	unpcklpd	%xmm6, %xmm6
	unpcklpd	%xmm2, %xmm2
	unpcklpd	%xmm11, %xmm11
	unpcklpd	%xmm13, %xmm13
	mulpd	%xmm0, %xmm2
	unpcklpd	%xmm10, %xmm10
	unpcklpd	%xmm4, %xmm4
	mulpd	-120(%rsp), %xmm14
	unpcklpd	%xmm3, %xmm3
	addpd	72(%rsp), %xmm14
	mulpd	%xmm1, %xmm15
	mulpd	%xmm0, %xmm7
	mulpd	%xmm1, %xmm6
	addpd	%xmm2, %xmm14
	mulpd	%xmm0, %xmm13
	movsd	192(%rsp), %xmm2
	mulpd	%xmm0, %xmm10
	mulpd	%xmm1, %xmm4
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm1, %xmm2
	mulpd	%xmm0, %xmm3
	addpd	%xmm2, %xmm14
	movsd	-384(%r12), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	96(%rbx), %xmm2
	addpd	%xmm2, %xmm14
	movsd	424(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	-120(%rsp), %xmm2
	movaps	%xmm14, 72(%rsp)
	addpd	88(%rsp), %xmm2
	movapd	%xmm2, %xmm14
	movsd	320(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm0, %xmm2
	movaps	%xmm2, 88(%rsp)
	movapd	88(%rsp), %xmm2
	addpd	%xmm14, %xmm2
	addpd	%xmm2, %xmm15
	movsd	-376(%r12), %xmm2
	unpcklpd	%xmm2, %xmm2
	movapd	%xmm15, %xmm14
	movapd	96(%rbx), %xmm15
	mulpd	%xmm2, %xmm15
	movsd	-352(%r12), %xmm2
	unpcklpd	%xmm2, %xmm2
	addpd	%xmm14, %xmm15
	movsd	432(%rsp), %xmm14
	unpcklpd	%xmm14, %xmm14
	movaps	%xmm15, 88(%rsp)
	movapd	-120(%rsp), %xmm15
	mulpd	%xmm15, %xmm8
	addpd	-72(%rsp), %xmm8
	mulpd	%xmm15, %xmm14
	mulpd	%xmm15, %xmm11
	addpd	%xmm8, %xmm7
	addpd	%xmm12, %xmm14
	movsd	200(%rsp), %xmm12
	addpd	%xmm6, %xmm7
	movapd	96(%rbx), %xmm6
	addpd	%xmm9, %xmm11
	movsd	208(%rsp), %xmm9
	unpcklpd	%xmm12, %xmm12
	mulpd	%xmm2, %xmm6
	unpcklpd	%xmm9, %xmm9
	addpd	%xmm14, %xmm13
	mulpd	%xmm1, %xmm12
	movapd	%xmm7, %xmm2
	addpd	%xmm11, %xmm10
	mulpd	%xmm1, %xmm9
	addpd	%xmm6, %xmm2
	addpd	%xmm12, %xmm13
	movsd	-368(%r12), %xmm12
	addpd	%xmm9, %xmm10
	movsd	-360(%r12), %xmm9
	movaps	%xmm2, -72(%rsp)
	unpcklpd	%xmm12, %xmm12
	movsd	440(%rsp), %xmm2
	mulpd	96(%rbx), %xmm12
	unpcklpd	%xmm9, %xmm9
	mulpd	96(%rbx), %xmm9
	unpcklpd	%xmm2, %xmm2
	movapd	%xmm2, %xmm6
	movsd	-344(%r12), %xmm2
	mulpd	%xmm15, %xmm6
	unpcklpd	%xmm2, %xmm2
	addpd	%xmm13, %xmm12
	addpd	%xmm10, %xmm9
	addpd	%xmm5, %xmm6
	movsd	328(%rsp), %xmm5
	unpcklpd	%xmm5, %xmm5
	mulpd	%xmm0, %xmm5
	addpd	%xmm6, %xmm5
	addpd	%xmm4, %xmm5
	movapd	96(%rbx), %xmm4
	mulpd	%xmm2, %xmm4
	movsd	448(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	addpd	%xmm4, %xmm5
	movapd	%xmm2, %xmm4
	movsd	216(%rsp), %xmm2
	mulpd	%xmm15, %xmm4
	addpd	-40(%rsp), %xmm4
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm1, %xmm2
	addpd	%xmm4, %xmm3
	addpd	%xmm2, %xmm3
	movsd	-336(%r12), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	96(%rbx), %xmm2
	addpd	%xmm2, %xmm3
	movsd	488(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm15, %xmm2
	movaps	%xmm3, -40(%rsp)
	movapd	-56(%rsp), %xmm3
	addpd	%xmm2, %xmm3
	movsd	336(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm0, %xmm2
	movsd	224(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm0, %xmm1
	addpd	%xmm3, %xmm2
	addpd	%xmm1, %xmm2
	movsd	-328(%r12), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	96(%rbx), %xmm1
	addpd	%xmm1, %xmm2
	movaps	%xmm2, -56(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	360(%rsp), %xmm1
	movapd	-24(%rsp), %xmm2
	movapd	344(%rsp), %xmm4
	movapd	104(%rsp), %xmm7
	unpcklpd	%xmm1, %xmm1
	movapd	112(%rbx), %xmm3
	movapd	232(%rsp), %xmm6
	movsd	128(%rsp), %xmm0
	mulpd	%xmm4, %xmm1
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm7, %xmm0
	addpd	%xmm1, %xmm2
	movsd	248(%rsp), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm6, %xmm1
	addpd	%xmm2, %xmm1
	addpd	%xmm0, %xmm1
	movsd	-320(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	368(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm4, %xmm0
	movaps	%xmm1, -24(%rsp)
	movapd	-8(%rsp), %xmm1
	addpd	%xmm0, %xmm1
	movsd	256(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm6, %xmm0
	addpd	%xmm0, %xmm1
	movsd	136(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm7, %xmm0
	addpd	%xmm0, %xmm1
	movsd	-312(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	376(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm4, %xmm0
	movaps	%xmm1, -8(%rsp)
	movapd	-88(%rsp), %xmm1
	addpd	%xmm0, %xmm1
	movsd	264(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm6, %xmm0
	addpd	%xmm0, %xmm1
	movsd	144(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm7, %xmm0
	addpd	%xmm0, %xmm1
	movsd	-304(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	384(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	movaps	%xmm1, -88(%rsp)
	mulpd	%xmm4, %xmm0
	movsd	272(%rsp), %xmm1
	addpd	-104(%rsp), %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm6, %xmm1
	addpd	%xmm1, %xmm0
	movsd	152(%rsp), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm7, %xmm1
	addpd	%xmm0, %xmm1
	movsd	-296(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	392(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	movaps	%xmm1, -104(%rsp)
	mulpd	%xmm4, %xmm0
	movsd	280(%rsp), %xmm1
	addpd	8(%rsp), %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm6, %xmm1
	addpd	%xmm1, %xmm0
	movsd	160(%rsp), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm7, %xmm1
	addpd	%xmm0, %xmm1
	movsd	-288(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	400(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	movaps	%xmm1, 8(%rsp)
	mulpd	%xmm4, %xmm0
	movsd	288(%rsp), %xmm1
	addpd	24(%rsp), %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm6, %xmm1
	addpd	%xmm1, %xmm0
	movsd	168(%rsp), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm7, %xmm1
	addpd	%xmm0, %xmm1
	movsd	-280(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	408(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	movaps	%xmm1, 24(%rsp)
	mulpd	%xmm4, %xmm0
	movsd	296(%rsp), %xmm1
	addpd	40(%rsp), %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm6, %xmm1
	addpd	%xmm1, %xmm0
	movsd	176(%rsp), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm7, %xmm1
	addpd	%xmm0, %xmm1
	movsd	-272(%r12), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm1
	movsd	416(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	movaps	%xmm1, 40(%rsp)
	mulpd	%xmm4, %xmm0
	movsd	304(%rsp), %xmm1
	addpd	56(%rsp), %xmm0
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm6, %xmm1
	addpd	%xmm0, %xmm1
	movsd	184(%rsp), %xmm0
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm7, %xmm0
	addpd	%xmm1, %xmm0
	movsd	-264(%r12), %xmm1
	unpcklpd	%xmm1, %xmm1
	mulpd	%xmm3, %xmm1
	addpd	%xmm1, %xmm0
	movaps	%xmm0, 56(%rsp)
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
#NO_APP
	subq	$-128, %rbx
#APP
# 1282 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end gebp micro kernel 1/half/quarterX4
# 0 "" 2
#NO_APP
	addq	$8, %rdi
	addq	$512, %r12
	cmpq	%r13, %rdi
	jl	.L299
	movq	%r15, %rbx
	movapd	%xmm5, %xmm6
	movapd	-88(%rsp), %xmm3
	movq	%rbp, %r15
	movapd	-104(%rsp), %xmm5
	movapd	%xmm12, %xmm15
	movq	%rbx, %rbp
	movapd	%xmm9, %xmm12
	movapd	-24(%rsp), %xmm0
	movapd	-8(%rsp), %xmm1
	movq	%r8, %rdi
	addpd	%xmm15, %xmm3
	movq	520(%rsp), %rbx
	movq	%r11, %r8
	movq	%r14, %r11
	movapd	8(%rsp), %xmm7
	movapd	40(%rsp), %xmm8
	movapd	56(%rsp), %xmm4
	addpd	%xmm12, %xmm5
	movq	%r10, %r14
	addpd	72(%rsp), %xmm0
	addpd	88(%rsp), %xmm1
	movq	%r9, %r10
	addpd	-72(%rsp), %xmm7
	addpd	24(%rsp), %xmm6
	addpd	-40(%rsp), %xmm8
	addpd	-56(%rsp), %xmm4
	movapd	456(%rsp), %xmm9
	leaq	(%r8,%rbx), %r9
	movdqa	472(%rsp), %xmm10
	movq	512(%rsp), %rbx
.L298:
	cmpq	%r11, %r13
	jge	.L300
	movq	%r13, %r12
	.p2align 4,,10
	.p2align 3
.L301:
#APP
# 1196 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1X4
# 0 "" 2
# 1197 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	(%r9), %xmm11
	movapd	(%rbx), %xmm2
	unpcklpd	%xmm11, %xmm11
	mulpd	%xmm2, %xmm11
	addpd	%xmm11, %xmm0
	movsd	8(%r9), %xmm11
	unpcklpd	%xmm11, %xmm11
	mulpd	%xmm2, %xmm11
	addpd	%xmm11, %xmm1
	movsd	16(%r9), %xmm11
	unpcklpd	%xmm11, %xmm11
	mulpd	%xmm2, %xmm11
	addpd	%xmm11, %xmm3
	movsd	24(%r9), %xmm11
	unpcklpd	%xmm11, %xmm11
	mulpd	%xmm2, %xmm11
	addpd	%xmm11, %xmm5
	movsd	32(%r9), %xmm11
	unpcklpd	%xmm11, %xmm11
	mulpd	%xmm2, %xmm11
	addpd	%xmm11, %xmm7
	movsd	40(%r9), %xmm11
	unpcklpd	%xmm11, %xmm11
	mulpd	%xmm2, %xmm11
	addpd	%xmm11, %xmm6
	movsd	48(%r9), %xmm11
	unpcklpd	%xmm11, %xmm11
	mulpd	%xmm2, %xmm11
	addpd	%xmm11, %xmm8
	movsd	56(%r9), %xmm11
	unpcklpd	%xmm11, %xmm11
	mulpd	%xmm11, %xmm2
	addpd	%xmm2, %xmm4
#APP
# 1207 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1X4
# 0 "" 2
#NO_APP
	addq	$1, %r12
	addq	$64, %r9
	addq	$16, %rbx
	cmpq	%r11, %r12
	jne	.L301
.L300:
	xorpd	.LC7(%rip), %xmm5
	xorpd	.LC7(%rip), %xmm1
	addq	$4, %rdi
	xorpd	.LC7(%rip), %xmm4
	xorpd	.LC7(%rip), %xmm6
	shufpd	$1, %xmm5, %xmm5
	addpd	%xmm3, %xmm5
	shufpd	$1, %xmm1, %xmm1
	addpd	%xmm1, %xmm0
	movupd	(%rdx), %xmm3
	shufpd	$1, %xmm4, %xmm4
	addpd	%xmm8, %xmm4
	shufpd	$1, %xmm6, %xmm6
	addpd	%xmm7, %xmm6
	pshufd	$238, %xmm5, %xmm1
	pshufd	$68, %xmm5, %xmm5
	pshufd	$238, %xmm0, %xmm2
	mulpd	%xmm9, %xmm5
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm10, %xmm1
	mulpd	%xmm9, %xmm0
	xorpd	.LC8(%rip), %xmm1
	mulpd	%xmm10, %xmm2
	xorpd	.LC8(%rip), %xmm2
	addpd	%xmm1, %xmm5
	pshufd	$238, %xmm4, %xmm1
	pshufd	$68, %xmm4, %xmm4
	mulpd	%xmm9, %xmm4
	addpd	%xmm2, %xmm0
	mulpd	%xmm10, %xmm1
	addpd	%xmm3, %xmm5
	movupd	(%rax), %xmm3
	xorpd	.LC8(%rip), %xmm1
	addpd	%xmm3, %xmm0
	addpd	%xmm1, %xmm4
	movups	%xmm0, (%rax)
	pshufd	$238, %xmm6, %xmm0
	pshufd	$68, %xmm6, %xmm6
	movq	504(%rsp), %rax
	mulpd	%xmm9, %xmm6
	movups	%xmm5, (%rdx)
	movupd	(%rsi), %xmm5
	mulpd	%xmm10, %xmm0
	addq	%rax, %r8
	xorpd	.LC8(%rip), %xmm0
	addpd	%xmm5, %xmm4
	movupd	(%rcx), %xmm5
	addpd	%xmm0, %xmm6
	addpd	%xmm5, %xmm6
	movups	%xmm6, (%rcx)
	movups	%xmm4, (%rsi)
	cmpq	%rbp, %rdi
	jl	.L302
	movq	%r13, %rbx
	movq	528(%rsp), %r13
.L311:
	movq	496(%rsp), %rax
	movq	536(%rsp), %r8
	movq	%rbp, %r9
	leaq	(%r14,%rax), %r12
	cmpq	672(%rsp), %rbp
	jge	.L309
	movq	%rbp, 688(%rsp)
	movq	672(%rsp), %rbp
	.p2align 4,,10
	.p2align 3
.L308:
	movq	8(%r15), %rax
	movq	%r14, %rdx
	prefetcht0	(%r14)
	imulq	%r9, %rax
	addq	%r10, %rax
	salq	$4, %rax
	addq	(%r15), %rax
	movq	%rax, %rsi
	movq	%r8, %rax
	testq	%rbx, %rbx
	jle	.L313
	pxor	%xmm0, %xmm0
	movq	%r12, %rcx
	xorl	%edi, %edi
	movq	%r10, %r12
	movapd	%xmm0, %xmm1
	movq	%r9, %r10
	movaps	%xmm9, -8(%rsp)
	movq	%rsi, %r9
	movaps	%xmm10, 8(%rsp)
	movq	%r8, %rsi
	movq	%r13, %r8
	.p2align 4,,10
	.p2align 3
.L305:
#APP
# 1336 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1350 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1350 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	8(%rax), %xmm5
	movupd	(%rdx), %xmm11
	movq	(%rax), %r13
	movsd	%xmm5, -120(%rsp)
#APP
# 1350 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1351 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1351 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	16(%rax), %xmm5
	movsd	24(%rax), %xmm3
	movupd	16(%rdx), %xmm8
	movsd	%xmm5, -104(%rsp)
	movsd	%xmm3, -88(%rsp)
#APP
# 1351 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1352 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1352 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	32(%rax), %xmm6
	movsd	40(%rax), %xmm4
	movupd	32(%rdx), %xmm7
	movsd	%xmm6, -72(%rsp)
	movsd	%xmm4, -56(%rsp)
#APP
# 1352 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1353 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1353 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	48(%rax), %xmm2
	movupd	48(%rdx), %xmm6
	movsd	56(%rax), %xmm15
	movsd	%xmm2, -40(%rsp)
#APP
# 1353 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1354 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1354 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	64(%rax), %xmm10
	movupd	64(%rdx), %xmm5
	movsd	72(%rax), %xmm14
	movsd	%xmm10, -24(%rsp)
#APP
# 1354 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1355 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1355 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movupd	80(%rdx), %xmm4
	movsd	80(%rax), %xmm10
	movsd	88(%rax), %xmm13
#APP
# 1355 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1356 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1356 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movupd	96(%rdx), %xmm3
	movsd	96(%rax), %xmm9
	movsd	104(%rax), %xmm12
#APP
# 1356 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1357 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1357 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movq	%r13, %xmm2
	unpcklpd	%xmm10, %xmm10
	unpcklpd	%xmm9, %xmm9
	unpcklpd	%xmm2, %xmm2
	unpcklpd	%xmm15, %xmm15
	unpcklpd	%xmm14, %xmm14
	mulpd	%xmm11, %xmm2
	unpcklpd	%xmm13, %xmm13
	unpcklpd	%xmm12, %xmm12
	mulpd	%xmm4, %xmm10
	mulpd	%xmm3, %xmm9
	mulpd	%xmm6, %xmm15
	mulpd	%xmm5, %xmm14
	addpd	%xmm2, %xmm1
	movsd	-104(%rsp), %xmm2
	mulpd	%xmm4, %xmm13
	mulpd	%xmm3, %xmm12
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm8, %xmm2
	addpd	%xmm2, %xmm1
	movsd	-72(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm7, %xmm2
	addpd	%xmm2, %xmm1
	movsd	-40(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm6, %xmm2
	addpd	%xmm2, %xmm1
	movsd	-24(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm5, %xmm2
	movupd	112(%rdx), %xmm5
	addpd	%xmm2, %xmm1
	movupd	112(%rdx), %xmm2
	addpd	%xmm10, %xmm1
	addpd	%xmm9, %xmm1
	movsd	112(%rax), %xmm9
	unpcklpd	%xmm9, %xmm9
	mulpd	%xmm2, %xmm9
	movsd	-120(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm11, %xmm2
	addpd	%xmm9, %xmm1
	addpd	%xmm2, %xmm0
	movsd	-88(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm8, %xmm2
	addpd	%xmm2, %xmm0
	movsd	-56(%rsp), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm7, %xmm2
	addpd	%xmm2, %xmm0
	movsd	120(%rax), %xmm2
	unpcklpd	%xmm2, %xmm2
	addpd	%xmm15, %xmm0
	mulpd	%xmm5, %xmm2
	addpd	%xmm14, %xmm0
	addpd	%xmm13, %xmm0
	addpd	%xmm12, %xmm0
	addpd	%xmm2, %xmm0
#APP
# 1357 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
#NO_APP
	subq	$-128, %rax
	subq	$-128, %rdx
#APP
# 1362 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end gebp micro kernel 1/half/quarterX1
# 0 "" 2
#NO_APP
	addq	$8, %rdi
	cmpq	%rdi, %rbx
	jg	.L305
	movq	496(%rsp), %rax
	movq	%r8, %r13
	movq	%rsi, %r8
	movq	%rcx, %rdx
	movapd	-8(%rsp), %xmm9
	movdqa	8(%rsp), %xmm10
	movq	%r9, %rsi
	movq	%r10, %r9
	addq	%r8, %rax
	movq	%r12, %r10
	movq	%rcx, %r12
.L304:
	cmpq	%r11, %rbx
	jge	.L306
	movq	%rbx, %rcx
	.p2align 4,,10
	.p2align 3
.L307:
#APP
# 1369 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#begin step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
# 1369 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#Note: these asm comments work around bug 935!
# 0 "" 2
#NO_APP
	movsd	(%rax), %xmm2
	movupd	(%rdx), %xmm5
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm5, %xmm2
	addpd	%xmm2, %xmm1
	movsd	8(%rax), %xmm2
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm5, %xmm2
	addpd	%xmm2, %xmm0
#APP
# 1369 "/usr/local/include/Eigen/src/Core/products/GeneralBlockPanelKernel.h" 1
	#end step of gebp micro kernel 1/half/quarterX1
# 0 "" 2
#NO_APP
	addq	$1, %rcx
	addq	$16, %rax
	addq	$16, %rdx
	cmpq	%rcx, %r11
	jne	.L307
.L306:
	xorpd	.LC7(%rip), %xmm0
	movupd	(%rsi), %xmm5
	addq	$1, %r9
	addq	%r13, %r8
	shufpd	$1, %xmm0, %xmm0
	addpd	%xmm1, %xmm0
	pshufd	$238, %xmm0, %xmm1
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm9, %xmm0
	mulpd	%xmm10, %xmm1
	xorpd	.LC8(%rip), %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm5, %xmm0
	movups	%xmm0, (%rsi)
	cmpq	%rbp, %r9
	jne	.L308
	movq	688(%rsp), %rbp
.L309:
	movq	552(%rsp), %rax
	addq	$1, %r10
	addq	%rax, %r14
	movq	560(%rsp), %rax
	cmpq	%rax, %r10
	jne	.L297
.L295:
	addq	$592, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L312:
	.cfi_restore_state
	pxor	%xmm4, %xmm4
	movq	%r14, %rbx
	movapd	%xmm4, %xmm8
	movapd	%xmm4, %xmm6
	movapd	%xmm4, %xmm7
	movapd	%xmm4, %xmm5
	movapd	%xmm4, %xmm3
	movapd	%xmm4, %xmm1
	movapd	%xmm4, %xmm0
	jmp	.L298
	.p2align 4,,10
	.p2align 3
.L313:
	pxor	%xmm0, %xmm0
	movapd	%xmm0, %xmm1
	jmp	.L304
	.cfi_endproc
.LFE14870:
	.size	_ZN5Eigen8internal22lhs_process_one_packetILi4ELl1ELl1ESt7complexIdES3_S3_NS0_12DoublePacketIDv2_dEES5_S6_NS0_9Packet1cdENS0_11gebp_traitsIS3_S3_Lb0ELb0ELi1ELi0EEENS0_16BlasLinearMapperIS3_lLi0ELi1EEENS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEEEclERKSD_PKS3_SI_S3_llllllilllll.constprop.0.isra.0, .-_ZN5Eigen8internal22lhs_process_one_packetILi4ELl1ELl1ESt7complexIdES3_S3_NS0_12DoublePacketIDv2_dEES5_S6_NS0_9Packet1cdENS0_11gebp_traitsIS3_S3_Lb0ELb0ELi1ELi0EEENS0_16BlasLinearMapperIS3_lLi0ELi1EEENS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEEEclERKSD_PKS3_SI_S3_llllllilllll.constprop.0.isra.0
	.align 2
	.p2align 4
	.type	_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb0ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0, @function
_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb0ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0:
.LFB14879:
	.cfi_startproc
	movq	%rdx, %r10
	movq	%r9, %rdx
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	sarq	$63, %rdx
	shrq	$62, %rdx
	leaq	(%r9,%rdx), %rax
	andl	$3, %eax
	subq	%rdx, %rax
	movq	%r9, %rdx
	subq	%rax, %rdx
	movq	%r8, %rax
	andq	$-8, %rax
	pushq	%rdx
	.cfi_def_cfa_offset 24
	movq	%r10, %rdx
	pushq	%r8
	.cfi_def_cfa_offset 32
	pushq	%r9
	.cfi_def_cfa_offset 40
	movq	%r8, %r9
	pushq	%rax
	.cfi_def_cfa_offset 48
	pushq	$0
	.cfi_def_cfa_offset 56
	pushq	$0
	.cfi_def_cfa_offset 64
	call	_ZN5Eigen8internal22lhs_process_one_packetILi4ELl1ELl1ESt7complexIdES3_S3_NS0_12DoublePacketIDv2_dEES5_S6_NS0_9Packet1cdENS0_11gebp_traitsIS3_S3_Lb0ELb0ELi1ELi0EEENS0_16BlasLinearMapperIS3_lLi0ELi1EEENS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEEEclERKSD_PKS3_SI_S3_llllllilllll.constprop.0.isra.0
	addq	$56, %rsp
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE14879:
	.size	_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb0ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0, .-_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb0ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0
	.p2align 4
	.type	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0, @function
_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0:
.LFB14881:
	.cfi_startproc
	testq	%rdi, %rdi
	je	.L330
	pushq	%rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	movq	%rdi, %rbx
	call	strlen@PLT
	movq	%rbx, %rsi
	leaq	_ZSt4cout(%rip), %rdi
	popq	%rbx
	.cfi_restore 3
	.cfi_def_cfa_offset 8
	movq	%rax, %rdx
	jmp	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	.p2align 4,,10
	.p2align 3
.L330:
	movq	_ZSt4cout(%rip), %rax
	leaq	_ZSt4cout(%rip), %rdi
	addq	-24(%rax), %rdi
	movl	32(%rdi), %esi
	orl	$1, %esi
	jmp	_ZNSt9basic_iosIcSt11char_traitsIcEE5clearESt12_Ios_Iostate@PLT
	.cfi_endproc
.LFE14881:
	.size	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0, .-_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	.p2align 4
	.type	_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE.isra.0, @function
_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE.isra.0:
.LFB14882:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movq	(%rdi), %rax
	movq	-24(%rax), %rax
	movq	240(%rdi,%rax), %rbp
	testq	%rbp, %rbp
	je	.L337
	cmpb	$0, 56(%rbp)
	movq	%rdi, %rbx
	movq	%rsi, %r12
	je	.L333
	movsbl	67(%rbp), %edx
.L334:
	movq	%r12, %rsi
	movq	%rbx, %rdi
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	jmp	_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EES4_@PLT
	.p2align 4,,10
	.p2align 3
.L333:
	.cfi_restore_state
	movq	%rbp, %rdi
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	0(%rbp), %rax
	movl	$10, %edx
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rcx
	movq	48(%rax), %rax
	cmpq	%rcx, %rax
	je	.L334
	movl	$10, %esi
	movq	%rbp, %rdi
	call	*%rax
	movsbl	%al, %edx
	jmp	.L334
.L337:
	call	_ZSt16__throw_bad_castv@PLT
	.cfi_endproc
.LFE14882:
	.size	_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE.isra.0, .-_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE.isra.0
	.p2align 4
	.type	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0, @function
_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0:
.LFB14884:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$8, %rsp
	.cfi_def_cfa_offset 32
	movq	(%rdi), %rax
	movq	-24(%rax), %rax
	movq	240(%rdi,%rax), %rbp
	testq	%rbp, %rbp
	je	.L344
	cmpb	$0, 56(%rbp)
	movq	%rdi, %rbx
	je	.L340
	movsbl	67(%rbp), %esi
.L341:
	movq	%rbx, %rdi
	call	_ZNSo3putEc@PLT
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	movq	%rax, %rdi
	popq	%rbp
	.cfi_def_cfa_offset 8
	jmp	_ZNSo5flushEv@PLT
	.p2align 4,,10
	.p2align 3
.L340:
	.cfi_restore_state
	movq	%rbp, %rdi
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	0(%rbp), %rax
	movl	$10, %esi
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rdx
	movq	48(%rax), %rax
	cmpq	%rdx, %rax
	je	.L341
	movq	%rbp, %rdi
	call	*%rax
	movsbl	%al, %esi
	jmp	.L341
.L344:
	call	_ZSt16__throw_bad_castv@PLT
	.cfi_endproc
.LFE14884:
	.size	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0, .-_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0
	.section	.text._ZNSt7__cxx119to_stringEi,"axG",@progbits,_ZNSt7__cxx119to_stringEi,comdat
	.p2align 4
	.weak	_ZNSt7__cxx119to_stringEi
	.type	_ZNSt7__cxx119to_stringEi, @function
_ZNSt7__cxx119to_stringEi:
.LFB1191:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA1191
	pushq	%r14
	.cfi_def_cfa_offset 16
	.cfi_offset 14, -16
	movl	%esi, %r14d
	pushq	%r13
	.cfi_def_cfa_offset 24
	.cfi_offset 13, -24
	shrl	$31, %r14d
	movq	%rdi, %r13
	pushq	%r12
	.cfi_def_cfa_offset 32
	.cfi_offset 12, -32
	pushq	%rbp
	.cfi_def_cfa_offset 40
	.cfi_offset 6, -40
	movl	%esi, %ebp
	negl	%ebp
	pushq	%rbx
	.cfi_def_cfa_offset 48
	.cfi_offset 3, -48
	cmovs	%esi, %ebp
	cmpl	$9, %ebp
	jbe	.L347
	cmpl	$99, %ebp
	jbe	.L348
	cmpl	$999, %ebp
	jbe	.L363
	movl	%ebp, %ebx
	cmpl	$9999, %ebp
	jbe	.L373
	movl	$5, %esi
	cmpq	$99999, %rbx
	jbe	.L356
	cmpq	$999999, %rbx
	jbe	.L374
	movl	$6, %r12d
	movl	$7, %esi
	cmpq	$9999999, %rbx
	jbe	.L355
	cmpq	$99999999, %rbx
	jbe	.L365
	cmpq	$999999999, %rbx
	jbe	.L366
	movl	$9, %r12d
.L357:
	leaq	16(%r13), %rax
	leal	1(%r12,%r14), %esi
	movq	%rax, 0(%r13)
.L372:
	movl	$45, %edx
	movq	%r13, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructEmc@PLT
	movzbl	%r14b, %edx
	leaq	_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits(%rip), %rcx
	addq	0(%r13), %rdx
	jmp	.L360
	.p2align 4,,10
	.p2align 3
.L375:
	movl	%ebx, %ebx
.L360:
	imulq	$1374389535, %rbx, %rbx
	movl	%ebp, %eax
	movl	%r12d, %edi
	shrq	$37, %rbx
	imull	$100, %ebx, %esi
	subl	%esi, %eax
	movl	%ebp, %esi
	movl	%ebx, %ebp
	addl	%eax, %eax
	leal	1(%rax), %r8d
	movzbl	(%rcx,%rax), %eax
	movzbl	(%rcx,%r8), %r8d
	movb	%r8b, (%rdx,%rdi)
	leal	-1(%r12), %edi
	subl	$2, %r12d
	movb	%al, (%rdx,%rdi)
	cmpl	$9999, %esi
	ja	.L375
	cmpl	$999, %esi
	jbe	.L359
.L362:
	addl	%ebp, %ebp
	leal	1(%rbp), %eax
	movzbl	(%rcx,%rbp), %ebp
	movzbl	(%rcx,%rax), %eax
	movb	%bpl, (%rdx)
	movb	%al, 1(%rdx)
	movq	%r13, %rax
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	popq	%rbp
	.cfi_def_cfa_offset 32
	popq	%r12
	.cfi_def_cfa_offset 24
	popq	%r13
	.cfi_def_cfa_offset 16
	popq	%r14
	.cfi_def_cfa_offset 8
	ret
.L347:
	.cfi_restore_state
	leaq	16(%rdi), %rax
	movl	$45, %edx
	leal	1(%r14), %esi
	movq	%rax, (%rdi)
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructEmc@PLT
	movzbl	%r14b, %edx
	addq	0(%r13), %rdx
.L359:
	addl	$48, %ebp
	movq	%r13, %rax
	movb	%bpl, (%rdx)
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	popq	%rbp
	.cfi_def_cfa_offset 32
	popq	%r12
	.cfi_def_cfa_offset 24
	popq	%r13
	.cfi_def_cfa_offset 16
	popq	%r14
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L366:
	.cfi_restore_state
	movl	$9, %esi
.L356:
	leal	-1(%rsi), %r12d
.L355:
	leaq	16(%r13), %rax
	addl	%r14d, %esi
	movq	%rax, 0(%r13)
	jmp	.L372
	.p2align 4,,10
	.p2align 3
.L365:
	movl	$8, %esi
	leaq	16(%r13), %rax
	movl	$7, %r12d
	movq	%rax, 0(%r13)
	addl	%r14d, %esi
	jmp	.L372
.L348:
	leaq	16(%rdi), %rax
	movl	$45, %edx
	leal	2(%r14), %esi
	movq	%rax, (%rdi)
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructEmc@PLT
	movzbl	%r14b, %edx
	leaq	_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits(%rip), %rcx
	addq	0(%r13), %rdx
	jmp	.L362
.L373:
	movl	$3, %r12d
	movl	$4, %esi
	jmp	.L355
.L363:
	movl	$2, %r12d
	movl	$3, %esi
	movl	%ebp, %ebx
	jmp	.L355
.L374:
	movl	$5, %r12d
	jmp	.L357
	.cfi_endproc
.LFE1191:
	.globl	__gxx_personality_v0
	.section	.gcc_except_table._ZNSt7__cxx119to_stringEi,"aG",@progbits,_ZNSt7__cxx119to_stringEi,comdat
.LLSDA1191:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE1191-.LLSDACSB1191
.LLSDACSB1191:
.LLSDACSE1191:
	.section	.text._ZNSt7__cxx119to_stringEi,"axG",@progbits,_ZNSt7__cxx119to_stringEi,comdat
	.size	_ZNSt7__cxx119to_stringEi, .-_ZNSt7__cxx119to_stringEi
	.section	.text.unlikely._ZN5Eigen8internal19throw_std_bad_allocEv,"axG",@progbits,_ZN5Eigen8internal19throw_std_bad_allocEv,comdat
	.weak	_ZN5Eigen8internal19throw_std_bad_allocEv
	.type	_ZN5Eigen8internal19throw_std_bad_allocEv, @function
_ZN5Eigen8internal19throw_std_bad_allocEv:
.LFB4700:
	.cfi_startproc
	pushq	%rax
	.cfi_def_cfa_offset 16
	movl	$8, %edi
	call	__cxa_allocate_exception@PLT
	movq	_ZNSt9bad_allocD1Ev@GOTPCREL(%rip), %rdx
	leaq	_ZTISt9bad_alloc(%rip), %rsi
	movq	%rax, %rdi
	leaq	16+_ZTVSt9bad_alloc(%rip), %rax
	movq	%rax, (%rdi)
	call	__cxa_throw@PLT
	.cfi_endproc
.LFE4700:
	.size	_ZN5Eigen8internal19throw_std_bad_allocEv, .-_ZN5Eigen8internal19throw_std_bad_allocEv
	.section	.rodata.str1.8
	.align 8
.LC25:
	.string	"void* Eigen::internal::aligned_malloc(std::size_t)"
	.align 8
.LC26:
	.string	"/usr/local/include/Eigen/src/Core/util/Memory.h"
	.align 8
.LC27:
	.string	"(size<16 || (std::size_t(result)%16)==0) && \"System's malloc returned an unaligned pointer. Compile with EIGEN_MALLOC_ALREADY_ALIGNED=0 to fallback to handmade aligned memory allocator.\""
	.align 8
.LC28:
	.string	"Eigen::MapBase<Derived, 0>::MapBase(PointerType, Eigen::Index, Eigen::Index) [with Derived = Eigen::Block<Eigen::Block<Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true>, -1, 1, true>; PointerType = std::complex<double>*; Eigen::Index = long int]"
	.section	.text.unlikely,"ax",@progbits
.LCOLDB29:
	.text
.LHOTB29:
	.p2align 4
	.type	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0, @function
_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0:
.LFB14885:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movapd	%xmm1, %xmm2
	movapd	%xmm0, %xmm4
	movapd	%xmm0, %xmm5
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r14
	pushq	%r13
	pushq	%r12
	.cfi_offset 14, -24
	.cfi_offset 13, -32
	.cfi_offset 12, -40
	movq	%rdx, %r12
	pushq	%rbx
	subq	$64, %rsp
	.cfi_offset 3, -48
	movsd	.LC23(%rip), %xmm3
	movq	(%rdi), %r13
	movq	%fs:40, %rax
	movq	%rax, -40(%rbp)
	xorl	%eax, %eax
	mulsd	%xmm3, %xmm2
	mulsd	%xmm3, %xmm4
	subsd	%xmm2, %xmm5
	addsd	%xmm1, %xmm4
	ucomisd	%xmm5, %xmm4
	jp	.L394
.L379:
	pxor	%xmm1, %xmm1
	movapd	%xmm4, %xmm2
	movapd	%xmm5, %xmm0
	mulsd	%xmm1, %xmm2
	mulsd	%xmm5, %xmm1
	subsd	%xmm2, %xmm0
	addsd	%xmm4, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L395
.L380:
	movq	8(%rsi), %rdx
	movq	%rdx, %rax
	shrq	$60, %rax
	jne	.L385
	movq	(%rsi), %rax
	movq	%rdx, %rbx
	xorl	%r14d, %r14d
	salq	$4, %rbx
	testq	%rax, %rax
	je	.L396
.L382:
	cmpq	$0, 8(%r12)
	movq	(%r12), %r8
	jns	.L386
	testq	%r8, %r8
	jne	.L397
.L386:
	movq	%rax, -64(%rbp)
	movq	8(%r13), %rsi
	leaq	-64(%rbp), %rcx
	leaq	-80(%rbp), %rdx
	movq	0(%r13), %rax
	movq	16(%r13), %rdi
	movq	$1, -56(%rbp)
	movq	%rsi, -72(%rbp)
	movq	%rax, -80(%rbp)
	call	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb1ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_.constprop.0
	cmpq	$131072, %rbx
	ja	.L398
.L378:
	movq	-40(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L399
	leaq	-32(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L398:
	.cfi_restore_state
	movq	%r14, %rdi
	call	free@PLT
	jmp	.L378
	.p2align 4,,10
	.p2align 3
.L396:
	cmpq	$131072, %rbx
	ja	.L383
	leaq	32(%rbx), %rax
	subq	%rax, %rsp
	leaq	15(%rsp), %r14
	andq	$-16, %r14
	movq	%r14, %rax
	jmp	.L382
	.p2align 4,,10
	.p2align 3
.L383:
	movq	%rbx, %rdi
	movsd	%xmm0, -96(%rbp)
	movsd	%xmm1, -88(%rbp)
	call	malloc@PLT
	movsd	-88(%rbp), %xmm1
	movsd	-96(%rbp), %xmm0
	testb	$15, %al
	movq	%rax, %r14
	je	.L384
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L397:
	leaq	.LC28(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
	.p2align 4,,10
	.p2align 3
.L384:
	testq	%rax, %rax
	jne	.L382
	jmp	.L385
.L395:
	pxor	%xmm3, %xmm3
	movapd	%xmm4, %xmm1
	movapd	%xmm5, %xmm0
	movq	%rsi, -88(%rbp)
	movsd	.LC24(%rip), %xmm2
	call	__muldc3@PLT
	movq	-88(%rbp), %rsi
	jmp	.L380
.L394:
	movsd	.LC24(%rip), %xmm2
	movq	%rsi, -88(%rbp)
	call	__muldc3@PLT
	movq	-88(%rbp), %rsi
	movapd	%xmm0, %xmm5
	movapd	%xmm1, %xmm4
	jmp	.L379
.L399:
	call	__stack_chk_fail@PLT
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.type	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0.cold, @function
_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0.cold:
.LFSB14885:
.L385:
	.cfi_def_cfa 6, 16
	.cfi_offset 3, -48
	.cfi_offset 6, -16
	.cfi_offset 12, -40
	.cfi_offset 13, -32
	.cfi_offset 14, -24
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
	.cfi_endproc
.LFE14885:
	.text
	.size	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0, .-_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0
	.section	.text.unlikely
	.size	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0.cold, .-_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0.cold
.LCOLDE29:
	.text
.LHOTE29:
	.section	.text._ZN5Eigen8internal14aligned_mallocEm,"axG",@progbits,_ZN5Eigen8internal14aligned_mallocEm,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal14aligned_mallocEm
	.type	_ZN5Eigen8internal14aligned_mallocEm, @function
_ZN5Eigen8internal14aligned_mallocEm:
.LFB4705:
	.cfi_startproc
	pushq	%rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	movq	%rdi, %rbx
	call	malloc@PLT
	cmpq	$15, %rbx
	jbe	.L401
	testb	$15, %al
	je	.L401
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
	.p2align 4,,10
	.p2align 3
.L401:
	testq	%rax, %rax
	je	.L413
.L400:
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L413:
	.cfi_restore_state
	testq	%rbx, %rbx
	je	.L400
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
	.cfi_endproc
.LFE4705:
	.size	_ZN5Eigen8internal14aligned_mallocEm, .-_ZN5Eigen8internal14aligned_mallocEm
	.section	.rodata.str1.8
	.align 8
.LC30:
	.string	"Eigen::internal::blas_data_mapper<Scalar, Index, StorageOrder, AlignmentType, 1>::blas_data_mapper(Scalar*, Index, Index) [with Scalar = std::complex<double>; Index = long int; int StorageOrder = 0; int AlignmentType = 0]"
	.align 8
.LC31:
	.string	"/usr/local/include/Eigen/src/Core/util/BlasUtil.h"
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC32:
	.string	"incr==1"
	.section	.text.unlikely
.LCOLDB33:
	.text
.LHOTB33:
	.p2align 4
	.type	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0, @function
_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0:
.LFB14888:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA14888
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$264, %rsp
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	24(%rbp), %rax
	movq	%rdx, -168(%rbp)
	movq	%rcx, -96(%rbp)
	movq	48(%rbp), %r12
	movq	%r8, -192(%rbp)
	movq	%r9, -120(%rbp)
	movq	%rax, -136(%rbp)
	movsd	%xmm0, -128(%rbp)
	movsd	%xmm1, -152(%rbp)
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	xorl	%eax, %eax
	cmpq	$1, 32(%rbp)
	jne	.L474
	movq	16(%r12), %r15
	movq	%rdi, %rcx
	movq	%rdi, %r14
	movq	32(%r12), %rax
	movabsq	$1152921504606846975, %r13
	movq	%rsi, %rbx
	cmpq	%rdi, %r15
	movq	24(%r12), %rdi
	movq	%rax, %rsi
	movq	%rax, -184(%rbp)
	cmovle	%r15, %rcx
	movq	%rsi, %rdx
	cmpq	%rbx, %rdi
	movq	%rdi, -88(%rbp)
	cmovg	%rbx, %rdi
	imulq	%rcx, %rax
	movq	%rcx, -248(%rbp)
	movq	%rdi, -112(%rbp)
	imulq	%rdi, %rdx
	cmpq	%rax, %r13
	jb	.L467
	movq	(%r12), %rsi
	salq	$4, %rax
	movq	%rax, -280(%rbp)
	movq	%rsi, -104(%rbp)
	testq	%rsi, %rsi
	je	.L475
	cmpq	%rdx, %r13
	jb	.L470
	movq	$0, -272(%rbp)
.L419:
	movq	8(%r12), %r13
	salq	$4, %rdx
	movq	%rdx, -288(%rbp)
	testq	%r13, %r13
	je	.L476
	xorl	%edi, %edi
	xorl	%eax, %eax
	cmpq	%r14, %r15
	jge	.L431
.L425:
	movq	-168(%rbp), %rsi
	cmpq	%rsi, -184(%rbp)
	sete	%al
	cmpq	%rbx, -88(%rbp)
	setge	%dl
	andl	%edx, %eax
.L431:
	testq	%r14, %r14
	jle	.L434
.L428:
	cmpq	$0, -168(%rbp)
	jle	.L434
	movq	-248(%rbp), %rsi
	movq	-96(%rbp), %rcx
	xorl	$1, %eax
	leaq	-80(%rbp), %r12
	movq	-192(%rbp), %rdx
	movb	%al, -202(%rbp)
	movq	%rdi, -296(%rbp)
	imulq	%rsi, %rdx
	salq	$4, %rsi
	movq	%rsi, -264(%rbp)
	movq	-136(%rbp), %rsi
	salq	$4, %rdx
	movq	%rsi, -200(%rbp)
	movq	16(%rbp), %rsi
	movq	%rdx, -256(%rbp)
	movq	-184(%rbp), %rdx
	salq	$4, %rdx
	subq	%rdx, %rcx
	movq	%rcx, -176(%rbp)
	movq	-112(%rbp), %rcx
	imulq	%rcx, %rsi
	imulq	40(%rbp), %rcx
	salq	$4, %rsi
	movq	%rsi, -160(%rbp)
	movq	-120(%rbp), %rsi
	subq	%rdx, %rsi
	movq	%rcx, %rdx
	movq	%r14, %rcx
	salq	$4, %rdx
	movq	%rsi, -224(%rbp)
	xorl	%esi, %esi
	movq	%rdx, -120(%rbp)
	movq	%rsi, %r14
	.p2align 4,,10
	.p2align 3
.L436:
	movq	-248(%rbp), %rax
	movq	%r14, %rdx
	movq	%rcx, %r15
	movq	%rcx, -240(%rbp)
	addq	%rax, %r14
	cmpq	%rcx, %r14
	movq	%r14, -232(%rbp)
	cmovle	%r14, %r15
	movq	%rbx, %r14
	subq	%rdx, %r15
	testq	%rdx, %rdx
	sete	%dl
	orb	-202(%rbp), %dl
	movq	%r15, -96(%rbp)
	xorl	%r8d, %r8d
	movb	%dl, -201(%rbp)
	.p2align 4,,10
	.p2align 3
.L433:
	movq	-184(%rbp), %rax
	movq	%r8, %rcx
	movq	-104(%rbp), %rdi
	movq	%r12, %rsi
	addq	%rax, %r8
	movq	-168(%rbp), %rax
	movq	%r8, %rbx
	movq	%r8, -88(%rbp)
	cmpq	%rax, %r8
	movq	%rax, %rdx
	movq	-176(%rbp), %rax
	cmovle	%r8, %rdx
	salq	$4, %rbx
	subq	%rcx, %rdx
	leaq	(%rbx,%rax), %rcx
	movq	-192(%rbp), %rax
	movq	%rcx, -80(%rbp)
	movq	-96(%rbp), %rcx
	movq	%rax, -72(%rbp)
	call	_ZN5Eigen8internal13gemm_pack_lhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi1EEELi1ELi1ENS0_9Packet1cdELi1ELb0ELb0EEclEPS3_RKS5_llll.constprop.0
	testq	%r14, %r14
	movq	-88(%rbp), %r8
	jle	.L441
	cmpb	$0, -201(%rbp)
	jne	.L437
	movq	-200(%rbp), %r15
	movq	%r8, -144(%rbp)
	xorl	%ebx, %ebx
	movq	%r12, %rax
	movq	%rdx, -136(%rbp)
	movsd	-152(%rbp), %xmm1
	movq	%rbx, %r12
	movq	%r15, %rbx
	movq	%r14, %r15
	movq	%rax, %r14
	.p2align 4,,10
	.p2align 3
.L438:
	movq	-112(%rbp), %rax
	movq	%r12, %rdx
	movsd	-128(%rbp), %xmm0
	movq	%r14, %rdi
	movq	-136(%rbp), %r8
	movq	-96(%rbp), %rcx
	movq	%rbx, -80(%rbp)
	addq	%rax, %r12
	movq	-104(%rbp), %rsi
	movq	40(%rbp), %rax
	movsd	%xmm1, -88(%rbp)
	cmpq	%r12, %r15
	movq	%r12, %r9
	cmovle	%r15, %r9
	movq	%rax, -72(%rbp)
	subq	%rdx, %r9
	movq	%r13, %rdx
	call	_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb1ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0
	movq	-120(%rbp), %rax
	movsd	-88(%rbp), %xmm1
	addq	%rax, %rbx
	cmpq	%r12, %r15
	jg	.L438
	movq	-144(%rbp), %r8
	movq	%r14, %r12
	movq	%r15, %r14
.L441:
	movq	-168(%rbp), %rax
	cmpq	%rax, %r8
	jl	.L433
	movq	%r14, %rbx
	movq	-240(%rbp), %rcx
	movq	-232(%rbp), %r14
	movq	-256(%rbp), %rsi
	movq	-264(%rbp), %rdi
	addq	%rsi, -176(%rbp)
	addq	%rdi, -200(%rbp)
	cmpq	%rcx, %r14
	jl	.L436
	movq	-296(%rbp), %rdi
.L434:
	cmpq	$131072, -288(%rbp)
	ja	.L429
.L430:
	cmpq	$131072, -280(%rbp)
	ja	.L477
.L414:
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L478
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L437:
	.cfi_restore_state
	movq	-224(%rbp), %rax
	movq	%r8, -216(%rbp)
	xorl	%r10d, %r10d
	movq	%rdx, -136(%rbp)
	leaq	(%rax,%rbx), %r11
	movq	-200(%rbp), %rbx
	movq	%r12, -88(%rbp)
	movq	%r14, %r12
	movq	%r11, %r15
	movq	%r10, %r14
	.p2align 4,,10
	.p2align 3
.L440:
	movq	-112(%rbp), %rax
	movq	%r14, %rdx
	movq	-88(%rbp), %rsi
	movq	%r13, %rdi
	movq	%r15, -80(%rbp)
	addq	%rax, %r14
	movq	16(%rbp), %rax
	cmpq	%r14, %r12
	movq	%r14, %r9
	cmovle	%r12, %r9
	movq	%rax, -72(%rbp)
	subq	%rdx, %r9
	movq	-136(%rbp), %rdx
	movq	%r9, %rcx
	movq	%r9, -144(%rbp)
	call	_ZN5Eigen8internal13gemm_pack_rhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi4ELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0
	movq	40(%rbp), %rax
	movq	-96(%rbp), %rcx
	movq	%r13, %rdx
	movsd	-128(%rbp), %xmm0
	movq	-104(%rbp), %rsi
	movq	%rbx, -80(%rbp)
	movsd	-152(%rbp), %xmm1
	movq	-144(%rbp), %r9
	movq	%rax, -72(%rbp)
	movq	-136(%rbp), %r8
	movq	-88(%rbp), %rdi
	call	_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb1ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0
	movq	-160(%rbp), %rax
	addq	%rax, %r15
	movq	-120(%rbp), %rax
	addq	%rax, %rbx
	cmpq	%r14, %r12
	jg	.L440
	movq	%r12, %r14
	movq	-216(%rbp), %r8
	movq	-88(%rbp), %r12
	jmp	.L441
.L424:
	movq	%rdx, %rdi
.LEHB0:
	call	_ZN5Eigen8internal14aligned_mallocEm
.LEHE0:
	cmpq	$0, 8(%r12)
	movq	%rax, %rdi
	jne	.L427
	movq	%rax, %r13
	cmpq	%r14, %r15
	jl	.L425
	xorl	%eax, %eax
	testq	%r14, %r14
	jg	.L428
.L429:
	call	free@PLT
	cmpq	$131072, -280(%rbp)
	jbe	.L414
.L477:
	movq	-272(%rbp), %rdi
	call	free@PLT
	jmp	.L414
.L475:
	cmpq	$131072, %rax
	ja	.L418
	addq	$32, %rax
	subq	%rax, %rsp
	leaq	15(%rsp), %rax
	andq	$-16, %rax
	movq	%rax, -272(%rbp)
	cmpq	%rdx, %r13
	jb	.L467
.L473:
	movq	%rax, -104(%rbp)
	jmp	.L419
.L476:
	movq	%rdx, %rax
	cmpq	$131072, %rdx
	ja	.L424
	addq	$32, %rax
	subq	%rax, %rsp
	leaq	15(%rsp), %rax
	andq	$-16, %rax
	movq	%rax, %rdi
	movq	%rax, %r13
	cmpq	%r14, %r15
	jl	.L425
	xorl	%eax, %eax
	testq	%r14, %r14
	jg	.L428
	jmp	.L430
	.p2align 4,,10
	.p2align 3
.L418:
	movq	%rax, %rdi
	movq	%rdx, -144(%rbp)
.LEHB1:
	call	_ZN5Eigen8internal14aligned_mallocEm
.LEHE1:
	cmpq	$0, (%r12)
	movq	-144(%rbp), %rdx
	movq	%rax, -272(%rbp)
	jne	.L479
	cmpq	%rdx, %r13
	jnb	.L473
	jmp	.L469
.L427:
	cmpq	%r14, %r15
	jge	.L480
	xorl	%edi, %edi
	movq	%rax, %r13
	jmp	.L425
.L478:
	call	__stack_chk_fail@PLT
.L479:
	cmpq	%rdx, %r13
	jb	.L468
	xorl	%ecx, %ecx
	movq	-272(%rbp), %rax
	movq	%rcx, -272(%rbp)
	jmp	.L473
.L474:
	leaq	.LC30(%rip), %rcx
	movl	$174, %edx
	leaq	.LC31(%rip), %rsi
	leaq	.LC32(%rip), %rdi
	call	__assert_fail@PLT
.L480:
	testq	%r14, %r14
	jle	.L481
	movq	%rax, %r13
	xorl	%edi, %edi
	xorl	%eax, %eax
	jmp	.L428
.L481:
	xorl	%edi, %edi
	jmp	.L429
.L458:
	movq	%rax, %rbx
	jmp	.L444
	.section	.gcc_except_table,"a",@progbits
.LLSDA14888:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE14888-.LLSDACSB14888
.LLSDACSB14888:
	.uleb128 .LEHB0-.LFB14888
	.uleb128 .LEHE0-.LEHB0
	.uleb128 .L458-.LFB14888
	.uleb128 0
	.uleb128 .LEHB1-.LFB14888
	.uleb128 .LEHE1-.LEHB1
	.uleb128 0
	.uleb128 0
.LLSDACSE14888:
	.text
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDAC14888
	.type	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0.cold, @function
_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0.cold:
.LFSB14888:
.L468:
	.cfi_def_cfa 6, 16
	.cfi_offset 3, -56
	.cfi_offset 6, -16
	.cfi_offset 12, -48
	.cfi_offset 13, -40
	.cfi_offset 14, -32
	.cfi_offset 15, -24
.LEHB2:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE2:
.L470:
.LEHB3:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE3:
.L467:
.LEHB4:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE4:
.L469:
.LEHB5:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE5:
.L460:
	movq	%rax, %rbx
	xorl	%eax, %eax
	movq	%rax, -272(%rbp)
.L444:
	cmpq	$131072, -280(%rbp)
	jbe	.L445
	movq	-272(%rbp), %rax
	movq	%rax, -104(%rbp)
.L422:
	movq	-104(%rbp), %rdi
	call	free@PLT
.L445:
	movq	%rbx, %rdi
.LEHB6:
	call	_Unwind_Resume@PLT
.LEHE6:
.L461:
	movq	%rax, %rbx
	movq	-272(%rbp), %rax
	movq	%rax, -104(%rbp)
	jmp	.L422
.L459:
	movq	%rax, %rbx
	jmp	.L422
	.cfi_endproc
.LFE14888:
	.section	.gcc_except_table
.LLSDAC14888:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSEC14888-.LLSDACSBC14888
.LLSDACSBC14888:
	.uleb128 .LEHB2-.LCOLDB33
	.uleb128 .LEHE2-.LEHB2
	.uleb128 .L459-.LCOLDB33
	.uleb128 0
	.uleb128 .LEHB3-.LCOLDB33
	.uleb128 .LEHE3-.LEHB3
	.uleb128 .L460-.LCOLDB33
	.uleb128 0
	.uleb128 .LEHB4-.LCOLDB33
	.uleb128 .LEHE4-.LEHB4
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB5-.LCOLDB33
	.uleb128 .LEHE5-.LEHB5
	.uleb128 .L461-.LCOLDB33
	.uleb128 0
	.uleb128 .LEHB6-.LCOLDB33
	.uleb128 .LEHE6-.LEHB6
	.uleb128 0
	.uleb128 0
.LLSDACSEC14888:
	.section	.text.unlikely
	.text
	.size	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0, .-_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0
	.section	.text.unlikely
	.size	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0.cold, .-_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0.cold
.LCOLDE33:
	.text
.LHOTE33:
	.section	.text.unlikely
.LCOLDB34:
	.text
.LHOTB34:
	.p2align 4
	.type	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi0ELb0ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0, @function
_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi0ELb0ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0:
.LFB14889:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA14889
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$264, %rsp
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	24(%rbp), %rax
	movq	%rdx, -176(%rbp)
	movq	%rcx, -120(%rbp)
	movq	48(%rbp), %rbx
	movq	%r8, -192(%rbp)
	movq	%r9, -136(%rbp)
	movq	%rax, -144(%rbp)
	movsd	%xmm0, -128(%rbp)
	movsd	%xmm1, -152(%rbp)
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	xorl	%eax, %eax
	cmpq	$1, 32(%rbp)
	jne	.L542
	movq	%rsi, %r12
	movq	16(%rbx), %rsi
	movq	%rdi, %r13
	movq	24(%rbx), %rcx
	movabsq	$1152921504606846975, %r15
	movq	32(%rbx), %rax
	cmpq	%r13, %rsi
	movq	%rsi, -88(%rbp)
	cmovg	%r13, %rsi
	cmpq	%r12, %rcx
	movq	%rax, %rdi
	movq	%rcx, -96(%rbp)
	movq	%rax, -184(%rbp)
	cmovg	%r12, %rcx
	movq	%rdi, %r14
	imulq	%rsi, %rax
	movq	%rsi, -272(%rbp)
	movq	%rcx, -112(%rbp)
	imulq	%rcx, %r14
	cmpq	%rax, %r15
	jb	.L535
	movq	(%rbx), %rdi
	salq	$4, %rax
	movq	%rax, -288(%rbp)
	movq	%rdi, -104(%rbp)
	testq	%rdi, %rdi
	je	.L543
	cmpq	%r14, %r15
	jb	.L538
	movq	$0, -280(%rbp)
.L487:
	movq	8(%rbx), %r15
	salq	$4, %r14
	movq	%r14, -296(%rbp)
	testq	%r15, %r15
	je	.L544
	xorl	%r14d, %r14d
	xorl	%eax, %eax
	cmpq	%r13, -88(%rbp)
	jge	.L499
.L493:
	movq	-176(%rbp), %rdi
	cmpq	%rdi, -184(%rbp)
	sete	%al
	cmpq	%r12, -96(%rbp)
	setge	%dl
	andl	%edx, %eax
.L499:
	testq	%r13, %r13
	jle	.L502
.L496:
	cmpq	$0, -176(%rbp)
	jle	.L502
	movq	-144(%rbp), %rsi
	xorl	$1, %eax
	movq	-120(%rbp), %rdi
	movq	%r14, -304(%rbp)
	movq	-192(%rbp), %rcx
	movq	-272(%rbp), %r9
	movb	%al, -210(%rbp)
	leaq	-80(%rbp), %rbx
	movq	%rsi, -208(%rbp)
	movq	-184(%rbp), %rsi
	salq	$4, %r9
	movq	%rsi, %rdx
	salq	$4, %rsi
	imulq	%rcx, %rdx
	movq	-112(%rbp), %rcx
	salq	$4, %rdx
	movq	%rdx, -200(%rbp)
	movq	16(%rbp), %rdx
	imulq	%rcx, %rdx
	imulq	40(%rbp), %rcx
	salq	$4, %rdx
	movq	%rdx, -160(%rbp)
	movq	%rsi, %rdx
	movq	-136(%rbp), %rsi
	subq	%rdx, %rsi
	movq	%rcx, %rdx
	movq	%r13, %rcx
	salq	$4, %rdx
	movq	%rsi, -232(%rbp)
	xorl	%esi, %esi
	movq	%rdx, -120(%rbp)
	movq	%rsi, %r13
	.p2align 4,,10
	.p2align 3
.L504:
	movq	-272(%rbp), %rax
	movq	%r13, %rdx
	movq	%rcx, %r14
	movq	%rdi, -168(%rbp)
	movq	%rdi, -240(%rbp)
	addq	%rax, %r13
	movq	%r9, -256(%rbp)
	cmpq	%rcx, %r13
	movq	%rcx, -264(%rbp)
	cmovle	%r13, %r14
	movq	%r13, -248(%rbp)
	movq	%r12, %r13
	subq	%rdx, %r14
	testq	%rdx, %rdx
	sete	%dl
	orb	-210(%rbp), %dl
	movq	%r14, -96(%rbp)
	xorl	%r8d, %r8d
	movb	%dl, -209(%rbp)
	.p2align 4,,10
	.p2align 3
.L501:
	movq	-184(%rbp), %rax
	movq	%r8, %rcx
	movq	-104(%rbp), %rdi
	movq	%rbx, %rsi
	movq	-176(%rbp), %rdx
	addq	%rax, %r8
	movq	-168(%rbp), %rax
	cmpq	%rdx, %r8
	movq	%r8, -88(%rbp)
	cmovle	%r8, %rdx
	movq	%rax, -80(%rbp)
	movq	-192(%rbp), %rax
	subq	%rcx, %rdx
	movq	-96(%rbp), %rcx
	movq	%rax, -72(%rbp)
	call	_ZN5Eigen8internal13gemm_pack_lhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi1ELi1ENS0_9Packet1cdELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0
	testq	%r13, %r13
	movq	-88(%rbp), %r8
	jle	.L509
	cmpb	$0, -209(%rbp)
	jne	.L505
	movq	-208(%rbp), %r14
	movq	%r8, -144(%rbp)
	xorl	%r12d, %r12d
	movq	%rbx, %rax
	movq	%rdx, -136(%rbp)
	movsd	-152(%rbp), %xmm1
	movq	%r12, %rbx
	movq	%r14, %r12
	movq	%r13, %r14
	movq	%rax, %r13
	.p2align 4,,10
	.p2align 3
.L506:
	movq	-112(%rbp), %rax
	movq	%rbx, %rdx
	movsd	-128(%rbp), %xmm0
	movq	%r13, %rdi
	movq	-136(%rbp), %r8
	movq	-96(%rbp), %rcx
	movq	%r12, -80(%rbp)
	addq	%rax, %rbx
	movq	-104(%rbp), %rsi
	movq	40(%rbp), %rax
	movsd	%xmm1, -88(%rbp)
	cmpq	%rbx, %r14
	movq	%rbx, %r9
	cmovle	%r14, %r9
	movq	%rax, -72(%rbp)
	subq	%rdx, %r9
	movq	%r15, %rdx
	call	_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb0ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0
	movq	-120(%rbp), %rax
	movsd	-88(%rbp), %xmm1
	addq	%rax, %r12
	cmpq	%rbx, %r14
	jg	.L506
	movq	-144(%rbp), %r8
	movq	%r13, %rbx
	movq	%r14, %r13
.L509:
	movq	-176(%rbp), %rax
	movq	-200(%rbp), %rsi
	addq	%rsi, -168(%rbp)
	cmpq	%rax, %r8
	jl	.L501
	movq	-256(%rbp), %r9
	movq	-240(%rbp), %rdi
	movq	%r13, %r12
	movq	-264(%rbp), %rcx
	movq	-248(%rbp), %r13
	addq	%r9, -208(%rbp)
	addq	%r9, %rdi
	cmpq	%rcx, %r13
	jl	.L504
	movq	-304(%rbp), %r14
.L502:
	cmpq	$131072, -296(%rbp)
	ja	.L497
.L498:
	cmpq	$131072, -288(%rbp)
	ja	.L545
.L482:
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L546
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L505:
	.cfi_restore_state
	movq	-232(%rbp), %rax
	movq	%r8, %r11
	xorl	%r10d, %r10d
	movq	%r8, -224(%rbp)
	salq	$4, %r11
	movq	-208(%rbp), %r12
	movq	%rdx, -136(%rbp)
	addq	%rax, %r11
	movq	%rbx, -88(%rbp)
	movq	%r13, %rbx
	movq	%r10, %r13
	movq	%r11, %r14
	.p2align 4,,10
	.p2align 3
.L508:
	movq	-112(%rbp), %rax
	movq	%r13, %rdx
	movq	-88(%rbp), %rsi
	movq	%r15, %rdi
	movq	%r14, -80(%rbp)
	addq	%rax, %r13
	movq	16(%rbp), %rax
	cmpq	%r13, %rbx
	movq	%r13, %r9
	cmovle	%rbx, %r9
	movq	%rax, -72(%rbp)
	subq	%rdx, %r9
	movq	-136(%rbp), %rdx
	movq	%r9, %rcx
	movq	%r9, -144(%rbp)
	call	_ZN5Eigen8internal13gemm_pack_rhsISt7complexIdElNS0_22const_blas_data_mapperIS3_lLi0EEELi4ELi0ELb0ELb0EEclEPS3_RKS5_llll.constprop.0
	movq	40(%rbp), %rax
	movq	-96(%rbp), %rcx
	movq	%r15, %rdx
	movsd	-128(%rbp), %xmm0
	movq	-104(%rbp), %rsi
	movq	%r12, -80(%rbp)
	movsd	-152(%rbp), %xmm1
	movq	-144(%rbp), %r9
	movq	%rax, -72(%rbp)
	movq	-136(%rbp), %r8
	movq	-88(%rbp), %rdi
	call	_ZN5Eigen8internal11gebp_kernelISt7complexIdES3_lNS0_16blas_data_mapperIS3_lLi0ELi0ELi1EEELi1ELi4ELb0ELb0EEclERKS5_PKS3_SA_lllS3_llll.constprop.0
	movq	-160(%rbp), %rax
	addq	%rax, %r14
	movq	-120(%rbp), %rax
	addq	%rax, %r12
	cmpq	%r13, %rbx
	jg	.L508
	movq	%rbx, %r13
	movq	-224(%rbp), %r8
	movq	-88(%rbp), %rbx
	jmp	.L509
.L492:
	movq	%r14, %rdi
.LEHB7:
	call	_ZN5Eigen8internal14aligned_mallocEm
.LEHE7:
	cmpq	$0, 8(%rbx)
	movq	%rax, %r14
	jne	.L495
	movq	%rax, %r15
	cmpq	%r13, -88(%rbp)
	jl	.L493
	xorl	%eax, %eax
	testq	%r13, %r13
	jg	.L496
.L497:
	movq	%r14, %rdi
	call	free@PLT
	cmpq	$131072, -288(%rbp)
	jbe	.L482
.L545:
	movq	-280(%rbp), %rdi
	call	free@PLT
	jmp	.L482
.L543:
	cmpq	$131072, %rax
	ja	.L486
	addq	$32, %rax
	subq	%rax, %rsp
	leaq	15(%rsp), %rax
	andq	$-16, %rax
	movq	%rax, -280(%rbp)
	cmpq	%r14, %r15
	jb	.L535
.L541:
	movq	%rax, -104(%rbp)
	jmp	.L487
.L544:
	movq	%r14, %rax
	cmpq	$131072, %r14
	ja	.L492
	addq	$32, %rax
	subq	%rax, %rsp
	leaq	15(%rsp), %rax
	andq	$-16, %rax
	movq	%rax, %r14
	movq	%rax, %r15
	cmpq	%r13, -88(%rbp)
	jl	.L493
	xorl	%eax, %eax
	testq	%r13, %r13
	jg	.L496
	jmp	.L498
	.p2align 4,,10
	.p2align 3
.L486:
	movq	%rax, %rdi
.LEHB8:
	call	_ZN5Eigen8internal14aligned_mallocEm
.LEHE8:
	cmpq	$0, (%rbx)
	movq	%rax, -280(%rbp)
	jne	.L547
	cmpq	%r14, %r15
	jnb	.L541
	jmp	.L537
.L495:
	cmpq	%r13, -88(%rbp)
	jge	.L548
	xorl	%r14d, %r14d
	movq	%rax, %r15
	jmp	.L493
.L546:
	call	__stack_chk_fail@PLT
.L547:
	cmpq	%r14, %r15
	jb	.L536
	xorl	%edx, %edx
	movq	-280(%rbp), %rax
	movq	%rdx, -280(%rbp)
	jmp	.L541
.L542:
	leaq	.LC30(%rip), %rcx
	movl	$174, %edx
	leaq	.LC31(%rip), %rsi
	leaq	.LC32(%rip), %rdi
	call	__assert_fail@PLT
.L548:
	testq	%r13, %r13
	jle	.L549
	movq	%rax, %r15
	xorl	%r14d, %r14d
	xorl	%eax, %eax
	jmp	.L496
.L549:
	xorl	%r14d, %r14d
	jmp	.L497
.L526:
	movq	%rax, %rbx
	jmp	.L512
	.section	.gcc_except_table
.LLSDA14889:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE14889-.LLSDACSB14889
.LLSDACSB14889:
	.uleb128 .LEHB7-.LFB14889
	.uleb128 .LEHE7-.LEHB7
	.uleb128 .L526-.LFB14889
	.uleb128 0
	.uleb128 .LEHB8-.LFB14889
	.uleb128 .LEHE8-.LEHB8
	.uleb128 0
	.uleb128 0
.LLSDACSE14889:
	.text
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDAC14889
	.type	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi0ELb0ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0.cold, @function
_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi0ELb0ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0.cold:
.LFSB14889:
.L536:
	.cfi_def_cfa 6, 16
	.cfi_offset 3, -56
	.cfi_offset 6, -16
	.cfi_offset 12, -48
	.cfi_offset 13, -40
	.cfi_offset 14, -32
	.cfi_offset 15, -24
.LEHB9:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE9:
.L538:
.LEHB10:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE10:
.L535:
.LEHB11:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE11:
.L537:
.LEHB12:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE12:
.L528:
	movq	%rax, %rbx
	xorl	%eax, %eax
	movq	%rax, -280(%rbp)
.L512:
	cmpq	$131072, -288(%rbp)
	jbe	.L513
	movq	-280(%rbp), %rax
	movq	%rax, -104(%rbp)
.L490:
	movq	-104(%rbp), %rdi
	call	free@PLT
.L513:
	movq	%rbx, %rdi
.LEHB13:
	call	_Unwind_Resume@PLT
.LEHE13:
.L529:
	movq	%rax, %rbx
	movq	-280(%rbp), %rax
	movq	%rax, -104(%rbp)
	jmp	.L490
.L527:
	movq	%rax, %rbx
	jmp	.L490
	.cfi_endproc
.LFE14889:
	.section	.gcc_except_table
.LLSDAC14889:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSEC14889-.LLSDACSBC14889
.LLSDACSBC14889:
	.uleb128 .LEHB9-.LCOLDB34
	.uleb128 .LEHE9-.LEHB9
	.uleb128 .L527-.LCOLDB34
	.uleb128 0
	.uleb128 .LEHB10-.LCOLDB34
	.uleb128 .LEHE10-.LEHB10
	.uleb128 .L528-.LCOLDB34
	.uleb128 0
	.uleb128 .LEHB11-.LCOLDB34
	.uleb128 .LEHE11-.LEHB11
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB12-.LCOLDB34
	.uleb128 .LEHE12-.LEHB12
	.uleb128 .L529-.LCOLDB34
	.uleb128 0
	.uleb128 .LEHB13-.LCOLDB34
	.uleb128 .LEHE13-.LEHB13
	.uleb128 0
	.uleb128 0
.LLSDACSEC14889:
	.section	.text.unlikely
	.text
	.size	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi0ELb0ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0, .-_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi0ELb0ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0
	.section	.text.unlikely
	.size	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi0ELb0ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0.cold, .-_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi0ELb0ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0.cold
.LCOLDE34:
	.text
.LHOTE34:
	.section	.text._ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i,"axG",@progbits,_ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i
	.type	_ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i, @function
_ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i:
.LFB4758:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rdx, %r8
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$40, %rsp
	.cfi_def_cfa_offset 64
	movq	%fs:40, %rax
	movq	%rax, 24(%rsp)
	xorl	%eax, %eax
	movl	$0, (%rdx)
	movl	$0, (%rsi)
	movl	$0, (%rdi)
	cmpl	$3, %ecx
	jg	.L608
	cmpl	$1, %ecx
	jg	.L609
.L550:
	movq	24(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L610
	addq	$40, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L609:
	.cfi_restore_state
	movl	$2, %eax
	xorl	%ecx, %ecx
#APP
# 996 "/usr/local/include/Eigen/src/Core/util/Memory.h" 1
	xchgq	%rbx, %r9; cpuid; xchgq	%rbx, %r9
# 0 "" 2
#NO_APP
	movl	%eax, (%rsp)
	xorl	%r10d, %r10d
	movl	%r9d, 4(%rsp)
	leaq	16(%rsp), %r9
	movl	%ecx, 8(%rsp)
	leaq	.L561(%rip), %rcx
	movl	%edx, 12(%rsp)
	leaq	2(%rsp), %rdx
	.p2align 4,,10
	.p2align 3
.L596:
	movzbl	(%rdx), %eax
	subl	$10, %eax
	cmpb	$-125, %al
	ja	.L594
	movzbl	%al, %eax
	movslq	(%rcx,%rax,4), %rax
	addq	%rcx, %rax
	jmp	*%rax
	.section	.rodata._ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i,"aG",@progbits,_ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i,comdat
	.align 4
	.align 4
.L561:
	.long	.L572-.L561
	.long	.L594-.L561
	.long	.L571-.L561
	.long	.L594-.L561
	.long	.L593-.L561
	.long	.L594-.L561
	.long	.L571-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L571-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L592-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L591-.L561
	.long	.L590-.L561
	.long	.L594-.L561
	.long	.L564-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L563-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L570-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L570-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L569-.L561
	.long	.L587-.L561
	.long	.L569-.L561
	.long	.L568-.L561
	.long	.L586-.L561
	.long	.L566-.L561
	.long	.L594-.L561
	.long	.L585-.L561
	.long	.L569-.L561
	.long	.L568-.L561
	.long	.L566-.L561
	.long	.L565-.L561
	.long	.L567-.L561
	.long	.L563-.L561
	.long	.L562-.L561
	.long	.L579-.L561
	.long	.L578-.L561
	.long	.L577-.L561
	.long	.L562-.L561
	.long	.L575-.L561
	.long	.L574-.L561
	.long	.L573-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L571-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L572-.L561
	.long	.L571-.L561
	.long	.L570-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L565-.L561
	.long	.L569-.L561
	.long	.L568-.L561
	.long	.L566-.L561
	.long	.L565-.L561
	.long	.L567-.L561
	.long	.L568-.L561
	.long	.L566-.L561
	.long	.L566-.L561
	.long	.L569-.L561
	.long	.L568-.L561
	.long	.L566-.L561
	.long	.L565-.L561
	.long	.L567-.L561
	.long	.L566-.L561
	.long	.L565-.L561
	.long	.L564-.L561
	.long	.L563-.L561
	.long	.L562-.L561
	.long	.L594-.L561
	.long	.L594-.L561
	.long	.L560-.L561
	.section	.text._ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i,"axG",@progbits,_ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i,comdat
	.p2align 4,,10
	.p2align 3
.L566:
	movl	$512, (%rsi)
	.p2align 4,,10
	.p2align 3
.L594:
	addq	$1, %rdx
	cmpq	%r9, %rdx
	jne	.L596
.L612:
	testb	%r10b, %r10b
	je	.L597
	movl	(%r8), %eax
	cmpl	%eax, (%rsi)
	je	.L611
.L597:
	sall	$10, (%rdi)
	sall	$10, (%rsi)
	sall	$10, (%r8)
	jmp	.L550
	.p2align 4,,10
	.p2align 3
.L569:
	addq	$1, %rdx
	movl	$128, (%rsi)
	cmpq	%r9, %rdx
	jne	.L596
	jmp	.L612
	.p2align 4,,10
	.p2align 3
.L571:
	addq	$1, %rdx
	movl	$16, (%rdi)
	cmpq	%r9, %rdx
	jne	.L596
	jmp	.L612
	.p2align 4,,10
	.p2align 3
.L565:
	addq	$1, %rdx
	movl	$1024, (%rsi)
	cmpq	%r9, %rdx
	jne	.L596
	jmp	.L612
	.p2align 4,,10
	.p2align 3
.L568:
	addq	$1, %rdx
	movl	$256, (%rsi)
	cmpq	%r9, %rdx
	jne	.L596
	jmp	.L612
.L578:
	movl	(%rsi), %eax
	testl	%eax, %eax
	je	.L595
.L563:
	addq	$1, %rdx
	movl	$4096, (%r8)
	cmpq	%r9, %rdx
	jne	.L596
	jmp	.L612
	.p2align 4,,10
	.p2align 3
.L570:
	addq	$1, %rdx
	movl	$32, (%rdi)
	cmpq	%r9, %rdx
	jne	.L596
	jmp	.L612
	.p2align 4,,10
	.p2align 3
.L562:
	addq	$1, %rdx
	movl	$8192, (%r8)
	cmpq	%r9, %rdx
	jne	.L596
	jmp	.L612
	.p2align 4,,10
	.p2align 3
.L567:
	addq	$1, %rdx
	movl	$2048, (%rsi)
	cmpq	%r9, %rdx
	jne	.L596
	jmp	.L612
	.p2align 4,,10
	.p2align 3
.L608:
	xorl	%r9d, %r9d
	movl	$4, %r10d
	jmp	.L555
	.p2align 4,,10
	.p2align 3
.L614:
	cmpl	$1, %eax
	je	.L613
.L552:
	addl	$1, %r9d
	testl	%ebx, %ebx
	je	.L550
	cmpl	$15, %r9d
	jg	.L550
.L555:
	movl	%r10d, %eax
	movl	%r9d, %ecx
#APP
# 967 "/usr/local/include/Eigen/src/Core/util/Memory.h" 1
	xchgq	%rbx, %r11; cpuid; xchgq	%rbx, %r11
# 0 "" 2
#NO_APP
	movl	%eax, %edx
	movl	%eax, %ebx
	andl	$13, %edx
	andl	$15, %ebx
	cmpl	$1, %edx
	jne	.L552
	movl	%r11d, %edx
	movl	%r11d, %ebp
	andl	$4095, %r11d
	sarl	$5, %eax
	sarl	$12, %edx
	shrl	$22, %ebp
	addl	$1, %r11d
	addl	$1, %ecx
	andl	$1023, %edx
	addl	$1, %ebp
	andl	$7, %eax
	addl	$1, %edx
	imull	%ebp, %edx
	imull	%r11d, %edx
	imull	%ecx, %edx
	cmpl	$2, %eax
	je	.L553
	cmpl	$3, %eax
	jne	.L614
	movl	%edx, (%r8)
	jmp	.L552
	.p2align 4,,10
	.p2align 3
.L553:
	movl	%edx, (%rsi)
	jmp	.L552
	.p2align 4,,10
	.p2align 3
.L613:
	movl	%edx, (%rdi)
	jmp	.L552
	.p2align 4,,10
	.p2align 3
.L564:
	addq	$1, %rdx
	movl	$2048, (%r8)
	cmpq	%r9, %rdx
	jne	.L596
	jmp	.L612
	.p2align 4,,10
	.p2align 3
.L572:
	addq	$1, %rdx
	movl	$8, (%rdi)
	cmpq	%r9, %rdx
	jne	.L596
	jmp	.L612
.L592:
	movl	$96, (%rsi)
	jmp	.L594
.L593:
	movl	$24, (%rdi)
	jmp	.L594
.L587:
	movl	$192, (%rsi)
	jmp	.L594
.L560:
	movl	$3072, (%r8)
	jmp	.L594
.L590:
	movl	$1024, (%r8)
	jmp	.L594
.L591:
	movl	$512, (%r8)
	jmp	.L594
.L573:
	movl	$6144, (%rsi)
	jmp	.L594
.L586:
	movl	$384, (%rsi)
	jmp	.L594
.L579:
	movl	$3072, (%rsi)
	jmp	.L594
.L585:
	movl	$0, (%rsi)
	jmp	.L594
.L574:
	movl	$16384, (%r8)
	jmp	.L594
.L575:
	movl	$12288, (%r8)
	jmp	.L594
.L577:
	movl	$6144, (%r8)
	jmp	.L594
	.p2align 4,,10
	.p2align 3
.L611:
	movl	$0, (%r8)
	jmp	.L597
.L595:
	movl	$4096, (%rsi)
	movl	$1, %r10d
	movl	$4096, (%r8)
	jmp	.L594
.L610:
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE4758:
	.size	_ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i, .-_ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i
	.section	.text._ZN5Eigen8internal15queryCacheSizesERiS1_S1_,"axG",@progbits,_ZN5Eigen8internal15queryCacheSizesERiS1_S1_,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal15queryCacheSizesERiS1_S1_
	.type	_ZN5Eigen8internal15queryCacheSizesERiS1_S1_, @function
_ZN5Eigen8internal15queryCacheSizesERiS1_S1_:
.LFB4760:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rdx, %r8
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$40, %rsp
	.cfi_def_cfa_offset 64
	movq	%fs:40, %rax
	movq	%rax, 24(%rsp)
	xorl	%eax, %eax
	movl	%eax, %ecx
#APP
# 1115 "/usr/local/include/Eigen/src/Core/util/Memory.h" 1
	xchgq	%rbx, %r9; cpuid; xchgq	%rbx, %r9
# 0 "" 2
#NO_APP
	cmpl	$1970169159, %r9d
	je	.L687
	cmpl	$1752462657, %r9d
	je	.L688
	cmpl	$1952801395, %edx
	sete	%r10b
	cmpl	$1766083905, %r9d
	sete	%dl
	testb	%dl, %r10b
	je	.L617
	cmpl	$561145204, %ecx
	je	.L620
	.p2align 4,,10
	.p2align 3
.L617:
	cmpl	$3, %eax
	jg	.L689
	cmpl	$1, %eax
	jg	.L690
.L630:
	movl	$0, (%r8)
	movl	$0, (%rsi)
	movl	$0, (%rdi)
.L615:
	movq	24(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L686
	addq	$40, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L688:
	.cfi_restore_state
	cmpl	$1145913699, %ecx
	jne	.L617
	cmpl	$1769238117, %edx
	jne	.L617
.L620:
	xorl	%r9d, %r9d
	movl	$-2147483648, %eax
	movl	%r9d, %ecx
#APP
# 1087 "/usr/local/include/Eigen/src/Core/util/Memory.h" 1
	xchgq	%rbx, %r10; cpuid; xchgq	%rbx, %r10
# 0 "" 2
#NO_APP
	cmpl	$-2147483643, %eax
	jbe	.L630
	movl	$-2147483643, %eax
	movl	%r9d, %ecx
#APP
# 1090 "/usr/local/include/Eigen/src/Core/util/Memory.h" 1
	xchgq	%rbx, %r10; cpuid; xchgq	%rbx, %r10
# 0 "" 2
#NO_APP
	sarl	$24, %ecx
	movl	$-2147483642, %eax
	sall	$10, %ecx
	movl	%ecx, (%rdi)
	movl	%r9d, %ecx
#APP
# 1093 "/usr/local/include/Eigen/src/Core/util/Memory.h" 1
	xchgq	%rbx, %rdi; cpuid; xchgq	%rbx, %rdi
# 0 "" 2
#NO_APP
	sarl	$16, %ecx
	leal	(%rdx,%rdx), %eax
	sall	$10, %ecx
	andl	$536346624, %eax
	movl	%ecx, (%rsi)
	movl	%eax, (%r8)
	jmp	.L615
	.p2align 4,,10
	.p2align 3
.L687:
	cmpl	$1231384169, %edx
	jne	.L617
	cmpl	$1818588270, %ecx
	jne	.L617
	movq	24(%rsp), %rdx
	subq	%fs:40, %rdx
	jne	.L686
	addq	$40, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	movl	%eax, %ecx
	movq	%r8, %rdx
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	jmp	_ZN5Eigen8internal21queryCacheSizes_intelERiS1_S1_i
	.p2align 4,,10
	.p2align 3
.L689:
	.cfi_restore_state
	movl	$0, (%r8)
	xorl	%r9d, %r9d
	movl	$4, %r10d
	movl	$0, (%rsi)
	movl	$0, (%rdi)
	jmp	.L627
	.p2align 4,,10
	.p2align 3
.L692:
	cmpl	$1, %eax
	je	.L691
.L624:
	addl	$1, %r9d
	testl	%ebx, %ebx
	je	.L615
	cmpl	$15, %r9d
	jg	.L615
.L627:
	movl	%r10d, %eax
	movl	%r9d, %ecx
#APP
# 967 "/usr/local/include/Eigen/src/Core/util/Memory.h" 1
	xchgq	%rbx, %r11; cpuid; xchgq	%rbx, %r11
# 0 "" 2
#NO_APP
	movl	%eax, %edx
	movl	%eax, %ebx
	andl	$13, %edx
	andl	$15, %ebx
	cmpl	$1, %edx
	jne	.L624
	movl	%r11d, %edx
	movl	%r11d, %ebp
	andl	$4095, %r11d
	sarl	$5, %eax
	sarl	$12, %edx
	shrl	$22, %ebp
	addl	$1, %r11d
	addl	$1, %ecx
	andl	$1023, %edx
	addl	$1, %ebp
	andl	$7, %eax
	addl	$1, %edx
	imull	%ebp, %edx
	imull	%r11d, %edx
	imull	%ecx, %edx
	cmpl	$2, %eax
	je	.L625
	cmpl	$3, %eax
	jne	.L692
	movl	%edx, (%r8)
	jmp	.L624
	.p2align 4,,10
	.p2align 3
.L625:
	movl	%edx, (%rsi)
	jmp	.L624
	.p2align 4,,10
	.p2align 3
.L691:
	movl	%edx, (%rdi)
	jmp	.L624
	.p2align 4,,10
	.p2align 3
.L690:
	movl	$0, (%r8)
	movl	$2, %eax
	xorl	%ecx, %ecx
	movl	$0, (%rsi)
	movl	$0, (%rdi)
#APP
# 996 "/usr/local/include/Eigen/src/Core/util/Memory.h" 1
	xchgq	%rbx, %r9; cpuid; xchgq	%rbx, %r9
# 0 "" 2
#NO_APP
	movl	%eax, (%rsp)
	xorl	%r10d, %r10d
	movl	%r9d, 4(%rsp)
	leaq	16(%rsp), %r9
	movl	%ecx, 8(%rsp)
	leaq	.L633(%rip), %rcx
	movl	%edx, 12(%rsp)
	leaq	2(%rsp), %rdx
	.p2align 4,,10
	.p2align 3
.L668:
	movzbl	(%rdx), %eax
	subl	$10, %eax
	cmpb	$-125, %al
	ja	.L666
	movzbl	%al, %eax
	movslq	(%rcx,%rax,4), %rax
	addq	%rcx, %rax
	jmp	*%rax
	.section	.rodata._ZN5Eigen8internal15queryCacheSizesERiS1_S1_,"aG",@progbits,_ZN5Eigen8internal15queryCacheSizesERiS1_S1_,comdat
	.align 4
	.align 4
.L633:
	.long	.L644-.L633
	.long	.L666-.L633
	.long	.L643-.L633
	.long	.L666-.L633
	.long	.L665-.L633
	.long	.L666-.L633
	.long	.L643-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L643-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L664-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L663-.L633
	.long	.L662-.L633
	.long	.L666-.L633
	.long	.L636-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L635-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L642-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L642-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L641-.L633
	.long	.L659-.L633
	.long	.L641-.L633
	.long	.L640-.L633
	.long	.L658-.L633
	.long	.L638-.L633
	.long	.L666-.L633
	.long	.L657-.L633
	.long	.L641-.L633
	.long	.L640-.L633
	.long	.L638-.L633
	.long	.L637-.L633
	.long	.L639-.L633
	.long	.L635-.L633
	.long	.L634-.L633
	.long	.L651-.L633
	.long	.L650-.L633
	.long	.L649-.L633
	.long	.L634-.L633
	.long	.L647-.L633
	.long	.L646-.L633
	.long	.L645-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L643-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L644-.L633
	.long	.L643-.L633
	.long	.L642-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L637-.L633
	.long	.L641-.L633
	.long	.L640-.L633
	.long	.L638-.L633
	.long	.L637-.L633
	.long	.L639-.L633
	.long	.L640-.L633
	.long	.L638-.L633
	.long	.L638-.L633
	.long	.L641-.L633
	.long	.L640-.L633
	.long	.L638-.L633
	.long	.L637-.L633
	.long	.L639-.L633
	.long	.L638-.L633
	.long	.L637-.L633
	.long	.L636-.L633
	.long	.L635-.L633
	.long	.L634-.L633
	.long	.L666-.L633
	.long	.L666-.L633
	.long	.L632-.L633
	.section	.text._ZN5Eigen8internal15queryCacheSizesERiS1_S1_,"axG",@progbits,_ZN5Eigen8internal15queryCacheSizesERiS1_S1_,comdat
	.p2align 4,,10
	.p2align 3
.L638:
	movl	$512, (%rsi)
	.p2align 4,,10
	.p2align 3
.L666:
	addq	$1, %rdx
	cmpq	%r9, %rdx
	jne	.L668
.L694:
	testb	%r10b, %r10b
	je	.L669
	movl	(%r8), %eax
	cmpl	%eax, (%rsi)
	je	.L693
.L669:
	sall	$10, (%rdi)
	sall	$10, (%rsi)
	sall	$10, (%r8)
	jmp	.L615
	.p2align 4,,10
	.p2align 3
.L641:
	addq	$1, %rdx
	movl	$128, (%rsi)
	cmpq	%r9, %rdx
	jne	.L668
	jmp	.L694
	.p2align 4,,10
	.p2align 3
.L640:
	addq	$1, %rdx
	movl	$256, (%rsi)
	cmpq	%r9, %rdx
	jne	.L668
	jmp	.L694
	.p2align 4,,10
	.p2align 3
.L643:
	addq	$1, %rdx
	movl	$16, (%rdi)
	cmpq	%r9, %rdx
	jne	.L668
	jmp	.L694
	.p2align 4,,10
	.p2align 3
.L637:
	addq	$1, %rdx
	movl	$1024, (%rsi)
	cmpq	%r9, %rdx
	jne	.L668
	jmp	.L694
.L650:
	movl	(%rsi), %eax
	testl	%eax, %eax
	je	.L667
.L635:
	addq	$1, %rdx
	movl	$4096, (%r8)
	cmpq	%r9, %rdx
	jne	.L668
	jmp	.L694
	.p2align 4,,10
	.p2align 3
.L639:
	addq	$1, %rdx
	movl	$2048, (%rsi)
	cmpq	%r9, %rdx
	jne	.L668
	jmp	.L694
	.p2align 4,,10
	.p2align 3
.L642:
	addq	$1, %rdx
	movl	$32, (%rdi)
	cmpq	%r9, %rdx
	jne	.L668
	jmp	.L694
	.p2align 4,,10
	.p2align 3
.L634:
	addq	$1, %rdx
	movl	$8192, (%r8)
	cmpq	%r9, %rdx
	jne	.L668
	jmp	.L694
	.p2align 4,,10
	.p2align 3
.L636:
	addq	$1, %rdx
	movl	$2048, (%r8)
	cmpq	%r9, %rdx
	jne	.L668
	jmp	.L694
	.p2align 4,,10
	.p2align 3
.L644:
	addq	$1, %rdx
	movl	$8, (%rdi)
	cmpq	%r9, %rdx
	jne	.L668
	jmp	.L694
.L665:
	movl	$24, (%rdi)
	jmp	.L666
.L645:
	movl	$6144, (%rsi)
	jmp	.L666
.L659:
	movl	$192, (%rsi)
	jmp	.L666
.L632:
	movl	$3072, (%r8)
	jmp	.L666
.L662:
	movl	$1024, (%r8)
	jmp	.L666
.L663:
	movl	$512, (%r8)
	jmp	.L666
.L664:
	movl	$96, (%rsi)
	jmp	.L666
.L646:
	movl	$16384, (%r8)
	jmp	.L666
.L647:
	movl	$12288, (%r8)
	jmp	.L666
.L649:
	movl	$6144, (%r8)
	jmp	.L666
.L658:
	movl	$384, (%rsi)
	jmp	.L666
.L651:
	movl	$3072, (%rsi)
	jmp	.L666
.L657:
	movl	$0, (%rsi)
	jmp	.L666
	.p2align 4,,10
	.p2align 3
.L693:
	movl	$0, (%r8)
	jmp	.L669
.L667:
	movl	$4096, (%rsi)
	movl	$1, %r10d
	movl	$4096, (%r8)
	jmp	.L666
.L686:
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE4760:
	.size	_ZN5Eigen8internal15queryCacheSizesERiS1_S1_, .-_ZN5Eigen8internal15queryCacheSizesERiS1_S1_
	.section	.text._ZN5Eigen8IOFormatC2EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c,"axG",@progbits,_ZN5Eigen8IOFormatC5EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen8IOFormatC2EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c
	.type	_ZN5Eigen8IOFormatC2EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c, @function
_ZN5Eigen8IOFormatC2EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c:
.LFB6075:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA6075
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	movq	%r9, %r14
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	movq	%r8, %r13
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	movq	%rcx, %r12
	leaq	16(%rdi), %rcx
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movl	%edx, %ebp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	%rdi, %rbx
	subq	$72, %rsp
	.cfi_def_cfa_offset 128
	movq	136(%rsp), %rax
	movl	%esi, 8(%rsp)
	movl	152(%rsp), %r15d
	movq	%rcx, (%rdi)
	movq	8(%rax), %rdx
	movq	(%rax), %rsi
	movq	%rcx, 16(%rsp)
	addq	%rsi, %rdx
.LEHB14:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
.LEHE14:
	leaq	48(%rbx), %rax
	leaq	32(%rbx), %rdi
	movq	%rax, 32(%rbx)
	movq	%rax, 24(%rsp)
	movq	144(%rsp), %rax
	movq	(%rax), %rsi
	movq	8(%rax), %rdx
	addq	%rsi, %rdx
.LEHB15:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
.LEHE15:
	leaq	80(%rbx), %rax
	movq	8(%r14), %rdx
	leaq	64(%rbx), %rdi
	movq	%rax, 64(%rbx)
	movq	(%r14), %rsi
	movq	%rax, 32(%rsp)
	addq	%rsi, %rdx
.LEHB16:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
.LEHE16:
	leaq	112(%rbx), %rax
	leaq	96(%rbx), %rdi
	movq	%rax, 96(%rbx)
	movq	%rax, 40(%rsp)
	movq	128(%rsp), %rax
	movq	(%rax), %rsi
	movq	8(%rax), %rdx
	addq	%rsi, %rdx
.LEHB17:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
.LEHE17:
	leaq	144(%rbx), %rax
	movq	8(%r13), %rdx
	leaq	128(%rbx), %rdi
	movq	%rax, 128(%rbx)
	movq	0(%r13), %rsi
	movq	%rax, 48(%rsp)
	addq	%rsi, %rdx
.LEHB18:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
.LEHE18:
	leaq	208(%rbx), %rax
	leaq	176(%rbx), %r13
	movb	$0, 176(%rbx)
	movq	%r13, 160(%rbx)
	leaq	192(%rbx), %rdi
	movq	$0, 168(%rbx)
	movq	8(%r12), %rdx
	movq	%rax, 192(%rbx)
	movq	(%r12), %rsi
	movq	%rax, 56(%rsp)
	addq	%rsi, %rdx
.LEHB19:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
.LEHE19:
	movl	8(%rsp), %eax
	movl	%ebp, 232(%rbx)
	andl	$1, %ebp
	movb	%r15b, 224(%rbx)
	movl	%eax, 228(%rbx)
	jne	.L695
	movq	40(%rbx), %r15
	movl	%r15d, %eax
	subl	$1, %eax
	js	.L695
	movslq	%r15d, %r15
	leaq	160(%rbx), %rcx
	movslq	%eax, %r12
	movl	%eax, %eax
	subq	$2, %r15
	movq	%rcx, 8(%rsp)
	subq	%rax, %r15
	jmp	.L700
	.p2align 4,,10
	.p2align 3
.L699:
	movb	$32, (%rax,%rbp)
	movq	160(%rbx), %rax
	subq	$1, %r12
	movq	%r14, 168(%rbx)
	movb	$0, 1(%rax,%rbp)
	cmpq	%r12, %r15
	je	.L695
.L700:
	movq	32(%rbx), %rax
	cmpb	$10, (%rax,%r12)
	je	.L695
	movq	168(%rbx), %rbp
	movq	160(%rbx), %rax
	leaq	1(%rbp), %r14
	cmpq	%rax, %r13
	je	.L714
	movq	176(%rbx), %rdx
.L698:
	cmpq	%r14, %rdx
	jnb	.L699
	movq	8(%rsp), %rdi
	movl	$1, %r8d
	xorl	%ecx, %ecx
	xorl	%edx, %edx
	movq	%rbp, %rsi
.LEHB20:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_mutateEmmPKcm@PLT
.LEHE20:
	movq	160(%rbx), %rax
	jmp	.L699
	.p2align 4,,10
	.p2align 3
.L714:
	movl	$15, %edx
	jmp	.L698
	.p2align 4,,10
	.p2align 3
.L695:
	addq	$72, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
.L715:
	.cfi_restore_state
	movq	%rax, %rbp
	jmp	.L712
.L720:
	movq	%rax, %rbp
.L701:
	movq	192(%rbx), %rdi
	cmpq	%rdi, 56(%rsp)
	je	.L703
	movq	208(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L703:
	movq	160(%rbx), %rdi
	cmpq	%r13, %rdi
	jne	.L723
.L704:
	movq	128(%rbx), %rdi
	cmpq	%rdi, 48(%rsp)
	je	.L706
	movq	144(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L706:
	movq	96(%rbx), %rdi
	cmpq	%rdi, 40(%rsp)
	je	.L708
	movq	112(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L708:
	movq	64(%rbx), %rdi
	cmpq	%rdi, 32(%rsp)
	je	.L710
	movq	80(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L710:
	movq	32(%rbx), %rdi
	cmpq	%rdi, 24(%rsp)
	je	.L712
	movq	48(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L712:
	movq	(%rbx), %rdi
	cmpq	%rdi, 16(%rsp)
	je	.L713
	movq	16(%rbx), %rsi
	addq	$1, %rsi
	call	_ZdlPvm@PLT
.L713:
	movq	%rbp, %rdi
.LEHB21:
	call	_Unwind_Resume@PLT
.LEHE21:
.L723:
	movq	176(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L704
.L719:
	movq	%rax, %rbp
	jmp	.L703
.L718:
	movq	%rax, %rbp
	jmp	.L706
.L717:
	movq	%rax, %rbp
	jmp	.L708
.L716:
	movq	%rax, %rbp
	jmp	.L710
	.cfi_endproc
.LFE6075:
	.section	.gcc_except_table
.LLSDA6075:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE6075-.LLSDACSB6075
.LLSDACSB6075:
	.uleb128 .LEHB14-.LFB6075
	.uleb128 .LEHE14-.LEHB14
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB15-.LFB6075
	.uleb128 .LEHE15-.LEHB15
	.uleb128 .L715-.LFB6075
	.uleb128 0
	.uleb128 .LEHB16-.LFB6075
	.uleb128 .LEHE16-.LEHB16
	.uleb128 .L716-.LFB6075
	.uleb128 0
	.uleb128 .LEHB17-.LFB6075
	.uleb128 .LEHE17-.LEHB17
	.uleb128 .L717-.LFB6075
	.uleb128 0
	.uleb128 .LEHB18-.LFB6075
	.uleb128 .LEHE18-.LEHB18
	.uleb128 .L718-.LFB6075
	.uleb128 0
	.uleb128 .LEHB19-.LFB6075
	.uleb128 .LEHE19-.LEHB19
	.uleb128 .L719-.LFB6075
	.uleb128 0
	.uleb128 .LEHB20-.LFB6075
	.uleb128 .LEHE20-.LEHB20
	.uleb128 .L720-.LFB6075
	.uleb128 0
	.uleb128 .LEHB21-.LFB6075
	.uleb128 .LEHE21-.LEHB21
	.uleb128 0
	.uleb128 0
.LLSDACSE6075:
	.section	.text._ZN5Eigen8IOFormatC2EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c,"axG",@progbits,_ZN5Eigen8IOFormatC5EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c,comdat
	.size	_ZN5Eigen8IOFormatC2EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c, .-_ZN5Eigen8IOFormatC2EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c
	.weak	_ZN5Eigen8IOFormatC1EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c
	.set	_ZN5Eigen8IOFormatC1EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c,_ZN5Eigen8IOFormatC2EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c
	.section	.rodata._ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_.str1.1,"aMS",@progbits,1
.LC35:
	.string	"basic_string::append"
	.section	.text._ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_,"axG",@progbits,_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_,comdat
	.p2align 4
	.weak	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_
	.type	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_, @function
_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_:
.LFB10269:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	movq	%rdx, %rax
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movq	%rdi, %rbx
	movq	8(%rsi), %r8
	movq	%rsi, %rdi
	movq	8(%rdx), %rdx
	movq	(%rsi), %r9
	leaq	16(%rdi), %r10
	movq	(%rax), %rsi
	leaq	(%rdx,%r8), %rcx
	cmpq	%r10, %r9
	je	.L745
	cmpq	%rcx, 16(%rdi)
	jnb	.L726
	leaq	16(%rax), %r10
	cmpq	%r10, %rsi
	je	.L735
.L727:
	movq	16(%rax), %r10
.L728:
	cmpq	%rcx, %r10
	jnb	.L746
.L726:
	movabsq	$4611686018427387903, %rax
	subq	%r8, %rax
	cmpq	%rdx, %rax
	jb	.L747
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_appendEPKcm@PLT
.L744:
	leaq	16(%rbx), %rdi
	movq	%rax, %rbp
	movq	%rdi, (%rbx)
	movq	(%rax), %rax
	leaq	16(%rbp), %r12
	cmpq	%r12, %rax
	je	.L748
	movq	%rax, (%rbx)
	movq	16(%rbp), %rax
	movq	%rax, 16(%rbx)
	movq	8(%rbp), %rax
.L734:
	movq	%rax, 8(%rbx)
	movq	%rbx, %rax
	movq	%r12, 0(%rbp)
	movq	$0, 8(%rbp)
	movb	$0, 16(%rbp)
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L746:
	.cfi_restore_state
	movq	%r9, %rcx
	xorl	%edx, %edx
	xorl	%esi, %esi
	movq	%rax, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE10_M_replaceEmmPKcm@PLT
	jmp	.L744
	.p2align 4,,10
	.p2align 3
.L748:
	movq	8(%rbp), %rdx
	movq	$-1, %rax
	addq	$1, %rdx
	je	.L734
	movq	%r12, %rsi
	call	memcpy@PLT
	movq	8(%rbp), %rax
	jmp	.L734
	.p2align 4,,10
	.p2align 3
.L735:
	movl	$15, %r10d
	jmp	.L728
	.p2align 4,,10
	.p2align 3
.L745:
	cmpq	$15, %rcx
	jbe	.L726
	leaq	16(%rax), %r10
	cmpq	%r10, %rsi
	jne	.L727
	jmp	.L726
.L747:
	leaq	.LC35(%rip), %rdi
	call	_ZSt20__throw_length_errorPKc@PLT
	.cfi_endproc
.LFE10269:
	.size	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_, .-_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_
	.section	.rodata._ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC2IiiEERKT_RKT0_.str1.8,"aMS",@progbits,1
	.align 8
.LC36:
	.string	"void Eigen::PlainObjectBase<Derived>::resize(Eigen::Index, Eigen::Index) [with Derived = Eigen::Matrix<std::complex<double>, -1, -1>; Eigen::Index = long int]"
	.align 8
.LC37:
	.string	"/usr/local/include/Eigen/src/Core/PlainObjectBase.h"
	.align 8
.LC38:
	.ascii	"(!(RowsAtCompileTime!=Dynamic) || (rows==RowsAtCompileTime))"
	.ascii	" && (!(ColsAtCompileTime!=Dynamic) || (cols==ColsAtCompileTi"
	.ascii	"me)) && (!(RowsAtCompi"
	.string	"leTime==Dynamic && MaxRowsAtCompileTime!=Dynamic) || (rows<=MaxRowsAtCompileTime)) && (!(ColsAtCompileTime==Dynamic && MaxColsAtCompileTime!=Dynamic) || (cols<=MaxColsAtCompileTime)) && rows>=0 && cols>=0 && \"Invalid sizes when resizing a matrix or array.\""
	.section	.text._ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC2IiiEERKT_RKT0_,"axG",@progbits,_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC5IiiEERKT_RKT0_,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC2IiiEERKT_RKT0_
	.type	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC2IiiEERKT_RKT0_, @function
_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC2IiiEERKT_RKT0_:
.LFB10561:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA10561
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	pxor	%xmm0, %xmm0
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movslq	(%rdx), %rbp
	movslq	(%rsi), %r12
	movq	$0, (%rdi)
	movq	%rbp, %rax
	movups	%xmm0, 8(%rdi)
	orl	%r12d, %eax
	js	.L762
	movq	%rdi, %rbx
	movq	%rbp, %rdi
	imulq	%r12, %rdi
	testq	%rdi, %rdi
	jne	.L763
	movq	%r12, 8(%rbx)
	movq	%rbp, 16(%rbx)
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L763:
	.cfi_restore_state
	movabsq	$1152921504606846975, %rax
	cmpq	%rax, %rdi
	jg	.L764
	salq	$4, %rdi
	call	malloc@PLT
	testb	$15, %al
	je	.L753
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L753:
	testq	%rax, %rax
	je	.L765
	movq	%rax, (%rbx)
	movq	%r12, 8(%rbx)
	movq	%rbp, 16(%rbx)
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
.L762:
	.cfi_restore_state
	leaq	.LC36(%rip), %rcx
	movl	$273, %edx
	leaq	.LC37(%rip), %rsi
	leaq	.LC38(%rip), %rdi
	call	__assert_fail@PLT
.L765:
.LEHB22:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.L764:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE22:
.L756:
	movq	%rax, %rbp
.L755:
	movq	(%rbx), %rdi
	call	free@PLT
	movq	%rbp, %rdi
.LEHB23:
	call	_Unwind_Resume@PLT
.LEHE23:
	.cfi_endproc
.LFE10561:
	.section	.gcc_except_table
.LLSDA10561:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE10561-.LLSDACSB10561
.LLSDACSB10561:
	.uleb128 .LEHB22-.LFB10561
	.uleb128 .LEHE22-.LEHB22
	.uleb128 .L756-.LFB10561
	.uleb128 0
	.uleb128 .LEHB23-.LFB10561
	.uleb128 .LEHE23-.LEHB23
	.uleb128 0
	.uleb128 0
.LLSDACSE10561:
	.section	.text._ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC2IiiEERKT_RKT0_,"axG",@progbits,_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC5IiiEERKT_RKT0_,comdat
	.size	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC2IiiEERKT_RKT0_, .-_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC2IiiEERKT_RKT0_
	.weak	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
	.set	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_,_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC2IiiEERKT_RKT0_
	.section	.rodata._ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE.str1.8,"aMS",@progbits,1
	.align 8
.LC39:
	.ascii	"Eigen::CommaInitializer<MatrixType>& Eigen::CommaInitializer"
	.ascii	"<MatrixType"
	.string	">::operator,(const Eigen::DenseBase<OtherDerived>&) [with OtherDerived = Eigen::CwiseNullaryOp<Eigen::internal::scalar_constant_op<std::complex<double> >, Eigen::Matrix<std::complex<double>, -1, -1> >; XprType = Eigen::Matrix<std::complex<double>, -1, -1>]"
	.section	.text._ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	.type	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE, @function
_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE:
.LFB10570:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsi, %r10
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	subq	$24, %rsp
	.cfi_def_cfa_offset 48
	movq	(%rdi), %rbx
	movq	16(%rdi), %rcx
	movq	8(%rsi), %r8
	movq	16(%rbx), %rax
	cmpq	%rax, %rcx
	je	.L798
	leaq	(%rcx,%r8), %rdx
	cmpq	%rax, %rdx
	jg	.L787
	movq	(%rsi), %rsi
	cmpq	24(%rdi), %rsi
	jne	.L799
.L769:
	movq	8(%rbx), %r9
	movq	%rcx, %rdx
	movq	8(%rdi), %r11
	imulq	%r9, %rdx
.L786:
	movq	%rsi, %rbp
	addq	%r11, %rdx
	notq	%rbp
	salq	$4, %rdx
	shrq	$63, %rbp
	addq	(%rbx), %rdx
	je	.L772
	movq	%rsi, %rbx
	orq	%r8, %rbx
	js	.L800
.L772:
	movq	%r11, %rbx
	orq	%rcx, %rbx
	orq	%r8, %rbx
	js	.L773
	testb	%bpl, %bpl
	je	.L773
	movq	%r9, %rbx
	subq	%rsi, %rbx
	cmpq	%rbx, %r11
	jg	.L773
	subq	%r8, %rax
	cmpq	%rcx, %rax
	jl	.L773
	movupd	16(%r10), %xmm3
	movaps	%xmm3, (%rsp)
	testb	$15, %dl
	jne	.L775
	testq	%r8, %r8
	jle	.L782
	testq	%rsi, %rsi
	jle	.L785
	salq	$4, %rsi
	salq	$4, %r9
	xorl	%ecx, %ecx
	addq	%rsi, %rdx
	.p2align 4,,10
	.p2align 3
.L784:
	movq	%rdx, %rax
	subq	%rsi, %rax
	.p2align 4,,10
	.p2align 3
.L783:
	movapd	(%rsp), %xmm2
	addq	$16, %rax
	movups	%xmm2, -16(%rax)
	cmpq	%rax, %rdx
	jne	.L783
	addq	$1, %rcx
	addq	%r9, %rdx
	cmpq	%r8, %rcx
	jne	.L784
	movq	8(%r10), %r8
.L785:
	movq	16(%rdi), %rcx
.L782:
	leaq	(%rcx,%r8), %rax
	movq	%rax, 16(%rdi)
	addq	$24, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	movq	%rdi, %rax
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
.L798:
	.cfi_restore_state
	movq	24(%rdi), %r11
	movq	(%rsi), %rsi
	testq	%r8, %r8
	je	.L801
.L768:
	addq	8(%rdi), %r11
	movq	8(%rbx), %r9
	movq	$0, 16(%rdi)
	leaq	(%r11,%rsi), %rdx
	movq	%r11, 8(%rdi)
	movq	%rsi, 24(%rdi)
	cmpq	%r9, %rdx
	jg	.L802
	cmpq	%r8, %rax
	jl	.L787
	xorl	%edx, %edx
	xorl	%ecx, %ecx
	jmp	.L786
.L775:
	testq	%r8, %r8
	jle	.L782
	testq	%rsi, %rsi
	jle	.L782
	salq	$4, %rsi
	movsd	(%rsp), %xmm1
	salq	$4, %r9
	xorl	%r10d, %r10d
	movsd	8(%rsp), %xmm0
	addq	%rsi, %rdx
	.p2align 4,,10
	.p2align 3
.L780:
	movq	%rdx, %rax
	subq	%rsi, %rax
	.p2align 4,,10
	.p2align 3
.L781:
	movsd	%xmm1, (%rax)
	addq	$16, %rax
	movsd	%xmm0, -8(%rax)
	cmpq	%rdx, %rax
	jne	.L781
	addq	$1, %r10
	addq	%r9, %rdx
	cmpq	%r8, %r10
	jne	.L780
	jmp	.L782
.L801:
	cmpq	%rsi, %r11
	jne	.L768
	jmp	.L769
.L802:
	leaq	.LC39(%rip), %rcx
	movl	$92, %edx
	leaq	.LC10(%rip), %rsi
	leaq	.LC13(%rip), %rdi
	call	__assert_fail@PLT
	.p2align 4,,10
	.p2align 3
.L773:
	leaq	.LC19(%rip), %rcx
	movl	$146, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC21(%rip), %rdi
	call	__assert_fail@PLT
.L800:
	leaq	.LC16(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L799:
	leaq	.LC39(%rip), %rcx
	movl	$97, %edx
	leaq	.LC10(%rip), %rsi
	leaq	.LC15(%rip), %rdi
	call	__assert_fail@PLT
.L787:
	leaq	.LC39(%rip), %rcx
	movl	$95, %edx
	leaq	.LC10(%rip), %rsi
	leaq	.LC14(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE10570:
	.size	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE, .-_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	.section	.rodata._ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll.str1.8,"aMS",@progbits,1
	.align 8
.LC40:
	.string	"Eigen::DenseCoeffsBase<Derived, 1>::Scalar& Eigen::DenseCoeffsBase<Derived, 1>::operator()(Eigen::Index, Eigen::Index) [with Derived = Eigen::Matrix<std::complex<double>, -1, -1>; Scalar = std::complex<double>; Eigen::Index = long int]"
	.section	.text._ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll,"axG",@progbits,_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	.type	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll, @function
_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll:
.LFB10575:
	.cfi_startproc
	testq	%rsi, %rsi
	js	.L804
	movq	8(%rdi), %rax
	cmpq	%rax, %rsi
	jge	.L804
	testq	%rdx, %rdx
	js	.L804
	cmpq	16(%rdi), %rdx
	jge	.L804
	imulq	%rax, %rdx
	leaq	(%rdx,%rsi), %rax
	salq	$4, %rax
	addq	(%rdi), %rax
	ret
.L804:
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	leaq	.LC40(%rip), %rcx
	movl	$366, %edx
	leaq	.LC4(%rip), %rsi
	leaq	.LC5(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE10575:
	.size	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll, .-_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	.section	.text._ZN5Eigen8IOFormatD2Ev,"axG",@progbits,_ZN5Eigen8IOFormatD5Ev,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen8IOFormatD2Ev
	.type	_ZN5Eigen8IOFormatD2Ev, @function
_ZN5Eigen8IOFormatD2Ev:
.LFB10588:
	.cfi_startproc
	pushq	%rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	movq	%rdi, %rbx
	movq	192(%rdi), %rdi
	leaq	208(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L812
	movq	208(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L812:
	movq	160(%rbx), %rdi
	leaq	176(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L813
	movq	176(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L813:
	movq	128(%rbx), %rdi
	leaq	144(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L814
	movq	144(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L814:
	movq	96(%rbx), %rdi
	leaq	112(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L815
	movq	112(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L815:
	movq	64(%rbx), %rdi
	leaq	80(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L816
	movq	80(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L816:
	movq	32(%rbx), %rdi
	leaq	48(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L817
	movq	48(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L817:
	movq	(%rbx), %rdi
	leaq	16(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L811
	movq	16(%rbx), %rsi
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	addq	$1, %rsi
	jmp	_ZdlPvm@PLT
	.p2align 4,,10
	.p2align 3
.L811:
	.cfi_restore_state
	popq	%rbx
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE10588:
	.size	_ZN5Eigen8IOFormatD2Ev, .-_ZN5Eigen8IOFormatD2Ev
	.weak	_ZN5Eigen8IOFormatD1Ev
	.set	_ZN5Eigen8IOFormatD1Ev,_ZN5Eigen8IOFormatD2Ev
	.section	.rodata._ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll.str1.8,"aMS",@progbits,1
	.align 8
.LC41:
	.string	"Eigen::DenseCoeffsBase<Derived, 1>::Scalar& Eigen::DenseCoeffsBase<Derived, 1>::operator()(Eigen::Index, Eigen::Index) [with Derived = Eigen::Matrix<double, -1, -1>; Scalar = double; Eigen::Index = long int]"
	.section	.text._ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll,"axG",@progbits,_ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	.type	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll, @function
_ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll:
.LFB10605:
	.cfi_startproc
	testq	%rsi, %rsi
	js	.L821
	movq	8(%rdi), %rax
	cmpq	%rax, %rsi
	jge	.L821
	testq	%rdx, %rdx
	js	.L821
	cmpq	16(%rdi), %rdx
	jge	.L821
	imulq	%rax, %rdx
	movq	(%rdi), %rax
	addq	%rsi, %rdx
	leaq	(%rax,%rdx,8), %rax
	ret
.L821:
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	leaq	.LC41(%rip), %rcx
	movl	$366, %edx
	leaq	.LC4(%rip), %rsi
	leaq	.LC5(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE10605:
	.size	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll, .-_ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	.section	.text._ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev,"axG",@progbits,_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED5Ev,comdat
	.align 2
	.p2align 4
	.weak	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev
	.type	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev, @function
_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev:
.LFB10628:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	movq	%rdi, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movq	8(%rdi), %rbp
	movq	(%rdi), %rbx
	cmpq	%rbx, %rbp
	je	.L829
	.p2align 4,,10
	.p2align 3
.L833:
	movq	(%rbx), %rdi
	leaq	16(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L830
	movq	16(%rbx), %rax
	addq	$32, %rbx
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	cmpq	%rbx, %rbp
	jne	.L833
.L832:
	movq	(%r12), %rbx
.L829:
	testq	%rbx, %rbx
	je	.L828
	movq	16(%r12), %rsi
	movq	%rbx, %rdi
	subq	%rbx, %rsi
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	jmp	_ZdlPvm@PLT
	.p2align 4,,10
	.p2align 3
.L830:
	.cfi_restore_state
	addq	$32, %rbx
	cmpq	%rbx, %rbp
	jne	.L833
	jmp	.L832
	.p2align 4,,10
	.p2align 3
.L828:
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE10628:
	.size	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev, .-_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev
	.weak	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED1Ev
	.set	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED1Ev,_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev
	.section	.text._ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED2Ev,"axG",@progbits,_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED5Ev,comdat
	.align 2
	.p2align 4
	.weak	_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED2Ev
	.type	_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED2Ev, @function
_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED2Ev:
.LFB10633:
	.cfi_startproc
	leaq	16+_ZTVNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEEE(%rip), %rax
	pushq	%rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	movq	%rdi, %rbx
	movq	%rax, (%rdi)
	movq	72(%rdi), %rdi
	leaq	88(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L837
	movq	88(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L837:
	leaq	16+_ZTVSt15basic_streambufIcSt11char_traitsIcEE(%rip), %rax
	leaq	56(%rbx), %rdi
	movq	%rax, (%rbx)
	popq	%rbx
	.cfi_def_cfa_offset 8
	jmp	_ZNSt6localeD1Ev@PLT
	.cfi_endproc
.LFE10633:
	.size	_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED2Ev, .-_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED2Ev
	.weak	_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED1Ev
	.set	_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED1Ev,_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED2Ev
	.section	.text._ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED0Ev,"axG",@progbits,_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED5Ev,comdat
	.align 2
	.p2align 4
	.weak	_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED0Ev
	.type	_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED0Ev, @function
_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED0Ev:
.LFB10635:
	.cfi_startproc
	leaq	16+_ZTVNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEEE(%rip), %rax
	pushq	%rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	movq	%rdi, %rbx
	movq	%rax, (%rdi)
	movq	72(%rdi), %rdi
	leaq	88(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L840
	movq	88(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L840:
	leaq	16+_ZTVSt15basic_streambufIcSt11char_traitsIcEE(%rip), %rax
	leaq	56(%rbx), %rdi
	movq	%rax, (%rbx)
	call	_ZNSt6localeD1Ev@PLT
	movq	%rbx, %rdi
	movl	$104, %esi
	popq	%rbx
	.cfi_def_cfa_offset 8
	jmp	_ZdlPvm@PLT
	.cfi_endproc
.LFE10635:
	.size	_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED0Ev, .-_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED0Ev
	.section	.rodata._ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC2IS4_EERS4_RKNS_9DenseBaseIT_EE.str1.8,"aMS",@progbits,1
	.align 8
.LC42:
	.string	"Eigen::CommaInitializer<MatrixType>::CommaInitializer(XprType&, const Eigen::DenseBase<OtherDerived>&) [with OtherDerived = Eigen::Matrix<std::complex<double>, -1, -1>; XprType = Eigen::Matrix<std::complex<double>, -1, -1>]"
	.align 8
.LC43:
	.string	"m_xpr.rows() >= other.rows() && m_xpr.cols() >= other.cols() && \"Cannot comma-initialize a 0x0 matrix (operator<<)\""
	.section	.text._ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC2IS4_EERS4_RKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC5IS4_EERS4_RKNS_9DenseBaseIT_EE,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC2IS4_EERS4_RKNS_9DenseBaseIT_EE
	.type	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC2IS4_EERS4_RKNS_9DenseBaseIT_EE, @function
_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC2IS4_EERS4_RKNS_9DenseBaseIT_EE:
.LFB10901:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	movq	16(%rdx), %r8
	movq	8(%rdx), %rcx
	movq	%rsi, (%rdi)
	movq	8(%rsi), %r9
	movq	$0, 8(%rdi)
	movq	%r8, 16(%rdi)
	movq	%rcx, 24(%rdi)
	cmpq	%r9, %rcx
	jg	.L843
	cmpq	%r8, 16(%rsi)
	jl	.L843
	movq	(%rsi), %rdi
	movq	%r8, %rax
	testq	%rdi, %rdi
	je	.L845
	orq	%rcx, %rax
	js	.L859
	movq	(%rdx), %rdx
	movq	%rdi, %r11
	testb	$15, %dil
	jne	.L860
.L848:
	testq	%r8, %r8
	jle	.L842
	testq	%rcx, %rcx
	jle	.L842
	salq	$4, %rcx
	xorl	%r10d, %r10d
	xorl	%edi, %edi
	.p2align 4,,10
	.p2align 3
.L851:
	movq	%r10, %rsi
	xorl	%eax, %eax
	salq	$4, %rsi
	addq	%r11, %rsi
	.p2align 4,,10
	.p2align 3
.L854:
	movupd	(%rdx,%rax), %xmm1
	movups	%xmm1, (%rsi,%rax)
	addq	$16, %rax
	cmpq	%rax, %rcx
	jne	.L854
	addq	$1, %rdi
	addq	%rcx, %rdx
	addq	%r9, %r10
	cmpq	%rdi, %r8
	jne	.L851
.L842:
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
.L860:
	.cfi_restore_state
	testq	%r8, %r8
	jle	.L842
	testq	%rcx, %rcx
	jle	.L842
	salq	$4, %rcx
	xorl	%r11d, %r11d
	xorl	%r10d, %r10d
	.p2align 4,,10
	.p2align 3
.L852:
	movq	%r11, %rsi
	xorl	%eax, %eax
	salq	$4, %rsi
	addq	%rdi, %rsi
	.p2align 4,,10
	.p2align 3
.L853:
	movsd	(%rdx,%rax), %xmm0
	movsd	%xmm0, (%rsi,%rax)
	movsd	8(%rdx,%rax), %xmm0
	movsd	%xmm0, 8(%rsi,%rax)
	addq	$16, %rax
	cmpq	%rcx, %rax
	jne	.L853
	addq	$1, %r10
	addq	%r9, %r11
	addq	%rcx, %rdx
	cmpq	%r10, %r8
	jne	.L852
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
.L845:
	.cfi_restore_state
	orq	%rcx, %rax
	js	.L849
	movq	(%rdx), %rdx
	xorl	%r11d, %r11d
	jmp	.L848
.L843:
	leaq	.LC42(%rip), %rcx
	movl	$46, %edx
	leaq	.LC10(%rip), %rsi
	leaq	.LC43(%rip), %rdi
	call	__assert_fail@PLT
.L849:
	leaq	.LC19(%rip), %rcx
	movl	$146, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC21(%rip), %rdi
	call	__assert_fail@PLT
.L859:
	leaq	.LC16(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE10901:
	.size	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC2IS4_EERS4_RKNS_9DenseBaseIT_EE, .-_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC2IS4_EERS4_RKNS_9DenseBaseIT_EE
	.weak	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC1IS4_EERS4_RKNS_9DenseBaseIT_EE
	.set	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC1IS4_EERS4_RKNS_9DenseBaseIT_EE,_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC2IS4_EERS4_RKNS_9DenseBaseIT_EE
	.section	.rodata._ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE.str1.8,"aMS",@progbits,1
	.align 8
.LC44:
	.string	"void Eigen::PlainObjectBase<Derived>::resize(Eigen::Index, Eigen::Index) [with Derived = Eigen::Matrix<double, -1, -1>; Eigen::Index = long int]"
	.section	.text._ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC5INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE
	.type	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE, @function
_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE:
.LFB11006:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA11006
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	pxor	%xmm0, %xmm0
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	movq	%rsi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movq	%rdi, %rbx
	subq	$16, %rsp
	.cfi_def_cfa_offset 48
	movups	%xmm0, 8(%rdi)
	movdqu	(%rsi), %xmm0
	movq	8(%rsi), %rsi
	movq	$0, (%rdi)
	movq	0(%rbp), %rdi
	movhlps	%xmm0, %xmm2
	movq	%xmm0, %rcx
	movq	%xmm2, %r8
	testq	%rdi, %rdi
	je	.L862
	testq	%rsi, %rsi
	je	.L862
	movabsq	$9223372036854775807, %rax
	cqto
	idivq	%r8
	cmpq	%rcx, %rax
	jl	.L896
	orq	%rdi, %rsi
	js	.L865
	imulq	%r8, %rcx
.L867:
	movabsq	$2305843009213693951, %rax
	cmpq	%rax, %rcx
	jg	.L897
	leaq	0(,%rcx,8), %r12
	movaps	%xmm0, (%rsp)
	movq	%r12, %rdi
	call	malloc@PLT
	cmpq	$15, %r12
	movdqa	(%rsp), %xmm0
	movq	%rax, %rdx
	ja	.L898
.L870:
	testq	%rax, %rax
	je	.L899
	leaq	-8(%r12), %rsi
	movups	%xmm0, 8(%rbx)
	movsd	16(%rbp), %xmm0
	shrq	$3, %rsi
	movq	%rax, (%rbx)
	addq	$1, %rsi
	cmpq	$8, %r12
	je	.L872
	movq	%rsi, %rcx
	movapd	%xmm0, %xmm1
	movq	%rax, %rdx
	shrq	%rcx
	unpcklpd	%xmm1, %xmm1
	salq	$4, %rcx
	addq	%rax, %rcx
	.p2align 4,,10
	.p2align 3
.L874:
	movups	%xmm1, (%rdx)
	addq	$16, %rdx
	cmpq	%rcx, %rdx
	jne	.L874
	testb	$1, %sil
	je	.L861
	andq	$-2, %rsi
	leaq	(%rax,%rsi,8), %rdx
.L872:
	movsd	%xmm0, (%rdx)
	addq	$16, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L862:
	.cfi_restore_state
	orq	%rdi, %rsi
	js	.L865
	imulq	%r8, %rcx
	testq	%rcx, %rcx
	jne	.L867
	movups	%xmm0, 8(%rbx)
.L861:
	addq	$16, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L898:
	.cfi_restore_state
	testb	$15, %al
	je	.L870
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L865:
	leaq	.LC44(%rip), %rcx
	movl	$273, %edx
	leaq	.LC37(%rip), %rsi
	leaq	.LC38(%rip), %rdi
	call	__assert_fail@PLT
.L899:
.LEHB24:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.L896:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.L897:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE24:
.L878:
	movq	%rax, %rbp
.L877:
	movq	(%rbx), %rdi
	call	free@PLT
	movq	%rbp, %rdi
.LEHB25:
	call	_Unwind_Resume@PLT
.LEHE25:
	.cfi_endproc
.LFE11006:
	.section	.gcc_except_table
.LLSDA11006:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE11006-.LLSDACSB11006
.LLSDACSB11006:
	.uleb128 .LEHB24-.LFB11006
	.uleb128 .LEHE24-.LEHB24
	.uleb128 .L878-.LFB11006
	.uleb128 0
	.uleb128 .LEHB25-.LFB11006
	.uleb128 .LEHE25-.LEHB25
	.uleb128 0
	.uleb128 0
.LLSDACSE11006:
	.section	.text._ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC5INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE,comdat
	.size	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE, .-_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE
	.weak	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC1INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE
	.set	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC1INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE,_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE
	.section	.text._ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC5INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE
	.type	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE, @function
_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE:
.LFB11024:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA11024
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	pxor	%xmm0, %xmm0
	movq	%rsi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 24
	.cfi_offset 3, -24
	movq	%rdi, %rbx
	subq	$8, %rsp
	.cfi_def_cfa_offset 32
	movq	$0, (%rdi)
	movups	%xmm0, 8(%rdi)
	movl	$72, %edi
	call	malloc@PLT
	testb	$15, %al
	je	.L901
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L901:
	testq	%rax, %rax
	je	.L907
	movq	0(%rbp), %rdx
	movdqa	.LC45(%rip), %xmm0
	movq	%rax, (%rbx)
	movq	16(%rbp), %rcx
	movups	%xmm0, 8(%rbx)
	movsd	8(%rdx), %xmm0
	movq	8(%rcx), %rcx
	movsd	%xmm0, (%rax)
	movsd	24(%rdx), %xmm0
	movq	%rcx, %rsi
	movsd	%xmm0, 8(%rax)
	movsd	40(%rdx), %xmm0
	salq	$4, %rsi
	movsd	%xmm0, 16(%rax)
	movsd	8(%rdx,%rsi), %xmm0
	movsd	%xmm0, 24(%rax)
	movsd	24(%rdx,%rsi), %xmm0
	movsd	%xmm0, 32(%rax)
	movsd	40(%rdx,%rsi), %xmm0
	movq	%rcx, %rsi
	salq	$5, %rsi
	movsd	%xmm0, 40(%rax)
	movsd	8(%rdx,%rsi), %xmm0
	movsd	%xmm0, 48(%rax)
	movsd	24(%rsi,%rdx), %xmm0
	movsd	%xmm0, 56(%rax)
	movsd	40(%rdx,%rsi), %xmm0
	movsd	%xmm0, 64(%rax)
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	ret
.L907:
	.cfi_restore_state
.LEHB26:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE26:
.L904:
	movq	%rax, %rbp
.L903:
	movq	(%rbx), %rdi
	call	free@PLT
	movq	%rbp, %rdi
.LEHB27:
	call	_Unwind_Resume@PLT
.LEHE27:
	.cfi_endproc
.LFE11024:
	.section	.gcc_except_table
.LLSDA11024:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE11024-.LLSDACSB11024
.LLSDACSB11024:
	.uleb128 .LEHB26-.LFB11024
	.uleb128 .LEHE26-.LEHB26
	.uleb128 .L904-.LFB11024
	.uleb128 0
	.uleb128 .LEHB27-.LFB11024
	.uleb128 .LEHE27-.LEHB27
	.uleb128 0
	.uleb128 0
.LLSDACSE11024:
	.section	.text._ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC5INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE,comdat
	.size	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE, .-_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE
	.weak	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC1INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE
	.set	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC1INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE,_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE
	.section	.rodata._ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_.str1.1,"aMS",@progbits,1
.LC47:
	.string	"vector::_M_realloc_insert"
	.section	.text._ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_,"axG",@progbits,_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_,comdat
	.align 2
	.p2align 4
	.weak	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_
	.type	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_, @function
_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_:
.LFB11172:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA11172
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movq	%rdx, %rcx
	movabsq	$288230376151711743, %rdx
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$56, %rsp
	.cfi_def_cfa_offset 112
	movq	8(%rdi), %r14
	movq	(%rdi), %r12
	movq	%r14, %rax
	subq	%r12, %rax
	sarq	$5, %rax
	cmpq	%rdx, %rax
	je	.L947
	cmpq	%r14, %r12
	movl	$1, %edx
	movq	%rsi, %r15
	movq	%rdi, %r13
	cmovne	%rax, %rdx
	movq	%rsi, %rbx
	addq	%rdx, %rax
	setc	%dl
	movq	%rax, 8(%rsp)
	subq	%r12, %r15
	movzbl	%dl, %edx
	testq	%rdx, %rdx
	jne	.L932
	testq	%rax, %rax
	jne	.L948
	xorl	%ebp, %ebp
.L914:
	addq	%rbp, %r15
	movq	(%rcx), %rsi
	movq	8(%rcx), %rdx
	leaq	16(%r15), %rax
	movq	%r15, %rdi
	movq	%rax, (%r15)
	addq	%rsi, %rdx
	movq	%rax, 16(%rsp)
.LEHB28:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
.LEHE28:
	cmpq	%r12, %rbx
	je	.L934
	leaq	16(%r12), %r15
	leaq	16(%rbx), %r8
	movq	%rbp, %rcx
	jmp	.L920
	.p2align 4,,10
	.p2align 3
.L916:
	movq	%rdx, (%rcx)
	movq	(%r15), %rdx
	movq	%rdx, 16(%rcx)
.L918:
	addq	$32, %r15
	movq	%rax, 8(%rcx)
	addq	$32, %rcx
	cmpq	%r15, %r8
	je	.L949
.L920:
	leaq	16(%rcx), %rdi
	movq	-8(%r15), %rax
	movq	%rdi, (%rcx)
	movq	-16(%r15), %rdx
	cmpq	%rdx, %r15
	jne	.L916
	movq	%rax, %rdx
	addq	$1, %rdx
	je	.L918
	movq	%r15, %rsi
	movq	%r8, 24(%rsp)
	movq	%rcx, 16(%rsp)
	call	memcpy@PLT
	movq	-8(%r15), %rax
	movq	16(%rsp), %rcx
	movq	24(%rsp), %r8
	jmp	.L918
	.p2align 4,,10
	.p2align 3
.L949:
	movq	%rbx, %r10
	subq	%r12, %r10
	addq	%rbp, %r10
.L915:
	addq	$32, %r10
	cmpq	%r14, %rbx
	je	.L921
	leaq	16(%rbx), %r15
	leaq	16(%r14), %r9
	movq	%r10, %rcx
	jmp	.L924
	.p2align 4,,10
	.p2align 3
.L922:
	movq	%rax, (%rcx)
	movq	(%r15), %rax
	movq	%rax, 16(%rcx)
.L923:
	addq	$32, %r15
	movq	%r8, 8(%rcx)
	addq	$32, %rcx
	cmpq	%r9, %r15
	je	.L950
.L924:
	movq	-16(%r15), %rax
	leaq	16(%rcx), %rdi
	movq	-8(%r15), %r8
	movq	%rdi, (%rcx)
	cmpq	%rax, %r15
	jne	.L922
	movq	%r8, %rdx
	addq	$1, %rdx
	je	.L923
	movq	%r15, %rsi
	movq	%r8, 40(%rsp)
	movq	%r9, 32(%rsp)
	movq	%r10, 24(%rsp)
	movq	%rcx, 16(%rsp)
	call	memcpy@PLT
	movq	16(%rsp), %rcx
	movq	24(%rsp), %r10
	movq	32(%rsp), %r9
	movq	40(%rsp), %r8
	jmp	.L923
	.p2align 4,,10
	.p2align 3
.L950:
	subq	%rbx, %r14
	addq	%r14, %r10
.L921:
	testq	%r12, %r12
	je	.L925
	movq	16(%r13), %rsi
	movq	%r12, %rdi
	movq	%r10, 16(%rsp)
	subq	%r12, %rsi
	call	_ZdlPvm@PLT
	movq	16(%rsp), %r10
.L925:
	movq	8(%rsp), %rax
	movq	%rbp, 0(%r13)
	movq	%r10, 8(%r13)
	addq	%rax, %rbp
	movq	%rbp, 16(%r13)
	addq	$56, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L932:
	.cfi_restore_state
	movabsq	$9223372036854775776, %rax
	movq	%rax, 8(%rsp)
	movq	%rax, %rdi
.L913:
	movq	%rcx, 16(%rsp)
.LEHB29:
	call	_Znwm@PLT
	movq	16(%rsp), %rcx
	movq	%rax, %rbp
	jmp	.L914
	.p2align 4,,10
	.p2align 3
.L934:
	movq	%rbp, %r10
	jmp	.L915
.L948:
	movq	%rax, %rsi
	movabsq	$288230376151711743, %rax
	cmpq	%rax, %rsi
	cmovbe	%rsi, %rax
	salq	$5, %rax
	movq	%rax, 8(%rsp)
	movq	%rax, %rdi
	jmp	.L913
.L947:
	leaq	.LC47(%rip), %rdi
	call	_ZSt20__throw_length_errorPKc@PLT
.LEHE29:
.L935:
	movq	%rax, %rdi
.L926:
	call	__cxa_begin_catch@PLT
	testq	%rbp, %rbp
	je	.L927
	movq	8(%rsp), %rsi
	movq	%rbp, %rdi
	call	_ZdlPvm@PLT
.L928:
.LEHB30:
	call	__cxa_rethrow@PLT
.LEHE30:
.L927:
	movq	(%r15), %rdi
	cmpq	%rdi, 16(%rsp)
	je	.L928
	movq	16(%r15), %rsi
	addq	$1, %rsi
	call	_ZdlPvm@PLT
	jmp	.L928
.L936:
	movq	%rax, %rbx
.L930:
	call	__cxa_end_catch@PLT
	movq	%rbx, %rdi
.LEHB31:
	call	_Unwind_Resume@PLT
.LEHE31:
	.cfi_endproc
.LFE11172:
	.section	.gcc_except_table
	.align 4
.LLSDA11172:
	.byte	0xff
	.byte	0x9b
	.uleb128 .LLSDATT11172-.LLSDATTD11172
.LLSDATTD11172:
	.byte	0x1
	.uleb128 .LLSDACSE11172-.LLSDACSB11172
.LLSDACSB11172:
	.uleb128 .LEHB28-.LFB11172
	.uleb128 .LEHE28-.LEHB28
	.uleb128 .L935-.LFB11172
	.uleb128 0x1
	.uleb128 .LEHB29-.LFB11172
	.uleb128 .LEHE29-.LEHB29
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB30-.LFB11172
	.uleb128 .LEHE30-.LEHB30
	.uleb128 .L936-.LFB11172
	.uleb128 0
	.uleb128 .LEHB31-.LFB11172
	.uleb128 .LEHE31-.LEHB31
	.uleb128 0
	.uleb128 0
.LLSDACSE11172:
	.byte	0x1
	.byte	0
	.align 4
	.long	0

.LLSDATT11172:
	.section	.text._ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_,"axG",@progbits,_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_,comdat
	.size	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_, .-_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_
	.section	.text._ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE9push_backERKS5_,"axG",@progbits,_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE9push_backERKS5_,comdat
	.align 2
	.p2align 4
	.weak	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE9push_backERKS5_
	.type	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE9push_backERKS5_, @function
_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE9push_backERKS5_:
.LFB10642:
	.cfi_startproc
	pushq	%rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	movq	%rdi, %rbx
	movq	%rsi, %rdx
	movq	8(%rdi), %rdi
	cmpq	16(%rbx), %rdi
	je	.L952
	leaq	16(%rdi), %rax
	movq	%rax, (%rdi)
	movq	8(%rdx), %rax
	movq	(%rsi), %rsi
	addq	%rsi, %rax
	movq	%rax, %rdx
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
	addq	$32, 8(%rbx)
	popq	%rbx
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L952:
	.cfi_restore_state
	movq	%rdi, %rsi
	movq	%rbx, %rdi
	popq	%rbx
	.cfi_def_cfa_offset 8
	jmp	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_M_realloc_insertIJRKS5_EEEvN9__gnu_cxx17__normal_iteratorIPS5_S7_EEDpOT_
	.cfi_endproc
.LFE10642:
	.size	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE9push_backERKS5_, .-_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE9push_backERKS5_
	.section	.text._ZNSt6vectorIdSaIdEE17_M_realloc_insertIJRKdEEEvN9__gnu_cxx17__normal_iteratorIPdS1_EEDpOT_,"axG",@progbits,_ZNSt6vectorIdSaIdEE17_M_realloc_insertIJRKdEEEvN9__gnu_cxx17__normal_iteratorIPdS1_EEDpOT_,comdat
	.align 2
	.p2align 4
	.weak	_ZNSt6vectorIdSaIdEE17_M_realloc_insertIJRKdEEEvN9__gnu_cxx17__normal_iteratorIPdS1_EEDpOT_
	.type	_ZNSt6vectorIdSaIdEE17_M_realloc_insertIJRKdEEEvN9__gnu_cxx17__normal_iteratorIPdS1_EEDpOT_, @function
_ZNSt6vectorIdSaIdEE17_M_realloc_insertIJRKdEEEvN9__gnu_cxx17__normal_iteratorIPdS1_EEDpOT_:
.LFB11178:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	movq	%rdx, %r15
	movabsq	$1152921504606846975, %rdx
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$24, %rsp
	.cfi_def_cfa_offset 80
	movq	8(%rdi), %r12
	movq	(%rdi), %r13
	movq	%r12, %rax
	subq	%r13, %rax
	sarq	$3, %rax
	cmpq	%rdx, %rax
	je	.L978
	cmpq	%r12, %r13
	movl	$1, %edx
	movq	%rdi, %rbp
	movq	%rsi, %r14
	cmovne	%rax, %rdx
	xorl	%ecx, %ecx
	addq	%rdx, %rax
	movq	%rsi, %rdx
	setc	%cl
	subq	%r13, %rdx
	testq	%rcx, %rcx
	jne	.L971
	testq	%rax, %rax
	jne	.L979
	xorl	%ebx, %ebx
	xorl	%ecx, %ecx
.L962:
	movsd	(%r15), %xmm0
	leaq	8(%rcx,%rdx), %r8
	subq	%r14, %r12
	leaq	(%r8,%r12), %r15
	movsd	%xmm0, (%rcx,%rdx)
	testq	%rdx, %rdx
	jg	.L980
	testq	%r12, %r12
	jg	.L967
	testq	%r13, %r13
	jne	.L965
.L968:
	movq	%rcx, 0(%rbp)
	movq	%r15, 8(%rbp)
	movq	%rbx, 16(%rbp)
	addq	$24, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L980:
	.cfi_restore_state
	movq	%rcx, %rdi
	movq	%r13, %rsi
	movq	%r8, (%rsp)
	call	memmove@PLT
	movq	%rax, %rcx
	testq	%r12, %r12
	jg	.L981
.L965:
	movq	16(%rbp), %rsi
	movq	%r13, %rdi
	movq	%rcx, (%rsp)
	subq	%r13, %rsi
	call	_ZdlPvm@PLT
	movq	(%rsp), %rcx
	jmp	.L968
	.p2align 4,,10
	.p2align 3
.L967:
	movq	%r12, %rdx
	movq	%r14, %rsi
	movq	%r8, %rdi
	movq	%rcx, (%rsp)
	call	memcpy@PLT
	movq	(%rsp), %rcx
	testq	%r13, %r13
	je	.L968
	jmp	.L965
	.p2align 4,,10
	.p2align 3
.L971:
	movabsq	$9223372036854775800, %rbx
.L961:
	movq	%rbx, %rdi
	movq	%rdx, (%rsp)
	call	_Znwm@PLT
	movq	(%rsp), %rdx
	movq	%rax, %rcx
	addq	%rax, %rbx
	jmp	.L962
	.p2align 4,,10
	.p2align 3
.L981:
	movq	(%rsp), %rdi
	movq	%r12, %rdx
	movq	%r14, %rsi
	movq	%rax, 8(%rsp)
	call	memcpy@PLT
	movq	8(%rsp), %rcx
	jmp	.L965
	.p2align 4,,10
	.p2align 3
.L979:
	movabsq	$1152921504606846975, %rcx
	cmpq	%rcx, %rax
	cmova	%rcx, %rax
	leaq	0(,%rax,8), %rbx
	jmp	.L961
.L978:
	leaq	.LC47(%rip), %rdi
	call	_ZSt20__throw_length_errorPKc@PLT
	.cfi_endproc
.LFE11178:
	.size	_ZNSt6vectorIdSaIdEE17_M_realloc_insertIJRKdEEEvN9__gnu_cxx17__normal_iteratorIPdS1_EEDpOT_, .-_ZNSt6vectorIdSaIdEE17_M_realloc_insertIJRKdEEEvN9__gnu_cxx17__normal_iteratorIPdS1_EEDpOT_
	.section	.text._ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC2ERKS3_,"axG",@progbits,_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC5ERKS3_,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC2ERKS3_
	.type	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC2ERKS3_, @function
_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC2ERKS3_:
.LFB11309:
	.cfi_startproc
	pushq	%r13
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	pushq	%r12
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	pushq	%rbp
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	movq	%rdi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
	subq	$8, %rsp
	.cfi_def_cfa_offset 48
	movq	8(%rsi), %r12
	testq	%r12, %r12
	je	.L983
	movq	%r12, %rax
	shrq	$60, %rax
	jne	.L986
	movq	%r12, %r13
	movq	%rsi, %rbx
	salq	$4, %r13
	movq	%r13, %rdi
	call	malloc@PLT
	testb	$15, %al
	je	.L985
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L985:
	testq	%rax, %rax
	je	.L986
	movq	%rax, 0(%rbp)
	movq	(%rbx), %rsi
	movq	%r13, %rdx
	movq	%rax, %rdi
	movq	%r12, 8(%rbp)
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	popq	%rbx
	.cfi_def_cfa_offset 32
	popq	%rbp
	.cfi_def_cfa_offset 24
	popq	%r12
	.cfi_def_cfa_offset 16
	popq	%r13
	.cfi_def_cfa_offset 8
	jmp	memcpy@PLT
	.p2align 4,,10
	.p2align 3
.L983:
	.cfi_restore_state
	movq	$0, (%rdi)
	movq	$0, 8(%rdi)
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	popq	%rbx
	.cfi_def_cfa_offset 32
	popq	%rbp
	.cfi_def_cfa_offset 24
	popq	%r12
	.cfi_def_cfa_offset 16
	popq	%r13
	.cfi_def_cfa_offset 8
	ret
.L986:
	.cfi_restore_state
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
	.cfi_endproc
.LFE11309:
	.size	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC2ERKS3_, .-_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC2ERKS3_
	.weak	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC1ERKS3_
	.set	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC1ERKS3_,_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC2ERKS3_
	.section	.text._ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC2ERKS3_,"axG",@progbits,_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC5ERKS3_,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC2ERKS3_
	.type	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC2ERKS3_, @function
_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC2ERKS3_:
.LFB11312:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	movq	%rdi, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	subq	$16, %rsp
	.cfi_def_cfa_offset 48
	movdqu	8(%rsi), %xmm0
	movhlps	%xmm0, %xmm1
	movq	%xmm0, %rax
	movq	%xmm1, %rdx
	imulq	%rax, %rdx
	testq	%rdx, %rdx
	je	.L992
	movq	%rdx, %rax
	shrq	$60, %rax
	jne	.L995
	movq	%rdx, %rbx
	movaps	%xmm0, (%rsp)
	movq	%rsi, %rbp
	salq	$4, %rbx
	movq	%rbx, %rdi
	call	malloc@PLT
	movdqa	(%rsp), %xmm0
	testb	$15, %al
	je	.L994
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L994:
	testq	%rax, %rax
	je	.L995
	movq	%rax, (%r12)
	movq	0(%rbp), %rsi
	movq	%rbx, %rdx
	movq	%rax, %rdi
	movups	%xmm0, 8(%r12)
	addq	$16, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	jmp	memcpy@PLT
	.p2align 4,,10
	.p2align 3
.L992:
	.cfi_restore_state
	movq	$0, (%rdi)
	movups	%xmm0, 8(%rdi)
	addq	$16, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
.L995:
	.cfi_restore_state
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
	.cfi_endproc
.LFE11312:
	.size	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC2ERKS3_, .-_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC2ERKS3_
	.weak	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC1ERKS3_
	.set	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC1ERKS3_,_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC2ERKS3_
	.section	.text._ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll,"axG",@progbits,_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
	.type	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll, @function
_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll:
.LFB11314:
	.cfi_startproc
	pushq	%r13
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	movq	%rsi, %rax
	pushq	%r12
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	pushq	%rbp
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	pushq	%rbx
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
	subq	$8, %rsp
	.cfi_def_cfa_offset 48
	orq	%rdx, %rax
	js	.L1021
	movq	%rdi, %rbx
	movq	%rsi, %r13
	movq	%rdx, %rbp
	testq	%rsi, %rsi
	je	.L1002
	testq	%rdx, %rdx
	je	.L1002
	movabsq	$9223372036854775807, %rax
	cqto
	idivq	%rbp
	cmpq	%rax, %rsi
	jg	.L1006
	movq	%rsi, %r12
	movq	16(%rdi), %rax
	imulq	8(%rdi), %rax
	imulq	%rbp, %r12
	cmpq	%r12, %rax
	jne	.L1022
.L1004:
	movq	%r13, 8(%rbx)
	movq	%rbp, 16(%rbx)
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	popq	%rbx
	.cfi_def_cfa_offset 32
	popq	%rbp
	.cfi_def_cfa_offset 24
	popq	%r12
	.cfi_def_cfa_offset 16
	popq	%r13
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L1002:
	.cfi_restore_state
	movq	%r13, %r12
	movq	8(%rbx), %rax
	imulq	16(%rbx), %rax
	imulq	%rbp, %r12
	cmpq	%rax, %r12
	je	.L1004
	movq	(%rbx), %rdi
	call	free@PLT
	testq	%r12, %r12
	jne	.L1008
	movq	$0, (%rbx)
	jmp	.L1004
	.p2align 4,,10
	.p2align 3
.L1022:
	movq	(%rbx), %rdi
	call	free@PLT
.L1008:
	movabsq	$1152921504606846975, %rax
	cmpq	%rax, %r12
	jg	.L1006
	movq	%r12, %rdi
	salq	$4, %rdi
	call	malloc@PLT
	testb	$15, %al
	je	.L1007
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L1007:
	testq	%rax, %rax
	je	.L1006
	movq	%rax, (%rbx)
	jmp	.L1004
.L1021:
	leaq	.LC36(%rip), %rcx
	movl	$273, %edx
	leaq	.LC37(%rip), %rsi
	leaq	.LC38(%rip), %rdi
	call	__assert_fail@PLT
.L1006:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
	.cfi_endproc
.LFE11314:
	.size	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll, .-_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
	.section	.text._ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE,"axG",@progbits,_ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE
	.type	_ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE, @function
_ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE:
.LFB10999:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA10999
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	movq	%rdx, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%rdi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$552, %rsp
	.cfi_def_cfa_offset 608
	movq	16(%rsi), %rcx
	movq	%fs:40, %rax
	movq	%rax, 536(%rsp)
	movq	8(%rsi), %rax
	movq	%rax, %rdx
	imulq	%rcx, %rdx
	testq	%rdx, %rdx
	je	.L1130
	movl	228(%r12), %edx
	movq	%rsi, %rbx
	cmpl	$-1, %edx
	je	.L1075
	cmpl	$-2, %edx
	je	.L1076
	movq	$0, 104(%rsp)
	movslq	%edx, %rsi
	movq	%rsi, 96(%rsp)
	testq	%rsi, %rsi
	jne	.L1027
.L1026:
	testb	$1, 232(%r12)
	jne	.L1078
.L1133:
	testq	%rcx, %rcx
	jle	.L1078
	leaq	64+_ZTVNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movq	16+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %r13
	xorl	%r14d, %r14d
	movq	$0, 40(%rsp)
	movq	%rcx, %xmm0
	movq	32+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %r15
	movdqa	%xmm0, %xmm3
	movdqa	%xmm0, %xmm4
	movhps	.LC48(%rip), %xmm3
	movhps	.LC49(%rip), %xmm4
	movaps	%xmm3, 64(%rsp)
	movaps	%xmm4, 48(%rsp)
	.p2align 4,,10
	.p2align 3
.L1045:
	testq	%rax, %rax
	jle	.L1028
	leaq	144(%rsp), %rax
	movq	$0, 8(%rsp)
	movq	%rax, 80(%rsp)
	leaq	272(%rsp), %rax
	movq	%rax, (%rsp)
	jmp	.L1044
	.p2align 4,,10
	.p2align 3
.L1132:
	movq	192(%rsp), %r8
	testq	%r8, %r8
	je	.L1093
	cmpq	%rax, %r8
	jb	.L1093
.L1035:
	movq	200(%rsp), %rcx
	xorl	%edx, %edx
	xorl	%esi, %esi
	subq	%rcx, %r8
.LEHB32:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE10_M_replaceEmmPKcm@PLT
.LEHE32:
.L1037:
	movq	120(%rsp), %rax
	movq	112(%rsp), %rdi
	cmpq	%rax, %r14
	cmovl	%rax, %r14
	movq	16(%rsp), %rax
	cmpq	%rax, %rdi
	je	.L1039
	movq	128(%rsp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1039:
	leaq	24+_ZTVNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movdqa	64(%rsp), %xmm2
	movq	240(%rsp), %rdi
	movq	%rax, 144(%rsp)
	addq	$80, %rax
	movq	%rax, 272(%rsp)
	movq	32(%rsp), %rax
	movaps	%xmm2, 160(%rsp)
	cmpq	%rax, %rdi
	je	.L1043
	movq	256(%rsp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1043:
	movq	24(%rsp), %rdi
	leaq	16+_ZTVSt15basic_streambufIcSt11char_traitsIcEE(%rip), %rax
	movq	%rax, 168(%rsp)
	call	_ZNSt6localeD1Ev@PLT
	movq	8+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	48+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movq	(%rsp), %rdi
	movq	-24(%rax), %rax
	movq	%rcx, 144(%rsp,%rax)
	movq	-24(%r15), %rax
	movq	40+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movq	%r15, 160(%rsp)
	movq	%rcx, 160(%rsp,%rax)
	movq	-24(%r13), %rax
	movq	24+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movq	%r13, 144(%rsp)
	movq	%rcx, 144(%rsp,%rax)
	leaq	16+_ZTVSt9basic_iosIcSt11char_traitsIcEE(%rip), %rax
	movq	%rax, 272(%rsp)
	movq	$0, 152(%rsp)
	call	_ZNSt8ios_baseD2Ev@PLT
	addq	$1, 8(%rsp)
	movq	8(%rbx), %rax
	movq	8(%rsp), %rcx
	cmpq	%rax, %rcx
	jge	.L1131
.L1044:
	movq	(%rsp), %rdi
	call	_ZNSt8ios_baseC2Ev@PLT
	leaq	16+_ZTVSt9basic_iosIcSt11char_traitsIcEE(%rip), %rax
	xorl	%edx, %edx
	xorl	%esi, %esi
	pxor	%xmm0, %xmm0
	movw	%dx, 496(%rsp)
	movq	24+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movups	%xmm0, 504(%rsp)
	movups	%xmm0, 520(%rsp)
	movq	%rax, 272(%rsp)
	movq	-24(%r13), %rax
	movq	$0, 488(%rsp)
	movq	%r13, 144(%rsp)
	movq	%rcx, 144(%rsp,%rax)
	movq	80(%rsp), %rax
	movq	$0, 152(%rsp)
	addq	-24(%r13), %rax
	movq	%rax, %rdi
.LEHB33:
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE4initEPSt15basic_streambufIcS1_E@PLT
.LEHE33:
	leaq	160(%rsp), %rax
	movq	%r15, 160(%rsp)
	xorl	%esi, %esi
	movq	%rax, 16(%rsp)
	addq	-24(%r15), %rax
	movq	%rax, %rdi
	movq	40+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	%rax, (%rdi)
.LEHB34:
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE4initEPSt15basic_streambufIcS1_E@PLT
.LEHE34:
	movq	8+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	48+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	pxor	%xmm0, %xmm0
	movdqa	48(%rsp), %xmm1
	movq	-24(%rax), %rax
	movq	%rcx, 144(%rsp,%rax)
	leaq	24+_ZTVNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	%rax, 144(%rsp)
	addq	$80, %rax
	movq	%rax, 272(%rsp)
	leaq	224(%rsp), %rax
	movq	%rax, %rdi
	movq	%rax, 24(%rsp)
	movaps	%xmm1, 160(%rsp)
	movaps	%xmm0, 176(%rsp)
	movaps	%xmm0, 192(%rsp)
	movaps	%xmm0, 208(%rsp)
	call	_ZNSt6localeC1Ev@PLT
	leaq	16+_ZTVNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	(%rsp), %rdi
	movl	$24, 232(%rsp)
	movq	%rax, 168(%rsp)
	leaq	256(%rsp), %rax
	movq	%rax, 32(%rsp)
	movq	%rax, 240(%rsp)
	leaq	168(%rsp), %rax
	movq	%rax, %rsi
	movb	$0, 256(%rsp)
	movq	$0, 248(%rsp)
	movq	%rax, 88(%rsp)
.LEHB35:
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE4initEPSt15basic_streambufIcS1_E@PLT
.LEHE35:
	movq	0(%rbp), %rax
	movq	(%rsp), %rdi
	movq	-24(%rax), %rsi
	addq	%rbp, %rsi
.LEHB36:
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE7copyfmtERKS2_@PLT
	movq	40(%rsp), %rax
	imulq	8(%rbx), %rax
	movq	8(%rsp), %rcx
	movq	(%rbx), %rdx
	movq	16(%rsp), %rdi
	addq	%rcx, %rax
	movsd	(%rdx,%rax,8), %xmm0
	call	_ZNSo9_M_insertIdEERSoT_@PLT
.LEHE36:
	leaq	128(%rsp), %rax
	movb	$0, 128(%rsp)
	leaq	112(%rsp), %rdi
	movq	%rax, 16(%rsp)
	movq	%rax, 112(%rsp)
	movq	208(%rsp), %rax
	movq	$0, 120(%rsp)
	testq	%rax, %rax
	jne	.L1132
	leaq	240(%rsp), %rsi
.LEHB37:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_assignERKS4_@PLT
.LEHE37:
	jmp	.L1037
.L1076:
	movq	$15, 96(%rsp)
.L1027:
	movq	0(%rbp), %rdx
	movq	96(%rsp), %rdi
	movq	-24(%rdx), %rsi
	movq	%rdi, 104(%rsp)
	addq	%rbp, %rsi
	movq	%rsi, %rdx
	movq	8(%rsi), %rsi
	movq	%rdi, 8(%rdx)
	movq	%rsi, 96(%rsp)
	testb	$1, 232(%r12)
	je	.L1133
.L1078:
	xorl	%r14d, %r14d
.L1028:
	movq	0(%rbp), %rax
	movq	-24(%rax), %r13
	addq	%rbp, %r13
	movq	16(%r13), %rax
	cmpb	$0, 225(%r13)
	movq	%rax, 32(%rsp)
	je	.L1048
	movzbl	224(%r13), %r15d
.L1049:
	movq	8(%r12), %rdx
	movq	(%r12), %rsi
	movq	%rbp, %rdi
	xorl	%r13d, %r13d
.LEHB38:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	cmpq	$0, 8(%rbx)
	jle	.L1070
	movb	%r15b, 24(%rsp)
	movq	%r12, %r15
	movq	%r13, %r12
	movq	%rbx, %r13
	.p2align 4,,10
	.p2align 3
.L1053:
	movq	72(%r15), %rdx
	movq	64(%r15), %rsi
	movq	%rbp, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	testq	%r14, %r14
	je	.L1056
	movq	0(%rbp), %rax
	movzbl	224(%r15), %edx
	movq	-24(%rax), %rbx
	addq	%rbp, %rbx
	cmpb	$0, 225(%rbx)
	movq	%rbx, %rax
	je	.L1134
.L1057:
	movb	%dl, 224(%rbx)
	movq	%r14, 16(%rax)
.L1056:
	movq	0(%r13), %rax
	movq	%rbp, %rdi
	movl	$1, %ebx
	movsd	(%rax,%r12,8), %xmm0
	call	_ZNSo9_M_insertIdEERSoT_@PLT
	cmpq	$1, 16(%r13)
	jg	.L1061
	jmp	.L1068
	.p2align 4,,10
	.p2align 3
.L1065:
	movb	%cl, 224(%rdx)
	movq	%r14, 16(%rax)
.L1064:
	movq	8(%r13), %rax
	movq	0(%r13), %rdx
	movq	%rbp, %rdi
	imulq	%rbx, %rax
	addq	$1, %rbx
	addq	%r12, %rax
	movsd	(%rdx,%rax,8), %xmm0
	call	_ZNSo9_M_insertIdEERSoT_@PLT
	cmpq	16(%r13), %rbx
	jge	.L1068
.L1061:
	movq	200(%r15), %rdx
	movq	192(%r15), %rsi
	movq	%rbp, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	testq	%r14, %r14
	je	.L1064
	movq	0(%rbp), %rax
	movzbl	224(%r15), %ecx
	movq	-24(%rax), %rdx
	addq	%rbp, %rdx
	cmpb	$0, 225(%rdx)
	movq	%rdx, %rax
	jne	.L1065
	movq	240(%rdx), %rdi
	testq	%rdi, %rdi
	je	.L1058
	cmpb	$0, 56(%rdi)
	jne	.L1066
	movb	%cl, 16(%rsp)
	movq	%rdx, 8(%rsp)
	movq	%rdi, (%rsp)
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	(%rsp), %rdi
	movq	8(%rsp), %rdx
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rsi
	movzbl	16(%rsp), %ecx
	movq	(%rdi), %rax
	movq	48(%rax), %rax
	cmpq	%rsi, %rax
	jne	.L1067
	movq	0(%rbp), %rax
	movq	-24(%rax), %rsi
	addq	%rbp, %rsi
	movq	%rsi, %rax
.L1066:
	movb	$1, 225(%rdx)
	jmp	.L1065
	.p2align 4,,10
	.p2align 3
.L1093:
	movq	%rax, %r8
	jmp	.L1035
	.p2align 4,,10
	.p2align 3
.L1131:
	addq	$1, 40(%rsp)
	movq	40(%rsp), %rcx
	cmpq	16(%rbx), %rcx
	jl	.L1045
	jmp	.L1028
	.p2align 4,,10
	.p2align 3
.L1068:
	movq	104(%r15), %rdx
	movq	96(%r15), %rsi
	movq	%rbp, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	8(%r13), %rax
	leaq	-1(%rax), %rdx
	cmpq	%r12, %rdx
	jg	.L1135
	addq	$1, %r12
	cmpq	%rax, %r12
	jge	.L1126
.L1069:
	movq	168(%r15), %rdx
	movq	160(%r15), %rsi
	movq	%rbp, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	jmp	.L1053
	.p2align 4,,10
	.p2align 3
.L1135:
	movq	136(%r15), %rdx
	movq	128(%r15), %rsi
	movq	%rbp, %rdi
	addq	$1, %r12
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	cmpq	8(%r13), %r12
	jl	.L1069
.L1126:
	movq	%r15, %r12
	movzbl	24(%rsp), %r15d
.L1070:
	movq	40(%r12), %rdx
	movq	32(%r12), %rsi
	movq	%rbp, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	cmpq	$0, 104(%rsp)
	je	.L1055
	movq	0(%rbp), %rax
	movq	96(%rsp), %rbx
	movq	-24(%rax), %rax
	movq	%rbx, 8(%rbp,%rax)
.L1055:
	testq	%r14, %r14
	je	.L1025
	movq	0(%rbp), %rax
	movq	-24(%rax), %rbx
	addq	%rbp, %rbx
	cmpb	$0, 225(%rbx)
	movq	%rbx, %rax
	je	.L1136
.L1071:
	movb	%r15b, 224(%rbx)
	movq	32(%rsp), %rbx
	movq	%rbx, 16(%rax)
.L1025:
	movq	536(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L1137
	addq	$552, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	movq	%rbp, %rax
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L1134:
	.cfi_restore_state
	movq	240(%rbx), %rdi
	testq	%rdi, %rdi
	je	.L1058
	cmpb	$0, 56(%rdi)
	jne	.L1059
	movb	%dl, 8(%rsp)
	movq	%rdi, (%rsp)
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	(%rsp), %rdi
	movzbl	8(%rsp), %edx
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rsi
	movq	(%rdi), %rax
	movq	48(%rax), %rax
	cmpq	%rsi, %rax
	jne	.L1060
	movq	0(%rbp), %rax
	movq	-24(%rax), %rcx
	addq	%rbp, %rcx
	movq	%rcx, %rax
.L1059:
	movb	$1, 225(%rbx)
	jmp	.L1057
	.p2align 4,,10
	.p2align 3
.L1067:
	movb	%cl, 8(%rsp)
	movl	$32, %esi
	movq	%rdx, (%rsp)
	call	*%rax
	movq	0(%rbp), %rax
	movzbl	8(%rsp), %ecx
	movq	(%rsp), %rdx
	movq	-24(%rax), %rsi
	addq	%rbp, %rsi
	movq	%rsi, %rax
	jmp	.L1066
.L1075:
	movq	$0, 104(%rsp)
	movq	$0, 96(%rsp)
	jmp	.L1026
.L1048:
	movq	240(%r13), %rdi
	testq	%rdi, %rdi
	je	.L1058
	cmpb	$0, 56(%rdi)
	je	.L1051
	movzbl	89(%rdi), %r15d
.L1052:
	movb	%r15b, 224(%r13)
	movb	$1, 225(%r13)
	jmp	.L1049
.L1130:
	movq	8(%r12), %rdx
	movq	(%r12), %rsi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	40(%r12), %rdx
	movq	32(%r12), %rsi
	movq	%rax, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	jmp	.L1025
.L1051:
	movq	%rdi, (%rsp)
	movl	$32, %r15d
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	(%rsp), %rdi
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rdx
	movq	(%rdi), %rax
	movq	48(%rax), %rax
	cmpq	%rdx, %rax
	je	.L1052
	movl	$32, %esi
	call	*%rax
	movl	%eax, %r15d
	jmp	.L1052
.L1136:
	movq	240(%rbx), %r12
	testq	%r12, %r12
	je	.L1058
	cmpb	$0, 56(%r12)
	jne	.L1072
	movq	%r12, %rdi
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	(%r12), %rax
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rdx
	movq	48(%rax), %rax
	cmpq	%rdx, %rax
	jne	.L1073
.L1129:
	movq	0(%rbp), %rax
	movq	-24(%rax), %rcx
	addq	%rbp, %rcx
	movq	%rcx, %rax
.L1072:
	movb	$1, 225(%rbx)
	jmp	.L1071
.L1060:
	movb	%dl, (%rsp)
	movl	$32, %esi
	call	*%rax
	movq	0(%rbp), %rax
	movzbl	(%rsp), %edx
	movq	-24(%rax), %rcx
	addq	%rbp, %rcx
	movq	%rcx, %rax
	jmp	.L1059
.L1073:
	movl	$32, %esi
	movq	%r12, %rdi
	call	*%rax
	jmp	.L1129
.L1058:
	call	_ZSt16__throw_bad_castv@PLT
.L1137:
	call	__stack_chk_fail@PLT
.L1092:
	movq	%rax, %rbx
	jmp	.L1040
.L1089:
	movq	%rax, %rbx
	jmp	.L1032
.L1090:
	movq	%rax, %rbx
	jmp	.L1033
.L1091:
	movq	%rax, %rbx
	jmp	.L1031
.L1088:
	movq	%rax, %rbx
	jmp	.L1042
.L1040:
	movq	112(%rsp), %rdi
	movq	16(%rsp), %rax
	cmpq	%rax, %rdi
	je	.L1042
	movq	128(%rsp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1042:
	movq	80(%rsp), %rdi
	call	_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEED1Ev@PLT
	movq	%rbx, %rdi
	call	_Unwind_Resume@PLT
.L1033:
	movq	88(%rsp), %rdi
	call	_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED1Ev
	movq	8+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	48+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movq	-24(%rax), %rax
	movq	%rcx, 144(%rsp,%rax)
	movq	-24(%r15), %rax
	movq	40+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movq	%r15, 160(%rsp)
	movq	%rcx, 160(%rsp,%rax)
.L1128:
	movq	-24(%r13), %rax
	movq	24+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movq	%r13, 144(%rsp)
	movq	%rcx, 144(%rsp,%rax)
	xorl	%eax, %eax
	movq	%rax, 152(%rsp)
.L1032:
	movq	(%rsp), %rdi
	leaq	16+_ZTVSt9basic_iosIcSt11char_traitsIcEE(%rip), %rax
	movq	%rax, 272(%rsp)
	call	_ZNSt8ios_baseD2Ev@PLT
	movq	%rbx, %rdi
	call	_Unwind_Resume@PLT
.LEHE38:
.L1031:
	jmp	.L1128
	.cfi_endproc
.LFE10999:
	.section	.gcc_except_table
.LLSDA10999:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE10999-.LLSDACSB10999
.LLSDACSB10999:
	.uleb128 .LEHB32-.LFB10999
	.uleb128 .LEHE32-.LEHB32
	.uleb128 .L1092-.LFB10999
	.uleb128 0
	.uleb128 .LEHB33-.LFB10999
	.uleb128 .LEHE33-.LEHB33
	.uleb128 .L1089-.LFB10999
	.uleb128 0
	.uleb128 .LEHB34-.LFB10999
	.uleb128 .LEHE34-.LEHB34
	.uleb128 .L1091-.LFB10999
	.uleb128 0
	.uleb128 .LEHB35-.LFB10999
	.uleb128 .LEHE35-.LEHB35
	.uleb128 .L1090-.LFB10999
	.uleb128 0
	.uleb128 .LEHB36-.LFB10999
	.uleb128 .LEHE36-.LEHB36
	.uleb128 .L1088-.LFB10999
	.uleb128 0
	.uleb128 .LEHB37-.LFB10999
	.uleb128 .LEHE37-.LEHB37
	.uleb128 .L1092-.LFB10999
	.uleb128 0
	.uleb128 .LEHB38-.LFB10999
	.uleb128 .LEHE38-.LEHB38
	.uleb128 0
	.uleb128 0
.LLSDACSE10999:
	.section	.text._ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE,"axG",@progbits,_ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE,comdat
	.size	_ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE, .-_ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE
	.section	.text._ZN5EigenlsINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS3_RKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5EigenlsINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS3_RKNS_9DenseBaseIT_EE,comdat
	.p2align 4
	.weak	_ZN5EigenlsINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS3_RKNS_9DenseBaseIT_EE
	.type	_ZN5EigenlsINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS3_RKNS_9DenseBaseIT_EE, @function
_ZN5EigenlsINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS3_RKNS_9DenseBaseIT_EE:
.LFB10607:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA10607
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	leaq	-368(%rbp), %rdx
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	leaq	-448(%rbp), %r14
	pushq	%r13
	leaq	-480(%rbp), %r15
	.cfi_offset 13, -40
	leaq	-416(%rbp), %r13
	pushq	%r12
	.cfi_offset 12, -48
	leaq	-384(%rbp), %r12
	leaq	-304(%rbp), %r10
	pushq	%rbx
	.cfi_offset 3, -56
	leaq	-320(%rbp), %rbx
	leaq	-432(%rbp), %r9
	leaq	-464(%rbp), %r8
	subq	$488, %rsp
	movq	%rdi, -520(%rbp)
	movzwl	.LC50(%rip), %ecx
	movq	%rsi, -528(%rbp)
	movzwl	.LC51(%rip), %edi
	leaq	-336(%rbp), %rsi
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	xorl	%eax, %eax
	leaq	-352(%rbp), %rax
	movw	%cx, -448(%rbp)
	leaq	-496(%rbp), %rcx
	movq	%rax, -512(%rbp)
	movq	%rax, -368(%rbp)
	leaq	-400(%rbp), %rax
	movw	%di, -480(%rbp)
	movq	%r10, %rdi
	movq	%rbx, -336(%rbp)
	movq	$0, -328(%rbp)
	movb	$0, -320(%rbp)
	movq	$0, -360(%rbp)
	movb	$0, -352(%rbp)
	movq	%r12, -400(%rbp)
	movq	$0, -392(%rbp)
	movb	$0, -384(%rbp)
	movq	%r13, -432(%rbp)
	movq	$0, -424(%rbp)
	movb	$0, -416(%rbp)
	movq	%r14, -464(%rbp)
	movq	$1, -456(%rbp)
	movq	%r15, -496(%rbp)
	movq	$1, -488(%rbp)
	pushq	$32
	pushq	%rsi
	movl	$-1, %esi
	pushq	%rdx
	xorl	%edx, %edx
	pushq	%rax
	movq	%r10, -504(%rbp)
.LEHB39:
	.cfi_escape 0x2e,0x20
	call	_ZN5Eigen8IOFormatC1EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c
.LEHE39:
	movq	-504(%rbp), %rdx
	movq	-528(%rbp), %rsi
	addq	$32, %rsp
	movq	-520(%rbp), %rdi
.LEHB40:
	.cfi_escape 0x2e,0
	call	_ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE
.LEHE40:
	movq	-504(%rbp), %rdi
	movq	%rax, -520(%rbp)
	call	_ZN5Eigen8IOFormatD1Ev
	movq	-496(%rbp), %rdi
	cmpq	%r15, %rdi
	je	.L1139
	movq	-480(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1139:
	movq	-464(%rbp), %rdi
	cmpq	%r14, %rdi
	je	.L1140
	movq	-448(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1140:
	movq	-432(%rbp), %rdi
	cmpq	%r13, %rdi
	je	.L1141
	movq	-416(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1141:
	movq	-400(%rbp), %rdi
	cmpq	%r12, %rdi
	je	.L1142
	movq	-384(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1142:
	movq	-368(%rbp), %rdi
	movq	-512(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L1143
	movq	-352(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1143:
	movq	-336(%rbp), %rdi
	cmpq	%rbx, %rdi
	je	.L1138
	movq	-320(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1138:
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L1157
	movq	-520(%rbp), %rax
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
.L1157:
	.cfi_restore_state
	call	__stack_chk_fail@PLT
.L1155:
	movq	%rax, -520(%rbp)
	jmp	.L1145
.L1154:
	movq	%rax, -504(%rbp)
	jmp	.L1146
.L1145:
	movq	-504(%rbp), %rdi
	call	_ZN5Eigen8IOFormatD1Ev
	movq	-520(%rbp), %rax
	movq	%rax, -504(%rbp)
.L1146:
	movq	-496(%rbp), %rdi
	cmpq	%r15, %rdi
	je	.L1147
	movq	-480(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1147:
	movq	-464(%rbp), %rdi
	cmpq	%r14, %rdi
	je	.L1148
	movq	-448(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1148:
	movq	-432(%rbp), %rdi
	cmpq	%r13, %rdi
	je	.L1149
	movq	-416(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1149:
	movq	-400(%rbp), %rdi
	cmpq	%r12, %rdi
	je	.L1150
	movq	-384(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1150:
	movq	-368(%rbp), %rdi
	movq	-512(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L1151
	movq	-352(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1151:
	movq	-336(%rbp), %rdi
	cmpq	%rbx, %rdi
	je	.L1152
	movq	-320(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1152:
	movq	-504(%rbp), %rdi
.LEHB41:
	call	_Unwind_Resume@PLT
.LEHE41:
	.cfi_endproc
.LFE10607:
	.section	.gcc_except_table
.LLSDA10607:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE10607-.LLSDACSB10607
.LLSDACSB10607:
	.uleb128 .LEHB39-.LFB10607
	.uleb128 .LEHE39-.LEHB39
	.uleb128 .L1154-.LFB10607
	.uleb128 0
	.uleb128 .LEHB40-.LFB10607
	.uleb128 .LEHE40-.LEHB40
	.uleb128 .L1155-.LFB10607
	.uleb128 0
	.uleb128 .LEHB41-.LFB10607
	.uleb128 .LEHE41-.LEHB41
	.uleb128 0
	.uleb128 0
.LLSDACSE10607:
	.section	.text._ZN5EigenlsINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS3_RKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5EigenlsINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS3_RKNS_9DenseBaseIT_EE,comdat
	.size	_ZN5EigenlsINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS3_RKNS_9DenseBaseIT_EE, .-_ZN5EigenlsINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS3_RKNS_9DenseBaseIT_EE
	.section	.text._ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE,comdat
	.p2align 4
	.weak	_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE
	.type	_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE, @function
_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE:
.LFB10586:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA10586
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	leaq	-384(%rbp), %rcx
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	leaq	-320(%rbp), %r15
	pushq	%r13
	leaq	-448(%rbp), %r14
	.cfi_offset 13, -40
	leaq	-352(%rbp), %r13
	pushq	%r12
	leaq	-368(%rbp), %rdx
	leaq	-432(%rbp), %r9
	pushq	%rbx
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	%rsi, %rbx
	leaq	-336(%rbp), %rsi
	leaq	-464(%rbp), %r8
	subq	$536, %rsp
	movq	%rdi, -576(%rbp)
	leaq	-480(%rbp), %rdi
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	xorl	%eax, %eax
	movq	%rcx, -552(%rbp)
	leaq	-400(%rbp), %rax
	movq	%rcx, -400(%rbp)
	leaq	-416(%rbp), %rcx
	movq	%rcx, -544(%rbp)
	movq	%rcx, -432(%rbp)
	movzwl	.LC50(%rip), %ecx
	movq	%rdi, -560(%rbp)
	movq	%rdi, -496(%rbp)
	movzwl	.LC51(%rip), %edi
	movw	%cx, -448(%rbp)
	leaq	-496(%rbp), %rcx
	movq	%r15, -336(%rbp)
	movq	$0, -328(%rbp)
	movb	$0, -320(%rbp)
	movq	%r13, -368(%rbp)
	movq	$0, -360(%rbp)
	movb	$0, -352(%rbp)
	movq	$0, -392(%rbp)
	movb	$0, -384(%rbp)
	movq	$0, -424(%rbp)
	movb	$0, -416(%rbp)
	movq	%r14, -464(%rbp)
	movq	$1, -456(%rbp)
	movq	$1, -488(%rbp)
	movw	%di, -480(%rbp)
	leaq	-304(%rbp), %rdi
	movq	%rdi, -568(%rbp)
	pushq	$32
	pushq	%rsi
	movl	$-1, %esi
	pushq	%rdx
	xorl	%edx, %edx
	pushq	%rax
.LEHB42:
	.cfi_escape 0x2e,0x20
	call	_ZN5Eigen8IOFormatC1EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c
.LEHE42:
	movq	(%rbx), %rax
	addq	$32, %rsp
	pxor	%xmm0, %xmm0
	movq	$0, -528(%rbp)
	movups	%xmm0, -520(%rbp)
	movq	8(%rax), %r12
	movq	(%rax), %rbx
	movq	16(%rax), %rax
	movq	%rax, -536(%rbp)
	orq	%r12, %rax
	je	.L1191
	js	.L1211
	testq	%r12, %r12
	je	.L1161
	movq	-536(%rbp), %rsi
	testq	%rsi, %rsi
	je	.L1161
	movabsq	$9223372036854775807, %rax
	cqto
	idivq	%rsi
	cmpq	%rax, %r12
	jg	.L1212
	movq	%rsi, %rax
	imulq	%r12, %rax
.L1163:
	movabsq	$2305843009213693951, %rdx
	cmpq	%rdx, %rax
	jg	.L1213
	leaq	0(,%rax,8), %rdi
.LEHB43:
	.cfi_escape 0x2e,0
	call	_ZN5Eigen8internal14aligned_mallocEm
.LEHE43:
	movq	-536(%rbp), %rsi
	movq	%r12, %xmm0
	movq	%rax, -528(%rbp)
	movq	%rax, -536(%rbp)
	movq	%rsi, %xmm1
	imulq	%r12, %rsi
	punpcklqdq	%xmm1, %xmm0
	movups	%xmm0, -520(%rbp)
	movq	%rsi, %r11
	testq	%rsi, %rsi
	jle	.L1159
	subq	$1, %rsi
	leaq	8(%rbx), %rdx
	cmpq	$3, %rsi
	jbe	.L1214
	leaq	0(,%r11,8), %rcx
	leaq	(%rax,%rcx), %rdi
	cmpq	%rdi, %rdx
	jnb	.L1195
	movq	%r11, %rdi
	salq	$4, %rdi
	addq	%rbx, %rdi
	cmpq	%rdi, %rax
	jnb	.L1195
.L1168:
	addq	%rax, %rcx
	.p2align 4,,10
	.p2align 3
.L1174:
	movsd	(%rdx), %xmm0
	addq	$8, %rax
	addq	$16, %rdx
	movsd	%xmm0, -8(%rax)
	cmpq	%rax, %rcx
	jne	.L1174
.L1159:
	movq	-568(%rbp), %r12
	movq	-576(%rbp), %rdi
	leaq	-528(%rbp), %rsi
	movq	%r12, %rdx
.LEHB44:
	call	_ZN5Eigen8internal12print_matrixINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS4_RKT_RKNS_8IOFormatE
.LEHE44:
	movq	-536(%rbp), %rdi
	movq	%rax, %rbx
	call	free@PLT
	movq	%r12, %rdi
	call	_ZN5Eigen8IOFormatD1Ev
	movq	-496(%rbp), %rdi
	movq	-560(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L1175
	movq	-480(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1175:
	movq	-464(%rbp), %rdi
	cmpq	%r14, %rdi
	je	.L1176
	movq	-448(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1176:
	movq	-432(%rbp), %rdi
	movq	-544(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L1177
	movq	-416(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1177:
	movq	-400(%rbp), %rdi
	movq	-552(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L1178
	movq	-384(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1178:
	movq	-368(%rbp), %rdi
	cmpq	%r13, %rdi
	je	.L1179
	movq	-352(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1179:
	movq	-336(%rbp), %rdi
	cmpq	%r15, %rdi
	je	.L1158
	movq	-320(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1158:
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L1215
	leaq	-40(%rbp), %rsp
	movq	%rbx, %rax
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L1161:
	.cfi_restore_state
	movq	-536(%rbp), %rax
	imulq	%r12, %rax
	testq	%rax, %rax
	jne	.L1163
	movq	-528(%rbp), %rax
	movq	%r12, %xmm0
	movhps	-536(%rbp), %xmm0
	movups	%xmm0, -520(%rbp)
	movq	%rax, -536(%rbp)
	jmp	.L1159
	.p2align 4,,10
	.p2align 3
.L1195:
	movq	%rsi, %rcx
	xorl	%edx, %edx
	movq	%r11, %rdi
	shrq	%rcx
	salq	$4, %rcx
	.p2align 4,,10
	.p2align 3
.L1170:
	movupd	8(%rbx,%rdx,2), %xmm0
	movhpd	24(%rbx,%rdx,2), %xmm0
	movups	%xmm0, (%rax,%rdx)
	addq	$16, %rdx
	cmpq	%rcx, %rdx
	jne	.L1170
	andq	$-2, %rsi
	movq	%rsi, %rdx
	salq	$4, %rdx
	movsd	8(%rbx,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rsi,8)
	addq	$1, %rsi
	cmpq	%rsi, %rdi
	jle	.L1159
	movq	%rsi, %rdx
	salq	$4, %rdx
	movsd	8(%rbx,%rdx), %xmm0
	movsd	%xmm0, (%rax,%rsi,8)
	jmp	.L1159
	.p2align 4,,10
	.p2align 3
.L1214:
	leaq	0(,%r11,8), %rcx
	jmp	.L1168
	.p2align 4,,10
	.p2align 3
.L1191:
	movq	$0, -536(%rbp)
	jmp	.L1159
.L1211:
	leaq	.LC44(%rip), %rcx
	movl	$273, %edx
	leaq	.LC37(%rip), %rsi
	leaq	.LC38(%rip), %rdi
	call	__assert_fail@PLT
.L1215:
	call	__stack_chk_fail@PLT
.L1213:
.LEHB45:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE45:
.L1193:
	movq	%rax, %rbx
	jmp	.L1182
.L1192:
	movq	%rax, %rbx
	jmp	.L1183
.L1181:
	movq	-536(%rbp), %rdi
	call	free@PLT
.L1182:
	movq	-568(%rbp), %rdi
	call	_ZN5Eigen8IOFormatD1Ev
.L1183:
	movq	-496(%rbp), %rdi
	movq	-560(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L1184
	movq	-480(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1184:
	movq	-464(%rbp), %rdi
	cmpq	%r14, %rdi
	je	.L1185
	movq	-448(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1185:
	movq	-432(%rbp), %rdi
	movq	-544(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L1186
	movq	-416(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1186:
	movq	-400(%rbp), %rdi
	movq	-552(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L1187
	movq	-384(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1187:
	movq	-368(%rbp), %rdi
	cmpq	%r13, %rdi
	je	.L1188
	movq	-352(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1188:
	movq	-336(%rbp), %rdi
	cmpq	%r15, %rdi
	je	.L1189
	movq	-320(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1189:
	movq	%rbx, %rdi
.LEHB46:
	call	_Unwind_Resume@PLT
.LEHE46:
.L1212:
.LEHB47:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE47:
.L1194:
	movq	%rax, %rbx
	jmp	.L1181
	.cfi_endproc
.LFE10586:
	.section	.gcc_except_table
.LLSDA10586:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE10586-.LLSDACSB10586
.LLSDACSB10586:
	.uleb128 .LEHB42-.LFB10586
	.uleb128 .LEHE42-.LEHB42
	.uleb128 .L1192-.LFB10586
	.uleb128 0
	.uleb128 .LEHB43-.LFB10586
	.uleb128 .LEHE43-.LEHB43
	.uleb128 .L1193-.LFB10586
	.uleb128 0
	.uleb128 .LEHB44-.LFB10586
	.uleb128 .LEHE44-.LEHB44
	.uleb128 .L1194-.LFB10586
	.uleb128 0
	.uleb128 .LEHB45-.LFB10586
	.uleb128 .LEHE45-.LEHB45
	.uleb128 .L1193-.LFB10586
	.uleb128 0
	.uleb128 .LEHB46-.LFB10586
	.uleb128 .LEHE46-.LEHB46
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB47-.LFB10586
	.uleb128 .LEHE47-.LEHB47
	.uleb128 .L1193-.LFB10586
	.uleb128 0
.LLSDACSE10586:
	.section	.text._ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE,comdat
	.size	_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE, .-_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE
	.section	.text._ZN5Eigen8internal12print_matrixINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS6_RKT_RKNS_8IOFormatE,"axG",@progbits,_ZN5Eigen8internal12print_matrixINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS6_RKT_RKNS_8IOFormatE,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal12print_matrixINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS6_RKT_RKNS_8IOFormatE
	.type	_ZN5Eigen8internal12print_matrixINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS6_RKT_RKNS_8IOFormatE, @function
_ZN5Eigen8internal12print_matrixINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS6_RKT_RKNS_8IOFormatE:
.LFB11184:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA11184
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	movq	%rdx, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%rdi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$552, %rsp
	.cfi_def_cfa_offset 608
	movq	16(%rsi), %rcx
	movq	%fs:40, %rax
	movq	%rax, 536(%rsp)
	movq	8(%rsi), %rax
	movq	%rax, %rdx
	imulq	%rcx, %rdx
	testq	%rdx, %rdx
	je	.L1323
	movslq	228(%r12), %rdx
	movq	%rsi, %rbx
	cmpl	$-1, %edx
	je	.L1268
	cmpl	$-2, %edx
	je	.L1269
	movq	%rdx, 96(%rsp)
	movq	$0, 104(%rsp)
	testq	%rdx, %rdx
	jne	.L1220
.L1219:
	testb	$1, 232(%r12)
	jne	.L1271
.L1326:
	testq	%rcx, %rcx
	jle	.L1271
	leaq	64+_ZTVNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rdx
	movq	16+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %r13
	xorl	%r14d, %r14d
	movq	$0, 40(%rsp)
	movq	%rdx, %xmm0
	movq	32+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %r15
	movdqa	%xmm0, %xmm3
	movdqa	%xmm0, %xmm4
	movhps	.LC48(%rip), %xmm3
	movhps	.LC49(%rip), %xmm4
	movaps	%xmm3, 48(%rsp)
	movaps	%xmm4, 64(%rsp)
	.p2align 4,,10
	.p2align 3
.L1238:
	testq	%rax, %rax
	jle	.L1221
	leaq	144(%rsp), %rax
	movq	$0, 8(%rsp)
	movq	%rax, 80(%rsp)
	leaq	272(%rsp), %rax
	movq	%rax, (%rsp)
	jmp	.L1237
	.p2align 4,,10
	.p2align 3
.L1325:
	movq	192(%rsp), %r8
	testq	%r8, %r8
	je	.L1286
	cmpq	%rax, %r8
	jb	.L1286
.L1228:
	movq	200(%rsp), %rcx
	xorl	%edx, %edx
	xorl	%esi, %esi
	subq	%rcx, %r8
.LEHB48:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE10_M_replaceEmmPKcm@PLT
.LEHE48:
.L1230:
	movq	120(%rsp), %rax
	movq	112(%rsp), %rdi
	cmpq	%rax, %r14
	cmovl	%rax, %r14
	movq	16(%rsp), %rax
	cmpq	%rax, %rdi
	je	.L1232
	movq	128(%rsp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1232:
	leaq	24+_ZTVNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movdqa	48(%rsp), %xmm2
	movq	240(%rsp), %rdi
	movq	%rax, 144(%rsp)
	addq	$80, %rax
	movq	%rax, 272(%rsp)
	movq	32(%rsp), %rax
	movaps	%xmm2, 160(%rsp)
	cmpq	%rax, %rdi
	je	.L1236
	movq	256(%rsp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1236:
	movq	24(%rsp), %rdi
	leaq	16+_ZTVSt15basic_streambufIcSt11char_traitsIcEE(%rip), %rax
	movq	%rax, 168(%rsp)
	call	_ZNSt6localeD1Ev@PLT
	movq	8+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	48+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movq	(%rsp), %rdi
	movq	-24(%rax), %rax
	movq	%rcx, 144(%rsp,%rax)
	movq	-24(%r15), %rax
	movq	40+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movq	%r15, 160(%rsp)
	movq	%rcx, 160(%rsp,%rax)
	movq	-24(%r13), %rax
	movq	24+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movq	%r13, 144(%rsp)
	movq	%rcx, 144(%rsp,%rax)
	leaq	16+_ZTVSt9basic_iosIcSt11char_traitsIcEE(%rip), %rax
	movq	%rax, 272(%rsp)
	movq	$0, 152(%rsp)
	call	_ZNSt8ios_baseD2Ev@PLT
	addq	$1, 8(%rsp)
	movq	8(%rbx), %rax
	movq	8(%rsp), %rcx
	cmpq	%rax, %rcx
	jge	.L1324
.L1237:
	movq	(%rsp), %rdi
	call	_ZNSt8ios_baseC2Ev@PLT
	leaq	16+_ZTVSt9basic_iosIcSt11char_traitsIcEE(%rip), %rax
	xorl	%edx, %edx
	xorl	%esi, %esi
	pxor	%xmm0, %xmm0
	movw	%dx, 496(%rsp)
	movq	24+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	movups	%xmm0, 504(%rsp)
	movups	%xmm0, 520(%rsp)
	movq	%rax, 272(%rsp)
	movq	-24(%r13), %rax
	movq	$0, 488(%rsp)
	movq	%r13, 144(%rsp)
	movq	%rcx, 144(%rsp,%rax)
	movq	80(%rsp), %rax
	movq	$0, 152(%rsp)
	addq	-24(%r13), %rax
	movq	%rax, %rdi
.LEHB49:
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE4initEPSt15basic_streambufIcS1_E@PLT
.LEHE49:
	leaq	160(%rsp), %rax
	movq	%r15, 160(%rsp)
	xorl	%esi, %esi
	movq	%rax, 16(%rsp)
	addq	-24(%r15), %rax
	movq	%rax, %rdi
	movq	40+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	%rax, (%rdi)
.LEHB50:
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE4initEPSt15basic_streambufIcS1_E@PLT
.LEHE50:
	movq	8+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	48+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rcx
	pxor	%xmm0, %xmm0
	movdqa	64(%rsp), %xmm1
	movq	-24(%rax), %rax
	movq	%rcx, 144(%rsp,%rax)
	leaq	24+_ZTVNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	%rax, 144(%rsp)
	addq	$80, %rax
	movq	%rax, 272(%rsp)
	leaq	224(%rsp), %rax
	movq	%rax, %rdi
	movq	%rax, 24(%rsp)
	movaps	%xmm1, 160(%rsp)
	movaps	%xmm0, 176(%rsp)
	movaps	%xmm0, 192(%rsp)
	movaps	%xmm0, 208(%rsp)
	call	_ZNSt6localeC1Ev@PLT
	leaq	16+_ZTVNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	(%rsp), %rdi
	movl	$24, 232(%rsp)
	movq	%rax, 168(%rsp)
	leaq	256(%rsp), %rax
	movq	%rax, 32(%rsp)
	movq	%rax, 240(%rsp)
	leaq	168(%rsp), %rax
	movq	%rax, %rsi
	movb	$0, 256(%rsp)
	movq	$0, 248(%rsp)
	movq	%rax, 88(%rsp)
.LEHB51:
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE4initEPSt15basic_streambufIcS1_E@PLT
.LEHE51:
	movq	0(%rbp), %rax
	movq	(%rsp), %rdi
	movq	-24(%rax), %rsi
	addq	%rbp, %rsi
.LEHB52:
	call	_ZNSt9basic_iosIcSt11char_traitsIcEE7copyfmtERKS2_@PLT
	movq	40(%rsp), %rsi
	imulq	8(%rbx), %rsi
	movq	8(%rsp), %rax
	movq	16(%rsp), %rdi
	addq	%rax, %rsi
	salq	$4, %rsi
	addq	(%rbx), %rsi
	call	_ZStlsIdcSt11char_traitsIcEERSt13basic_ostreamIT0_T1_ES6_RKSt7complexIT_E@PLT
.LEHE52:
	leaq	128(%rsp), %rax
	movb	$0, 128(%rsp)
	leaq	112(%rsp), %rdi
	movq	%rax, 16(%rsp)
	movq	%rax, 112(%rsp)
	movq	208(%rsp), %rax
	movq	$0, 120(%rsp)
	testq	%rax, %rax
	jne	.L1325
	leaq	240(%rsp), %rsi
.LEHB53:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE9_M_assignERKS4_@PLT
.LEHE53:
	jmp	.L1230
.L1269:
	movq	$15, 96(%rsp)
.L1220:
	movq	0(%rbp), %rdx
	movq	-24(%rdx), %rdi
	addq	%rbp, %rdi
	movq	8(%rdi), %rsi
	movq	%rdi, %rdx
	movq	96(%rsp), %rdi
	movq	%rsi, 96(%rsp)
	movq	%rdi, 8(%rdx)
	movq	%rdi, 104(%rsp)
	testb	$1, 232(%r12)
	je	.L1326
.L1271:
	xorl	%r14d, %r14d
.L1221:
	movq	0(%rbp), %rax
	movq	-24(%rax), %r13
	addq	%rbp, %r13
	movq	16(%r13), %rax
	cmpb	$0, 225(%r13)
	movq	%rax, 32(%rsp)
	je	.L1241
	movzbl	224(%r13), %r15d
.L1242:
	movq	8(%r12), %rdx
	movq	(%r12), %rsi
	movq	%rbp, %rdi
	xorl	%r13d, %r13d
.LEHB54:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	cmpq	$0, 8(%rbx)
	jle	.L1263
	movb	%r15b, 24(%rsp)
	movq	%r12, %r15
	.p2align 4,,10
	.p2align 3
.L1246:
	movq	72(%r15), %rdx
	movq	64(%r15), %rsi
	movq	%rbp, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	testq	%r14, %r14
	je	.L1249
	movq	0(%rbp), %rax
	movzbl	224(%r15), %edx
	movq	-24(%rax), %r12
	addq	%rbp, %r12
	cmpb	$0, 225(%r12)
	movq	%r12, %rax
	je	.L1327
.L1250:
	movb	%dl, 224(%r12)
	movq	%r14, 16(%rax)
.L1249:
	movq	%r13, %rsi
	movq	%rbp, %rdi
	movl	$1, %r12d
	salq	$4, %rsi
	addq	(%rbx), %rsi
	call	_ZStlsIdcSt11char_traitsIcEERSt13basic_ostreamIT0_T1_ES6_RKSt7complexIT_E@PLT
	cmpq	$1, 16(%rbx)
	jg	.L1254
	jmp	.L1261
	.p2align 4,,10
	.p2align 3
.L1258:
	movb	%cl, 224(%rdx)
	movq	%r14, 16(%rax)
.L1257:
	movq	8(%rbx), %rsi
	movq	%rbp, %rdi
	imulq	%r12, %rsi
	addq	$1, %r12
	addq	%r13, %rsi
	salq	$4, %rsi
	addq	(%rbx), %rsi
	call	_ZStlsIdcSt11char_traitsIcEERSt13basic_ostreamIT0_T1_ES6_RKSt7complexIT_E@PLT
	cmpq	16(%rbx), %r12
	jge	.L1261
.L1254:
	movq	200(%r15), %rdx
	movq	192(%r15), %rsi
	movq	%rbp, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	testq	%r14, %r14
	je	.L1257
	movq	0(%rbp), %rax
	movzbl	224(%r15), %ecx
	movq	-24(%rax), %rdx
	addq	%rbp, %rdx
	cmpb	$0, 225(%rdx)
	movq	%rdx, %rax
	jne	.L1258
	movq	240(%rdx), %rdi
	testq	%rdi, %rdi
	je	.L1251
	cmpb	$0, 56(%rdi)
	jne	.L1259
	movb	%cl, 16(%rsp)
	movq	%rdx, 8(%rsp)
	movq	%rdi, (%rsp)
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	(%rsp), %rdi
	movzbl	16(%rsp), %ecx
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rdx
	movq	(%rdi), %rax
	movq	48(%rax), %rax
	cmpq	%rdx, %rax
	movq	8(%rsp), %rdx
	jne	.L1260
	movq	0(%rbp), %rax
	movq	-24(%rax), %rdi
	addq	%rbp, %rdi
	movq	%rdi, %rax
.L1259:
	movb	$1, 225(%rdx)
	jmp	.L1258
	.p2align 4,,10
	.p2align 3
.L1286:
	movq	%rax, %r8
	jmp	.L1228
	.p2align 4,,10
	.p2align 3
.L1324:
	addq	$1, 40(%rsp)
	movq	40(%rsp), %rdx
	cmpq	16(%rbx), %rdx
	jl	.L1238
	jmp	.L1221
	.p2align 4,,10
	.p2align 3
.L1261:
	movq	104(%r15), %rdx
	movq	96(%r15), %rsi
	movq	%rbp, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	8(%rbx), %rax
	leaq	-1(%rax), %rdx
	cmpq	%r13, %rdx
	jg	.L1328
	addq	$1, %r13
	cmpq	%rax, %r13
	jge	.L1319
.L1262:
	movq	168(%r15), %rdx
	movq	160(%r15), %rsi
	movq	%rbp, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	jmp	.L1246
	.p2align 4,,10
	.p2align 3
.L1328:
	movq	136(%r15), %rdx
	movq	128(%r15), %rsi
	movq	%rbp, %rdi
	addq	$1, %r13
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	cmpq	%r13, 8(%rbx)
	jg	.L1262
.L1319:
	movq	%r15, %r12
	movzbl	24(%rsp), %r15d
.L1263:
	movq	40(%r12), %rdx
	movq	32(%r12), %rsi
	movq	%rbp, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	cmpq	$0, 104(%rsp)
	je	.L1248
	movq	0(%rbp), %rax
	movq	96(%rsp), %rbx
	movq	-24(%rax), %rax
	movq	%rbx, 8(%rbp,%rax)
.L1248:
	testq	%r14, %r14
	je	.L1218
	movq	0(%rbp), %rax
	movq	-24(%rax), %rbx
	addq	%rbp, %rbx
	cmpb	$0, 225(%rbx)
	movq	%rbx, %rax
	je	.L1329
.L1264:
	movb	%r15b, 224(%rbx)
	movq	32(%rsp), %rbx
	movq	%rbx, 16(%rax)
.L1218:
	movq	536(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L1330
	addq	$552, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	movq	%rbp, %rax
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L1327:
	.cfi_restore_state
	movq	240(%r12), %rdi
	testq	%rdi, %rdi
	je	.L1251
	cmpb	$0, 56(%rdi)
	jne	.L1252
	movb	%dl, 8(%rsp)
	movq	%rdi, (%rsp)
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	(%rsp), %rdi
	movzbl	8(%rsp), %edx
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rsi
	movq	(%rdi), %rax
	movq	48(%rax), %rax
	cmpq	%rsi, %rax
	jne	.L1253
	movq	0(%rbp), %rax
	movq	-24(%rax), %rcx
	addq	%rbp, %rcx
	movq	%rcx, %rax
.L1252:
	movb	$1, 225(%r12)
	jmp	.L1250
	.p2align 4,,10
	.p2align 3
.L1260:
	movb	%cl, 8(%rsp)
	movl	$32, %esi
	movq	%rdx, (%rsp)
	call	*%rax
	movq	0(%rbp), %rax
	movzbl	8(%rsp), %ecx
	movq	-24(%rax), %rdx
	addq	%rbp, %rdx
	movq	%rdx, %rax
	movq	(%rsp), %rdx
	jmp	.L1259
.L1268:
	movq	$0, 104(%rsp)
	movq	$0, 96(%rsp)
	jmp	.L1219
.L1241:
	movq	240(%r13), %rdi
	testq	%rdi, %rdi
	je	.L1251
	cmpb	$0, 56(%rdi)
	je	.L1244
	movzbl	89(%rdi), %r15d
.L1245:
	movb	%r15b, 224(%r13)
	movb	$1, 225(%r13)
	jmp	.L1242
.L1323:
	movq	8(%r12), %rdx
	movq	(%r12), %rsi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	40(%r12), %rdx
	movq	32(%r12), %rsi
	movq	%rax, %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	jmp	.L1218
.L1244:
	movq	%rdi, (%rsp)
	movl	$32, %r15d
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	(%rsp), %rdi
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rdx
	movq	(%rdi), %rax
	movq	48(%rax), %rax
	cmpq	%rdx, %rax
	je	.L1245
	movl	$32, %esi
	call	*%rax
	movl	%eax, %r15d
	jmp	.L1245
.L1329:
	movq	240(%rbx), %r12
	testq	%r12, %r12
	je	.L1251
	cmpb	$0, 56(%r12)
	jne	.L1265
	movq	%r12, %rdi
	call	_ZNKSt5ctypeIcE13_M_widen_initEv@PLT
	movq	(%r12), %rax
	leaq	_ZNKSt5ctypeIcE8do_widenEc(%rip), %rdx
	movq	48(%rax), %rax
	cmpq	%rdx, %rax
	jne	.L1266
.L1322:
	movq	0(%rbp), %rax
	movq	-24(%rax), %rdx
	addq	%rbp, %rdx
	movq	%rdx, %rax
.L1265:
	movb	$1, 225(%rbx)
	jmp	.L1264
.L1253:
	movb	%dl, (%rsp)
	movl	$32, %esi
	call	*%rax
	movq	0(%rbp), %rax
	movq	-24(%rax), %rdx
	addq	%rbp, %rdx
	movq	%rdx, %rax
	movzbl	(%rsp), %edx
	jmp	.L1252
.L1266:
	movl	$32, %esi
	movq	%r12, %rdi
	call	*%rax
	jmp	.L1322
.L1251:
	call	_ZSt16__throw_bad_castv@PLT
.L1330:
	call	__stack_chk_fail@PLT
.L1285:
	movq	%rax, %rbx
	jmp	.L1233
.L1282:
	movq	%rax, %rbx
	jmp	.L1225
.L1283:
	movq	%rax, %rbx
	jmp	.L1226
.L1284:
	movq	%rax, %rbx
	jmp	.L1224
.L1281:
	movq	%rax, %rbx
	jmp	.L1235
.L1233:
	movq	112(%rsp), %rdi
	movq	16(%rsp), %rax
	cmpq	%rax, %rdi
	je	.L1235
	movq	128(%rsp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1235:
	movq	80(%rsp), %rdi
	call	_ZNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEED1Ev@PLT
	movq	%rbx, %rdi
	call	_Unwind_Resume@PLT
.L1226:
	movq	88(%rsp), %rdi
	call	_ZNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEED1Ev
	movq	8+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rax
	movq	48+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rdx
	movq	-24(%rax), %rax
	movq	%rdx, 144(%rsp,%rax)
	movq	-24(%r15), %rax
	movq	40+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rdx
	movq	%r15, 160(%rsp)
	movq	%rdx, 160(%rsp,%rax)
.L1321:
	movq	-24(%r13), %rax
	movq	24+_ZTTNSt7__cxx1118basic_stringstreamIcSt11char_traitsIcESaIcEEE(%rip), %rdx
	movq	%r13, 144(%rsp)
	movq	%rdx, 144(%rsp,%rax)
	xorl	%eax, %eax
	movq	%rax, 152(%rsp)
.L1225:
	movq	(%rsp), %rdi
	leaq	16+_ZTVSt9basic_iosIcSt11char_traitsIcEE(%rip), %rax
	movq	%rax, 272(%rsp)
	call	_ZNSt8ios_baseD2Ev@PLT
	movq	%rbx, %rdi
	call	_Unwind_Resume@PLT
.LEHE54:
.L1224:
	jmp	.L1321
	.cfi_endproc
.LFE11184:
	.section	.gcc_except_table
.LLSDA11184:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE11184-.LLSDACSB11184
.LLSDACSB11184:
	.uleb128 .LEHB48-.LFB11184
	.uleb128 .LEHE48-.LEHB48
	.uleb128 .L1285-.LFB11184
	.uleb128 0
	.uleb128 .LEHB49-.LFB11184
	.uleb128 .LEHE49-.LEHB49
	.uleb128 .L1282-.LFB11184
	.uleb128 0
	.uleb128 .LEHB50-.LFB11184
	.uleb128 .LEHE50-.LEHB50
	.uleb128 .L1284-.LFB11184
	.uleb128 0
	.uleb128 .LEHB51-.LFB11184
	.uleb128 .LEHE51-.LEHB51
	.uleb128 .L1283-.LFB11184
	.uleb128 0
	.uleb128 .LEHB52-.LFB11184
	.uleb128 .LEHE52-.LEHB52
	.uleb128 .L1281-.LFB11184
	.uleb128 0
	.uleb128 .LEHB53-.LFB11184
	.uleb128 .LEHE53-.LEHB53
	.uleb128 .L1285-.LFB11184
	.uleb128 0
	.uleb128 .LEHB54-.LFB11184
	.uleb128 .LEHE54-.LEHB54
	.uleb128 0
	.uleb128 0
.LLSDACSE11184:
	.section	.text._ZN5Eigen8internal12print_matrixINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS6_RKT_RKNS_8IOFormatE,"axG",@progbits,_ZN5Eigen8internal12print_matrixINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS6_RKT_RKNS_8IOFormatE,comdat
	.size	_ZN5Eigen8internal12print_matrixINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS6_RKT_RKNS_8IOFormatE, .-_ZN5Eigen8internal12print_matrixINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS6_RKT_RKNS_8IOFormatE
	.section	.text._ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE,comdat
	.p2align 4
	.weak	_ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE
	.type	_ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE, @function
_ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE:
.LFB10650:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA10650
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	leaq	-368(%rbp), %rdx
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	leaq	-448(%rbp), %r14
	pushq	%r13
	leaq	-480(%rbp), %r15
	.cfi_offset 13, -40
	leaq	-416(%rbp), %r13
	pushq	%r12
	.cfi_offset 12, -48
	leaq	-384(%rbp), %r12
	leaq	-304(%rbp), %r10
	pushq	%rbx
	.cfi_offset 3, -56
	leaq	-320(%rbp), %rbx
	leaq	-432(%rbp), %r9
	leaq	-464(%rbp), %r8
	subq	$488, %rsp
	movq	%rdi, -520(%rbp)
	movzwl	.LC50(%rip), %ecx
	movq	%rsi, -528(%rbp)
	movzwl	.LC51(%rip), %edi
	leaq	-336(%rbp), %rsi
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	xorl	%eax, %eax
	leaq	-352(%rbp), %rax
	movw	%cx, -448(%rbp)
	leaq	-496(%rbp), %rcx
	movq	%rax, -512(%rbp)
	movq	%rax, -368(%rbp)
	leaq	-400(%rbp), %rax
	movw	%di, -480(%rbp)
	movq	%r10, %rdi
	movq	%rbx, -336(%rbp)
	movq	$0, -328(%rbp)
	movb	$0, -320(%rbp)
	movq	$0, -360(%rbp)
	movb	$0, -352(%rbp)
	movq	%r12, -400(%rbp)
	movq	$0, -392(%rbp)
	movb	$0, -384(%rbp)
	movq	%r13, -432(%rbp)
	movq	$0, -424(%rbp)
	movb	$0, -416(%rbp)
	movq	%r14, -464(%rbp)
	movq	$1, -456(%rbp)
	movq	%r15, -496(%rbp)
	movq	$1, -488(%rbp)
	pushq	$32
	pushq	%rsi
	movl	$-1, %esi
	pushq	%rdx
	xorl	%edx, %edx
	pushq	%rax
	movq	%r10, -504(%rbp)
.LEHB55:
	.cfi_escape 0x2e,0x20
	call	_ZN5Eigen8IOFormatC1EiiRKNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEES8_S8_S8_S8_S8_c
.LEHE55:
	movq	-504(%rbp), %rdx
	movq	-528(%rbp), %rsi
	addq	$32, %rsp
	movq	-520(%rbp), %rdi
.LEHB56:
	.cfi_escape 0x2e,0
	call	_ZN5Eigen8internal12print_matrixINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS6_RKT_RKNS_8IOFormatE
.LEHE56:
	movq	-504(%rbp), %rdi
	movq	%rax, -520(%rbp)
	call	_ZN5Eigen8IOFormatD1Ev
	movq	-496(%rbp), %rdi
	cmpq	%r15, %rdi
	je	.L1332
	movq	-480(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1332:
	movq	-464(%rbp), %rdi
	cmpq	%r14, %rdi
	je	.L1333
	movq	-448(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1333:
	movq	-432(%rbp), %rdi
	cmpq	%r13, %rdi
	je	.L1334
	movq	-416(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1334:
	movq	-400(%rbp), %rdi
	cmpq	%r12, %rdi
	je	.L1335
	movq	-384(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1335:
	movq	-368(%rbp), %rdi
	movq	-512(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L1336
	movq	-352(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1336:
	movq	-336(%rbp), %rdi
	cmpq	%rbx, %rdi
	je	.L1331
	movq	-320(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1331:
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L1350
	movq	-520(%rbp), %rax
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
.L1350:
	.cfi_restore_state
	call	__stack_chk_fail@PLT
.L1348:
	movq	%rax, -520(%rbp)
	jmp	.L1338
.L1347:
	movq	%rax, -504(%rbp)
	jmp	.L1339
.L1338:
	movq	-504(%rbp), %rdi
	call	_ZN5Eigen8IOFormatD1Ev
	movq	-520(%rbp), %rax
	movq	%rax, -504(%rbp)
.L1339:
	movq	-496(%rbp), %rdi
	cmpq	%r15, %rdi
	je	.L1340
	movq	-480(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1340:
	movq	-464(%rbp), %rdi
	cmpq	%r14, %rdi
	je	.L1341
	movq	-448(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1341:
	movq	-432(%rbp), %rdi
	cmpq	%r13, %rdi
	je	.L1342
	movq	-416(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1342:
	movq	-400(%rbp), %rdi
	cmpq	%r12, %rdi
	je	.L1343
	movq	-384(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1343:
	movq	-368(%rbp), %rdi
	movq	-512(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L1344
	movq	-352(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1344:
	movq	-336(%rbp), %rdi
	cmpq	%rbx, %rdi
	je	.L1345
	movq	-320(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L1345:
	movq	-504(%rbp), %rdi
.LEHB57:
	call	_Unwind_Resume@PLT
.LEHE57:
	.cfi_endproc
.LFE10650:
	.section	.gcc_except_table
.LLSDA10650:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE10650-.LLSDACSB10650
.LLSDACSB10650:
	.uleb128 .LEHB55-.LFB10650
	.uleb128 .LEHE55-.LEHB55
	.uleb128 .L1347-.LFB10650
	.uleb128 0
	.uleb128 .LEHB56-.LFB10650
	.uleb128 .LEHE56-.LEHB56
	.uleb128 .L1348-.LFB10650
	.uleb128 0
	.uleb128 .LEHB57-.LFB10650
	.uleb128 .LEHE57-.LEHB57
	.uleb128 0
	.uleb128 0
.LLSDACSE10650:
	.section	.text._ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE,"axG",@progbits,_ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE,comdat
	.size	_ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE, .-_ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE
	.section	.text._ZN5Eigen8internal28conditional_aligned_new_autoISt7complexIdELb1EEEPT_m,"axG",@progbits,_ZN5Eigen8internal28conditional_aligned_new_autoISt7complexIdELb1EEEPT_m,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal28conditional_aligned_new_autoISt7complexIdELb1EEEPT_m
	.type	_ZN5Eigen8internal28conditional_aligned_new_autoISt7complexIdELb1EEEPT_m, @function
_ZN5Eigen8internal28conditional_aligned_new_autoISt7complexIdELb1EEEPT_m:
.LFB11598:
	.cfi_startproc
	testq	%rdi, %rdi
	je	.L1356
	movq	%rdi, %rax
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	shrq	$60, %rax
	jne	.L1355
	salq	$4, %rdi
	call	malloc@PLT
	testb	$15, %al
	je	.L1354
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L1354:
	testq	%rax, %rax
	je	.L1355
	addq	$8, %rsp
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L1356:
	xorl	%eax, %eax
	ret
.L1355:
	.cfi_def_cfa_offset 16
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
	.cfi_endproc
.LFE11598:
	.size	_ZN5Eigen8internal28conditional_aligned_new_autoISt7complexIdELb1EEEPT_m, .-_ZN5Eigen8internal28conditional_aligned_new_autoISt7complexIdELb1EEEPT_m
	.section	.rodata._ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS4_l.str1.8,"aMS",@progbits,1
	.align 8
.LC52:
	.string	"Eigen::MapBase<Derived, 0>::MapBase(PointerType, Eigen::Index, Eigen::Index) [with Derived = Eigen::Block<Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true>; PointerType = std::complex<double>*; Eigen::Index = long int]"
	.align 8
.LC53:
	.string	"Eigen::Block<XprType, BlockRows, BlockCols, InnerPanel>::Block(XprType&, Eigen::Index) [with XprType = Eigen::Matrix<std::complex<double>, -1, -1>; int BlockRows = -1; int BlockCols = 1; bool InnerPanel = true; Eigen::Index = long int]"
	.align 8
.LC54:
	.string	"(i>=0) && ( ((BlockRows==1) && (BlockCols==XprType::ColsAtCompileTime) && i<xpr.rows()) ||((BlockRows==XprType::RowsAtCompileTime) && (BlockCols==1) && i<xpr.cols()))"
	.section	.text._ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS4_l,"axG",@progbits,_ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC5ERS4_l,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS4_l
	.type	_ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS4_l, @function
_ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS4_l:
.LFB12070:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	movq	8(%rsi), %rcx
	movq	%rdx, %rax
	imulq	%rcx, %rax
	salq	$4, %rax
	addq	(%rsi), %rax
	movq	%rcx, 8(%rdi)
	movq	%rax, (%rdi)
	je	.L1365
	testq	%rcx, %rcx
	js	.L1369
.L1365:
	movq	%rsi, 24(%rdi)
	movq	$0, 32(%rdi)
	movq	%rdx, 40(%rdi)
	movq	%rcx, 48(%rdi)
	testq	%rdx, %rdx
	js	.L1366
	cmpq	16(%rsi), %rdx
	jge	.L1366
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
.L1369:
	.cfi_restore_state
	leaq	.LC52(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L1366:
	leaq	.LC53(%rip), %rcx
	movl	$120, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC54(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE12070:
	.size	_ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS4_l, .-_ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS4_l
	.weak	_ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS4_l
	.set	_ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS4_l,_ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS4_l
	.section	.rodata._ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS5_l.str1.8,"aMS",@progbits,1
	.align 8
.LC55:
	.string	"Eigen::MapBase<Derived, 0>::MapBase(PointerType, Eigen::Index, Eigen::Index) [with Derived = Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true>; PointerType = const std::complex<double>*; Eigen::Index = long int]"
	.align 8
.LC56:
	.string	"Eigen::Block<XprType, BlockRows, BlockCols, InnerPanel>::Block(XprType&, Eigen::Index) [with XprType = const Eigen::Matrix<std::complex<double>, -1, -1>; int BlockRows = -1; int BlockCols = 1; bool InnerPanel = true; Eigen::Index = long int]"
	.section	.text._ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS5_l,"axG",@progbits,_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC5ERS5_l,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS5_l
	.type	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS5_l, @function
_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS5_l:
.LFB12073:
	.cfi_startproc
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	movq	8(%rsi), %rcx
	movq	%rdx, %rax
	imulq	%rcx, %rax
	salq	$4, %rax
	addq	(%rsi), %rax
	movq	%rcx, 8(%rdi)
	movq	%rax, (%rdi)
	testq	%rcx, %rcx
	jns	.L1371
	testq	%rax, %rax
	jne	.L1375
.L1371:
	movq	%rsi, 24(%rdi)
	movq	$0, 32(%rdi)
	movq	%rdx, 40(%rdi)
	movq	%rcx, 48(%rdi)
	testq	%rdx, %rdx
	js	.L1372
	cmpq	16(%rsi), %rdx
	jge	.L1372
	addq	$8, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 8
	ret
.L1375:
	.cfi_restore_state
	leaq	.LC55(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L1372:
	leaq	.LC56(%rip), %rcx
	movl	$120, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC54(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE12073:
	.size	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS5_l, .-_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS5_l
	.weak	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS5_l
	.set	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS5_l,_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC2ERS5_l
	.section	.text._ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_,"axG",@progbits,_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_
	.type	_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_, @function
_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_:
.LFB12323:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	movq	%rcx, %r13
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	movq	%rsi, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%rdi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	movq	%rdx, %rbx
	subq	$56, %rsp
	.cfi_def_cfa_offset 112
	movq	%fs:40, %rax
	movq	%rax, 40(%rsp)
	xorl	%eax, %eax
	movzbl	_ZGVZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes(%rip), %eax
	testb	%al, %al
	je	.L1436
.L1378:
	movq	_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes(%rip), %r8
	movq	8+_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes(%rip), %r10
	movq	16+_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes(%rip), %rdi
	cmpq	$1, %r13
	jle	.L1383
	movl	$320, %edx
	cmpq	$25743, %r8
	jg	.L1384
	movabsq	$7378697629483820647, %rax
	leaq	-64(%r8), %rcx
	imulq	%rcx
	sarq	$63, %rcx
	movl	$8, %eax
	sarq	$5, %rdx
	subq	%rcx, %rdx
	cmpq	%rax, %rdx
	cmovl	%rax, %rdx
.L1384:
	movq	0(%rbp), %rsi
	cmpq	%rdx, %rsi
	jle	.L1385
	movq	%rdx, %rsi
	andq	$-8, %rsi
	movq	%rsi, 0(%rbp)
.L1385:
	movq	%r10, %rax
	salq	$6, %rsi
	xorl	%edx, %edx
	subq	%r8, %rax
	divq	%rsi
	movq	(%rbx), %rsi
	movq	%rax, %rcx
	leaq	-1(%rsi,%r13), %rax
	cqto
	idivq	%r13
	cmpq	%rax, %rcx
	jle	.L1437
	addq	$3, %rax
	movq	%rax, %rcx
	sarq	$63, %rcx
	shrq	$62, %rcx
	leaq	(%rax,%rcx), %rdx
	andl	$3, %edx
	subq	%rcx, %rdx
	subq	%rdx, %rax
	cmpq	%rsi, %rax
	cmovle	%rax, %rsi
	movq	%rsi, %rcx
.L1387:
	movq	%rcx, (%rbx)
	cmpq	%rdi, %r10
	jge	.L1376
	movq	0(%rbp), %rcx
	subq	%r10, %rdi
	xorl	%edx, %edx
	movq	%rdi, %rax
	imulq	%r13, %rcx
	salq	$4, %rcx
	divq	%rcx
	movq	(%r12), %rcx
	movq	%rax, %rdi
	leaq	-1(%rcx,%r13), %rax
	cqto
	idivq	%r13
	cmpq	%rax, %rdi
	jge	.L1389
	testq	%rdi, %rdi
	jle	.L1389
	movq	%rdi, (%r12)
	.p2align 4,,10
	.p2align 3
.L1376:
	movq	40(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L1438
	addq	$56, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L1437:
	.cfi_restore_state
	movq	%rcx, %rdx
	sarq	$63, %rdx
	shrq	$62, %rdx
	leaq	(%rcx,%rdx), %rax
	andl	$3, %eax
	subq	%rdx, %rax
	subq	%rax, %rcx
	jmp	.L1387
	.p2align 4,,10
	.p2align 3
.L1383:
	movq	(%r12), %r13
	movq	(%rbx), %r11
	movq	0(%rbp), %r9
	cmpq	%r13, %r11
	movq	%r13, %rax
	cmovge	%r11, %rax
	cmpq	%r9, %rax
	cmovl	%r9, %rax
	cmpq	$47, %rax
	jle	.L1376
	movabsq	$7378697629483820647, %rax
	leaq	-64(%r8), %rcx
	imulq	%rcx
	movq	%rcx, %rax
	sarq	$63, %rax
	movq	%rdx, %rsi
	sarq	$5, %rsi
	subq	%rax, %rsi
	andq	$-8, %rsi
	jle	.L1439
	movq	%r9, %rax
	cmpq	%r9, %rsi
	jl	.L1395
	salq	$6, %rax
	movq	%r9, %rbp
	xorl	%edx, %edx
	movq	%r9, 8(%rsp)
	movq	%rax, (%rsp)
	salq	$5, %rbp
	movl	$1572864, %eax
	movq	%r9, %r15
	divq	%rbp
	movq	%rax, %r14
.L1394:
	movq	%r13, %rdx
	movq	%rcx, %rax
	movq	(%rsp), %rcx
	imulq	%r15, %rdx
	salq	$4, %rdx
	subq	%rdx, %rax
	cmpq	%rcx, %rax
	jl	.L1397
	movq	8(%rsp), %rcx
	xorl	%edx, %edx
	salq	$4, %rcx
	divq	%rcx
.L1398:
	cmpq	%r14, %rax
	movq	%r14, %rcx
	cmovle	%rax, %rcx
	andq	$-4, %rcx
	cmpq	%r11, %rcx
	jge	.L1399
	movq	%r11, %rax
	cqto
	idivq	%rcx
	testq	%rdx, %rdx
	je	.L1400
	movq	%rcx, %rsi
	leaq	4(,%rax,4), %rdi
	subq	%rdx, %rsi
	movq	%rsi, %rax
	cqto
	idivq	%rdi
	salq	$2, %rax
	subq	%rax, %rcx
.L1400:
	movq	%rcx, (%rbx)
	jmp	.L1376
	.p2align 4,,10
	.p2align 3
.L1436:
	leaq	_ZGVZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes(%rip), %r14
	movq	%r14, %rdi
	call	__cxa_guard_acquire@PLT
	testl	%eax, %eax
	je	.L1378
	pcmpeqd	%xmm0, %xmm0
	leaq	36(%rsp), %rdx
	leaq	28(%rsp), %rdi
	movq	$-1, 16+_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes(%rip)
	leaq	32(%rsp), %rsi
	movups	%xmm0, _ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes(%rip)
	call	_ZN5Eigen8internal15queryCacheSizesERiS1_S1_
	movslq	28(%rsp), %rax
	movl	$32768, %edx
	movq	%r14, %rdi
	testq	%rax, %rax
	cmovle	%rdx, %rax
	movl	$262144, %edx
	movq	%rax, _ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes(%rip)
	movslq	32(%rsp), %rax
	testq	%rax, %rax
	cmovle	%rdx, %rax
	movl	$2097152, %edx
	movq	%rax, 8+_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes(%rip)
	movslq	36(%rsp), %rax
	testq	%rax, %rax
	cmovle	%rdx, %rax
	movq	%rax, 16+_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes(%rip)
	call	__cxa_guard_release@PLT
	jmp	.L1378
	.p2align 4,,10
	.p2align 3
.L1439:
	cmpq	$1, %r9
	jle	.L1440
	movq	$1, 8(%rsp)
	movl	$49152, %r14d
	movl	$1, %esi
	movl	$1, %r15d
	movq	$64, (%rsp)
	jmp	.L1393
	.p2align 4,,10
	.p2align 3
.L1389:
	cmpq	%rcx, %rax
	cmovg	%rcx, %rax
	movq	%rax, (%r12)
	jmp	.L1376
	.p2align 4,,10
	.p2align 3
.L1399:
	cmpq	%r15, %r9
	jne	.L1376
	imulq	%r11, %r9
	salq	$4, %r9
	cmpq	$1024, %r9
	jg	.L1401
	movq	%r8, %rax
	movq	%r13, %rsi
.L1402:
	movq	8(%rsp), %rbx
	xorl	%edx, %edx
	leaq	(%rbx,%rbx,2), %rcx
	salq	$4, %rcx
	divq	%rcx
	cmpq	%rsi, %rax
	cmovle	%rax, %rsi
	movq	%rsi, %rcx
	testq	%rsi, %rsi
	je	.L1376
	movq	%r13, %rax
	cqto
	idivq	%rsi
	movq	%rax, %r13
	testq	%rdx, %rdx
	je	.L1403
	movq	%rsi, %rax
	addq	$1, %r13
	subq	%rdx, %rax
	cqto
	idivq	%r13
	subq	%rax, %rcx
.L1403:
	movq	%rcx, (%r12)
	jmp	.L1376
	.p2align 4,,10
	.p2align 3
.L1397:
	salq	$6, %rsi
	movl	$4718592, %eax
	xorl	%edx, %edx
	divq	%rsi
	jmp	.L1398
	.p2align 4,,10
	.p2align 3
.L1395:
	cqto
	idivq	%rsi
	testq	%rdx, %rdx
	je	.L1441
	leaq	-1(%rsi), %r11
	leaq	8(,%rax,8), %r13
	movq	%rsi, %r15
	subq	%rdx, %r11
	movq	%r11, %rax
	cqto
	idivq	%r13
	xorl	%edx, %edx
	salq	$3, %rax
	subq	%rax, %r15
	movq	%r15, %rax
	movq	%r15, %r11
	movq	%r15, 8(%rsp)
	salq	$6, %rax
	salq	$5, %r11
	movq	%rax, (%rsp)
	movl	$1572864, %eax
	divq	%r11
	movq	%rax, %r14
.L1393:
	movq	%r15, 0(%rbp)
	movq	(%r12), %r13
	movq	(%rbx), %r11
	jmp	.L1394
	.p2align 4,,10
	.p2align 3
.L1440:
	movq	%r9, %rax
	movq	%r9, %rsi
	xorl	%edx, %edx
	movq	%r9, 8(%rsp)
	salq	$6, %rax
	salq	$5, %rsi
	movq	%r9, %r15
	movq	%rax, (%rsp)
	movl	$1572864, %eax
	divq	%rsi
	movl	$1, %esi
	movq	%rax, %r14
	jmp	.L1394
	.p2align 4,,10
	.p2align 3
.L1401:
	testq	%rdi, %rdi
	je	.L1410
	cmpq	$32768, %r9
	jle	.L1442
.L1410:
	movq	%r13, %rsi
	movl	$1572864, %eax
	jmp	.L1402
	.p2align 4,,10
	.p2align 3
.L1441:
	movq	%rsi, %rax
	movq	%rsi, %r11
	xorl	%edx, %edx
	movq	%rsi, 8(%rsp)
	salq	$6, %rax
	salq	$5, %r11
	movq	%rsi, %r15
	movq	%rax, (%rsp)
	movl	$1572864, %eax
	divq	%r11
	movq	%rax, %r14
	jmp	.L1393
.L1442:
	movl	$576, %esi
	movq	%r10, %rax
	cmpq	%rsi, %r13
	cmovle	%r13, %rsi
	jmp	.L1402
.L1438:
	call	__stack_chk_fail@PLT
	.cfi_endproc
.LFE12323:
	.size	_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_, .-_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_
	.section	.text._ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb0ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_,"axG",@progbits,_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb0ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb0ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_
	.type	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb0ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_, @function
_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb0ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_:
.LFB12542:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	leaq	-3(%rdi), %rax
	movapd	%xmm0, %xmm13
	movapd	%xmm1, %xmm14
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%rsi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$344, %rsp
	.cfi_def_cfa_offset 400
	movq	8(%rdx), %rbx
	movq	(%rdx), %r11
	movq	%rax, 88(%rsp)
	leaq	-1(%rdi), %rax
	movq	%rbx, %rsi
	movq	%rax, 200(%rsp)
	xorl	%eax, %eax
	salq	$4, %rsi
	movq	%rdi, 160(%rsp)
	movq	%rcx, 16(%rsp)
	movq	%r8, 80(%rsp)
	movq	%r9, 184(%rsp)
	movq	%r11, 176(%rsp)
	movq	%rbx, 192(%rsp)
	movq	%rsi, 168(%rsp)
	cmpq	$32000, %rsi
	ja	.L1444
	subq	$7, %rdi
	movq	%rdi, 104(%rsp)
	testq	%rdi, %rdi
	jle	.L1444
	movq	%rbx, %rax
	movq	%r9, %rcx
	movq	%r8, %r15
	movq	%rbx, %r10
	salq	$5, %rax
	leaq	(%rbx,%rbx,2), %r13
	movq	%rbx, %rdi
	movq	%r11, 72(%rsp)
	movq	%rax, %r9
	movq	%rcx, %rax
	salq	$4, %rcx
	leaq	(%rbx,%rbx,4), %r14
	salq	$7, %rax
	movq	%r13, %r8
	salq	$6, %rdi
	movq	%rbp, %rdx
	leaq	0(,%rbx,8), %r12
	movq	%rax, 136(%rsp)
	leaq	(%r15,%rcx), %rax
	xorl	%r15d, %r15d
	movq	$0, 8(%rsp)
	subq	%rbx, %r12
	leaq	(%rax,%rcx), %rbx
	salq	$4, %r8
	movq	%rbx, 144(%rsp)
	addq	%rcx, %rbx
	movq	%r8, %rbp
	salq	$7, %r10
	movq	%rbx, 96(%rsp)
	addq	%rcx, %rbx
	salq	$4, %r14
	movq	%rdx, %r8
	movq	%rbx, 152(%rsp)
	salq	$5, %r13
	addq	%rcx, %rbx
	salq	$4, %r12
	movq	.LC8(%rip), %xmm15
	movq	%rbx, 112(%rsp)
	addq	%rcx, %rbx
	movq	%rbx, 120(%rsp)
	addq	%rcx, %rbx
	movq	%r15, %rcx
	movq	%rdi, %r15
	movq	%rbx, 128(%rsp)
	movq	%rcx, %rdi
	xorl	%ebx, %ebx
	.p2align 4,,10
	.p2align 3
.L1445:
	testq	%r8, %r8
	jle	.L1489
	movq	16(%rsp), %rdx
	movq	%rax, 24(%rsp)
	pxor	%xmm9, %xmm9
	movq	72(%rsp), %rcx
	movapd	%xmm9, %xmm3
	movapd	%xmm9, %xmm2
	movapd	%xmm9, %xmm4
	movq	(%rdx), %r11
	movapd	%xmm9, %xmm7
	xorl	%edx, %edx
	movapd	%xmm9, %xmm8
	movapd	%xmm9, %xmm5
	movapd	%xmm9, %xmm6
	.p2align 4,,10
	.p2align 3
.L1446:
	movq	%rdx, %rax
	movdqu	(%rcx), %xmm12
	addq	$1, %rdx
	salq	$4, %rax
	movupd	(%r11,%rax), %xmm0
	pshufd	$238, %xmm12, %xmm11
	pshufd	$68, %xmm12, %xmm10
	leaq	(%rdi,%rcx), %rax
	movdqu	(%rax,%rsi), %xmm12
	addq	$16, %rcx
	mulpd	%xmm0, %xmm10
	pshufd	$78, %xmm0, %xmm1
	mulpd	%xmm1, %xmm11
	xorpd	%xmm15, %xmm11
	addpd	%xmm11, %xmm10
	pshufd	$238, %xmm12, %xmm11
	mulpd	%xmm1, %xmm11
	addpd	%xmm10, %xmm2
	pshufd	$68, %xmm12, %xmm10
	movdqu	(%rax,%r9), %xmm12
	mulpd	%xmm0, %xmm10
	xorpd	%xmm15, %xmm11
	addpd	%xmm11, %xmm10
	pshufd	$238, %xmm12, %xmm11
	mulpd	%xmm1, %xmm11
	addpd	%xmm10, %xmm4
	pshufd	$68, %xmm12, %xmm10
	movdqu	(%rax,%rbp), %xmm12
	mulpd	%xmm0, %xmm10
	xorpd	%xmm15, %xmm11
	addpd	%xmm11, %xmm10
	pshufd	$238, %xmm12, %xmm11
	mulpd	%xmm1, %xmm11
	addpd	%xmm10, %xmm3
	pshufd	$68, %xmm12, %xmm10
	movdqu	(%rax,%r15), %xmm12
	mulpd	%xmm0, %xmm10
	xorpd	%xmm15, %xmm11
	addpd	%xmm11, %xmm10
	pshufd	$238, %xmm12, %xmm11
	mulpd	%xmm1, %xmm11
	addpd	%xmm10, %xmm5
	pshufd	$68, %xmm12, %xmm10
	movdqu	(%rax,%r14), %xmm12
	mulpd	%xmm0, %xmm10
	xorpd	%xmm15, %xmm11
	addpd	%xmm11, %xmm10
	pshufd	$238, %xmm12, %xmm11
	mulpd	%xmm1, %xmm11
	addpd	%xmm10, %xmm9
	pshufd	$68, %xmm12, %xmm10
	movdqu	(%rax,%r13), %xmm12
	mulpd	%xmm0, %xmm10
	xorpd	%xmm15, %xmm11
	addpd	%xmm11, %xmm10
	pshufd	$238, %xmm12, %xmm11
	mulpd	%xmm1, %xmm11
	addpd	%xmm10, %xmm8
	pshufd	$68, %xmm12, %xmm10
	mulpd	%xmm0, %xmm10
	xorpd	%xmm15, %xmm11
	addpd	%xmm11, %xmm10
	addpd	%xmm10, %xmm7
	movupd	(%rax,%r12), %xmm10
	pshufd	$238, %xmm10, %xmm11
	pshufd	$68, %xmm10, %xmm10
	mulpd	%xmm11, %xmm1
	mulpd	%xmm10, %xmm0
	xorpd	%xmm15, %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm0, %xmm6
	cmpq	%rdx, %r8
	jne	.L1446
	movapd	%xmm3, %xmm10
	unpckhpd	%xmm3, %xmm3
	movq	24(%rsp), %rax
	movapd	%xmm4, %xmm11
	movapd	%xmm3, %xmm12
	movapd	%xmm2, %xmm3
	unpckhpd	%xmm4, %xmm4
	movhpd	%xmm5, 24(%rsp)
	unpckhpd	%xmm3, %xmm3
	movhpd	%xmm6, 32(%rsp)
	movhpd	%xmm7, 40(%rsp)
	movhpd	%xmm8, 48(%rsp)
	movhpd	%xmm9, 56(%rsp)
.L1456:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	mulsd	%xmm3, %xmm1
	mulsd	%xmm2, %xmm0
	subsd	%xmm1, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm2, %xmm1
	movsd	%xmm1, 64(%rsp)
	movapd	%xmm13, %xmm1
	mulsd	%xmm3, %xmm1
	addsd	64(%rsp), %xmm1
	ucomisd	%xmm0, %xmm1
	jp	.L1490
.L1447:
	movq	80(%rsp), %rcx
	movapd	%xmm13, %xmm2
	mulsd	%xmm4, %xmm2
	addsd	8(%rcx,%rbx), %xmm1
	addsd	(%rcx,%rbx), %xmm0
	movsd	%xmm1, 8(%rcx,%rbx)
	movapd	%xmm14, %xmm1
	movsd	%xmm0, (%rcx,%rbx)
	mulsd	%xmm4, %xmm1
	movapd	%xmm13, %xmm0
	mulsd	%xmm11, %xmm0
	subsd	%xmm1, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm11, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm0, %xmm1
	jp	.L1491
.L1448:
	addsd	8(%rax,%rbx), %xmm1
	addsd	(%rax,%rbx), %xmm0
	movapd	%xmm13, %xmm2
	mulsd	%xmm12, %xmm2
	movsd	%xmm0, (%rax,%rbx)
	movapd	%xmm13, %xmm0
	movsd	%xmm1, 8(%rax,%rbx)
	mulsd	%xmm10, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm12, %xmm1
	subsd	%xmm1, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm10, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm0, %xmm1
	jp	.L1492
.L1449:
	movq	144(%rsp), %rdx
	movsd	24(%rsp), %xmm3
	movapd	%xmm14, %xmm2
	mulsd	%xmm5, %xmm2
	addsd	8(%rdx,%rbx), %xmm1
	addsd	(%rdx,%rbx), %xmm0
	movsd	%xmm1, 8(%rdx,%rbx)
	movapd	%xmm3, %xmm1
	movsd	%xmm0, (%rdx,%rbx)
	mulsd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	mulsd	%xmm5, %xmm0
	mulsd	%xmm13, %xmm3
	subsd	%xmm1, %xmm0
	movapd	%xmm3, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1493
.L1450:
	movq	96(%rsp), %rcx
	movsd	56(%rsp), %xmm5
	movapd	%xmm14, %xmm2
	mulsd	%xmm9, %xmm2
	addsd	8(%rcx,%rbx), %xmm1
	addsd	(%rcx,%rbx), %xmm0
	movsd	%xmm1, 8(%rcx,%rbx)
	movapd	%xmm5, %xmm1
	movsd	%xmm0, (%rcx,%rbx)
	mulsd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	mulsd	%xmm9, %xmm0
	mulsd	%xmm13, %xmm5
	subsd	%xmm1, %xmm0
	movapd	%xmm5, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1494
.L1451:
	movq	152(%rsp), %rdx
	movsd	48(%rsp), %xmm5
	movapd	%xmm14, %xmm2
	mulsd	%xmm8, %xmm2
	addsd	8(%rdx,%rbx), %xmm1
	addsd	(%rdx,%rbx), %xmm0
	movsd	%xmm1, 8(%rdx,%rbx)
	movapd	%xmm5, %xmm1
	movsd	%xmm0, (%rdx,%rbx)
	mulsd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	mulsd	%xmm8, %xmm0
	mulsd	%xmm13, %xmm5
	subsd	%xmm1, %xmm0
	movapd	%xmm5, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1495
.L1452:
	movq	112(%rsp), %rcx
	movsd	40(%rsp), %xmm3
	movapd	%xmm14, %xmm2
	mulsd	%xmm7, %xmm2
	addsd	8(%rcx,%rbx), %xmm1
	addsd	(%rcx,%rbx), %xmm0
	movsd	%xmm1, 8(%rcx,%rbx)
	movapd	%xmm3, %xmm1
	movsd	%xmm0, (%rcx,%rbx)
	mulsd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	mulsd	%xmm7, %xmm0
	mulsd	%xmm13, %xmm3
	subsd	%xmm1, %xmm0
	movapd	%xmm3, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1496
.L1453:
	movq	120(%rsp), %rdx
	movsd	32(%rsp), %xmm7
	movapd	%xmm14, %xmm2
	mulsd	%xmm6, %xmm2
	addsd	8(%rdx,%rbx), %xmm1
	addsd	(%rdx,%rbx), %xmm0
	movsd	%xmm1, 8(%rdx,%rbx)
	movapd	%xmm7, %xmm1
	movsd	%xmm0, (%rdx,%rbx)
	mulsd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	mulsd	%xmm6, %xmm0
	mulsd	%xmm13, %xmm7
	subsd	%xmm1, %xmm0
	movapd	%xmm7, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1497
.L1454:
	addq	$8, 8(%rsp)
	subq	%r10, %rdi
	addq	%r10, %rsi
	addq	%r10, %r9
	addq	%r10, 72(%rsp)
	movq	8(%rsp), %rdx
	addq	%r10, %rbp
	addq	%r10, %r15
	movq	128(%rsp), %rcx
	addq	%r10, %r14
	addq	%r10, %r13
	addq	%r10, %r12
	addsd	8(%rcx,%rbx), %xmm1
	addsd	(%rcx,%rbx), %xmm0
	movsd	%xmm1, 8(%rcx,%rbx)
	movsd	%xmm0, (%rcx,%rbx)
	movq	136(%rsp), %rcx
	addq	%rcx, %rbx
	cmpq	%rdx, 104(%rsp)
	jg	.L1445
	movq	160(%rsp), %rax
	movq	%r8, %rbp
	subq	$8, %rax
	shrq	$3, %rax
	leaq	8(,%rax,8), %rax
.L1444:
	movq	88(%rsp), %rdi
	cmpq	%rdi, %rax
	jge	.L1457
	movq	184(%rsp), %r10
	leaq	3(%rax), %rbx
	movq	80(%rsp), %rdx
	leaq	1(%rax), %r13
	leaq	2(%rax), %r12
	movq	176(%rsp), %rcx
	movq	192(%rsp), %r9
	movq	%rax, 24(%rsp)
	movq	%r10, %rsi
	movq	.LC8(%rip), %xmm4
	salq	$6, %rsi
	salq	$6, %r9
	movq	%rsi, 8(%rsp)
	movq	%r10, %rsi
	imulq	%rbx, %rsi
	movq	%rsi, %r11
	movq	%rsi, %r8
	movq	%rsi, %rdi
	subq	%r10, %rsi
	salq	$4, %r11
	salq	$4, %rsi
	leaq	(%rdx,%r11), %r15
	movq	168(%rsp), %rdx
	movq	%rdx, %r14
	imulq	%rdx, %r13
	imulq	%rdx, %r12
	imulq	%rdx, %rbx
	leaq	(%r10,%r10), %rdx
	imulq	%rax, %r14
	addq	%rcx, %r13
	subq	%rdx, %rdi
	addq	%rcx, %r12
	salq	$4, %rdi
	addq	%rcx, %rbx
	addq	%rcx, %r14
	leaq	(%rdx,%r10), %rcx
	movq	%rax, %r10
	subq	%rcx, %r8
	salq	$4, %r8
	.p2align 4,,10
	.p2align 3
.L1458:
	testq	%rbp, %rbp
	jle	.L1498
	movq	16(%rsp), %rax
	pxor	%xmm3, %xmm3
	xorl	%edx, %edx
	movapd	%xmm3, %xmm5
	movapd	%xmm3, %xmm6
	movapd	%xmm3, %xmm7
	movq	(%rax), %rcx
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L1459:
	movupd	(%rcx,%rax), %xmm0
	movdqu	(%r14,%rax), %xmm2
	addq	$1, %rdx
	pshufd	$78, %xmm0, %xmm1
	pshufd	$238, %xmm2, %xmm8
	pshufd	$68, %xmm2, %xmm2
	mulpd	%xmm1, %xmm8
	mulpd	%xmm0, %xmm2
	xorpd	%xmm4, %xmm8
	addpd	%xmm8, %xmm2
	addpd	%xmm2, %xmm7
	movdqu	0(%r13,%rax), %xmm2
	pshufd	$238, %xmm2, %xmm8
	pshufd	$68, %xmm2, %xmm2
	mulpd	%xmm1, %xmm8
	mulpd	%xmm0, %xmm2
	xorpd	%xmm4, %xmm8
	addpd	%xmm8, %xmm2
	addpd	%xmm2, %xmm6
	movdqu	(%r12,%rax), %xmm2
	pshufd	$238, %xmm2, %xmm8
	pshufd	$68, %xmm2, %xmm2
	mulpd	%xmm1, %xmm8
	mulpd	%xmm0, %xmm2
	xorpd	%xmm4, %xmm8
	addpd	%xmm8, %xmm2
	addpd	%xmm2, %xmm5
	movdqu	(%rbx,%rax), %xmm2
	pshufd	$238, %xmm2, %xmm2
	mulpd	%xmm2, %xmm1
	movdqu	(%rbx,%rax), %xmm2
	addq	$16, %rax
	pshufd	$68, %xmm2, %xmm2
	mulpd	%xmm2, %xmm0
	xorpd	%xmm4, %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm0, %xmm3
	cmpq	%rdx, %rbp
	jne	.L1459
	movapd	%xmm3, %xmm0
	movapd	%xmm3, %xmm9
	movapd	%xmm5, %xmm3
	unpckhpd	%xmm3, %xmm3
	unpckhpd	%xmm0, %xmm0
	movapd	%xmm3, %xmm10
	movapd	%xmm6, %xmm3
	movapd	%xmm0, %xmm8
	unpckhpd	%xmm3, %xmm3
	movapd	%xmm3, %xmm11
	movapd	%xmm7, %xmm3
	unpckhpd	%xmm3, %xmm3
.L1465:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm14, %xmm2
	mulsd	%xmm3, %xmm1
	mulsd	%xmm7, %xmm0
	mulsd	%xmm7, %xmm2
	subsd	%xmm1, %xmm0
	movapd	%xmm13, %xmm1
	mulsd	%xmm3, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1499
.L1460:
	movq	%r15, %rax
	movapd	%xmm14, %xmm2
	mulsd	%xmm6, %xmm2
	subq	%r11, %rax
	addsd	8(%rax,%r8), %xmm1
	addsd	(%rax,%r8), %xmm0
	movsd	%xmm1, 8(%rax,%r8)
	movapd	%xmm14, %xmm1
	movsd	%xmm0, (%rax,%r8)
	mulsd	%xmm11, %xmm1
	movapd	%xmm13, %xmm0
	mulsd	%xmm6, %xmm0
	subsd	%xmm1, %xmm0
	movapd	%xmm13, %xmm1
	mulsd	%xmm11, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1500
.L1461:
	addsd	8(%rax,%rdi), %xmm1
	addsd	(%rax,%rdi), %xmm0
	movapd	%xmm14, %xmm2
	mulsd	%xmm5, %xmm2
	movsd	%xmm0, (%rax,%rdi)
	movapd	%xmm13, %xmm0
	movsd	%xmm1, 8(%rax,%rdi)
	mulsd	%xmm5, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm10, %xmm1
	subsd	%xmm1, %xmm0
	movapd	%xmm13, %xmm1
	mulsd	%xmm10, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1501
.L1462:
	addsd	8(%rax,%rsi), %xmm1
	addsd	(%rax,%rsi), %xmm0
	movapd	%xmm14, %xmm2
	mulsd	%xmm9, %xmm2
	movsd	%xmm0, (%rax,%rsi)
	movapd	%xmm13, %xmm0
	movsd	%xmm1, 8(%rax,%rsi)
	mulsd	%xmm9, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm8, %xmm1
	subsd	%xmm1, %xmm0
	movapd	%xmm13, %xmm1
	mulsd	%xmm8, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1502
.L1463:
	addsd	8(%r15), %xmm1
	addq	$4, %r10
	addq	%r9, %r14
	addq	%r9, %r13
	addsd	(%r15), %xmm0
	movq	8(%rsp), %rax
	addq	%r9, %r12
	addq	%r9, %rbx
	movsd	%xmm1, 8(%r15)
	movsd	%xmm0, (%r15)
	addq	%rax, %r15
	movq	88(%rsp), %rax
	cmpq	%rax, %r10
	jl	.L1458
	movq	160(%rsp), %rdi
	movq	24(%rsp), %rax
	leaq	-4(%rdi), %rdx
	subq	%rax, %rdx
	andq	$-4, %rdx
	leaq	4(%rax,%rdx), %rax
.L1457:
	movq	200(%rsp), %r9
	cmpq	%r9, %rax
	jge	.L1466
	movq	184(%rsp), %rsi
	leaq	1(%rax), %rbx
	movq	80(%rsp), %rcx
	movq	%rax, %r15
	movq	192(%rsp), %r8
	movq	.LC8(%rip), %xmm4
	movq	%rsi, %rdx
	movq	%rsi, %rdi
	imulq	%rbx, %rdx
	salq	$5, %rdi
	salq	$5, %r8
	movq	%rdx, %r14
	subq	%rsi, %rdx
	movq	168(%rsp), %rsi
	salq	$4, %rdx
	salq	$4, %r14
	movq	%rsi, %r12
	imulq	%rsi, %rbx
	leaq	(%rcx,%rdx), %r13
	addq	%rcx, %r14
	imulq	%rax, %r12
	movq	176(%rsp), %rdx
	addq	%rdx, %rbx
	addq	%rdx, %r12
.L1467:
	testq	%rbp, %rbp
	jle	.L1503
	movq	16(%rsp), %rsi
	pxor	%xmm6, %xmm6
	xorl	%edx, %edx
	xorl	%ecx, %ecx
	movapd	%xmm6, %xmm5
	movq	(%rsi), %rsi
	.p2align 4,,10
	.p2align 3
.L1468:
	movupd	(%rsi,%rdx), %xmm1
	movdqu	(%r12,%rdx), %xmm7
	addq	$1, %rcx
	pshufd	$78, %xmm1, %xmm2
	pshufd	$238, %xmm7, %xmm3
	pshufd	$68, %xmm7, %xmm0
	movdqu	(%rbx,%rdx), %xmm7
	mulpd	%xmm2, %xmm3
	addq	$16, %rdx
	mulpd	%xmm1, %xmm0
	xorpd	%xmm4, %xmm3
	addpd	%xmm3, %xmm0
	addpd	%xmm0, %xmm6
	pshufd	$238, %xmm7, %xmm0
	mulpd	%xmm0, %xmm2
	pshufd	$68, %xmm7, %xmm0
	mulpd	%xmm0, %xmm1
	xorpd	%xmm4, %xmm2
	addpd	%xmm2, %xmm1
	addpd	%xmm1, %xmm5
	cmpq	%rcx, %rbp
	jne	.L1468
	movapd	%xmm5, %xmm7
	movapd	%xmm6, %xmm3
	movapd	%xmm6, %xmm2
	unpckhpd	%xmm7, %xmm7
	unpckhpd	%xmm3, %xmm3
.L1472:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm14, %xmm6
	mulsd	%xmm3, %xmm1
	mulsd	%xmm2, %xmm0
	mulsd	%xmm2, %xmm6
	subsd	%xmm1, %xmm0
	movapd	%xmm13, %xmm1
	mulsd	%xmm3, %xmm1
	addsd	%xmm6, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1504
.L1469:
	addsd	8(%r13), %xmm1
	addsd	0(%r13), %xmm0
	movapd	%xmm14, %xmm2
	mulsd	%xmm5, %xmm2
	movsd	%xmm0, 0(%r13)
	movapd	%xmm13, %xmm0
	movsd	%xmm1, 8(%r13)
	mulsd	%xmm5, %xmm0
	movapd	%xmm14, %xmm1
	mulsd	%xmm7, %xmm1
	subsd	%xmm1, %xmm0
	movapd	%xmm13, %xmm1
	mulsd	%xmm7, %xmm1
	addsd	%xmm2, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1505
.L1470:
	addsd	8(%r14), %xmm1
	addq	$2, %r15
	addq	%rdi, %r13
	addq	%r8, %r12
	addsd	(%r14), %xmm0
	addq	%r8, %rbx
	movsd	%xmm1, 8(%r14)
	movsd	%xmm0, (%r14)
	addq	%rdi, %r14
	cmpq	%r9, %r15
	jl	.L1467
	movq	160(%rsp), %rdi
	leaq	-2(%rdi), %rdx
	subq	%rax, %rdx
	andq	$-2, %rdx
	leaq	2(%rax,%rdx), %rax
.L1466:
	movq	160(%rsp), %r15
	cmpq	%rax, %r15
	jle	.L1443
	movq	184(%rsp), %rdi
	movq	168(%rsp), %r14
	movq	.LC8(%rip), %xmm4
	movq	%rdi, %r13
	imulq	%rax, %rdi
	movq	%r14, %r12
	imulq	%rax, %r12
	salq	$4, %r13
	movq	%rdi, %rbx
	movq	80(%rsp), %rdi
	salq	$4, %rbx
	addq	%rbx, %rdi
	movq	%rdi, %rbx
	movq	176(%rsp), %rdi
	addq	%r12, %rdi
	movq	%rdi, %r12
	movq	16(%rsp), %rdi
.L1474:
	pxor	%xmm2, %xmm2
	movapd	%xmm2, %xmm3
	testq	%rbp, %rbp
	jle	.L1478
	movq	(%rdi), %rsi
	xorl	%edx, %edx
	pxor	%xmm2, %xmm2
	xorl	%ecx, %ecx
	.p2align 4,,10
	.p2align 3
.L1475:
	movdqu	(%rsi,%rdx), %xmm6
	addq	$1, %rcx
	pshufd	$78, %xmm6, %xmm1
	movdqu	(%r12,%rdx), %xmm6
	pshufd	$238, %xmm6, %xmm0
	mulpd	%xmm0, %xmm1
	pshufd	$68, %xmm6, %xmm0
	movupd	(%rsi,%rdx), %xmm6
	addq	$16, %rdx
	mulpd	%xmm6, %xmm0
	xorpd	%xmm4, %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm0, %xmm2
	cmpq	%rcx, %rbp
	jne	.L1475
	movapd	%xmm2, %xmm6
	unpckhpd	%xmm6, %xmm6
	movapd	%xmm6, %xmm3
.L1478:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm14, %xmm5
	mulsd	%xmm3, %xmm1
	mulsd	%xmm2, %xmm0
	mulsd	%xmm2, %xmm5
	subsd	%xmm1, %xmm0
	movapd	%xmm13, %xmm1
	mulsd	%xmm3, %xmm1
	addsd	%xmm5, %xmm1
	ucomisd	%xmm1, %xmm0
	jp	.L1506
.L1476:
	addsd	8(%rbx), %xmm1
	addsd	(%rbx), %xmm0
	addq	$1, %rax
	addq	%r14, %r12
	movsd	%xmm0, (%rbx)
	movsd	%xmm1, 8(%rbx)
	addq	%r13, %rbx
	cmpq	%rax, %r15
	jne	.L1474
.L1443:
	addq	$344, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L1498:
	.cfi_restore_state
	pxor	%xmm7, %xmm7
	movapd	%xmm7, %xmm3
	movapd	%xmm7, %xmm6
	movapd	%xmm7, %xmm11
	movapd	%xmm7, %xmm5
	movapd	%xmm7, %xmm10
	movapd	%xmm7, %xmm9
	movapd	%xmm7, %xmm8
	jmp	.L1465
	.p2align 4,,10
	.p2align 3
.L1489:
	pxor	%xmm9, %xmm9
	movapd	%xmm9, %xmm8
	movapd	%xmm9, %xmm7
	movsd	%xmm9, 56(%rsp)
	movapd	%xmm9, %xmm6
	movapd	%xmm9, %xmm2
	movapd	%xmm9, %xmm3
	movapd	%xmm9, %xmm5
	movsd	%xmm9, 48(%rsp)
	movapd	%xmm9, %xmm4
	movapd	%xmm9, %xmm11
	movapd	%xmm9, %xmm12
	movsd	%xmm9, 40(%rsp)
	movsd	%xmm9, 32(%rsp)
	movapd	%xmm9, %xmm10
	movsd	%xmm9, 24(%rsp)
	jmp	.L1456
.L1503:
	pxor	%xmm2, %xmm2
	movapd	%xmm2, %xmm3
	movapd	%xmm2, %xmm5
	movapd	%xmm2, %xmm7
	jmp	.L1472
.L1497:
	movsd	32(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm6, %xmm2
	movapd	%xmm13, %xmm0
	movq	%r8, 216(%rsp)
	movq	%r10, 208(%rsp)
	movq	%rdi, 64(%rsp)
	movq	%rsi, 56(%rsp)
	movq	%r9, 48(%rsp)
	movq	%rax, 40(%rsp)
	movsd	%xmm14, 32(%rsp)
	movsd	%xmm13, 24(%rsp)
	call	__muldc3@PLT
	movq	64(%rsp), %rdi
	movq	56(%rsp), %rsi
	movq	216(%rsp), %r8
	movq	208(%rsp), %r10
	movq	48(%rsp), %r9
	movq	40(%rsp), %rax
	movq	.LC8(%rip), %xmm15
	movsd	32(%rsp), %xmm14
	movsd	24(%rsp), %xmm13
	jmp	.L1454
.L1502:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm8, %xmm3
	movq	%r11, 104(%rsp)
	movapd	%xmm9, %xmm2
	movq	%r9, 96(%rsp)
	movq	%rsi, 72(%rsp)
	movq	%rdi, 64(%rsp)
	movq	%r8, 56(%rsp)
	movq	%r10, 48(%rsp)
	movsd	%xmm14, 40(%rsp)
	movsd	%xmm13, 32(%rsp)
	call	__muldc3@PLT
	movq	104(%rsp), %r11
	movq	96(%rsp), %r9
	movq	72(%rsp), %rsi
	movq	64(%rsp), %rdi
	movq	56(%rsp), %r8
	movq	48(%rsp), %r10
	movq	.LC8(%rip), %xmm4
	movsd	40(%rsp), %xmm14
	movsd	32(%rsp), %xmm13
	jmp	.L1463
.L1501:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm5, %xmm2
	movq	%r11, 128(%rsp)
	movapd	%xmm10, %xmm3
	movq	%r9, 120(%rsp)
	movq	%rsi, 112(%rsp)
	movq	%rdi, 104(%rsp)
	movq	%r8, 96(%rsp)
	movq	%r10, 72(%rsp)
	movq	%rax, 48(%rsp)
	movsd	%xmm9, 64(%rsp)
	movsd	%xmm8, 56(%rsp)
	movsd	%xmm14, 40(%rsp)
	movsd	%xmm13, 32(%rsp)
	call	__muldc3@PLT
	movq	120(%rsp), %r9
	movq	112(%rsp), %rsi
	movq	128(%rsp), %r11
	movq	104(%rsp), %rdi
	movq	96(%rsp), %r8
	movq	72(%rsp), %r10
	movq	.LC8(%rip), %xmm4
	movsd	64(%rsp), %xmm9
	movsd	56(%rsp), %xmm8
	movq	48(%rsp), %rax
	movsd	40(%rsp), %xmm14
	movsd	32(%rsp), %xmm13
	jmp	.L1462
.L1500:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm6, %xmm2
	movq	%r11, 144(%rsp)
	movapd	%xmm11, %xmm3
	movq	%r9, 136(%rsp)
	movq	%rsi, 128(%rsp)
	movq	%rdi, 120(%rsp)
	movq	%r8, 112(%rsp)
	movq	%r10, 104(%rsp)
	movq	%rax, 48(%rsp)
	movsd	%xmm5, 96(%rsp)
	movsd	%xmm10, 72(%rsp)
	movsd	%xmm9, 64(%rsp)
	movsd	%xmm8, 56(%rsp)
	movsd	%xmm14, 40(%rsp)
	movsd	%xmm13, 32(%rsp)
	call	__muldc3@PLT
	movq	120(%rsp), %rdi
	movq	112(%rsp), %r8
	movq	144(%rsp), %r11
	movq	136(%rsp), %r9
	movq	128(%rsp), %rsi
	movq	104(%rsp), %r10
	movq	.LC8(%rip), %xmm4
	movsd	96(%rsp), %xmm5
	movsd	72(%rsp), %xmm10
	movsd	64(%rsp), %xmm9
	movsd	56(%rsp), %xmm8
	movq	48(%rsp), %rax
	movsd	40(%rsp), %xmm14
	movsd	32(%rsp), %xmm13
	jmp	.L1461
.L1499:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm7, %xmm2
	movq	%r11, 152(%rsp)
	movq	%r9, 144(%rsp)
	movq	%rsi, 136(%rsp)
	movq	%rdi, 128(%rsp)
	movq	%r8, 120(%rsp)
	movq	%r10, 112(%rsp)
	movsd	%xmm6, 104(%rsp)
	movsd	%xmm11, 96(%rsp)
	movsd	%xmm5, 72(%rsp)
	movsd	%xmm10, 64(%rsp)
	movsd	%xmm9, 56(%rsp)
	movsd	%xmm8, 48(%rsp)
	movsd	%xmm14, 40(%rsp)
	movsd	%xmm13, 32(%rsp)
	call	__muldc3@PLT
	movq	152(%rsp), %r11
	movq	144(%rsp), %r9
	movq	136(%rsp), %rsi
	movq	128(%rsp), %rdi
	movq	120(%rsp), %r8
	movq	112(%rsp), %r10
	movq	.LC8(%rip), %xmm4
	movsd	104(%rsp), %xmm6
	movsd	96(%rsp), %xmm11
	movsd	72(%rsp), %xmm5
	movsd	64(%rsp), %xmm10
	movsd	56(%rsp), %xmm9
	movsd	48(%rsp), %xmm8
	movsd	40(%rsp), %xmm14
	movsd	32(%rsp), %xmm13
	jmp	.L1460
.L1490:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movq	%r8, 328(%rsp)
	movq	%r10, 320(%rsp)
	movq	%rdi, 312(%rsp)
	movq	%rsi, 304(%rsp)
	movq	%r9, 296(%rsp)
	movq	%rax, 288(%rsp)
	movsd	%xmm9, 280(%rsp)
	movsd	%xmm8, 272(%rsp)
	movsd	%xmm7, 264(%rsp)
	movsd	%xmm6, 256(%rsp)
	movsd	%xmm5, 248(%rsp)
	movsd	%xmm4, 240(%rsp)
	movsd	%xmm11, 232(%rsp)
	movsd	%xmm12, 224(%rsp)
	movsd	%xmm10, 216(%rsp)
	movsd	%xmm14, 208(%rsp)
	movsd	%xmm13, 64(%rsp)
	call	__muldc3@PLT
	movq	328(%rsp), %r8
	movq	320(%rsp), %r10
	movq	312(%rsp), %rdi
	movq	304(%rsp), %rsi
	movq	296(%rsp), %r9
	movq	288(%rsp), %rax
	movq	.LC8(%rip), %xmm15
	movsd	64(%rsp), %xmm13
	movsd	280(%rsp), %xmm9
	movsd	272(%rsp), %xmm8
	movsd	264(%rsp), %xmm7
	movsd	256(%rsp), %xmm6
	movsd	248(%rsp), %xmm5
	movsd	240(%rsp), %xmm4
	movsd	232(%rsp), %xmm11
	movsd	224(%rsp), %xmm12
	movsd	216(%rsp), %xmm10
	movsd	208(%rsp), %xmm14
	jmp	.L1447
.L1491:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm4, %xmm3
	movq	%r8, 312(%rsp)
	movapd	%xmm11, %xmm2
	movq	%r10, 304(%rsp)
	movq	%rdi, 296(%rsp)
	movq	%rsi, 288(%rsp)
	movq	%r9, 280(%rsp)
	movq	%rax, 272(%rsp)
	movsd	%xmm9, 264(%rsp)
	movsd	%xmm8, 256(%rsp)
	movsd	%xmm7, 248(%rsp)
	movsd	%xmm6, 240(%rsp)
	movsd	%xmm5, 232(%rsp)
	movsd	%xmm12, 224(%rsp)
	movsd	%xmm10, 216(%rsp)
	movsd	%xmm14, 208(%rsp)
	movsd	%xmm13, 64(%rsp)
	call	__muldc3@PLT
	movq	312(%rsp), %r8
	movq	304(%rsp), %r10
	movq	296(%rsp), %rdi
	movq	288(%rsp), %rsi
	movq	280(%rsp), %r9
	movq	272(%rsp), %rax
	movq	.LC8(%rip), %xmm15
	movsd	64(%rsp), %xmm13
	movsd	264(%rsp), %xmm9
	movsd	256(%rsp), %xmm8
	movsd	248(%rsp), %xmm7
	movsd	240(%rsp), %xmm6
	movsd	232(%rsp), %xmm5
	movsd	224(%rsp), %xmm12
	movsd	216(%rsp), %xmm10
	movsd	208(%rsp), %xmm14
	jmp	.L1448
.L1492:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm12, %xmm3
	movq	%r8, 296(%rsp)
	movapd	%xmm10, %xmm2
	movq	%r10, 288(%rsp)
	movq	%rdi, 280(%rsp)
	movq	%rsi, 272(%rsp)
	movq	%r9, 264(%rsp)
	movq	%rax, 256(%rsp)
	movsd	%xmm9, 248(%rsp)
	movsd	%xmm8, 240(%rsp)
	movsd	%xmm7, 232(%rsp)
	movsd	%xmm6, 224(%rsp)
	movsd	%xmm5, 216(%rsp)
	movsd	%xmm14, 208(%rsp)
	movsd	%xmm13, 64(%rsp)
	call	__muldc3@PLT
	movq	296(%rsp), %r8
	movq	288(%rsp), %r10
	movq	280(%rsp), %rdi
	movq	272(%rsp), %rsi
	movq	264(%rsp), %r9
	movq	256(%rsp), %rax
	movq	.LC8(%rip), %xmm15
	movsd	64(%rsp), %xmm13
	movsd	248(%rsp), %xmm9
	movsd	240(%rsp), %xmm8
	movsd	232(%rsp), %xmm7
	movsd	224(%rsp), %xmm6
	movsd	216(%rsp), %xmm5
	movsd	208(%rsp), %xmm14
	jmp	.L1449
.L1493:
	movsd	24(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm5, %xmm2
	movapd	%xmm13, %xmm0
	movq	%r8, 280(%rsp)
	movq	%r10, 272(%rsp)
	movq	%rdi, 264(%rsp)
	movq	%rsi, 256(%rsp)
	movq	%r9, 248(%rsp)
	movq	%rax, 240(%rsp)
	movsd	%xmm9, 232(%rsp)
	movsd	%xmm8, 224(%rsp)
	movsd	%xmm7, 216(%rsp)
	movsd	%xmm6, 208(%rsp)
	movsd	%xmm14, 64(%rsp)
	movsd	%xmm13, 24(%rsp)
	call	__muldc3@PLT
	movq	280(%rsp), %r8
	movq	272(%rsp), %r10
	movq	264(%rsp), %rdi
	movq	256(%rsp), %rsi
	movq	248(%rsp), %r9
	movq	240(%rsp), %rax
	movq	.LC8(%rip), %xmm15
	movsd	64(%rsp), %xmm14
	movsd	216(%rsp), %xmm7
	movsd	24(%rsp), %xmm13
	movsd	232(%rsp), %xmm9
	movsd	224(%rsp), %xmm8
	movsd	208(%rsp), %xmm6
	jmp	.L1450
.L1494:
	movsd	56(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm9, %xmm2
	movq	%r8, 264(%rsp)
	movq	%r10, 256(%rsp)
	movq	%rdi, 248(%rsp)
	movq	%rsi, 240(%rsp)
	movq	%r9, 232(%rsp)
	movq	%rax, 224(%rsp)
	movsd	%xmm8, 216(%rsp)
	movsd	%xmm7, 208(%rsp)
	movsd	%xmm6, 64(%rsp)
	movsd	%xmm14, 56(%rsp)
	movsd	%xmm13, 24(%rsp)
	call	__muldc3@PLT
	movq	264(%rsp), %r8
	movq	256(%rsp), %r10
	movq	248(%rsp), %rdi
	movq	240(%rsp), %rsi
	movq	232(%rsp), %r9
	movsd	64(%rsp), %xmm6
	movq	224(%rsp), %rax
	movsd	56(%rsp), %xmm14
	movq	.LC8(%rip), %xmm15
	movsd	24(%rsp), %xmm13
	movsd	216(%rsp), %xmm8
	movsd	208(%rsp), %xmm7
	jmp	.L1451
.L1495:
	movsd	48(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm8, %xmm2
	movq	%r8, 248(%rsp)
	movq	%r10, 240(%rsp)
	movq	%rdi, 232(%rsp)
	movq	%rsi, 224(%rsp)
	movq	%r9, 216(%rsp)
	movq	%rax, 208(%rsp)
	movsd	%xmm7, 64(%rsp)
	movsd	%xmm6, 56(%rsp)
	movsd	%xmm14, 48(%rsp)
	movsd	%xmm13, 24(%rsp)
	call	__muldc3@PLT
	movq	248(%rsp), %r8
	movq	240(%rsp), %r10
	movq	232(%rsp), %rdi
	movq	224(%rsp), %rsi
	movq	216(%rsp), %r9
	movq	208(%rsp), %rax
	movq	.LC8(%rip), %xmm15
	movsd	64(%rsp), %xmm7
	movsd	56(%rsp), %xmm6
	movsd	48(%rsp), %xmm14
	movsd	24(%rsp), %xmm13
	jmp	.L1452
.L1496:
	movsd	40(%rsp), %xmm3
	movapd	%xmm14, %xmm1
	movapd	%xmm7, %xmm2
	movapd	%xmm13, %xmm0
	movq	%r8, 232(%rsp)
	movq	%r10, 224(%rsp)
	movq	%rdi, 216(%rsp)
	movq	%rsi, 208(%rsp)
	movq	%r9, 64(%rsp)
	movq	%rax, 56(%rsp)
	movsd	%xmm6, 48(%rsp)
	movsd	%xmm14, 40(%rsp)
	movsd	%xmm13, 24(%rsp)
	call	__muldc3@PLT
	movq	232(%rsp), %r8
	movq	224(%rsp), %r10
	movq	216(%rsp), %rdi
	movq	208(%rsp), %rsi
	movq	64(%rsp), %r9
	movq	56(%rsp), %rax
	movq	.LC8(%rip), %xmm15
	movsd	48(%rsp), %xmm6
	movsd	40(%rsp), %xmm14
	movsd	24(%rsp), %xmm13
	jmp	.L1453
.L1505:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movapd	%xmm7, %xmm3
	movq	%r9, 56(%rsp)
	movapd	%xmm5, %xmm2
	movq	%rdi, 48(%rsp)
	movq	%r8, 40(%rsp)
	movq	%rax, 32(%rsp)
	movsd	%xmm14, 24(%rsp)
	movsd	%xmm13, 8(%rsp)
	call	__muldc3@PLT
	movq	56(%rsp), %r9
	movq	48(%rsp), %rdi
	movq	40(%rsp), %r8
	movq	32(%rsp), %rax
	movq	.LC8(%rip), %xmm4
	movsd	24(%rsp), %xmm14
	movsd	8(%rsp), %xmm13
	jmp	.L1470
.L1504:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movq	%r9, 72(%rsp)
	movq	%rdi, 64(%rsp)
	movq	%r8, 56(%rsp)
	movq	%rax, 48(%rsp)
	movsd	%xmm5, 40(%rsp)
	movsd	%xmm7, 32(%rsp)
	movsd	%xmm14, 24(%rsp)
	movsd	%xmm13, 8(%rsp)
	call	__muldc3@PLT
	movq	72(%rsp), %r9
	movq	64(%rsp), %rdi
	movq	56(%rsp), %r8
	movq	48(%rsp), %rax
	movq	.LC8(%rip), %xmm4
	movsd	40(%rsp), %xmm5
	movsd	32(%rsp), %xmm7
	movsd	24(%rsp), %xmm14
	movsd	8(%rsp), %xmm13
	jmp	.L1469
.L1506:
	movapd	%xmm14, %xmm1
	movapd	%xmm13, %xmm0
	movq	%rdi, 32(%rsp)
	movq	%rax, 24(%rsp)
	movsd	%xmm14, 16(%rsp)
	movsd	%xmm13, 8(%rsp)
	call	__muldc3@PLT
	movq	32(%rsp), %rdi
	movq	24(%rsp), %rax
	movq	.LC8(%rip), %xmm4
	movsd	16(%rsp), %xmm14
	movsd	8(%rsp), %xmm13
	jmp	.L1476
	.cfi_endproc
.LFE12542:
	.size	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb0ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_, .-_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb0ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_
	.section	.rodata.str1.8
	.align 8
.LC57:
	.string	"void Eigen::PlainObjectBase<Derived>::resize(Eigen::Index, Eigen::Index) [with Derived = Eigen::Matrix<std::complex<double>, -1, 1>; Eigen::Index = long int]"
	.align 8
.LC58:
	.ascii	"Eigen::MapBase<De"
	.string	"rived, 0>::MapBase(PointerType, Eigen::Index, Eigen::Index) [with Derived = Eigen::Block<Eigen::Transpose<Eigen::Block<Eigen::Matrix<std::complex<double>, -1, -1>, 1, -1, false> >, -1, 1, true>; PointerType = std::complex<double>*; Eigen::Index = long int]"
	.section	.text.unlikely
.LCOLDB59:
	.text
.LHOTB59:
	.p2align 4
	.type	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0, @function
_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0:
.LFB14910:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA14910
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rcx, %r8
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	movq	%rdx, %r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$104, %rsp
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	40(%rsi), %rbx
	movq	(%rsi), %rdx
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	xorl	%eax, %eax
	movq	16(%rsi), %r15
	movq	24(%rsi), %rcx
	movq	$0, -112(%rbp)
	testq	%rbx, %rbx
	js	.L1539
	movq	%rdi, %r13
	movl	$0, %r12d
	jne	.L1540
.L1509:
	movsd	8(%r8), %xmm5
	movsd	(%r8), %xmm0
	pxor	%xmm1, %xmm1
	movapd	%xmm5, %xmm2
	movapd	%xmm0, %xmm4
	mulsd	%xmm1, %xmm2
	mulsd	%xmm0, %xmm1
	subsd	%xmm2, %xmm4
	addsd	%xmm5, %xmm1
	ucomisd	%xmm4, %xmm1
	jp	.L1541
.L1515:
	pxor	%xmm0, %xmm0
	movapd	%xmm1, %xmm2
	movapd	%xmm4, %xmm6
	mulsd	%xmm0, %xmm2
	mulsd	%xmm4, %xmm0
	subsd	%xmm2, %xmm6
	addsd	%xmm1, %xmm0
	movsd	%xmm6, -128(%rbp)
	ucomisd	%xmm6, %xmm0
	movsd	%xmm0, -120(%rbp)
	jp	.L1542
.L1516:
	movabsq	$1152921504606846976, %rax
	cmpq	%rax, %rbx
	je	.L1543
	salq	$4, %rbx
	movq	%r12, %rax
	xorl	%r15d, %r15d
	testq	%r12, %r12
	je	.L1544
.L1520:
	cmpq	$0, 16(%r14)
	movq	(%r14), %r8
	jns	.L1522
	testq	%r8, %r8
	jne	.L1545
.L1522:
	movq	24(%r14), %r9
	movq	%rax, -80(%rbp)
	leaq	-80(%rbp), %rcx
	leaq	-96(%rbp), %rdx
	movq	8(%r13), %rsi
	movq	0(%r13), %rax
	movq	$1, -72(%rbp)
	movq	16(%r13), %rdi
	movsd	-128(%rbp), %xmm0
	movsd	-120(%rbp), %xmm1
	movq	8(%r9), %r9
	movq	%rax, -96(%rbp)
	movq	%rsi, -88(%rbp)
	call	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb0ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_
	cmpq	$131072, %rbx
	ja	.L1546
.L1523:
	movq	%r12, %rdi
	call	free@PLT
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L1547
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L1546:
	.cfi_restore_state
	movq	%r15, %rdi
	call	free@PLT
	jmp	.L1523
	.p2align 4,,10
	.p2align 3
.L1544:
	cmpq	$131072, %rbx
	ja	.L1521
	leaq	32(%rbx), %rax
	subq	%rax, %rsp
	leaq	15(%rsp), %rax
	andq	$-16, %rax
	movq	%rax, %r15
	jmp	.L1520
	.p2align 4,,10
	.p2align 3
.L1540:
	movabsq	$1152921504606846975, %rax
	cmpq	%rax, %rbx
	jg	.L1534
	movq	%rbx, %rdi
	movq	%r8, -136(%rbp)
	salq	$4, %rdi
	movq	%rcx, -128(%rbp)
	movq	%rdx, -120(%rbp)
	call	malloc@PLT
	movq	-120(%rbp), %rdx
	movq	-128(%rbp), %rcx
	movq	%rax, %r12
	andl	$15, %eax
	movq	-136(%rbp), %r8
	je	.L1511
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
	.p2align 4,,10
	.p2align 3
.L1521:
	movq	%rbx, %rdi
.LEHB58:
	call	_ZN5Eigen8internal14aligned_mallocEm
.LEHE58:
	movq	%rax, %r15
	jmp	.L1520
.L1511:
	testq	%r12, %r12
	je	.L1535
	imulq	8(%rdx), %r15
	movq	%r12, -112(%rbp)
	movapd	.LC7(%rip), %xmm1
	addq	%r15, %rcx
	salq	$4, %rcx
	addq	(%rdx), %rcx
	xorl	%edx, %edx
	.p2align 4,,10
	.p2align 3
.L1513:
	movupd	(%rcx,%rax), %xmm0
	addq	$1, %rdx
	xorpd	%xmm1, %xmm0
	movaps	%xmm0, (%r12,%rax)
	addq	$16, %rax
	cmpq	%rbx, %rdx
	jne	.L1513
	jmp	.L1509
.L1541:
	movsd	.LC24(%rip), %xmm2
	pxor	%xmm3, %xmm3
	movapd	%xmm5, %xmm1
	call	__muldc3@PLT
	movapd	%xmm0, %xmm4
	jmp	.L1515
.L1542:
	movsd	.LC24(%rip), %xmm2
	movapd	%xmm4, %xmm0
	pxor	%xmm3, %xmm3
	call	__muldc3@PLT
	movsd	%xmm0, -128(%rbp)
	movsd	%xmm1, -120(%rbp)
	jmp	.L1516
.L1545:
	leaq	.LC58(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L1539:
	leaq	.LC57(%rip), %rcx
	movl	$273, %edx
	leaq	.LC37(%rip), %rsi
	leaq	.LC38(%rip), %rdi
	call	__assert_fail@PLT
.L1547:
	call	__stack_chk_fail@PLT
.L1532:
	jmp	.L1533
.L1543:
	jmp	.L1517
	.section	.gcc_except_table
.LLSDA14910:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE14910-.LLSDACSB14910
.LLSDACSB14910:
	.uleb128 .LEHB58-.LFB14910
	.uleb128 .LEHE58-.LEHB58
	.uleb128 .L1532-.LFB14910
	.uleb128 0
.LLSDACSE14910:
	.text
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDAC14910
	.type	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0.cold, @function
_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0.cold:
.LFSB14910:
.L1535:
	.cfi_def_cfa 6, 16
	.cfi_offset 3, -56
	.cfi_offset 6, -16
	.cfi_offset 12, -48
	.cfi_offset 13, -40
	.cfi_offset 14, -32
	.cfi_offset 15, -24
.LEHB59:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE59:
.L1529:
	movq	-112(%rbp), %rdi
	movq	%rax, %rbx
	call	free@PLT
	movq	%rbx, %rdi
.LEHB60:
	call	_Unwind_Resume@PLT
.L1528:
.L1533:
	movq	%rax, %rbx
	movq	%r12, %rdi
	call	free@PLT
	movq	%rbx, %rdi
	call	_Unwind_Resume@PLT
.LEHE60:
.L1534:
.LEHB61:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE61:
.L1517:
.LEHB62:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE62:
	.cfi_endproc
.LFE14910:
	.section	.gcc_except_table
.LLSDAC14910:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSEC14910-.LLSDACSBC14910
.LLSDACSBC14910:
	.uleb128 .LEHB59-.LCOLDB59
	.uleb128 .LEHE59-.LEHB59
	.uleb128 .L1529-.LCOLDB59
	.uleb128 0
	.uleb128 .LEHB60-.LCOLDB59
	.uleb128 .LEHE60-.LEHB60
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB61-.LCOLDB59
	.uleb128 .LEHE61-.LEHB61
	.uleb128 .L1529-.LCOLDB59
	.uleb128 0
	.uleb128 .LEHB62-.LCOLDB59
	.uleb128 .LEHE62-.LEHB62
	.uleb128 .L1528-.LCOLDB59
	.uleb128 0
.LLSDACSEC14910:
	.section	.text.unlikely
	.text
	.size	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0, .-_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0
	.section	.text.unlikely
	.size	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0.cold, .-_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0.cold
.LCOLDE59:
	.text
.LHOTE59:
	.section	.text.unlikely
.LCOLDB60:
	.text
.LHOTB60:
	.p2align 4
	.type	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS5_IS7_Li1ELin1ELi1ELi1ELin1EEEEENS4_INS_5BlockIS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNSN_6ScalarE.isra.0, @function
_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS5_IS7_Li1ELin1ELi1ELi1ELin1EEEEENS4_INS_5BlockIS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNSN_6ScalarE.isra.0:
.LFB14914:
	.cfi_startproc
	pxor	%xmm4, %xmm4
	movapd	%xmm1, %xmm2
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movapd	%xmm0, %xmm5
	mulsd	%xmm4, %xmm2
	mulsd	%xmm0, %xmm4
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r14
	pushq	%r13
	pushq	%r12
	.cfi_offset 14, -24
	.cfi_offset 13, -32
	.cfi_offset 12, -40
	movq	%rdx, %r12
	pushq	%rbx
	subsd	%xmm2, %xmm5
	.cfi_offset 3, -48
	movq	%rdi, %rbx
	addsd	%xmm1, %xmm4
	subq	$64, %rsp
	movq	%fs:40, %rax
	movq	%rax, -40(%rbp)
	xorl	%eax, %eax
	ucomisd	%xmm4, %xmm5
	jp	.L1564
.L1549:
	pxor	%xmm1, %xmm1
	movapd	%xmm4, %xmm2
	movapd	%xmm5, %xmm0
	mulsd	%xmm1, %xmm2
	mulsd	%xmm5, %xmm1
	subsd	%xmm2, %xmm0
	addsd	%xmm4, %xmm1
	ucomisd	%xmm0, %xmm1
	jp	.L1565
.L1550:
	movq	8(%rsi), %rdx
	movq	%rdx, %rax
	shrq	$60, %rax
	jne	.L1555
	movq	(%rsi), %rax
	salq	$4, %rdx
	xorl	%r13d, %r13d
	movq	%rdx, %r14
	testq	%rax, %rax
	je	.L1566
.L1552:
	cmpq	$0, 16(%r12)
	movq	(%r12), %r8
	jns	.L1556
	testq	%r8, %r8
	jne	.L1567
.L1556:
	movq	%rax, -64(%rbp)
	movq	8(%rbx), %rsi
	leaq	-64(%rbp), %rcx
	leaq	-80(%rbp), %rdx
	movq	24(%r12), %r9
	movq	(%rbx), %rax
	movq	$1, -56(%rbp)
	movq	16(%rbx), %rdi
	movq	%rsi, -72(%rbp)
	movq	8(%r9), %r9
	movq	%rax, -80(%rbp)
	call	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi1EEELi1ELb0ES3_NS4_IS3_lLi0EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_
	cmpq	$131072, %r14
	ja	.L1568
.L1548:
	movq	-40(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L1569
	leaq	-32(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L1568:
	.cfi_restore_state
	movq	%r13, %rdi
	call	free@PLT
	jmp	.L1548
	.p2align 4,,10
	.p2align 3
.L1566:
	cmpq	$131072, %rdx
	ja	.L1553
	leaq	32(%rdx), %rax
	subq	%rax, %rsp
	leaq	15(%rsp), %r13
	andq	$-16, %r13
	movq	%r13, %rax
	jmp	.L1552
	.p2align 4,,10
	.p2align 3
.L1553:
	movq	%rdx, %rdi
	movsd	%xmm0, -96(%rbp)
	movsd	%xmm1, -88(%rbp)
	call	malloc@PLT
	movsd	-88(%rbp), %xmm1
	movsd	-96(%rbp), %xmm0
	testb	$15, %al
	movq	%rax, %r13
	je	.L1554
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L1567:
	leaq	.LC58(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
	.p2align 4,,10
	.p2align 3
.L1554:
	testq	%rax, %rax
	jne	.L1552
	jmp	.L1555
.L1565:
	pxor	%xmm3, %xmm3
	movapd	%xmm4, %xmm1
	movapd	%xmm5, %xmm0
	movq	%rsi, -88(%rbp)
	movsd	.LC24(%rip), %xmm2
	call	__muldc3@PLT
	movq	-88(%rbp), %rsi
	jmp	.L1550
.L1564:
	movsd	.LC24(%rip), %xmm2
	pxor	%xmm3, %xmm3
	movq	%rsi, -88(%rbp)
	call	__muldc3@PLT
	movq	-88(%rbp), %rsi
	movapd	%xmm0, %xmm5
	movapd	%xmm1, %xmm4
	jmp	.L1549
.L1569:
	call	__stack_chk_fail@PLT
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.type	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS5_IS7_Li1ELin1ELi1ELi1ELin1EEEEENS4_INS_5BlockIS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNSN_6ScalarE.isra.0.cold, @function
_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS5_IS7_Li1ELin1ELi1ELi1ELin1EEEEENS4_INS_5BlockIS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNSN_6ScalarE.isra.0.cold:
.LFSB14914:
.L1555:
	.cfi_def_cfa 6, 16
	.cfi_offset 3, -48
	.cfi_offset 6, -16
	.cfi_offset 12, -40
	.cfi_offset 13, -32
	.cfi_offset 14, -24
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
	.cfi_endproc
.LFE14914:
	.text
	.size	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS5_IS7_Li1ELin1ELi1ELi1ELin1EEEEENS4_INS_5BlockIS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNSN_6ScalarE.isra.0, .-_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS5_IS7_Li1ELin1ELi1ELi1ELin1EEEEENS4_INS_5BlockIS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNSN_6ScalarE.isra.0
	.section	.text.unlikely
	.size	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS5_IS7_Li1ELin1ELi1ELi1ELin1EEEEENS4_INS_5BlockIS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNSN_6ScalarE.isra.0.cold, .-_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS5_IS7_Li1ELin1ELi1ELi1ELin1EEEEENS4_INS_5BlockIS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNSN_6ScalarE.isra.0.cold
.LCOLDE60:
	.text
.LHOTE60:
	.section	.rodata._ZNK5Eigen8internal17product_evaluatorINS_7ProductINS2_INS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li1EEELi8ENS_10DenseShapeESH_S6_S6_E5coeffEll.str1.8,"aMS",@progbits,1
	.align 8
.LC61:
	.string	"Eigen::MapBase<Derived, 0>::MapBase(PointerType, Eigen::Index, Eigen::Index) [with Derived = Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, 1, -1, false>; PointerType = const std::complex<double>*; Eigen::Index = long int]"
	.align 8
.LC62:
	.string	"Eigen::Block<XprType, BlockRows, BlockCols, InnerPanel>::Block(XprType&, Eigen::Index) [with XprType = const Eigen::Matrix<std::complex<double>, -1, -1>; int BlockRows = 1; int BlockCols = -1; bool InnerPanel = false; Eigen::Index = long int]"
	.align 8
.LC63:
	.ascii	"Eigen::CwiseBinaryOp<BinaryOp, Lhs, Rhs>::CwiseBinaryOp(cons"
	.ascii	"t Lhs&, const Rhs&, const BinaryOp&) [with BinaryOp = Eigen:"
	.ascii	":internal::scalar_product_op<std::complex<double>, std::comp"
	.ascii	"lex<double> >; LhsType = const Eigen::Transpose<const Eigen:"
	.ascii	":Block<const Eigen::Matrix<std::complex<double>, -1, -1>, 1,"
	.ascii	" -1, false> >; RhsType = const Eigen::Block<co"
	.string	"nst Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true>; Lhs = Eigen::Transpose<const Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, 1, -1, false> >; Rhs = Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true>]"
	.align 8
.LC64:
	.string	"/usr/local/include/Eigen/src/Core/CwiseBinaryOp.h"
	.align 8
.LC65:
	.string	"aLhs.rows() == aRhs.rows() && aLhs.cols() == aRhs.cols()"
	.align 8
.LC66:
	.ascii	"typename Eigen::internal::traits<T>::Scalar Eigen::DenseBase"
	.ascii	"<Derived>::redux(const Func&) const [with BinaryOp = Eigen::"
	.ascii	"internal::scalar_sum_op<std::complex<double>, std::complex<d"
	.ascii	"ouble> >; Derived = Eigen::CwiseBinaryOp<Eigen::internal::sc"
	.ascii	"alar_product_op<std::complex<double>, std::complex<double> >"
	.ascii	", const E"
	.string	"igen::Transpose<const Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, 1, -1, false> >, const Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true> >; typename Eigen::internal::traits<T>::Scalar = std::complex<double>]"
	.align 8
.LC67:
	.string	"/usr/local/include/Eigen/src/Core/Redux.h"
	.align 8
.LC68:
	.string	"this->rows()>0 && this->cols()>0 && \"you are using an empty matrix\""
	.section	.text._ZNK5Eigen8internal17product_evaluatorINS_7ProductINS2_INS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li1EEELi8ENS_10DenseShapeESH_S6_S6_E5coeffEll,"axG",@progbits,_ZNK5Eigen8internal17product_evaluatorINS_7ProductINS2_INS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li1EEELi8ENS_10DenseShapeESH_S6_S6_E5coeffEll,comdat
	.align 2
	.p2align 4
	.weak	_ZNK5Eigen8internal17product_evaluatorINS_7ProductINS2_INS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li1EEELi8ENS_10DenseShapeESH_S6_S6_E5coeffEll
	.type	_ZNK5Eigen8internal17product_evaluatorINS_7ProductINS2_INS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li1EEELi8ENS_10DenseShapeESH_S6_S6_E5coeffEll, @function
_ZNK5Eigen8internal17product_evaluatorINS_7ProductINS2_INS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li1EEELi8ENS_10DenseShapeESH_S6_S6_E5coeffEll:
.LFB12964:
	.cfi_startproc
	pushq	%r13
	.cfi_def_cfa_offset 16
	.cfi_offset 13, -16
	pushq	%r12
	.cfi_def_cfa_offset 24
	.cfi_offset 12, -24
	pushq	%rbp
	.cfi_def_cfa_offset 32
	.cfi_offset 6, -32
	movq	%rsi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 40
	.cfi_offset 3, -40
	salq	$4, %rbp
	subq	$40, %rsp
	.cfi_def_cfa_offset 80
	movq	16(%rdi), %r8
	addq	(%rdi), %rbp
	je	.L1571
	testq	%r8, %r8
	js	.L1590
.L1571:
	testq	%rsi, %rsi
	js	.L1572
	movq	8(%rdi), %r13
	cmpq	%r13, %rsi
	jge	.L1572
	movq	24(%rdi), %rcx
	movq	%rdx, %rax
	movq	8(%rcx), %r12
	imulq	%r12, %rax
	salq	$4, %rax
	addq	(%rcx), %rax
	je	.L1574
	testq	%r12, %r12
	js	.L1591
.L1574:
	testq	%rdx, %rdx
	js	.L1575
	cmpq	16(%rcx), %rdx
	jge	.L1575
	cmpq	%r8, %r12
	jne	.L1592
	pxor	%xmm6, %xmm6
	movaps	%xmm6, (%rsp)
	testq	%r12, %r12
	jne	.L1593
.L1578:
	movsd	(%rsp), %xmm0
	movsd	8(%rsp), %xmm1
	addq	$40, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 40
	popq	%rbx
	.cfi_def_cfa_offset 32
	popq	%rbp
	.cfi_def_cfa_offset 24
	popq	%r12
	.cfi_def_cfa_offset 16
	popq	%r13
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L1593:
	.cfi_restore_state
	jle	.L1594
	movupd	0(%rbp), %xmm2
	movupd	(%rax), %xmm1
	movapd	%xmm2, %xmm0
	movapd	%xmm2, %xmm3
	movapd	%xmm1, %xmm4
	unpcklpd	%xmm2, %xmm0
	unpckhpd	%xmm2, %xmm3
	shufpd	$1, %xmm1, %xmm4
	mulpd	%xmm1, %xmm0
	mulpd	%xmm3, %xmm4
	movapd	%xmm0, %xmm3
	addpd	%xmm4, %xmm3
	subpd	%xmm4, %xmm0
	movapd	%xmm3, %xmm7
	unpckhpd	%xmm3, %xmm3
	ucomisd	%xmm3, %xmm0
	movsd	%xmm0, %xmm7
	movaps	%xmm7, (%rsp)
	jp	.L1595
.L1580:
	cmpq	$1, %r12
	je	.L1578
	salq	$4, %r13
	salq	$4, %r12
	leaq	16(%rax), %rbx
	addq	%r13, %rbp
	addq	%rax, %r12
	jmp	.L1583
	.p2align 4,,10
	.p2align 3
.L1581:
	addpd	(%rsp), %xmm4
	addq	$16, %rbx
	addq	%r13, %rbp
	movaps	%xmm4, (%rsp)
	cmpq	%r12, %rbx
	je	.L1578
.L1583:
	movupd	(%rbx), %xmm1
	movupd	0(%rbp), %xmm2
	movapd	%xmm1, %xmm0
	movapd	%xmm1, %xmm4
	movapd	%xmm2, %xmm3
	unpcklpd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm4
	shufpd	$1, %xmm2, %xmm3
	mulpd	%xmm2, %xmm0
	mulpd	%xmm3, %xmm4
	movapd	%xmm0, %xmm3
	addpd	%xmm4, %xmm3
	subpd	%xmm4, %xmm0
	movapd	%xmm3, %xmm4
	unpckhpd	%xmm3, %xmm3
	ucomisd	%xmm3, %xmm0
	movsd	%xmm0, %xmm4
	jnp	.L1581
	movapd	%xmm2, %xmm7
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm1
	addq	$16, %rbx
	unpckhpd	%xmm7, %xmm7
	addq	%r13, %rbp
	movapd	%xmm7, %xmm3
	call	__muldc3@PLT
	unpcklpd	%xmm1, %xmm0
	addpd	(%rsp), %xmm0
	movaps	%xmm0, (%rsp)
	cmpq	%r12, %rbx
	jne	.L1583
	jmp	.L1578
.L1591:
	leaq	.LC55(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L1590:
	leaq	.LC61(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L1572:
	leaq	.LC62(%rip), %rcx
	movl	$120, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC54(%rip), %rdi
	call	__assert_fail@PLT
.L1575:
	leaq	.LC56(%rip), %rcx
	movl	$120, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC54(%rip), %rdi
	call	__assert_fail@PLT
.L1595:
	movapd	%xmm2, %xmm7
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm1
	movq	%rax, 24(%rsp)
	unpckhpd	%xmm7, %xmm7
	movapd	%xmm7, %xmm3
	call	__muldc3@PLT
	movq	24(%rsp), %rax
	movapd	%xmm0, %xmm7
	unpcklpd	%xmm1, %xmm7
	movaps	%xmm7, (%rsp)
	jmp	.L1580
.L1592:
	leaq	.LC63(%rip), %rcx
	movl	$116, %edx
	leaq	.LC64(%rip), %rsi
	leaq	.LC65(%rip), %rdi
	call	__assert_fail@PLT
.L1594:
	leaq	.LC66(%rip), %rcx
	movl	$411, %edx
	leaq	.LC67(%rip), %rsi
	leaq	.LC68(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE12964:
	.size	_ZNK5Eigen8internal17product_evaluatorINS_7ProductINS2_INS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li1EEELi8ENS_10DenseShapeESH_S6_S6_E5coeffEll, .-_ZNK5Eigen8internal17product_evaluatorINS_7ProductINS2_INS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li1EEELi8ENS_10DenseShapeESH_S6_S6_E5coeffEll
	.section	.text._ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll,"axG",@progbits,_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll,comdat
	.align 2
	.p2align 4
	.weak	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll
	.type	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll, @function
_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll:
.LFB13330:
	.cfi_startproc
	movq	96(%rdi), %rcx
	imulq	%rdx, %rcx
	addq	56(%rdi), %rdx
	addq	%rsi, %rcx
	addq	64(%rdi), %rsi
	imulq	24(%rdi), %rsi
	salq	$4, %rcx
	addq	80(%rdi), %rcx
	movupd	(%rcx), %xmm0
	addq	%rsi, %rdx
	movapd	%xmm0, %xmm2
	salq	$4, %rdx
	addq	16(%rdi), %rdx
	unpcklpd	%xmm0, %xmm2
	movupd	(%rdx), %xmm1
	unpckhpd	%xmm0, %xmm0
	mulpd	%xmm1, %xmm2
	shufpd	$1, %xmm1, %xmm1
	mulpd	%xmm1, %xmm0
	movapd	%xmm2, %xmm1
	subpd	%xmm0, %xmm1
	addpd	%xmm2, %xmm0
	movsd	%xmm1, %xmm0
	movaps	%xmm0, -24(%rsp)
	movsd	-16(%rsp), %xmm1
	movsd	-24(%rsp), %xmm0
	ret
	.cfi_endproc
.LFE13330:
	.size	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll, .-_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll
	.section	.text._ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEELi1ELin1ELb1EEEEEEEKNSB_ISF_Lin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll,"axG",@progbits,_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEELi1ELin1ELb1EEEEEEEKNSB_ISF_Lin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll,comdat
	.align 2
	.p2align 4
	.weak	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEELi1ELin1ELb1EEEEEEEKNSB_ISF_Lin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll
	.type	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEELi1ELin1ELb1EEEEEEEKNSB_ISF_Lin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll, @function
_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEELi1ELin1ELb1EEEEEEEKNSB_ISF_Lin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll:
.LFB13368:
	.cfi_startproc
	movq	120(%rdi), %rcx
	imulq	%rdx, %rcx
	addq	56(%rdi), %rdx
	addq	%rsi, %rcx
	addq	88(%rdi), %rsi
	addq	64(%rdi), %rsi
	imulq	24(%rdi), %rsi
	salq	$4, %rcx
	addq	104(%rdi), %rcx
	movupd	(%rcx), %xmm0
	addq	%rsi, %rdx
	movapd	%xmm0, %xmm2
	salq	$4, %rdx
	addq	16(%rdi), %rdx
	unpcklpd	%xmm0, %xmm2
	movupd	(%rdx), %xmm1
	unpckhpd	%xmm0, %xmm0
	mulpd	%xmm1, %xmm2
	shufpd	$1, %xmm1, %xmm1
	mulpd	%xmm1, %xmm0
	movapd	%xmm2, %xmm1
	subpd	%xmm0, %xmm1
	addpd	%xmm2, %xmm0
	movsd	%xmm1, %xmm0
	movaps	%xmm0, -24(%rsp)
	movsd	-16(%rsp), %xmm1
	movsd	-24(%rsp), %xmm0
	ret
	.cfi_endproc
.LFE13368:
	.size	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEELi1ELin1ELb1EEEEEEEKNSB_ISF_Lin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll, .-_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEELi1ELin1ELb1EEEEEEEKNSB_ISF_Lin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll
	.section	.rodata._ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_.str1.8,"aMS",@progbits,1
	.align 8
.LC69:
	.ascii	"Eigen::Block<XprType, BlockRows, BlockCols, InnerPanel>::Blo"
	.ascii	"ck(XprType&, Eigen::Index) [with XprTyp"
	.string	"e = const Eigen::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std::complex<double> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >; int BlockRows = 1; int BlockCols = -1; bool InnerPanel = true; Eigen::Index = long int]"
	.align 8
.LC70:
	.ascii	"Eigen::CwiseBinaryOp<BinaryOp, Lhs, Rhs>::CwiseBinaryOp(cons"
	.ascii	"t Lhs&, const Rhs&, const BinaryOp&) [with BinaryOp = Eigen:"
	.ascii	":internal::scalar_product_op<std::complex<double>, std::comp"
	.ascii	"lex<double> >; LhsType = const Eigen::Transpose<const Eigen:"
	.ascii	":Block<const Eigen::CwiseUnaryOp<Eigen::internal::scalar_con"
	.ascii	"jugate_op<std::complex<double> >, const Eigen::Transpose<con"
	.ascii	"st Eigen::Matrix<std::complex<double>, -1, -1> > >, 1, -1, t"
	.ascii	"rue> >; RhsType = const Eigen::Block<const Eigen::Matrix<std"
	.ascii	"::complex<double>, -1, -1>, -1, 1, true>; Lhs = Eigen::Trans"
	.ascii	"pose<const Eigen::Block<const Ei"
	.string	"gen::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std::complex<double> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >, 1, -1, true> >; Rhs = Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true>]"
	.align 8
.LC71:
	.ascii	"typename Eigen::internal::traits<T>::Scalar Eigen::DenseBase"
	.ascii	"<Derived>::redux(const Func&) const [with BinaryOp = Eigen::"
	.ascii	"internal::scalar_sum_op<std::complex<double>, std::complex<d"
	.ascii	"ouble> >; Derived = Eigen::CwiseBinaryOp<Eigen::internal::sc"
	.ascii	"alar_product_op<std::complex<double>, std::complex<double> >"
	.ascii	", const Eigen::Transpose<const Eigen::Block<const Eigen::Cwi"
	.ascii	"seUnaryOp<Eigen::internal::scalar_conjugate_op<std::complex<"
	.ascii	"do"
	.string	"uble> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >, 1, -1, true> >, const Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true> >; typename Eigen::internal::traits<T>::Scalar = std::complex<double>]"
	.section	.text._ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_,"axG",@progbits,_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_
	.type	_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_, @function
_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_:
.LFB13533:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	movq	%rdi, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$72, %rsp
	.cfi_def_cfa_offset 128
	movq	(%rsi), %rax
	movq	16(%rsi), %rbx
	movq	8(%rdi), %rcx
	movq	16(%rax), %r12
	movq	16(%rbx), %r14
	movq	%rax, 32(%rsp)
	movq	%r12, %xmm0
	movq	%r14, %xmm3
	punpcklqdq	%xmm3, %xmm0
	cmpq	%rcx, %r12
	jne	.L1599
	cmpq	16(%rdi), %r14
	je	.L1600
.L1599:
	movq	%r12, %rax
	orq	%r14, %rax
	js	.L1645
	testq	%r12, %r12
	je	.L1602
	testq	%r14, %r14
	je	.L1602
	movabsq	$9223372036854775807, %rax
	cqto
	idivq	%r14
	cmpq	%rax, %r12
	jg	.L1644
	movq	%r12, %r13
	imulq	16(%rbp), %rcx
	imulq	%r14, %r13
	cmpq	%r13, %rcx
	jne	.L1646
.L1604:
	movups	%xmm0, 8(%rbp)
.L1600:
	movq	0(%rbp), %rax
	testq	%r14, %r14
	jle	.L1598
	testq	%r12, %r12
	jle	.L1598
	movq	32(%rsp), %rsi
	movq	%rax, 48(%rsp)
	xorl	%edx, %edx
	xorl	%r11d, %r11d
	movq	%r14, 56(%rsp)
	movapd	.LC7(%rip), %xmm4
	movq	16(%rsi), %rdi
	movq	8(%rsi), %r13
	movq	.LC8(%rip), %xmm3
	movq	%rdi, (%rsp)
	movq	%r13, %rdi
	salq	$4, %rdi
	movq	%rdi, 24(%rsp)
	.p2align 4,,10
	.p2align 3
.L1620:
	movq	48(%rsp), %rax
	movq	%rdx, %rcx
	movq	%rdx, 40(%rsp)
	xorl	%r9d, %r9d
	salq	$4, %rcx
	xorl	%ebp, %ebp
	xorl	%esi, %esi
	addq	%rax, %rcx
	jmp	.L1619
	.p2align 4,,10
	.p2align 3
.L1614:
	movq	24(%rsp), %rax
	addq	$1, %rsi
	movsd	%xmm1, (%rcx)
	addq	%r13, %rbp
	movsd	%xmm0, 8(%rcx)
	addq	$16, %rcx
	addq	%rax, %r9
	cmpq	%r12, %rsi
	je	.L1647
.L1619:
	cmpq	%rsi, (%rsp)
	jle	.L1648
	movq	8(%rbx), %rdx
	movq	%r11, %rax
	imulq	%rdx, %rax
	salq	$4, %rax
	addq	(%rbx), %rax
	testq	%rdx, %rdx
	jns	.L1611
	testq	%rax, %rax
	jne	.L1649
.L1611:
	cmpq	16(%rbx), %r11
	jge	.L1650
	cmpq	%r13, %rdx
	jne	.L1651
	pxor	%xmm1, %xmm1
	movapd	%xmm1, %xmm0
	testq	%rdx, %rdx
	je	.L1614
	jle	.L1652
	movq	32(%rsp), %rdi
	movdqu	(%rax), %xmm5
	movq	(%rdi), %r15
	pshufd	$78, %xmm5, %xmm1
	movupd	(%rax), %xmm5
	movupd	(%r15,%r9), %xmm0
	xorpd	%xmm4, %xmm0
	pshufd	$238, %xmm0, %xmm2
	pshufd	$68, %xmm0, %xmm0
	mulpd	%xmm2, %xmm1
	mulpd	%xmm5, %xmm0
	xorpd	%xmm3, %xmm1
	addpd	%xmm1, %xmm0
	cmpq	$1, %rdx
	je	.L1616
	movdqu	16(%rax), %xmm7
	movupd	16(%rax), %xmm6
	movq	%rdx, %rdi
	movq	%rdx, %r14
	movupd	16(%r15,%r9), %xmm5
	sarq	%rdi
	andq	$-2, %r14
	pshufd	$78, %xmm7, %xmm1
	xorpd	%xmm4, %xmm5
	pshufd	$238, %xmm5, %xmm2
	pshufd	$68, %xmm5, %xmm5
	mulpd	%xmm2, %xmm1
	mulpd	%xmm6, %xmm5
	xorpd	%xmm3, %xmm1
	addpd	%xmm1, %xmm5
	cmpq	$1, %rdi
	je	.L1617
	leaq	48(%rax), %r8
	leaq	48(%r15,%r9), %rdi
	movl	$2, %r10d
	.p2align 4,,10
	.p2align 3
.L1618:
	movupd	-16(%rdi), %xmm1
	addq	$2, %r10
	addq	$32, %r8
	addq	$32, %rdi
	movdqu	-48(%r8), %xmm7
	xorpd	%xmm4, %xmm1
	pshufd	$238, %xmm1, %xmm6
	pshufd	$78, %xmm7, %xmm2
	movupd	-48(%r8), %xmm7
	pshufd	$68, %xmm1, %xmm1
	mulpd	%xmm6, %xmm2
	mulpd	%xmm7, %xmm1
	movdqu	-32(%r8), %xmm7
	xorpd	%xmm3, %xmm2
	addpd	%xmm2, %xmm1
	pshufd	$78, %xmm7, %xmm2
	addpd	%xmm1, %xmm0
	movupd	-32(%rdi), %xmm1
	xorpd	%xmm4, %xmm1
	pshufd	$238, %xmm1, %xmm6
	pshufd	$68, %xmm1, %xmm1
	mulpd	%xmm6, %xmm2
	movupd	-32(%r8), %xmm6
	mulpd	%xmm6, %xmm1
	xorpd	%xmm3, %xmm2
	addpd	%xmm2, %xmm1
	addpd	%xmm1, %xmm5
	cmpq	%r10, %r14
	jg	.L1618
.L1617:
	addpd	%xmm5, %xmm0
	cmpq	%r14, %rdx
	jle	.L1616
	movq	%r14, %rdx
	addq	%rbp, %r14
	salq	$4, %r14
	salq	$4, %rdx
	movupd	(%r15,%r14), %xmm2
	movupd	(%rax,%rdx), %xmm1
	xorpd	%xmm4, %xmm2
	pshufd	$78, %xmm1, %xmm5
	pshufd	$238, %xmm2, %xmm6
	pshufd	$68, %xmm2, %xmm2
	mulpd	%xmm6, %xmm5
	mulpd	%xmm2, %xmm1
	xorpd	%xmm3, %xmm5
	addpd	%xmm5, %xmm1
	addpd	%xmm1, %xmm0
.L1616:
	movapd	%xmm0, %xmm1
	unpckhpd	%xmm0, %xmm0
	jmp	.L1614
	.p2align 4,,10
	.p2align 3
.L1647:
	movq	40(%rsp), %rdx
	movq	56(%rsp), %rax
	addq	$1, %r11
	addq	%r12, %rdx
	cmpq	%rax, %r11
	jne	.L1620
.L1598:
	addq	$72, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
.L1602:
	.cfi_restore_state
	movq	%r12, %r13
	imulq	16(%rbp), %rcx
	imulq	%r14, %r13
	cmpq	%rcx, %r13
	je	.L1604
	movq	0(%rbp), %rdi
	movaps	%xmm0, (%rsp)
	call	free@PLT
	testq	%r13, %r13
	movdqa	(%rsp), %xmm0
	jne	.L1621
	movq	$0, 0(%rbp)
	jmp	.L1604
.L1646:
	movq	0(%rbp), %rdi
	movaps	%xmm0, (%rsp)
	call	free@PLT
	movdqa	(%rsp), %xmm0
.L1621:
	movabsq	$1152921504606846975, %rax
	cmpq	%rax, %r13
	jg	.L1644
	movq	%r13, %rdi
	movaps	%xmm0, (%rsp)
	salq	$4, %rdi
	call	malloc@PLT
	movdqa	(%rsp), %xmm0
	testb	$15, %al
	je	.L1607
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L1607:
	testq	%rax, %rax
	je	.L1644
	movq	%rax, 0(%rbp)
	jmp	.L1604
.L1650:
	leaq	.LC56(%rip), %rcx
	movl	$120, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC54(%rip), %rdi
	call	__assert_fail@PLT
.L1648:
	leaq	.LC69(%rip), %rcx
	movl	$120, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC54(%rip), %rdi
	call	__assert_fail@PLT
.L1651:
	leaq	.LC70(%rip), %rcx
	movl	$116, %edx
	leaq	.LC64(%rip), %rsi
	leaq	.LC65(%rip), %rdi
	call	__assert_fail@PLT
.L1649:
	leaq	.LC55(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L1652:
	leaq	.LC71(%rip), %rcx
	movl	$411, %edx
	leaq	.LC67(%rip), %rsi
	leaq	.LC68(%rip), %rdi
	call	__assert_fail@PLT
.L1645:
	leaq	.LC36(%rip), %rcx
	movl	$273, %edx
	leaq	.LC37(%rip), %rsi
	leaq	.LC38(%rip), %rdi
	call	__assert_fail@PLT
.L1644:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
	.cfi_endproc
.LFE13533:
	.size	_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_, .-_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_
	.section	.text._ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC2ERSE_l,"axG",@progbits,_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC5ERSE_l,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC2ERSE_l
	.type	_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC2ERSE_l, @function
_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC2ERSE_l:
.LFB13565:
	.cfi_startproc
	movq	(%rsi), %rax
	movq	%rdx, 16(%rdi)
	movq	$0, 24(%rdi)
	movq	%rax, (%rdi)
	movq	(%rsi), %rax
	movq	8(%rax), %rcx
	movq	%rcx, 40(%rdi)
	testq	%rdx, %rdx
	js	.L1654
	cmpq	16(%rax), %rdx
	jge	.L1654
	ret
.L1654:
	pushq	%rax
	.cfi_def_cfa_offset 16
	leaq	.LC69(%rip), %rcx
	movl	$120, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC54(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE13565:
	.size	_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC2ERSE_l, .-_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC2ERSE_l
	.weak	_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC1ERSE_l
	.set	_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC1ERSE_l,_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC2ERSE_l
	.section	.text._ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l,"axG",@progbits,_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l,comdat
	.align 2
	.p2align 4
	.weak	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	.type	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l, @function
_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l:
.LFB14374:
	.cfi_startproc
	movq	%rsi, %rax
	addq	56(%rdi), %rsi
	salq	$4, %rsi
	addq	24(%rdi), %rsi
	salq	$4, %rax
	addq	64(%rdi), %rax
	movupd	(%rsi), %xmm1
	xorpd	.LC7(%rip), %xmm1
	movupd	(%rax), %xmm0
	pshufd	$238, %xmm1, %xmm3
	pshufd	$68, %xmm1, %xmm1
	pshufd	$78, %xmm0, %xmm2
	mulpd	%xmm1, %xmm0
	mulpd	%xmm3, %xmm2
	xorpd	.LC8(%rip), %xmm2
	addpd	%xmm2, %xmm0
	ret
	.cfi_endproc
.LFE14374:
	.size	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l, .-_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	.section	.text._ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l,"axG",@progbits,_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l,comdat
	.align 2
	.p2align 4
	.weak	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	.type	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l, @function
_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l:
.LFB14389:
	.cfi_startproc
	movq	%rsi, %rax
	salq	$4, %rax
	addq	88(%rdi), %rax
	movupd	(%rax), %xmm0
	movq	56(%rdi), %rax
	addq	80(%rdi), %rax
	addq	%rsi, %rax
	pshufd	$78, %xmm0, %xmm2
	salq	$4, %rax
	addq	24(%rdi), %rax
	movupd	(%rax), %xmm1
	xorpd	.LC7(%rip), %xmm1
	pshufd	$238, %xmm1, %xmm3
	pshufd	$68, %xmm1, %xmm1
	mulpd	%xmm1, %xmm0
	mulpd	%xmm3, %xmm2
	xorpd	.LC8(%rip), %xmm2
	addpd	%xmm2, %xmm0
	ret
	.cfi_endproc
.LFE14389:
	.size	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l, .-_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	.section	.rodata._ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_.str1.8,"aMS",@progbits,1
	.align 8
.LC72:
	.ascii	"static void Eigen::internal::generic_product_impl<Lhs, Rhs, "
	.ascii	"Eigen::DenseShape, Eigen::DenseShape, 8>::scaleAndAddTo(Dest"
	.ascii	"&, const Lhs&, const Rhs&, const Scalar&) [with Dest = Eigen"
	.ascii	"::Matrix<std::complex<double>,"
	.string	" -1, -1>; Lhs = Eigen::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std::complex<double> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >; Rhs = Eigen::Matrix<std::complex<double>, -1, -1>; Scalar = std::complex<double>]"
	.align 8
.LC73:
	.string	"/usr/local/include/Eigen/src/Core/products/GeneralMatrixMatrix.h"
	.align 8
.LC74:
	.string	"dst.rows()==a_lhs.rows() && dst.cols()==a_rhs.cols()"
	.align 8
.LC75:
	.ascii	"Eigen::MapBase<"
	.string	"Derived, 0>::MapBase(PointerType, Eigen::Index, Eigen::Index) [with Derived = Eigen::Block<const Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true>, -1, 1, true>; PointerType = const std::complex<double>*; Eigen::Index = long int]"
	.align 8
.LC76:
	.ascii	"typename Eigen::ScalarBinaryOpTraits<typename Eigen::interna"
	.ascii	"l::traits<T>::Scalar, typename Eigen::internal::traits<Other"
	.ascii	"Derived>::Scalar>::ReturnType Eigen::MatrixBase<Derived>::do"
	.ascii	"t(const Eigen::MatrixBase<OtherDerived>&) const [with OtherD"
	.ascii	"erived = Eigen::Block<const Eigen::Block<const Eigen::Matrix"
	.ascii	"<std::complex<double>, -1, -1>, -1, 1, true>, -1, 1, true>; "
	.ascii	"Derived = Eigen::CwiseUnaryOp<Eigen::internal::scalar_conjug"
	.ascii	"ate_op<std::complex<double> >, const Eigen::Block<const Eige"
	.ascii	"n::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std::co"
	.ascii	"mplex<double> >, const Eigen::Transpose<const Eigen::Matrix<"
	.ascii	"std::complex<double>, -1, -1> > >, 1, -1, true> >; typename "
	.ascii	"Eigen::ScalarBinaryOpTraits<typename Eigen::internal::t"
	.string	"raits<T>::Scalar, typename Eigen::internal::traits<OtherDerived>::Scalar>::ReturnType = std::complex<double>; typename Eigen::internal::traits<OtherDerived>::Scalar = std::complex<double>; typename Eigen::internal::traits<T>::Scalar = std::complex<double>]"
	.align 8
.LC77:
	.string	"/usr/local/include/Eigen/src/Core/Dot.h"
	.section	.rodata._ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_.str1.1,"aMS",@progbits,1
.LC78:
	.string	"size() == other.size()"
	.section	.rodata._ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_.str1.8
	.align 8
.LC79:
	.ascii	"typename Eigen::internal::traits<T>::Scalar Eigen::DenseBase"
	.ascii	"<Derived>::redux(const Func&) const [with BinaryOp = Eigen::"
	.ascii	"internal::scalar_sum_op<std::complex<double>, std::complex<d"
	.ascii	"ouble> >; Derived = Eigen::CwiseBinaryOp<Eigen::internal::sc"
	.ascii	"alar_conj_product_op<std::complex<double>, std::complex<doub"
	.ascii	"le> >, const Eigen::Transpose<const Eigen::CwiseUnaryOp<Eige"
	.ascii	"n::internal::scalar_conjugate_op<std::complex<double> >, con"
	.ascii	"st Eigen::Block<const Eigen::CwiseUnaryOp<Eigen::internal::s"
	.ascii	"calar_conjugate_op<std::complex<double> >, const Eigen::Tran"
	.ascii	"spose<con"
	.string	"st Eigen::Matrix<std::complex<double>, -1, -1> > >, 1, -1, true> > >, const Eigen::Block<const Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true>, -1, 1, true> >; typename Eigen::internal::traits<T>::Scalar = std::complex<double>]"
	.align 8
.LC80:
	.string	"Eigen::MapBase<Derived, 0>::MapBase(PointerType, Eigen::Index, Eigen::Index) [with Derived = Eigen::Block<Eigen::Matrix<std::complex<double>, -1, -1>, 1, -1, false>; PointerType = std::complex<double>*; Eigen::Index = long int]"
	.align 8
.LC81:
	.ascii	"typename Eigen::ScalarBinaryOpTraits<typename Eigen::interna"
	.ascii	"l::traits<T>::Scalar, typename Eigen::internal::traits<Other"
	.ascii	"Derived>::Scalar>::ReturnType Eigen::MatrixBase<Derived>::do"
	.ascii	"t(const Eigen::MatrixBase<OtherDerived>&) const [with OtherD"
	.ascii	"erived = Eigen::Block<const Eigen::Matrix<std::complex<doubl"
	.ascii	"e>, -1, -1>, -1, 1, true>; Derived = Eigen::CwiseUnaryOp<Eig"
	.ascii	"en::internal::scalar_conjugate_op<std::complex<double> >, co"
	.ascii	"nst Eigen::Block<const Eigen::Block<const Eigen::CwiseUnaryO"
	.ascii	"p<Eigen::internal::scalar_conjugate_op<std::complex<double> "
	.ascii	">, const Eigen::Transpose<const Eigen::Matrix<std::complex<d"
	.ascii	"ouble>, -1, -1> > >, 1, -1, true>, 1, -1, true> >; typename "
	.ascii	"Eigen::ScalarBinaryOpTraits<typename Eigen::internal::t"
	.string	"raits<T>::Scalar, typename Eigen::internal::traits<OtherDerived>::Scalar>::ReturnType = std::complex<double>; typename Eigen::internal::traits<OtherDerived>::Scalar = std::complex<double>; typename Eigen::internal::traits<T>::Scalar = std::complex<double>]"
	.align 8
.LC82:
	.ascii	"typename Eigen::internal::traits<T>::Scalar Eigen::DenseBase"
	.ascii	"<Derived>::redux(const Func&) const [with BinaryOp = Eigen::"
	.ascii	"internal::scalar_sum_op<std::complex<double>, std::complex<d"
	.ascii	"ouble> >; Derived = Eigen::CwiseBinaryOp<Eigen::internal::sc"
	.ascii	"alar_conj_product_op<std::complex<double>, std::complex<doub"
	.ascii	"le> >, const Eigen::Transpose<const Eigen::CwiseUnaryOp<Eige"
	.ascii	"n::internal::scalar_conjugate_op<std::complex<double> >, con"
	.ascii	"st Eigen::Block<const Eigen::Block<const Eigen::CwiseUnaryOp"
	.ascii	"<Eigen::internal::scalar_conjugate_op<std::complex<double> >"
	.ascii	", const E"
	.string	"igen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >, 1, -1, true>, 1, -1, true> > >, const Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true> >; typename Eigen::internal::traits<T>::Scalar = std::complex<double>]"
	.section	.text._ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_,"axG",@progbits,_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_
	.type	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_, @function
_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_:
.LFB13241:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA13241
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$360, %rsp
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	(%rsi), %r9
	movq	8(%rdi), %r13
	movq	%rdx, -368(%rbp)
	movq	16(%r9), %r8
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	xorl	%eax, %eax
	cmpq	%r8, %r13
	jne	.L1661
	movq	16(%rdi), %r14
	movq	16(%rdx), %rdx
	movq	%rdi, %rbx
	cmpq	%rdx, %r14
	jne	.L1661
	movq	8(%r9), %r15
	movq	%rsi, %r12
	testq	%r15, %r15
	sete	%al
	testq	%r13, %r13
	sete	%sil
	orb	%sil, %al
	jne	.L1660
	testq	%r14, %r14
	jne	.L1714
.L1660:
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L1713
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L1714:
	.cfi_restore_state
	cmpq	$1, %r14
	je	.L1715
	cmpq	$1, %r13
	je	.L1716
	movsd	(%rcx), %xmm5
	movsd	8(%rcx), %xmm1
	movsd	.LC23(%rip), %xmm3
	movapd	%xmm1, %xmm2
	movapd	%xmm5, %xmm4
	movapd	%xmm5, %xmm0
	mulsd	%xmm3, %xmm2
	mulsd	%xmm3, %xmm4
	subsd	%xmm2, %xmm0
	addsd	%xmm1, %xmm4
	ucomisd	%xmm4, %xmm0
	jp	.L1717
.L1686:
	pxor	%xmm1, %xmm1
	movapd	%xmm4, %xmm2
	movapd	%xmm0, %xmm3
	mulsd	%xmm1, %xmm2
	mulsd	%xmm0, %xmm1
	subsd	%xmm2, %xmm3
	addsd	%xmm4, %xmm1
	ucomisd	%xmm1, %xmm3
	jp	.L1718
.L1687:
	movq	%r13, %xmm2
	movq	%r14, %xmm6
	pxor	%xmm0, %xmm0
	movl	$1, %ecx
	punpcklqdq	%xmm6, %xmm2
	leaq	-152(%rbp), %rdx
	leaq	-160(%rbp), %rsi
	movq	%r9, -384(%rbp)
	leaq	-144(%rbp), %rdi
	movsd	%xmm1, -400(%rbp)
	leaq	-176(%rbp), %r13
	movsd	%xmm3, -392(%rbp)
	movaps	%xmm0, -176(%rbp)
	movq	%r15, -144(%rbp)
	movaps	%xmm2, -160(%rbp)
	call	_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_
	movq	-144(%rbp), %rax
	subq	$8, %rsp
	movq	-160(%rbp), %rdx
	movq	-384(%rbp), %r9
	movq	-368(%rbp), %r10
	imulq	%rax, %rdx
	movsd	-392(%rbp), %xmm3
	movsd	-400(%rbp), %xmm1
	imulq	-152(%rbp), %rax
	movq	(%r9), %rcx
	movq	16(%r10), %rsi
	movsd	%xmm3, -216(%rbp)
	movapd	%xmm3, %xmm0
	movq	%rdx, -136(%rbp)
	movq	8(%r9), %rdx
	movq	%rax, -128(%rbp)
	movq	(%r12), %rax
	movsd	%xmm1, -208(%rbp)
	movq	%rdx, %r8
	movq	16(%rax), %rdi
	pushq	%r13
	pushq	8(%rbx)
	pushq	$1
	pushq	(%rbx)
	pushq	8(%r10)
	movq	(%r10), %r9
.LEHB63:
	.cfi_escape 0x2e,0x30
	call	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0
.LEHE63:
	movq	-176(%rbp), %rdi
	addq	$48, %rsp
	call	free@PLT
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L1713
	movq	-168(%rbp), %rdi
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	jmp	free@PLT
	.p2align 4,,10
	.p2align 3
.L1716:
	.cfi_restore_state
	movq	(%rdi), %r13
	testq	%r13, %r13
	je	.L1678
	testq	%r14, %r14
	jns	.L1678
	leaq	.LC80(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
	.p2align 4,,10
	.p2align 3
.L1678:
	movq	(%r12), %rax
	movq	%rax, -352(%rbp)
	cmpq	$1, %rdx
	je	.L1719
	movq	-368(%rbp), %rdi
	pxor	%xmm0, %xmm0
	leaq	-176(%rbp), %rdx
	leaq	-240(%rbp), %rsi
	movq	%rax, -304(%rbp)
	movq	%r15, -264(%rbp)
	movq	%r13, -176(%rbp)
	movq	%r14, -160(%rbp)
	movq	%rbx, -152(%rbp)
	movq	$1, -128(%rbp)
	movq	%rax, -240(%rbp)
	movq	$0, -224(%rbp)
	movq	$0, -216(%rbp)
	movq	%r15, -200(%rbp)
	movaps	%xmm0, -288(%rbp)
	movaps	%xmm0, -144(%rbp)
.LEHB64:
	call	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0
	jmp	.L1660
	.p2align 4,,10
	.p2align 3
.L1719:
	movq	-368(%rbp), %rsi
	xorl	%edx, %edx
	leaq	-240(%rbp), %rdi
	movq	%rcx, -392(%rbp)
	movq	%r9, -384(%rbp)
	call	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS5_l
	movq	-232(%rbp), %r12
	movq	-384(%rbp), %r9
	movq	-392(%rbp), %rcx
	cmpq	%r12, %r15
	jne	.L1720
	movq	-240(%rbp), %rdx
	movq	-216(%rbp), %rax
	testq	%r15, %r15
	jle	.L1721
	movq	(%r9), %rsi
	movq	8(%rax), %rax
	pxor	%xmm0, %xmm0
	movq	%rdx, -88(%rbp)
	leaq	-176(%rbp), %r15
	movups	%xmm0, -104(%rbp)
	movq	%rsi, -152(%rbp)
	movq	%r15, %rdi
	xorl	%esi, %esi
	movq	%rcx, -368(%rbp)
	movq	%r12, -144(%rbp)
	movq	$0, -136(%rbp)
	movq	$0, -128(%rbp)
	movq	$0, -120(%rbp)
	movq	%rax, -72(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	cmpq	$1, %r12
	movq	-368(%rbp), %rcx
	movapd	%xmm0, %xmm2
	je	.L1681
	movq	%r12, %rbx
	movq	%r12, %r14
	movl	$1, %esi
	movq	%r15, %rdi
	sarq	%rbx
	movq	%rcx, -384(%rbp)
	andq	$-2, %r14
	movaps	%xmm0, -368(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	cmpq	$1, %rbx
	movapd	-368(%rbp), %xmm2
	movq	-384(%rbp), %rcx
	movapd	%xmm0, %xmm1
	je	.L1682
	movl	$2, %ebx
.L1683:
	movq	%rbx, %rsi
	movq	%r15, %rdi
	movq	%rcx, -392(%rbp)
	movaps	%xmm1, -384(%rbp)
	movaps	%xmm2, -368(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	-368(%rbp), %xmm2
	leaq	1(%rbx), %rsi
	movq	%r15, %rdi
	addq	$2, %rbx
	addpd	%xmm0, %xmm2
	movaps	%xmm2, -368(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	-384(%rbp), %xmm1
	cmpq	%rbx, %r14
	movapd	-368(%rbp), %xmm2
	movq	-392(%rbp), %rcx
	addpd	%xmm0, %xmm1
	jg	.L1683
.L1682:
	addpd	%xmm1, %xmm2
	cmpq	%r14, %r12
	jle	.L1681
	movq	%r14, %rsi
	movq	%r15, %rdi
	movq	%rcx, -384(%rbp)
	movaps	%xmm2, -368(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEELi1ELin1ELb1EEEEEEEKNSB_ISE_Lin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	-368(%rbp), %xmm2
	movq	-384(%rbp), %rcx
	addpd	%xmm0, %xmm2
.L1681:
	movupd	(%rcx), %xmm1
	movapd	%xmm2, %xmm3
	shufpd	$1, %xmm2, %xmm3
	movapd	%xmm1, %xmm4
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm4
	unpcklpd	%xmm1, %xmm0
	mulpd	%xmm3, %xmm4
	mulpd	%xmm2, %xmm0
	movapd	%xmm4, %xmm3
	addpd	%xmm0, %xmm3
	subpd	%xmm4, %xmm0
	movapd	%xmm3, %xmm4
	unpckhpd	%xmm3, %xmm3
	ucomisd	%xmm0, %xmm3
	movsd	%xmm0, %xmm4
	jp	.L1722
.L1684:
	movupd	0(%r13), %xmm0
	addpd	%xmm4, %xmm0
	movups	%xmm0, 0(%r13)
	jmp	.L1660
	.p2align 4,,10
	.p2align 3
.L1715:
	movq	(%rdi), %r14
	movq	%r13, -296(%rbp)
	movq	%r14, -304(%rbp)
	testq	%r14, %r14
	je	.L1665
	testq	%r13, %r13
	js	.L1723
.L1665:
	movq	-368(%rbp), %rsi
	pxor	%xmm0, %xmm0
	xorl	%edx, %edx
	leaq	-240(%rbp), %rdi
	movq	%r8, -256(%rbp)
	movq	%rcx, -400(%rbp)
	movq	%r9, -392(%rbp)
	movq	%r8, -384(%rbp)
	movq	%rbx, -280(%rbp)
	movaps	%xmm0, -272(%rbp)
	call	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS5_l
	movq	-384(%rbp), %r8
	movq	-392(%rbp), %r9
	movq	-400(%rbp), %rcx
	cmpq	$1, %r8
	je	.L1724
	movq	(%r12), %rax
	movsd	(%rcx), %xmm0
	leaq	-304(%rbp), %rdx
	leaq	-176(%rbp), %rsi
	movdqa	-240(%rbp), %xmm7
	movsd	8(%rcx), %xmm1
	leaq	-352(%rbp), %rdi
	movq	%rax, -352(%rbp)
	movq	-192(%rbp), %rax
	movaps	%xmm7, -176(%rbp)
	movdqa	-224(%rbp), %xmm7
	movq	%rax, -128(%rbp)
	movaps	%xmm7, -160(%rbp)
	movdqa	-208(%rbp), %xmm7
	movaps	%xmm7, -144(%rbp)
	call	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0
	jmp	.L1660
.L1724:
	movq	-216(%rbp), %rax
	movq	-232(%rbp), %rbx
	movq	-240(%rbp), %rdx
	movq	8(%rax), %rax
	testq	%rbx, %rbx
	jns	.L1669
	testq	%rdx, %rdx
	jne	.L1725
.L1669:
	cmpq	%rbx, %r15
	jne	.L1726
	testq	%r15, %r15
	jle	.L1727
	movq	(%r9), %rsi
	leaq	-176(%rbp), %r15
	movq	%rcx, -368(%rbp)
	movq	%r15, %rdi
	movq	%rbx, -144(%rbp)
	movq	%rsi, -152(%rbp)
	xorl	%esi, %esi
	movq	$0, -136(%rbp)
	movq	$0, -128(%rbp)
	movq	$0, -120(%rbp)
	movq	%rdx, -112(%rbp)
	movq	%rax, -96(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	cmpq	$1, %rbx
	movq	-368(%rbp), %rcx
	movapd	%xmm0, %xmm2
	je	.L1672
	movq	%rbx, %r12
	movq	%rbx, %r13
	movl	$1, %esi
	movq	%r15, %rdi
	sarq	%r12
	movq	%rcx, -384(%rbp)
	andq	$-2, %r13
	movaps	%xmm0, -368(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	subq	$1, %r12
	movapd	-368(%rbp), %xmm2
	movq	-384(%rbp), %rcx
	movapd	%xmm0, %xmm1
	je	.L1673
	movl	$2, %r12d
.L1674:
	movq	%r12, %rsi
	movq	%r15, %rdi
	movq	%rcx, -392(%rbp)
	movaps	%xmm1, -384(%rbp)
	movaps	%xmm2, -368(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	-368(%rbp), %xmm2
	leaq	1(%r12), %rsi
	movq	%r15, %rdi
	addq	$2, %r12
	addpd	%xmm0, %xmm2
	movaps	%xmm2, -368(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	-384(%rbp), %xmm1
	cmpq	%r12, %r13
	movapd	-368(%rbp), %xmm2
	movq	-392(%rbp), %rcx
	addpd	%xmm0, %xmm1
	jg	.L1674
.L1673:
	addpd	%xmm1, %xmm2
	cmpq	%r13, %rbx
	jle	.L1672
	movq	%r13, %rsi
	movq	%r15, %rdi
	movq	%rcx, -384(%rbp)
	movaps	%xmm2, -368(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	-368(%rbp), %xmm2
	movq	-384(%rbp), %rcx
	addpd	%xmm0, %xmm2
.L1672:
	movupd	(%rcx), %xmm1
	movapd	%xmm2, %xmm3
	shufpd	$1, %xmm2, %xmm3
	movapd	%xmm1, %xmm4
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm4
	unpcklpd	%xmm1, %xmm0
	mulpd	%xmm3, %xmm4
	mulpd	%xmm2, %xmm0
	movapd	%xmm4, %xmm3
	addpd	%xmm0, %xmm3
	subpd	%xmm4, %xmm0
	movapd	%xmm3, %xmm4
	unpckhpd	%xmm3, %xmm3
	ucomisd	%xmm0, %xmm3
	movsd	%xmm0, %xmm4
	jp	.L1728
.L1675:
	movupd	(%r14), %xmm0
	addpd	%xmm4, %xmm0
	movups	%xmm0, (%r14)
	jmp	.L1660
.L1661:
	leaq	.LC72(%rip), %rcx
	movl	$470, %edx
	leaq	.LC73(%rip), %rsi
	leaq	.LC74(%rip), %rdi
	call	__assert_fail@PLT
.L1713:
	call	__stack_chk_fail@PLT
.L1718:
	movsd	.LC24(%rip), %xmm2
	pxor	%xmm3, %xmm3
	movapd	%xmm4, %xmm1
	movq	%r9, -384(%rbp)
	call	__muldc3@PLT
	movq	-384(%rbp), %r9
	movapd	%xmm0, %xmm3
	jmp	.L1687
.L1717:
	movsd	.LC24(%rip), %xmm2
	movapd	%xmm5, %xmm0
	movq	%r9, -384(%rbp)
	call	__muldc3@PLT
	movq	-384(%rbp), %r9
	movapd	%xmm1, %xmm4
	jmp	.L1686
.L1722:
	movapd	%xmm2, %xmm7
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm1
	unpckhpd	%xmm7, %xmm7
	movapd	%xmm7, %xmm3
	call	__muldc3@PLT
	movapd	%xmm0, %xmm4
	unpcklpd	%xmm1, %xmm4
	jmp	.L1684
.L1728:
	movapd	%xmm2, %xmm7
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm1
	unpckhpd	%xmm7, %xmm7
	movapd	%xmm7, %xmm3
	call	__muldc3@PLT
	movapd	%xmm0, %xmm4
	unpcklpd	%xmm1, %xmm4
	jmp	.L1675
.L1721:
	leaq	.LC82(%rip), %rcx
	movl	$411, %edx
	leaq	.LC67(%rip), %rsi
	leaq	.LC68(%rip), %rdi
	call	__assert_fail@PLT
.L1720:
	leaq	.LC81(%rip), %rcx
	movl	$82, %edx
	leaq	.LC77(%rip), %rsi
	leaq	.LC78(%rip), %rdi
	call	__assert_fail@PLT
.L1727:
	leaq	.LC79(%rip), %rcx
	movl	$411, %edx
	leaq	.LC67(%rip), %rsi
	leaq	.LC68(%rip), %rdi
	call	__assert_fail@PLT
.L1726:
	leaq	.LC76(%rip), %rcx
	movl	$82, %edx
	leaq	.LC77(%rip), %rsi
	leaq	.LC78(%rip), %rdi
	call	__assert_fail@PLT
.L1723:
	leaq	.LC52(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L1725:
	leaq	.LC75(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L1693:
	movq	%rax, %rbx
.L1689:
	movq	-176(%rbp), %rdi
	call	free@PLT
	movq	-168(%rbp), %rdi
	call	free@PLT
	movq	%rbx, %rdi
	call	_Unwind_Resume@PLT
.LEHE64:
	.cfi_endproc
.LFE13241:
	.section	.gcc_except_table
.LLSDA13241:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE13241-.LLSDACSB13241
.LLSDACSB13241:
	.uleb128 .LEHB63-.LFB13241
	.uleb128 .LEHE63-.LEHB63
	.uleb128 .L1693-.LFB13241
	.uleb128 0
	.uleb128 .LEHB64-.LFB13241
	.uleb128 .LEHE64-.LEHB64
	.uleb128 0
	.uleb128 0
.LLSDACSE13241:
	.section	.text._ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_,"axG",@progbits,_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_,comdat
	.size	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_, .-_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_
	.section	.rodata._ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE6evalToIS9_EEvRT_RKSD_RSA_.str1.8,"aMS",@progbits,1
	.align 8
.LC83:
	.ascii	"Eigen::Product<Lhs, Rhs, Option>::Product(const Lhs&, const "
	.ascii	"Rhs&) [with _Lhs = Eigen::CwiseUnaryOp<Eigen::internal::scal"
	.ascii	"ar_conjugate_op<std::complex<double> >, const Eigen::Transpo"
	.ascii	"se<const Eigen::Matrix<std::complex<double>, -1, -1> > >; _R"
	.ascii	"hs = Eigen::Matrix<std::c"
	.string	"omplex<double>, -1, -1>; int Option = 1; Lhs = Eigen::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std::complex<double> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >; Rhs = Eigen::Matrix<std::complex<double>, -1, -1>]"
	.align 8
.LC84:
	.string	"/usr/local/include/Eigen/src/Core/Product.h"
	.align 8
.LC85:
	.string	"lhs.cols() == rhs.rows() && \"invalid matrix product\" && \"if you wanted a coeff-wise or a dot product use the respective explicit functions\""
	.section	.text._ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE6evalToIS9_EEvRT_RKSD_RSA_,"axG",@progbits,_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE6evalToIS9_EEvRT_RKSD_RSA_,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE6evalToIS9_EEvRT_RKSD_RSA_
	.type	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE6evalToIS9_EEvRT_RKSD_RSA_, @function
_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE6evalToIS9_EEvRT_RKSD_RSA_:
.LFB13091:
	.cfi_startproc
	pushq	%r12
	.cfi_def_cfa_offset 16
	.cfi_offset 12, -16
	movq	%rsi, %r12
	pushq	%rbp
	.cfi_def_cfa_offset 24
	.cfi_offset 6, -24
	movq	%rdx, %rbp
	pushq	%rbx
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	movq	%rdi, %rbx
	subq	$48, %rsp
	.cfi_def_cfa_offset 80
	movq	8(%rdx), %rsi
	movq	8(%rdi), %rdx
	movq	%fs:40, %rax
	movq	%rax, 40(%rsp)
	xorl	%eax, %eax
	movq	16(%rdi), %rcx
	leaq	(%rsi,%rdx), %rax
	addq	%rcx, %rax
	cmpq	$19, %rax
	jg	.L1730
	testq	%rsi, %rsi
	jg	.L1749
.L1730:
	movq	%rdx, %rax
	orq	%rcx, %rax
	js	.L1750
	imulq	%rcx, %rdx
	testq	%rdx, %rdx
	je	.L1734
	salq	$4, %rdx
	movq	(%rbx), %rdi
	je	.L1734
	xorl	%esi, %esi
	call	memset@PLT
	.p2align 4,,10
	.p2align 3
.L1734:
	leaq	16(%rsp), %rcx
	movq	%rbp, %rdx
	movq	%r12, %rsi
	movq	%rbx, %rdi
	movq	.LC86(%rip), %xmm0
	movaps	%xmm0, 16(%rsp)
	call	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_
.L1729:
	movq	40(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L1751
	addq	$48, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 32
	popq	%rbx
	.cfi_def_cfa_offset 24
	popq	%rbp
	.cfi_def_cfa_offset 16
	popq	%r12
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L1749:
	.cfi_restore_state
	movq	(%r12), %rax
	movq	%rbp, 32(%rsp)
	movq	%rax, 16(%rsp)
	cmpq	8(%rax), %rsi
	jne	.L1752
	leaq	15(%rsp), %rdx
	leaq	16(%rsp), %rsi
	call	_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_
	jmp	.L1729
.L1750:
	call	_ZN5Eigen9DenseBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8ConstantEllRKS3_.part.0
.L1751:
	call	__stack_chk_fail@PLT
.L1752:
	leaq	.LC83(%rip), %rcx
	movl	$96, %edx
	leaq	.LC84(%rip), %rsi
	leaq	.LC85(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE13091:
	.size	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE6evalToIS9_EEvRT_RKSD_RSA_, .-_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE6evalToIS9_EEvRT_RKSD_RSA_
	.section	.text._ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_,"axG",@progbits,_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC5ERKSF_,comdat
	.align 2
	.p2align 4
	.weak	_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_
	.type	_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_, @function
_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_:
.LFB14242:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA14242
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pxor	%xmm0, %xmm0
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$104, %rsp
	.cfi_def_cfa_offset 160
	movq	16(%rsi), %r14
	movq	(%rsi), %r15
	movq	%fs:40, %rax
	movq	%rax, 88(%rsp)
	xorl	%eax, %eax
	movq	$0, (%rdi)
	movq	16(%r14), %r12
	movq	16(%r15), %rbp
	movq	$-1, 8(%rdi)
	movq	$0, 16(%rdi)
	movq	%r12, %rax
	movq	%rbp, %xmm1
	movq	%r12, %xmm2
	movups	%xmm0, 24(%rdi)
	orq	%rbp, %rax
	punpcklqdq	%xmm2, %xmm1
	js	.L1788
	movq	%rdi, %rbx
	movq	%rsi, %r13
	leaq	16(%rdi), %r8
	testq	%rbp, %rbp
	je	.L1755
	testq	%r12, %r12
	je	.L1755
	movabsq	$9223372036854775807, %rax
	cqto
	idivq	%r12
	cmpq	%rax, %rbp
	jg	.L1789
	movq	%r12, %rdi
	imulq	%rbp, %rdi
.L1757:
	movabsq	$1152921504606846975, %rax
	cmpq	%rax, %rdi
	jg	.L1790
	salq	$4, %rdi
	movl	$1, %esi
	movq	%r8, 8(%rsp)
	movaps	%xmm1, 16(%rsp)
	call	calloc@PLT
	movq	8(%rsp), %r8
	movdqa	16(%rsp), %xmm1
	testb	$15, %al
	je	.L1760
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
	.p2align 4,,10
	.p2align 3
.L1755:
	movq	%r12, %rdi
	imulq	%rbp, %rdi
	testq	%rdi, %rdi
	jne	.L1757
	movups	%xmm1, 24(%rbx)
	movq	8(%r14), %rax
	movq	%rbp, 8(%rbx)
	addq	%rax, %rbp
	movq	$0, (%rbx)
	addq	%r12, %rbp
	cmpq	$19, %rbp
	jg	.L1764
	testq	%rax, %rax
	jle	.L1764
.L1762:
	movq	%r15, 64(%rsp)
	movq	%r14, 80(%rsp)
	cmpq	%rax, 8(%r15)
	jne	.L1791
	leaq	47(%rsp), %rdx
	leaq	64(%rsp), %rsi
	movq	%r8, %rdi
.LEHB65:
	call	_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_
.L1753:
	movq	88(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L1792
	addq	$104, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
.L1760:
	.cfi_restore_state
	testq	%rax, %rax
	je	.L1793
	movq	%rax, 16(%rbx)
	movq	%rax, (%rbx)
	movups	%xmm1, 24(%rbx)
	movq	8(%r14), %rax
	movq	%rbp, 8(%rbx)
	addq	%rax, %rbp
	addq	%r12, %rbp
	cmpq	$19, %rbp
	jg	.L1764
	testq	%rax, %rax
	jg	.L1762
.L1764:
	leaq	48(%rsp), %rcx
	movq	%r14, %rdx
	movq	%r13, %rsi
	movq	%r8, %rdi
	movq	.LC86(%rip), %xmm0
	movaps	%xmm0, 48(%rsp)
	call	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_
.LEHE65:
	jmp	.L1753
.L1788:
	leaq	.LC36(%rip), %rcx
	movl	$273, %edx
	leaq	.LC37(%rip), %rsi
	leaq	.LC38(%rip), %rdi
	call	__assert_fail@PLT
.L1792:
	call	__stack_chk_fail@PLT
.L1791:
	leaq	.LC83(%rip), %rcx
	movl	$96, %edx
	leaq	.LC84(%rip), %rsi
	leaq	.LC85(%rip), %rdi
	call	__assert_fail@PLT
.L1793:
.LEHB66:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.L1789:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.L1772:
.L1787:
	movq	%rax, %rbp
	jmp	.L1768
.L1790:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE66:
.L1771:
	jmp	.L1787
.L1768:
	movq	16(%rbx), %rdi
	call	free@PLT
	movq	%rbp, %rdi
.LEHB67:
	call	_Unwind_Resume@PLT
.LEHE67:
	.cfi_endproc
.LFE14242:
	.section	.gcc_except_table
.LLSDA14242:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE14242-.LLSDACSB14242
.LLSDACSB14242:
	.uleb128 .LEHB65-.LFB14242
	.uleb128 .LEHE65-.LEHB65
	.uleb128 .L1771-.LFB14242
	.uleb128 0
	.uleb128 .LEHB66-.LFB14242
	.uleb128 .LEHE66-.LEHB66
	.uleb128 .L1772-.LFB14242
	.uleb128 0
	.uleb128 .LEHB67-.LFB14242
	.uleb128 .LEHE67-.LEHB67
	.uleb128 0
	.uleb128 0
.LLSDACSE14242:
	.section	.text._ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_,"axG",@progbits,_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC5ERKSF_,comdat
	.size	_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_, .-_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_
	.weak	_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC1ERKSF_
	.set	_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC1ERKSF_,_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_
	.section	.rodata._ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_.str1.8,"aMS",@progbits,1
	.align 8
.LC87:
	.ascii	"static void Eigen::internal::generic_product_impl<Lhs, Rhs, "
	.ascii	"Eigen::DenseShape, Eigen::DenseShape, 8>::scaleAndAddTo(Dest"
	.ascii	"&, const Lhs&, const Rhs&, const Scalar&) [with Dest = Eigen"
	.ascii	"::Matrix<std::complex<double>, -1, -1>; Lhs = Eigen::Product"
	.ascii	"<Eigen::CwiseUnaryOp<Eigen::intern"
	.string	"al::scalar_conjugate_op<std::complex<double> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >, Eigen::Matrix<std::complex<double>, -1, -1>, 0>; Rhs = Eigen::Matrix<std::complex<double>, -1, -1>; Scalar = std::complex<double>]"
	.align 8
.LC88:
	.ascii	"typename Eigen::ScalarBinaryOpTraits<typename Eigen::interna"
	.ascii	"l::traits<T>::Scalar, typename Eigen::internal::traits<Other"
	.ascii	"Derived>::Scalar>::ReturnType Eigen::MatrixBase<Derived>::do"
	.ascii	"t(const Eigen::MatrixBase<OtherDerived>&) const [with OtherD"
	.ascii	"erived = Eigen::Block<const Eigen::Block<const Eigen::Matrix"
	.ascii	"<std::complex<double>, -1, -1>, -1, 1, true>, -1, 1, true>; "
	.ascii	"Derived = Eigen::CwiseUnaryOp<Eigen::internal::scalar_conjug"
	.ascii	"ate_op<std::complex<double> >, const Eigen::Block<const Eige"
	.ascii	"n::Product<Eigen::CwiseUnaryOp<Eigen::internal::scalar_conju"
	.ascii	"gate_op<std::complex<double> >, const Eigen::Transpose<const"
	.ascii	" Eigen::Matrix<std::complex<double>, -1, -1> > >, Eigen::Mat"
	.ascii	"rix<std::complex<double>, -1, -1>, 0>, 1, -1, false> >; type"
	.ascii	"name Eigen::ScalarBinaryOpTraits<typename Eigen::internal::t"
	.string	"raits<T>::Scalar, typename Eigen::internal::traits<OtherDerived>::Scalar>::ReturnType = std::complex<double>; typename Eigen::internal::traits<OtherDerived>::Scalar = std::complex<double>; typename Eigen::internal::traits<T>::Scalar = std::complex<double>]"
	.align 8
.LC89:
	.ascii	"typename Eigen::internal::traits<T>::Scalar Eigen::DenseBase"
	.ascii	"<Derived>::redux(const Func&) const [with BinaryOp = Eigen::"
	.ascii	"internal::scalar_sum_op<std::complex<double>, std::complex<d"
	.ascii	"ouble> >; Derived = Eigen::CwiseBinaryOp<Eigen::internal::sc"
	.ascii	"alar_conj_product_op<std::complex<double>, std::complex<doub"
	.ascii	"le> >, const Eigen::Transpose<const Eigen::CwiseUnaryOp<Eige"
	.ascii	"n::internal::scalar_conjugate_op<std::complex<double> >, con"
	.ascii	"st Eigen::Block<const Eigen::Product<Eigen::CwiseUnaryOp<Eig"
	.ascii	"en::internal::scalar_conjugate_op<std::complex<double> >, co"
	.ascii	"nst Eigen::Transpose<const Eigen::Matrix<std::complex<double"
	.ascii	">, -1, -1> > >"
	.string	", Eigen::Matrix<std::complex<double>, -1, -1>, 0>, 1, -1, false> > >, const Eigen::Block<const Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true>, -1, 1, true> >; typename Eigen::internal::traits<T>::Scalar = std::complex<double>]"
	.align 8
.LC90:
	.ascii	"typename Eigen::ScalarBinaryOpTraits<typename Eigen::interna"
	.ascii	"l::traits<T>::Scalar, typename Eigen::internal::traits<Other"
	.ascii	"Derived>::Scalar>::ReturnType Eigen::MatrixBase<Derived>::do"
	.ascii	"t(const Eigen::MatrixBase<OtherDerived>&) const [with OtherD"
	.ascii	"erived = Eigen::Block<const Eigen::Matrix<std::complex<doubl"
	.ascii	"e>, -1, -1>, -1, 1, true>; Derived = Eigen::CwiseUnaryOp<Eig"
	.ascii	"en::internal::scalar_conjugate_op<std::complex<double> >, co"
	.ascii	"nst Eigen::Block<const Eigen::Block<const Eigen::Product<Eig"
	.ascii	"en::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std::c"
	.ascii	"omplex<double> >, const Eigen::Transpose<const Eigen::Matrix"
	.ascii	"<std::complex<double>, -1, -1> > >, Eigen::Matrix<std::compl"
	.ascii	"ex<double>, -1, -1>, 0>, 1, -1, false>, 1, -1, true> >; type"
	.ascii	"name Eigen::ScalarBinaryOpTraits<typename Eigen::internal::t"
	.string	"raits<T>::Scalar, typename Eigen::internal::traits<OtherDerived>::Scalar>::ReturnType = std::complex<double>; typename Eigen::internal::traits<OtherDerived>::Scalar = std::complex<double>; typename Eigen::internal::traits<T>::Scalar = std::complex<double>]"
	.align 8
.LC91:
	.ascii	"typename Eigen::internal::traits<T>::Scalar Eigen::DenseBase"
	.ascii	"<Derived>::redux(const Func&) const [with BinaryOp = Eigen::"
	.ascii	"internal::scalar_sum_op<std::complex<double>, std::complex<d"
	.ascii	"ouble> >; Derived = Eigen::CwiseBinaryOp<Eigen::internal::sc"
	.ascii	"alar_conj_product_op<std::complex<double>, std::complex<doub"
	.ascii	"le> >, const Eigen::Transpose<const Eigen::CwiseUnaryOp<Eige"
	.ascii	"n::internal::scalar_conjugate_op<std::complex<double> >, con"
	.ascii	"st Eigen::Block<const Eigen::Block<const Eigen::Product<Eige"
	.ascii	"n::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std::co"
	.ascii	"mplex<double> >, const Eigen::Transpose<const Eigen::Matrix<"
	.ascii	"std::complex<d"
	.string	"ouble>, -1, -1> > >, Eigen::Matrix<std::complex<double>, -1, -1>, 0>, 1, -1, false>, 1, -1, true> > >, const Eigen::Block<const Eigen::Matrix<std::complex<double>, -1, -1>, -1, 1, true> >; typename Eigen::internal::traits<T>::Scalar = std::complex<double>]"
	.align 8
.LC92:
	.string	"void Eigen::PlainObjectBase<Derived>::resize(Eigen::Index, Eigen::Index) [with Derived = Eigen::Matrix<std::complex<double>, 1, -1>; Eigen::Index = long int]"
	.section	.text._ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_,"axG",@progbits,_ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_,comdat
	.p2align 4
	.weak	_ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_
	.type	_ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_, @function
_ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_:
.LFB11847:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA11847
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$856, %rsp
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	%rdi, -864(%rbp)
	movq	(%rsi), %r8
	movq	16(%r8), %r9
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	movq	8(%rdi), %rax
	cmpq	%r9, %rax
	jne	.L1795
	movq	%rsi, %rbx
	movq	16(%rdi), %r14
	movq	16(%rdx), %rsi
	movq	%rdx, %r12
	cmpq	%rsi, %r14
	jne	.L1795
	movq	%rcx, %r13
	movq	16(%rbx), %rcx
	movq	16(%rcx), %r15
	testq	%r15, %r15
	sete	%dl
	testq	%rax, %rax
	sete	%dil
	orb	%dil, %dl
	jne	.L1794
	testq	%r14, %r14
	jne	.L1893
.L1794:
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L1891
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
	.p2align 4,,10
	.p2align 3
.L1893:
	.cfi_restore_state
	cmpq	$1, %r14
	je	.L1894
	cmpq	$1, %rax
	je	.L1895
	movq	$0, -496(%rbp)
	pxor	%xmm0, %xmm0
	movups	%xmm0, -488(%rbp)
	movq	16(%r8), %rsi
	movq	16(%rcx), %rcx
	testq	%rsi, %rsi
	je	.L1833
	testq	%rcx, %rcx
	je	.L1833
	movabsq	$9223372036854775807, %rax
	cqto
	idivq	%rcx
	cmpq	%rax, %rsi
	jg	.L1896
.L1833:
	leaq	-496(%rbp), %r14
	movq	%rcx, %rdx
	movq	%r14, %rdi
.LEHB68:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
	movq	(%rbx), %rax
	movq	16(%rbx), %r15
	movq	16(%rax), %rsi
	movq	-488(%rbp), %rax
	movq	16(%r15), %rcx
	cmpq	%rax, %rsi
	jne	.L1834
	movq	-480(%rbp), %rdx
	cmpq	%rdx, %rcx
	je	.L1835
.L1834:
	movq	%rcx, %rdx
	movq	%r14, %rdi
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
	movq	16(%rbx), %r15
	movq	-488(%rbp), %rax
	movq	-480(%rbp), %rdx
.L1835:
	movq	8(%r15), %rsi
	leaq	(%rsi,%rax), %rcx
	addq	%rdx, %rcx
	cmpq	$19, %rcx
	jg	.L1836
	testq	%rsi, %rsi
	jg	.L1897
.L1836:
	movq	%rdx, %rcx
	orq	%rax, %rcx
	js	.L1898
	imulq	%rax, %rdx
	testq	%rdx, %rdx
	je	.L1842
	salq	$4, %rdx
	movq	-496(%rbp), %rdi
	je	.L1842
	xorl	%esi, %esi
	call	memset@PLT
.L1842:
	leaq	-384(%rbp), %rcx
	movq	%r15, %rdx
	movq	%rbx, %rsi
	movq	%r14, %rdi
	movq	.LC86(%rip), %xmm0
	movaps	%xmm0, -384(%rbp)
	call	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_
.LEHE68:
	leaq	-256(%rbp), %rax
	movq	%rax, -872(%rbp)
.L1838:
	movsd	8(%r13), %xmm1
	movsd	0(%r13), %xmm5
	pxor	%xmm0, %xmm0
	movapd	%xmm1, %xmm2
	movapd	%xmm5, %xmm4
	mulsd	%xmm0, %xmm2
	mulsd	%xmm5, %xmm0
	subsd	%xmm2, %xmm4
	addsd	%xmm1, %xmm0
	ucomisd	%xmm4, %xmm0
	jp	.L1899
.L1843:
	pxor	%xmm1, %xmm1
	movapd	%xmm0, %xmm3
	movapd	%xmm4, %xmm2
	mulsd	%xmm1, %xmm3
	mulsd	%xmm4, %xmm1
	subsd	%xmm3, %xmm2
	addsd	%xmm0, %xmm1
	ucomisd	%xmm2, %xmm1
	jp	.L1900
.L1844:
	movq	-864(%rbp), %r15
	pxor	%xmm0, %xmm0
	movq	-480(%rbp), %rax
	leaq	-232(%rbp), %rdx
	leaq	-240(%rbp), %rsi
	leaq	-224(%rbp), %rdi
	movl	$1, %ecx
	movsd	%xmm2, -888(%rbp)
	movdqu	8(%r15), %xmm7
	movq	%rax, -224(%rbp)
	movsd	%xmm1, -880(%rbp)
	movaps	%xmm0, -256(%rbp)
	movaps	%xmm7, -864(%rbp)
	movaps	%xmm7, -240(%rbp)
	call	_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_
	movq	-224(%rbp), %rax
	movq	-240(%rbp), %rdx
	subq	$8, %rsp
	movsd	-888(%rbp), %xmm2
	movsd	-880(%rbp), %xmm1
	imulq	%rax, %rdx
	movq	16(%r12), %rsi
	imulq	-232(%rbp), %rax
	movapd	%xmm2, %xmm0
	movq	%rdx, -216(%rbp)
	movq	%rax, -208(%rbp)
	movq	(%rbx), %rax
	movq	16(%rax), %rdi
	movsd	%xmm2, -360(%rbp)
	movsd	%xmm1, -352(%rbp)
	pushq	-872(%rbp)
	pushq	8(%r15)
	pushq	$1
	pushq	(%r15)
	pushq	8(%r12)
	movq	(%r12), %r9
	movq	-488(%rbp), %r8
	movq	-496(%rbp), %rcx
	movq	-480(%rbp), %rdx
.LEHB69:
	.cfi_escape 0x2e,0x30
	call	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi0ELb0ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0
.LEHE69:
	movq	-256(%rbp), %rdi
	addq	$48, %rsp
	call	free@PLT
	movq	-248(%rbp), %rdi
	call	free@PLT
	movq	-496(%rbp), %rdi
	call	free@PLT
	jmp	.L1794
	.p2align 4,,10
	.p2align 3
.L1897:
	movq	(%rbx), %rax
	movq	%r15, -240(%rbp)
	movq	%rax, -256(%rbp)
	cmpq	8(%rax), %rsi
	jne	.L1901
	leaq	-256(%rbp), %rax
	leaq	-608(%rbp), %rdx
	movq	%r14, %rdi
	movq	%rax, %rsi
	movq	%rax, -872(%rbp)
.LEHB70:
	.cfi_escape 0x2e,0
	call	_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_
.LEHE70:
	jmp	.L1838
	.p2align 4,,10
	.p2align 3
.L1895:
	movq	-864(%rbp), %rax
	movq	(%rax), %rax
	movq	%rax, -880(%rbp)
	testq	%rax, %rax
	je	.L1813
	testq	%r14, %r14
	jns	.L1813
	leaq	.LC80(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
	.p2align 4,,10
	.p2align 3
.L1813:
	movq	(%rbx), %rax
	pxor	%xmm0, %xmm0
	movq	%rcx, -816(%rbp)
	movq	%r15, -784(%rbp)
	movq	%rax, -832(%rbp)
	movups	%xmm0, -808(%rbp)
	cmpq	$1, %rsi
	je	.L1902
	leaq	-256(%rbp), %rax
	leaq	-832(%rbp), %rsi
	movq	$0, -384(%rbp)
	movq	%rax, %rdi
	movq	%rax, -872(%rbp)
	movq	$0, -376(%rbp)
.LEHB71:
	call	_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_
.LEHE71:
	pxor	%xmm0, %xmm0
	movups	%xmm0, -216(%rbp)
	testq	%r15, %r15
	js	.L1903
	movq	%r15, %rdi
.LEHB72:
	call	_ZN5Eigen8internal28conditional_aligned_new_autoISt7complexIdELb1EEEPT_m
.LEHE72:
	movq	-248(%rbp), %rdx
	movq	%rax, -384(%rbp)
	movq	%rax, %rbx
	movq	%r15, -376(%rbp)
	movq	-208(%rbp), %rax
	movq	-216(%rbp), %r8
	movq	-256(%rbp), %rcx
	cmpq	$1, %rdx
	jne	.L1904
	addq	%r8, %rax
	salq	$4, %rax
	addq	%rax, %rcx
	testq	%r15, %r15
	cmovg	%r15, %rdx
	xorl	%eax, %eax
	salq	$4, %rdx
	.p2align 4,,10
	.p2align 3
.L1828:
	movupd	(%rcx,%rax), %xmm6
	movups	%xmm6, (%rbx,%rax)
	addq	$16, %rax
	cmpq	%rdx, %rax
	jne	.L1828
.L1827:
	movq	-240(%rbp), %rdi
	call	free@PLT
	movq	-880(%rbp), %rax
	pxor	%xmm0, %xmm0
	movsd	8(%r13), %xmm1
	movaps	%xmm0, -224(%rbp)
	movq	-872(%rbp), %rdx
	movsd	0(%r13), %xmm0
	movq	%r12, %rdi
	movq	%rax, -256(%rbp)
	movq	-864(%rbp), %rax
	leaq	-384(%rbp), %rsi
	movq	%r14, -240(%rbp)
	movq	%rax, -232(%rbp)
	movq	$1, -208(%rbp)
.LEHB73:
	call	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS5_IS7_Li1ELin1ELi1ELi1ELin1EEEEENS4_INS_5BlockIS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNSN_6ScalarE.isra.0
.LEHE73:
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L1891
	leaq	-40(%rbp), %rsp
	movq	%rbx, %rdi
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	jmp	free@PLT
	.p2align 4,,10
	.p2align 3
.L1902:
	.cfi_restore_state
	movq	-808(%rbp), %rsi
	movq	-800(%rbp), %rdx
	leaq	-768(%rbp), %rdi
	movq	%rcx, -864(%rbp)
	movq	%rax, -704(%rbp)
	movq	%rsi, -680(%rbp)
	movq	%rdx, -672(%rbp)
	movq	%rsi, -576(%rbp)
	movq	%r12, %rsi
	movq	%rdx, -568(%rbp)
	xorl	%edx, %edx
	movq	%r15, -656(%rbp)
	movq	%rax, -600(%rbp)
	movq	%r15, -552(%rbp)
	call	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS5_l
	cmpq	-760(%rbp), %r15
	movq	-864(%rbp), %rcx
	jne	.L1905
	movq	-552(%rbp), %rax
	movdqa	-768(%rbp), %xmm5
	movq	%rcx, -224(%rbp)
	movq	-600(%rbp), %rdi
	movq	-576(%rbp), %rsi
	movq	$0, -176(%rbp)
	movq	-568(%rbp), %rdx
	movq	%rax, -440(%rbp)
	movq	%rax, -192(%rbp)
	movdqa	-752(%rbp), %xmm6
	movq	-720(%rbp), %rax
	movaps	%xmm5, -144(%rbp)
	movdqa	-736(%rbp), %xmm5
	movq	%rdi, -488(%rbp)
	movq	%rsi, -464(%rbp)
	movq	%rdx, -456(%rbp)
	movq	%rdi, -240(%rbp)
	movq	%rsi, -216(%rbp)
	movq	%rdx, -208(%rbp)
	movq	%r15, -160(%rbp)
	movq	%rax, -96(%rbp)
	movaps	%xmm6, -128(%rbp)
	movaps	%xmm5, -112(%rbp)
	testq	%r15, %r15
	jle	.L1906
	leaq	-240(%rbp), %rsi
	leaq	-368(%rbp), %rdi
.LEHB74:
	call	_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_
.LEHE74:
	movq	-120(%rbp), %rax
	pxor	%xmm0, %xmm0
	xorl	%edx, %edx
	movups	%xmm0, -328(%rbp)
	movq	-144(%rbp), %rbx
	leaq	-384(%rbp), %r14
	xorl	%esi, %esi
	movups	%xmm0, -296(%rbp)
	movq	8(%rax), %rax
	movq	%r14, %rdi
	movq	%rbx, -280(%rbp)
	movq	%rax, -264(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNSB_IKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEELi1ELin1ELb1EEEEEEEKNSB_ISF_Lin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll
	movsd	%xmm0, -848(%rbp)
	movsd	%xmm1, -840(%rbp)
	movapd	-848(%rbp), %xmm0
	movapd	%xmm0, %xmm2
	shufpd	$1, %xmm0, %xmm2
	cmpq	$1, %r15
	je	.L1907
	movq	-360(%rbp), %rcx
	movq	-368(%rbp), %rdx
	salq	$4, %r15
	movl	$16, %eax
	cmpq	$1, %rcx
	jne	.L1908
	.p2align 4,,10
	.p2align 3
.L1821:
	movupd	(%rdx,%rax), %xmm3
	movupd	(%rbx,%rax), %xmm1
	movapd	%xmm3, %xmm0
	unpcklpd	%xmm1, %xmm1
	shufpd	$1, %xmm3, %xmm0
	mulpd	%xmm0, %xmm1
	movupd	(%rbx,%rax), %xmm0
	addq	$16, %rax
	unpckhpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	movapd	%xmm0, %xmm3
	addpd	%xmm1, %xmm3
	subpd	%xmm0, %xmm1
	movsd	%xmm3, %xmm1
	addpd	%xmm1, %xmm2
	cmpq	%r15, %rax
	jne	.L1821
	movq	%xmm2, %r14
	unpckhpd	%xmm2, %xmm2
	movq	%xmm2, %rbx
.L1817:
	movq	-352(%rbp), %rdi
	call	free@PLT
	movupd	0(%r13), %xmm1
	movq	%r14, %xmm2
	movq	%rbx, %xmm0
	unpcklpd	%xmm2, %xmm2
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm1, %xmm2
	mulpd	%xmm1, %xmm0
	shufpd	$1, %xmm2, %xmm2
	movapd	%xmm2, %xmm3
	addpd	%xmm0, %xmm3
	subpd	%xmm2, %xmm0
	movapd	%xmm3, %xmm2
	unpckhpd	%xmm3, %xmm3
	ucomisd	%xmm0, %xmm3
	movsd	%xmm0, %xmm2
	jp	.L1909
.L1822:
	movq	-880(%rbp), %rax
	movupd	(%rax), %xmm6
	movapd	%xmm6, %xmm0
	movaps	%xmm6, -864(%rbp)
	addpd	%xmm2, %xmm0
	movups	%xmm0, (%rax)
	jmp	.L1794
	.p2align 4,,10
	.p2align 3
.L1894:
	movq	-864(%rbp), %rsi
	movq	(%rsi), %rsi
	movq	%rsi, -888(%rbp)
	testq	%rsi, %rsi
	je	.L1799
	testq	%rax, %rax
	js	.L1910
.L1799:
	xorl	%edx, %edx
	leaq	-832(%rbp), %rdi
	movq	%r12, %rsi
	movq	%r9, -880(%rbp)
	movq	%rcx, -872(%rbp)
	movq	%r8, -864(%rbp)
	call	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS5_l
	movq	-880(%rbp), %r9
	movq	-864(%rbp), %r8
	movq	-872(%rbp), %rcx
	cmpq	$1, %r9
	je	.L1911
	movq	$0, -256(%rbp)
	pxor	%xmm0, %xmm0
	movups	%xmm0, -248(%rbp)
	movq	16(%r8), %rsi
	movq	16(%rcx), %rdx
	movq	%rsi, %rax
	orq	%rdx, %rax
	leaq	-256(%rbp), %rax
	movq	%rax, -872(%rbp)
	je	.L1808
	movq	%rax, %rdi
.LEHB75:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
	movq	16(%rbx), %rcx
.L1808:
	movq	-872(%rbp), %rdi
	movq	%rcx, %rdx
	movq	%rbx, %rsi
	call	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE6evalToIS9_EEvRT_RKSD_RSA_
.LEHE75:
	movsd	8(%r13), %xmm1
	movsd	0(%r13), %xmm5
	pxor	%xmm0, %xmm0
	movq	-832(%rbp), %rbx
	movapd	%xmm1, %xmm2
	movapd	%xmm5, %xmm4
	mulsd	%xmm0, %xmm2
	mulsd	%xmm5, %xmm0
	subsd	%xmm2, %xmm4
	addsd	%xmm1, %xmm0
	ucomisd	%xmm4, %xmm0
	jp	.L1912
.L1809:
	pxor	%xmm1, %xmm1
	movapd	%xmm0, %xmm3
	movapd	%xmm4, %xmm2
	mulsd	%xmm1, %xmm3
	mulsd	%xmm4, %xmm1
	subsd	%xmm3, %xmm2
	addsd	%xmm0, %xmm1
	ucomisd	%xmm2, %xmm1
	jp	.L1913
.L1810:
	movq	-248(%rbp), %rdi
	movq	-256(%rbp), %rax
	leaq	-496(%rbp), %rcx
	leaq	-384(%rbp), %rdx
	movq	-888(%rbp), %r8
	movq	-240(%rbp), %rsi
	movapd	%xmm2, %xmm0
	movq	%rbx, -496(%rbp)
	movq	%rdi, -376(%rbp)
	movq	$1, -488(%rbp)
	movq	%rax, -384(%rbp)
	call	_ZN5Eigen8internal29general_matrix_vector_productIlSt7complexIdENS0_22const_blas_data_mapperIS3_lLi0EEELi0ELb0ES3_NS4_IS3_lLi1EEELb0ELi0EE3runEllRKS5_RKS6_PS3_lS3_.isra.0
	movq	-256(%rbp), %rdi
	call	free@PLT
	jmp	.L1794
.L1904:
	imulq	%rdx, %rax
	movq	%rdx, %rdi
	movl	$1, %edx
	movq	%rbx, %rsi
	salq	$4, %rdi
	addq	%r8, %rax
	salq	$4, %rax
	addq	%rcx, %rax
	testq	%r15, %r15
	cmovg	%r15, %rdx
	salq	$4, %rdx
	addq	%rbx, %rdx
	.p2align 4,,10
	.p2align 3
.L1826:
	movupd	(%rax), %xmm6
	addq	$16, %rsi
	addq	%rdi, %rax
	movups	%xmm6, -16(%rsi)
	cmpq	%rdx, %rsi
	jne	.L1826
	jmp	.L1827
.L1911:
	movq	(%rbx), %rax
	movq	-808(%rbp), %r8
	movq	-832(%rbp), %r12
	movq	-824(%rbp), %rbx
	movq	%rax, -768(%rbp)
	movq	%rax, -696(%rbp)
	movq	8(%r8), %rax
	testq	%r12, %r12
	je	.L1801
	testq	%rbx, %rbx
	js	.L1914
.L1801:
	movdqa	-832(%rbp), %xmm5
	movdqa	-816(%rbp), %xmm6
	movq	-784(%rbp), %rdx
	movups	%xmm5, -472(%rbp)
	movdqa	-800(%rbp), %xmm5
	movq	%rdx, -424(%rbp)
	movups	%xmm6, -456(%rbp)
	movups	%xmm5, -440(%rbp)
	cmpq	%rbx, %r15
	jne	.L1915
	movq	%r12, -496(%rbp)
	movdqa	-480(%rbp), %xmm5
	pxor	%xmm0, %xmm0
	movq	%r15, -488(%rbp)
	movdqa	-496(%rbp), %xmm6
	movq	-696(%rbp), %rdx
	movdqa	-448(%rbp), %xmm7
	movq	$0, -416(%rbp)
	movaps	%xmm6, -176(%rbp)
	movdqa	-464(%rbp), %xmm6
	movaps	%xmm5, -160(%rbp)
	movdqa	-432(%rbp), %xmm5
	movaps	%xmm6, -144(%rbp)
	movdqa	-416(%rbp), %xmm6
	movq	%rdx, -600(%rbp)
	movq	%rdx, -240(%rbp)
	movq	%rcx, -224(%rbp)
	movq	%r15, -192(%rbp)
	movq	%rax, -400(%rbp)
	movq	%rax, -80(%rbp)
	movups	%xmm0, -216(%rbp)
	movaps	%xmm7, -128(%rbp)
	movaps	%xmm5, -112(%rbp)
	movaps	%xmm6, -96(%rbp)
	testq	%r15, %r15
	jle	.L1916
	leaq	-240(%rbp), %rsi
	leaq	-368(%rbp), %rdi
	movq	%r9, -880(%rbp)
	movq	%r8, -864(%rbp)
	leaq	-384(%rbp), %r14
.LEHB76:
	call	_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_
	movq	-864(%rbp), %r8
	xorl	%edx, %edx
	xorl	%esi, %esi
	pxor	%xmm0, %xmm0
	movq	%r14, %rdi
	movq	%r12, -304(%rbp)
	movups	%xmm0, -328(%rbp)
	movq	8(%r8), %rax
	movq	%rax, -288(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll
	cmpq	$1, %rbx
	movq	-880(%rbp), %r9
	movsd	%xmm0, -864(%rbp)
	movsd	%xmm1, -872(%rbp)
	je	.L1804
	.p2align 4,,10
	.p2align 3
.L1805:
	movq	%r9, %rsi
	xorl	%edx, %edx
	movq	%r14, %rdi
	movq	%r9, -880(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll
	movq	-880(%rbp), %r9
	addsd	-864(%rbp), %xmm0
	addsd	-872(%rbp), %xmm1
	addq	$1, %r9
	movsd	%xmm0, -864(%rbp)
	movsd	%xmm1, -872(%rbp)
	cmpq	%r9, %rbx
	jne	.L1805
.L1804:
	movq	-352(%rbp), %rdi
	call	free@PLT
	movsd	-864(%rbp), %xmm0
	movupd	0(%r13), %xmm1
	movsd	-872(%rbp), %xmm2
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm1, %xmm0
	unpcklpd	%xmm2, %xmm2
	mulpd	%xmm1, %xmm2
	movapd	%xmm0, %xmm3
	shufpd	$1, %xmm2, %xmm2
	addpd	%xmm2, %xmm3
	subpd	%xmm2, %xmm0
	movapd	%xmm3, %xmm2
	unpckhpd	%xmm3, %xmm3
	ucomisd	%xmm0, %xmm3
	movsd	%xmm0, %xmm2
	jp	.L1917
.L1806:
	movq	-888(%rbp), %rax
	movupd	(%rax), %xmm7
	movapd	%xmm7, %xmm0
	movaps	%xmm7, -864(%rbp)
	addpd	%xmm2, %xmm0
	movups	%xmm0, (%rax)
	jmp	.L1794
.L1795:
	leaq	.LC87(%rip), %rcx
	movl	$470, %edx
	leaq	.LC73(%rip), %rsi
	leaq	.LC74(%rip), %rdi
	call	__assert_fail@PLT
.L1908:
	salq	$4, %rcx
	leaq	16(%rbx), %rax
	addq	%r15, %rbx
	addq	%rcx, %rdx
.L1819:
	movupd	(%rdx), %xmm3
	movupd	(%rax), %xmm1
	addq	$16, %rax
	addq	%rcx, %rdx
	movapd	%xmm3, %xmm0
	unpcklpd	%xmm1, %xmm1
	shufpd	$1, %xmm3, %xmm0
	mulpd	%xmm0, %xmm1
	movupd	-16(%rax), %xmm0
	unpckhpd	%xmm0, %xmm0
	mulpd	%xmm3, %xmm0
	movapd	%xmm0, %xmm3
	addpd	%xmm1, %xmm3
	subpd	%xmm0, %xmm1
	movsd	%xmm3, %xmm1
	addpd	%xmm1, %xmm2
	cmpq	%rax, %rbx
	jne	.L1819
	movq	%xmm2, %r14
	unpckhpd	%xmm2, %xmm2
	movq	%xmm2, %rbx
	jmp	.L1817
.L1907:
	movq	%xmm0, %rbx
	unpckhpd	%xmm0, %xmm0
	movq	%xmm0, %r14
	jmp	.L1817
.L1891:
	call	__stack_chk_fail@PLT
.L1900:
	movsd	.LC24(%rip), %xmm2
	movapd	%xmm0, %xmm1
	pxor	%xmm3, %xmm3
	movapd	%xmm4, %xmm0
	call	__muldc3@PLT
	movapd	%xmm0, %xmm2
	jmp	.L1844
.L1899:
	movsd	.LC24(%rip), %xmm2
	movapd	%xmm5, %xmm0
	pxor	%xmm3, %xmm3
	call	__muldc3@PLT
	movapd	%xmm0, %xmm4
	movapd	%xmm1, %xmm0
	jmp	.L1843
.L1898:
	call	_ZN5Eigen9DenseBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8ConstantEllRKS3_.part.0
.L1916:
	leaq	.LC89(%rip), %rcx
	movl	$411, %edx
	leaq	.LC67(%rip), %rsi
	leaq	.LC68(%rip), %rdi
	call	__assert_fail@PLT
.L1915:
	leaq	.LC88(%rip), %rcx
	movl	$82, %edx
	leaq	.LC77(%rip), %rsi
	leaq	.LC78(%rip), %rdi
	call	__assert_fail@PLT
.L1914:
	leaq	.LC75(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L1910:
	leaq	.LC52(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L1909:
	movapd	%xmm1, %xmm0
	movq	%rbx, %xmm2
	unpckhpd	%xmm1, %xmm1
	movq	%r14, %xmm3
	call	__muldc3@PLT
	movapd	%xmm0, %xmm2
	unpcklpd	%xmm1, %xmm2
	jmp	.L1822
.L1917:
	movsd	-864(%rbp), %xmm2
	movsd	-872(%rbp), %xmm3
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm1
	call	__muldc3@PLT
	movapd	%xmm0, %xmm2
	unpcklpd	%xmm1, %xmm2
	jmp	.L1806
.L1906:
	leaq	.LC91(%rip), %rcx
	movl	$411, %edx
	leaq	.LC67(%rip), %rsi
	leaq	.LC68(%rip), %rdi
	call	__assert_fail@PLT
.L1905:
	leaq	.LC90(%rip), %rcx
	movl	$82, %edx
	leaq	.LC77(%rip), %rsi
	leaq	.LC78(%rip), %rdi
	call	__assert_fail@PLT
.L1901:
	leaq	.LC83(%rip), %rcx
	movl	$96, %edx
	leaq	.LC84(%rip), %rsi
	leaq	.LC85(%rip), %rdi
	call	__assert_fail@PLT
.L1903:
	leaq	.LC92(%rip), %rcx
	movl	$273, %edx
	leaq	.LC37(%rip), %rsi
	leaq	.LC38(%rip), %rdi
	call	__assert_fail@PLT
.L1913:
	movsd	.LC24(%rip), %xmm2
	movapd	%xmm0, %xmm1
	pxor	%xmm3, %xmm3
	movapd	%xmm4, %xmm0
	call	__muldc3@PLT
	movapd	%xmm0, %xmm2
	jmp	.L1810
.L1912:
	movsd	.LC24(%rip), %xmm2
	movapd	%xmm5, %xmm0
	pxor	%xmm3, %xmm3
	call	__muldc3@PLT
	movapd	%xmm0, %xmm4
	movapd	%xmm1, %xmm0
	jmp	.L1809
.L1852:
	movq	%rax, %rbx
	jmp	.L1848
.L1855:
	movq	%rax, %rbx
	jmp	.L1830
.L1848:
	movq	-256(%rbp), %rdi
	call	free@PLT
	movq	-248(%rbp), %rdi
	call	free@PLT
.L1892:
	movq	-496(%rbp), %rdi
	call	free@PLT
	movq	%rbx, %rdi
	call	_Unwind_Resume@PLT
.L1830:
	movq	-240(%rbp), %rdi
	call	free@PLT
	movq	%rbx, %rdi
	call	_Unwind_Resume@PLT
.LEHE76:
.L1856:
	movq	%rax, %rbx
	jmp	.L1846
.L1896:
.LEHB77:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE77:
.L1846:
	jmp	.L1892
.L1854:
	movq	%rax, %r14
	jmp	.L1832
.L1853:
	movq	%rax, %rbx
	jmp	.L1811
.L1832:
	movq	%rbx, %rdi
	call	free@PLT
	movq	%r14, %rdi
.LEHB78:
	call	_Unwind_Resume@PLT
.L1811:
	movq	-256(%rbp), %rdi
	call	free@PLT
	movq	%rbx, %rdi
	call	_Unwind_Resume@PLT
.LEHE78:
	.cfi_endproc
.LFE11847:
	.section	.gcc_except_table
.LLSDA11847:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE11847-.LLSDACSB11847
.LLSDACSB11847:
	.uleb128 .LEHB68-.LFB11847
	.uleb128 .LEHE68-.LEHB68
	.uleb128 .L1856-.LFB11847
	.uleb128 0
	.uleb128 .LEHB69-.LFB11847
	.uleb128 .LEHE69-.LEHB69
	.uleb128 .L1852-.LFB11847
	.uleb128 0
	.uleb128 .LEHB70-.LFB11847
	.uleb128 .LEHE70-.LEHB70
	.uleb128 .L1856-.LFB11847
	.uleb128 0
	.uleb128 .LEHB71-.LFB11847
	.uleb128 .LEHE71-.LEHB71
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB72-.LFB11847
	.uleb128 .LEHE72-.LEHB72
	.uleb128 .L1855-.LFB11847
	.uleb128 0
	.uleb128 .LEHB73-.LFB11847
	.uleb128 .LEHE73-.LEHB73
	.uleb128 .L1854-.LFB11847
	.uleb128 0
	.uleb128 .LEHB74-.LFB11847
	.uleb128 .LEHE74-.LEHB74
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB75-.LFB11847
	.uleb128 .LEHE75-.LEHB75
	.uleb128 .L1853-.LFB11847
	.uleb128 0
	.uleb128 .LEHB76-.LFB11847
	.uleb128 .LEHE76-.LEHB76
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB77-.LFB11847
	.uleb128 .LEHE77-.LEHB77
	.uleb128 .L1856-.LFB11847
	.uleb128 0
	.uleb128 .LEHB78-.LFB11847
	.uleb128 .LEHE78-.LEHB78
	.uleb128 0
	.uleb128 0
.LLSDACSE11847:
	.section	.text._ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_,"axG",@progbits,_ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_,comdat
	.size	_ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_, .-_ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_
	.section	.rodata.str1.8
	.align 8
.LC93:
	.ascii	"Eigen::DenseCoeffsBase<Derived, 0>::CoeffReturnType Eigen::D"
	.ascii	"enseCoeffsBase<Derived, 0>::operator()(Eigen::Index, Eigen::"
	.ascii	"Index) const [with Derived = Eigen::Product<Eigen::Product<E"
	.ascii	"igen::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std:"
	.ascii	":complex<d"
	.string	"ouble> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >, Eigen::Matrix<std::complex<double>, -1, -1>, 0>, Eigen::Matrix<std::complex<double>, -1, -1>, 0>; CoeffReturnType = const std::complex<double>; Eigen::Index = long int]"
	.align 8
.LC94:
	.ascii	"Eigen::Product<Lhs, Rhs, Option>::Product(const Lhs&, const "
	.ascii	"Rhs&) [with _Lhs = Eigen::Product<Eigen::CwiseUnaryOp<Eigen:"
	.ascii	":internal::scalar_conjugate_op<std::complex<double> >, const"
	.ascii	" Eigen::Transpose<const Eigen::Matrix<std::complex<double>, "
	.ascii	"-1, -1> > >, Eigen::Matrix<std::complex<double>, -1, -1>, 0>"
	.ascii	"; _Rhs = Eigen::Matrix<std::complex<double>, -1, -1>; int Op"
	.ascii	"tion = 1; Lhs = Eigen::Product<Ei"
	.string	"gen::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std::complex<double> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >, Eigen::Matrix<std::complex<double>, -1, -1>, 0>; Rhs = Eigen::Matrix<std::complex<double>, -1, -1>]"
	.align 8
.LC95:
	.ascii	"void Eigen::internal::resize_if_allowed(DstXprType&, const S"
	.ascii	"rcXprType&, const assign_op<T1, T2>&) [with DstXprType = Eig"
	.ascii	"en::Matrix<std::complex<double>, -1, -1>; SrcXprType = Eigen"
	.ascii	"::Product<Eigen::Product<Eigen::CwiseUnaryOp<scalar_conjugat"
	.ascii	"e_"
	.string	"op<std::complex<double> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >, Eigen::Matrix<std::complex<double>, -1, -1>, 0>, Eigen::Matrix<std::complex<double>, -1, -1>, 1>; T1 = std::complex<double>; T2 = std::complex<double>]"
	.align 8
.LC96:
	.string	"/usr/local/include/Eigen/src/Core/AssignEvaluator.h"
	.align 8
.LC97:
	.string	"dst.rows() == dstRows && dst.cols() == dstCols"
	.section	.text.unlikely
	.align 2
.LCOLDB98:
	.text
.LHOTB98:
	.align 2
	.p2align 4
	.type	_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll, @function
_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll:
.LFB10576:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA10576
	pushq	%r15
	.cfi_def_cfa_offset 16
	.cfi_offset 15, -16
	pushq	%r14
	.cfi_def_cfa_offset 24
	.cfi_offset 14, -24
	pushq	%r13
	.cfi_def_cfa_offset 32
	.cfi_offset 13, -32
	pushq	%r12
	.cfi_def_cfa_offset 40
	.cfi_offset 12, -40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	.cfi_offset 6, -48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	.cfi_offset 3, -56
	subq	$264, %rsp
	.cfi_def_cfa_offset 320
	movq	%fs:40, %rax
	movq	%rax, 248(%rsp)
	xorl	%eax, %eax
	testq	%rsi, %rsi
	js	.L1919
	movq	(%rdi), %r14
	movq	%rdx, %rbx
	movq	16(%r14), %r12
	testq	%rdx, %rdx
	js	.L1919
	movq	%rsi, %rbp
	cmpq	%r12, %rsi
	jge	.L1919
	movq	24(%rdi), %r13
	movq	%r12, %xmm1
	movq	16(%r13), %r15
	movq	%r15, %xmm3
	punpcklqdq	%xmm3, %xmm1
	cmpq	%r15, %rdx
	jge	.L1919
	movabsq	$9223372036854775807, %rax
	pxor	%xmm0, %xmm0
	movq	$0, 128(%rsp)
	cqto
	movups	%xmm0, 152(%rsp)
	idivq	%r15
	movq	$-1, 136(%rsp)
	movq	$0, 144(%rsp)
	cmpq	%rax, %r12
	jg	.L1979
	movabsq	$1152921504606846975, %rdx
	movq	%r12, %rax
	imulq	%r15, %rax
	cmpq	%rdx, %rax
	jg	.L1980
	salq	$4, %rax
	movq	%rdi, 8(%rsp)
	movl	$1, %esi
	movq	%rax, %rdi
	movaps	%xmm1, 16(%rsp)
	call	calloc@PLT
	movq	8(%rsp), %r8
	movdqa	16(%rsp), %xmm1
	testb	$15, %al
	je	.L1923
	leaq	.LC25(%rip), %rcx
	movl	$185, %edx
	leaq	.LC26(%rip), %rsi
	leaq	.LC27(%rip), %rdi
	call	__assert_fail@PLT
.L1923:
	testq	%rax, %rax
	je	.L1981
	movups	%xmm1, 152(%rsp)
	movq	8(%r13), %rdx
	movq	%rax, 144(%rsp)
	movq	%rax, 128(%rsp)
	leaq	(%r12,%rdx), %rax
	addq	%r15, %rax
	movq	%r12, 136(%rsp)
	cmpq	$19, %rax
	jg	.L1927
	testq	%rdx, %rdx
	jg	.L1985
.L1927:
	movq	.LC86(%rip), %xmm0
	leaq	176(%rsp), %rcx
	movq	%r13, %rdx
	movq	%r8, %rsi
	leaq	144(%rsp), %rdi
	movaps	%xmm0, 176(%rsp)
.LEHB79:
	call	_ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_
.LEHE79:
.L1952:
	imulq	136(%rsp), %rbx
	movq	144(%rsp), %rdi
	leaq	(%rbx,%rbp), %rax
	salq	$4, %rax
	addq	128(%rsp), %rax
	movsd	(%rax), %xmm0
	movsd	8(%rax), %xmm1
	movsd	%xmm0, 16(%rsp)
	movsd	%xmm1, 8(%rsp)
	call	free@PLT
	movsd	8(%rsp), %xmm1
	movsd	16(%rsp), %xmm0
	movq	248(%rsp), %rax
	subq	%fs:40, %rax
	jne	.L1986
	addq	$264, %rsp
	.cfi_remember_state
	.cfi_def_cfa_offset 56
	popq	%rbx
	.cfi_def_cfa_offset 48
	popq	%rbp
	.cfi_def_cfa_offset 40
	popq	%r12
	.cfi_def_cfa_offset 32
	popq	%r13
	.cfi_def_cfa_offset 24
	popq	%r14
	.cfi_def_cfa_offset 16
	popq	%r15
	.cfi_def_cfa_offset 8
	ret
.L1985:
	.cfi_restore_state
	movq	16(%r8), %r12
	movq	(%r8), %rax
	movq	%r13, %xmm3
	movq	%r12, %xmm0
	movq	%rax, 96(%rsp)
	punpcklqdq	%xmm3, %xmm0
	movaps	%xmm0, 112(%rsp)
	cmpq	16(%r12), %rdx
	jne	.L1987
	movq	$0, 176(%rsp)
	pxor	%xmm0, %xmm0
	movups	%xmm0, 184(%rsp)
	movq	16(%r14), %rsi
	movq	16(%r12), %rdx
	movq	%rsi, %rax
	orq	%rdx, %rax
	je	.L1930
	leaq	176(%rsp), %r15
	movq	%r15, %rdi
.LEHB80:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
	movq	8(%r12), %rdx
	movq	184(%rsp), %rax
	movq	192(%rsp), %rcx
	leaq	(%rdx,%rax), %rsi
	addq	%rcx, %rsi
	cmpq	$19, %rsi
	jg	.L1931
	testq	%rdx, %rdx
	jg	.L1954
.L1931:
	movq	%rax, %rdi
	orq	%rcx, %rdi
	js	.L1988
	imulq	%rcx, %rax
	testq	%rax, %rax
	je	.L1937
	salq	$4, %rax
	movq	%rax, %rdx
	je	.L1937
	movq	176(%rsp), %rdi
	xorl	%esi, %esi
	call	memset@PLT
.L1937:
	leaq	48(%rsp), %rcx
	leaq	96(%rsp), %rsi
	movq	%r12, %rdx
	movq	%r15, %rdi
	movq	.LC86(%rip), %xmm0
	movaps	%xmm0, 48(%rsp)
	call	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE13scaleAndAddToIS9_EEvRT_RKSD_RSA_RKS5_
.LEHE80:
.L1933:
	movq	%r13, %xmm0
	movq	184(%rsp), %rax
	movq	152(%rsp), %r9
	movhps	176(%rsp), %xmm0
	movups	%xmm0, 200(%rsp)
	movq	8(%r13), %xmm0
	movq	%rax, 216(%rsp)
	movq	0(%r13), %rax
	movhps	16(%r12), %xmm0
	movups	%xmm0, 232(%rsp)
	movq	16(%r14), %r14
	movq	16(%r13), %r13
	movq	%rax, 224(%rsp)
	cmpq	%r9, %r14
	jne	.L1938
	movq	160(%rsp), %r12
	cmpq	%r12, %r13
	je	.L1942
.L1938:
	leaq	144(%rsp), %rdi
	movq	%r13, %rdx
	movq	%r14, %rsi
.LEHB81:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
.LEHE81:
	movq	152(%rsp), %r9
	cmpq	%r9, %r14
	jne	.L1944
	movq	160(%rsp), %r12
	cmpq	%r12, %r13
	jne	.L1944
.L1942:
	movq	144(%rsp), %r10
	testq	%r12, %r12
	jle	.L1945
	movq	%r9, %r13
	movq	.LC8(%rip), %xmm3
	xorl	%r8d, %r8d
	xorl	%r11d, %r11d
	salq	$4, %r13
	.p2align 4,,10
	.p2align 3
.L1946:
	movq	%r8, %rdi
	salq	$4, %rdi
	cmpq	%r8, %r9
	jle	.L1951
	.p2align 4,,10
	.p2align 3
.L1949:
	movq	240(%rsp), %rsi
	testq	%rsi, %rsi
	jle	.L1956
	movq	232(%rsp), %rax
	movq	216(%rsp), %r14
	pxor	%xmm2, %xmm2
	xorl	%edx, %edx
	movq	208(%rsp), %rcx
	imulq	%r11, %rax
	salq	$4, %r14
	addq	%rdi, %rcx
	salq	$4, %rax
	addq	224(%rsp), %rax
	.p2align 4,,10
	.p2align 3
.L1948:
	movdqu	(%rcx), %xmm5
	movdqu	(%rax), %xmm4
	addq	$1, %rdx
	addq	$16, %rax
	movupd	-16(%rax), %xmm7
	addq	%r14, %rcx
	pshufd	$238, %xmm5, %xmm0
	pshufd	$78, %xmm4, %xmm1
	mulpd	%xmm0, %xmm1
	pshufd	$68, %xmm5, %xmm0
	mulpd	%xmm7, %xmm0
	xorpd	%xmm3, %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm0, %xmm2
	cmpq	%rdx, %rsi
	jne	.L1948
.L1947:
	addq	$1, %r8
	movaps	%xmm2, (%r10,%rdi)
	addq	$16, %rdi
	cmpq	%r8, %r9
	jne	.L1949
.L1951:
	xorl	%r8d, %r8d
	testq	%r9, %r9
	cmovle	%r9, %r8
	addq	$1, %r11
	addq	%r13, %r10
	cmpq	%r12, %r11
	jne	.L1946
.L1945:
	movq	176(%rsp), %rdi
	call	free@PLT
	jmp	.L1952
	.p2align 4,,10
	.p2align 3
.L1956:
	pxor	%xmm2, %xmm2
	jmp	.L1947
.L1930:
	movq	8(%r12), %rdx
	leaq	176(%rsp), %r15
	leaq	-1(%rdx), %rax
	cmpq	$18, %rax
	ja	.L1937
.L1954:
	movq	%r14, 64(%rsp)
	movq	%r12, 80(%rsp)
	cmpq	%rdx, 8(%r14)
	jne	.L1989
	leaq	47(%rsp), %rdx
	leaq	64(%rsp), %rsi
	leaq	176(%rsp), %rdi
.LEHB82:
	call	_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_
.LEHE82:
	jmp	.L1933
.L1987:
	leaq	.LC94(%rip), %rcx
	movl	$96, %edx
	leaq	.LC84(%rip), %rsi
	leaq	.LC85(%rip), %rdi
	call	__assert_fail@PLT
.L1986:
	call	__stack_chk_fail@PLT
.L1919:
	leaq	.LC93(%rip), %rcx
	movl	$118, %edx
	leaq	.LC4(%rip), %rsi
	leaq	.LC5(%rip), %rdi
	call	__assert_fail@PLT
.L1989:
	leaq	.LC83(%rip), %rcx
	movl	$96, %edx
	leaq	.LC84(%rip), %rsi
	leaq	.LC85(%rip), %rdi
	call	__assert_fail@PLT
.L1944:
	leaq	.LC95(%rip), %rcx
	movl	$765, %edx
	leaq	.LC96(%rip), %rsi
	leaq	.LC97(%rip), %rdi
	call	__assert_fail@PLT
.L1988:
	call	_ZN5Eigen9DenseBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8ConstantEllRKS3_.part.0
.L1959:
	movq	%rax, %rbx
	jmp	.L1953
.L1960:
	movq	%rax, %rbx
	jmp	.L1940
.L1957:
	movq	%rax, %rbx
	jmp	.L1941
	.section	.gcc_except_table
.LLSDA10576:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE10576-.LLSDACSB10576
.LLSDACSB10576:
	.uleb128 .LEHB79-.LFB10576
	.uleb128 .LEHE79-.LEHB79
	.uleb128 .L1957-.LFB10576
	.uleb128 0
	.uleb128 .LEHB80-.LFB10576
	.uleb128 .LEHE80-.LEHB80
	.uleb128 .L1960-.LFB10576
	.uleb128 0
	.uleb128 .LEHB81-.LFB10576
	.uleb128 .LEHE81-.LEHB81
	.uleb128 .L1959-.LFB10576
	.uleb128 0
	.uleb128 .LEHB82-.LFB10576
	.uleb128 .LEHE82-.LEHB82
	.uleb128 .L1960-.LFB10576
	.uleb128 0
.LLSDACSE10576:
	.text
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDAC10576
	.type	_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll.cold, @function
_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll.cold:
.LFSB10576:
.L1953:
	.cfi_def_cfa_offset 320
	.cfi_offset 3, -56
	.cfi_offset 6, -48
	.cfi_offset 12, -40
	.cfi_offset 13, -32
	.cfi_offset 14, -24
	.cfi_offset 15, -16
	movq	176(%rsp), %rdi
	call	free@PLT
.L1941:
	movq	144(%rsp), %rdi
	call	free@PLT
	movq	%rbx, %rdi
.LEHB83:
	call	_Unwind_Resume@PLT
.LEHE83:
.L1940:
	movq	176(%rsp), %rdi
	call	free@PLT
	jmp	.L1941
.L1981:
.LEHB84:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.L1980:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.L1979:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE84:
.L1958:
	movq	144(%rsp), %rdi
	movq	%rax, %rbx
	call	free@PLT
	movq	%rbx, %rdi
.LEHB85:
	call	_Unwind_Resume@PLT
.LEHE85:
	.cfi_endproc
.LFE10576:
	.section	.gcc_except_table
.LLSDAC10576:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSEC10576-.LLSDACSBC10576
.LLSDACSBC10576:
	.uleb128 .LEHB83-.LCOLDB98
	.uleb128 .LEHE83-.LEHB83
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB84-.LCOLDB98
	.uleb128 .LEHE84-.LEHB84
	.uleb128 .L1958-.LCOLDB98
	.uleb128 0
	.uleb128 .LEHB85-.LCOLDB98
	.uleb128 .LEHE85-.LEHB85
	.uleb128 0
	.uleb128 0
.LLSDACSEC10576:
	.section	.text.unlikely
	.text
	.size	_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll, .-_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll
	.section	.text.unlikely
	.size	_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll.cold, .-_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll.cold
.LCOLDE98:
	.text
.LHOTE98:
	.section	.rodata.str1.1
.LC99:
	.string	"x"
.LC100:
	.string	"y"
.LC101:
	.string	"z"
.LC102:
	.string	":: Basic infos\n"
	.section	.rodata.str1.8
	.align 8
.LC103:
	.string	"      number of atoms:                "
	.section	.rodata.str1.1
.LC104:
	.string	"\n"
	.section	.rodata.str1.8
	.align 8
.LC105:
	.string	"      number of occupied spinors:     "
	.align 8
.LC106:
	.string	"      number of virtual spinors:      "
	.align 8
.LC107:
	.string	"      strength of magnetic field:     "
	.align 8
.LC108:
	.string	"      B field vector (unscaled):   "
	.section	.rodata.str1.1
.LC109:
	.string	"   "
	.section	.rodata.str1.8
	.align 8
.LC110:
	.string	"\n\ncalculating orbital rotation matrix\n\n"
	.section	.rodata.str1.1
.LC111:
	.string	" :: calculating rhs...   "
.LC112:
	.string	"done!\n"
	.section	.rodata.str1.8
	.align 8
.LC113:
	.string	" :: solving CPHF equation...   "
	.section	.rodata.str1.1
.LC114:
	.string	"      saving to disk...   "
.LC115:
	.string	"u"
.LC116:
	.string	"_"
	.section	.rodata.str1.8
	.align 8
.LC117:
	.string	"\n\ndone calculating orbital rotation matrix\n\n"
	.section	.rodata.str1.1
.LC118:
	.string	"\nsplitting sbraket files...\n"
.LC119:
	.string	"b"
.LC120:
	.string	"k"
	.section	.rodata.str1.8
	.align 8
.LC122:
	.ascii	"Eigen::Product<Lhs, Rhs, Option>::Product(const Lhs&, const "
	.ascii	"Rhs&) [with _Lhs = Eigen::CwiseUnaryOp<Eigen::internal::scal"
	.ascii	"ar_conjugate_op<std::complex<double> >, const Eigen::Transpo"
	.ascii	"se<const Eigen::Matrix<std::complex<double>, -1, -1> > >; _R"
	.ascii	"hs = Eigen::Matrix<std::c"
	.string	"omplex<double>, -1, -1>; int Option = 0; Lhs = Eigen::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std::complex<double> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >; Rhs = Eigen::Matrix<std::complex<double>, -1, -1>]"
	.align 8
.LC123:
	.ascii	"Eigen::Product<Lhs, Rhs, Option>::Product(const Lhs&, const "
	.ascii	"Rhs&) [with _Lhs = Eigen::Product<Eigen::CwiseUnaryOp<Eigen:"
	.ascii	":internal::scalar_conjugate_op<std::complex<double> >, const"
	.ascii	" Eigen::Transpose<const Eigen::Matrix<std::complex<double>, "
	.ascii	"-1, -1> > >, Eigen::Matrix<std::complex<double>, -1, -1>, 0>"
	.ascii	"; _Rhs = Eigen::Matrix<std::complex<double>, -1, -1>; int Op"
	.ascii	"tion = 0; Lhs = Eigen::Product<Ei"
	.string	"gen::CwiseUnaryOp<Eigen::internal::scalar_conjugate_op<std::complex<double> >, const Eigen::Transpose<const Eigen::Matrix<std::complex<double>, -1, -1> > >, Eigen::Matrix<std::complex<double>, -1, -1>, 0>; Rhs = Eigen::Matrix<std::complex<double>, -1, -1>]"
	.section	.rodata.str1.1
.LC124:
	.string	"Ia "
.LC125:
	.string	"   Jb "
.LC126:
	.string	"bk"
	.section	.rodata.str1.8
	.align 8
.LC129:
	.string	"Eigen::DenseCoeffsBase<Derived, 1>::Scalar& Eigen::DenseCoeffsBase<Derived, 1>::operator()(Eigen::Index) [with Derived = Eigen::Matrix<std::complex<double>, -1, 1>; Scalar = std::complex<double>; Eigen::Index = long int]"
	.section	.rodata.str1.1
.LC130:
	.string	"index >= 0 && index < size()"
.LC132:
	.string	"\n\n\nBerry-curvature term 1:\n"
.LC133:
	.string	"\n\n"
.LC134:
	.string	"\n\n\nBerry-curvature term 2:\n"
.LC135:
	.string	"\n\n\nBerry-curvature term 3:\n"
.LC136:
	.string	"\n\n\nBerry-curvature term 4:\n"
.LC137:
	.string	"\n\n\nBerry-curvature term 5:\n"
.LC138:
	.string	"\n\n\nBerry-curvature total:\n"
	.section	.rodata.str1.8
	.align 8
.LC139:
	.string	"================================================================================\n"
	.align 8
.LC140:
	.string	"Eigen::CwiseNullaryOp<NullaryOp, MatrixType>::CwiseNullaryOp(Eigen::Index, Eigen::Index, const NullaryOp&) [with NullaryOp = Eigen::internal::scalar_constant_op<double>; PlainObjectType = Eigen::Matrix<double, -1, -1>; Eigen::Index = long int]"
	.align 8
.LC141:
	.string	"Eigen::Block<XprType, BlockRows, BlockCols, InnerPanel>::Block(XprType&, Eigen::Index, Eigen::Index) [with XprType = Eigen::Matrix<std::complex<double>, -1, -1>; int BlockRows = 3; int BlockCols = 3; bool InnerPanel = false; Eigen::Index = long int]"
	.align 8
.LC142:
	.string	"startRow >= 0 && BlockRows >= 0 && startRow + BlockRows <= xpr.rows() && startCol >= 0 && BlockCols >= 0 && startCol + BlockCols <= xpr.cols()"
	.align 8
.LC143:
	.ascii	"Eigen::CwiseBinaryOp<BinaryOp, Lhs, Rhs>::CwiseBinaryOp(cons"
	.ascii	"t Lhs&, const Rhs&, const BinaryOp&) [with BinaryOp = Eigen:"
	.ascii	":internal:"
	.string	":scalar_difference_op<double, double>; LhsType = const Eigen::Matrix<double, -1, -1>; RhsType = const Eigen::Transpose<const Eigen::Matrix<double, -1, -1> >; Lhs = Eigen::Matrix<double, -1, -1>; Rhs = Eigen::Transpose<const Eigen::Matrix<double, -1, -1> >]"
	.align 8
.LC144:
	.string	"Eigen::CwiseNullaryOp<NullaryOp, MatrixType>::CwiseNullaryOp(Eigen::Index, Eigen::Index, const NullaryOp&) [with NullaryOp = Eigen::internal::scalar_constant_op<double>; PlainObjectType = const Eigen::Matrix<double, -1, -1>; Eigen::Index = long int]"
	.section	.rodata.str1.1
.LC146:
	.string	"\n\ncharge fluctuation:\n"
.LC147:
	.string	"coord"
.LC148:
	.string	"\npartial charges:\n"
	.section	.rodata.str1.8
	.align 8
.LC149:
	.string	"Atom\telec. charge \tnuc. charge\ttotal\n"
	.section	.rodata.str1.1
.LC150:
	.string	" %u%s\t%10.7f\t%10.7f\t%10.7f\n"
	.section	.rodata.str1.8
	.align 8
.LC151:
	.string	"-----------------------------------------------------\n"
	.section	.rodata.str1.1
.LC152:
	.string	" sum\t%10.7f\t%10.7f\t%10.7f\n"
	.section	.rodata.str1.8
	.align 8
.LC153:
	.string	"\n\n =================== time stats ===================\n"
	.align 8
.LC155:
	.string	"   Calculate electronic Hessian: %.3fs\n"
	.align 8
.LC156:
	.string	"   Calculate CPHF:               %.3fs\n"
	.align 8
.LC157:
	.string	"   Calculate Berry-Curvature:    %.3fs\n"
	.align 8
.LC158:
	.string	"   --------------------------------------------------------"
	.align 8
.LC159:
	.string	"   Total:                        %.3fs\n\n"
	.section	.rodata.str1.1
.LC160:
	.string	"\n\n\nf-Debug"
.LC161:
	.string	"smat"
.LC162:
	.string	"smatSAO:\n"
.LC163:
	.string	"smatcao"
.LC164:
	.string	"smat:\n"
	.section	.rodata
	.align 8
.LC121:
	.long	0
	.long	0
	.long	0
	.long	0
	.align 8
.LC131:
	.long	0
	.long	1072693248
	.long	0
	.long	0
	.section	.text.unlikely
.LCOLDB166:
	.section	.text.startup,"ax",@progbits
.LHOTB166:
	.p2align 4
	.globl	main
	.type	main, @function
main:
.LFB9961:
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDA9961
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	leaq	.LC99(%rip), %rsi
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	subq	$4744, %rsp
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	%fs:40, %rax
	movq	%rax, -56(%rbp)
	xorl	%eax, %eax
	leaq	-1056(%rbp), %rax
	movq	%rax, %rdi
	movq	%rax, -4600(%rbp)
.LEHB86:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC2IS3_EEPKcRKS3_.constprop.3
.LEHE86:
	leaq	-1024(%rbp), %rdi
	leaq	.LC100(%rip), %rsi
.LEHB87:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC2IS3_EEPKcRKS3_.constprop.3
.LEHE87:
	leaq	-992(%rbp), %rdi
	leaq	.LC101(%rip), %rsi
.LEHB88:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC2IS3_EEPKcRKS3_.constprop.3
.LEHE88:
	subq	$8, %rsp
	leaq	-4232(%rbp), %rax
	leaq	-4256(%rbp), %rcx
	pushq	%rax
	leaq	-4264(%rbp), %rdx
	leaq	-4268(%rbp), %rsi
	leaq	-4272(%rbp), %rdi
	leaq	-4240(%rbp), %r9
	leaq	-4248(%rbp), %r8
.LEHB89:
	.cfi_escape 0x2e,0x10
	call	_Z4infoRiS_S_RdS0_S0_S0_@PLT
.LEHE89:
	popq	%rax
	leaq	.LC102(%rip), %rdi
	popq	%rdx
.LEHB90:
	.cfi_escape 0x2e,0
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
.LEHE90:
	movl	$38, %edx
	leaq	.LC103(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
.LEHB91:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
.LEHE91:
	movl	-4272(%rbp), %esi
	leaq	_ZSt4cout(%rip), %rdi
.LEHB92:
	call	_ZNSolsEi@PLT
.LEHE92:
	movq	%rax, %rdi
	leaq	.LC104(%rip), %rsi
.LEHB93:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
.LEHE93:
	movl	$38, %edx
	leaq	.LC105(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
.LEHB94:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
.LEHE94:
	movl	-4268(%rbp), %esi
	leaq	_ZSt4cout(%rip), %rdi
.LEHB95:
	call	_ZNSolsEi@PLT
.LEHE95:
	movq	%rax, %rdi
	leaq	.LC104(%rip), %rsi
.LEHB96:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
.LEHE96:
	movl	$38, %edx
	leaq	.LC106(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
.LEHB97:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
.LEHE97:
	movl	-4264(%rbp), %esi
	leaq	_ZSt4cout(%rip), %rdi
.LEHB98:
	call	_ZNSolsEi@PLT
.LEHE98:
	movq	%rax, %rdi
	leaq	.LC104(%rip), %rsi
.LEHB99:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
.LEHE99:
	movl	$38, %edx
	leaq	.LC107(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
.LEHB100:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
.LEHE100:
	movsd	-4232(%rbp), %xmm0
	leaq	_ZSt4cout(%rip), %rdi
.LEHB101:
	call	_ZNSo9_M_insertIdEERSoT_@PLT
.LEHE101:
	movq	%rax, %rdi
	leaq	.LC104(%rip), %rsi
.LEHB102:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
.LEHE102:
	movl	$35, %edx
	leaq	.LC108(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
.LEHB103:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
.LEHE103:
	movsd	-4256(%rbp), %xmm0
	leaq	_ZSt4cout(%rip), %rdi
.LEHB104:
	call	_ZNSo9_M_insertIdEERSoT_@PLT
.LEHE104:
	movl	$3, %edx
	leaq	.LC109(%rip), %rsi
	movq	%rax, %rdi
	movq	%rax, %rbx
.LEHB105:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
.LEHE105:
	movsd	-4248(%rbp), %xmm0
	movq	%rbx, %rdi
.LEHB106:
	call	_ZNSo9_M_insertIdEERSoT_@PLT
.LEHE106:
	movl	$3, %edx
	leaq	.LC109(%rip), %rsi
	movq	%rax, %rdi
	movq	%rax, %rbx
.LEHB107:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
.LEHE107:
	movsd	-4240(%rbp), %xmm0
	movq	%rbx, %rdi
.LEHB108:
	call	_ZNSo9_M_insertIdEERSoT_@PLT
.LEHE108:
	movq	%rax, %rdi
	leaq	.LC104(%rip), %rsi
.LEHB109:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
.LEHE109:
	leaq	-4064(%rbp), %rax
	pxor	%xmm0, %xmm0
	leaq	-4096(%rbp), %rsi
	movq	$0, -4080(%rbp)
	movq	%rax, %rdi
	movq	%rax, -4392(%rbp)
	movaps	%xmm0, -4096(%rbp)
.LEHB110:
	call	_Z10readSpinorRSt6vectorIdSaIdEE@PLT
.LEHE110:
	movl	-4264(%rbp), %eax
	addl	-4268(%rbp), %eax
	movl	%eax, -4532(%rbp)
	movl	%eax, -4260(%rbp)
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	pxor	%xmm0, %xmm0
	leaq	-4000(%rbp), %rsi
	movq	$0, -4032(%rbp)
	movq	%rax, -4728(%rbp)
	leaq	-4032(%rbp), %rax
	movq	%rax, %rdi
	movq	%rsi, -4352(%rbp)
	movq	$0, -4000(%rbp)
	movq	%rax, -4336(%rbp)
	movups	%xmm0, -4024(%rbp)
	movups	%xmm0, -3992(%rbp)
.LEHB111:
	call	_Z11calcStabmatRN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEES4_@PLT
	movq	%rax, %rsi
	leaq	-4224(%rbp), %rax
	movq	%rax, %rdi
	movq	%rax, -4368(%rbp)
	call	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC1ERKS3_
.LEHE111:
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	leaq	-960(%rbp), %rbx
	movq	%rax, -4736(%rbp)
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	leaq	.LC110(%rip), %rdi
	movq	%rax, -4744(%rbp)
.LEHB112:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	movl	-4272(%rbp), %eax
	xorl	%r15d, %r15d
	testl	%eax, %eax
	jle	.L2040
.L2039:
	xorl	%r13d, %r13d
.L2046:
	leaq	.LC111(%rip), %rdi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	leaq	-1840(%rbp), %rax
	movq	-4368(%rbp), %rcx
	movl	%r13d, %edx
	movl	%r15d, %esi
	movq	%rax, %rdi
	movq	%rax, %r12
	call	_Z8berryRHSiiRKN5Eigen6MatrixISt7complexIdELin1ELi1ELi0ELin1ELi1EEE@PLT
.LEHE112:
	leaq	.LC112(%rip), %rdi
.LEHB113:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	leaq	.LC113(%rip), %rdi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	leaq	-1616(%rbp), %r14
	movq	%r12, %rsi
	movq	%r14, %rdi
	call	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC1ERKS3_
.LEHE113:
	leaq	-1376(%rbp), %r12
	movq	-4352(%rbp), %rsi
	movq	%r12, %rdi
.LEHB114:
	call	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC1ERKS3_
.LEHE114:
	leaq	-1504(%rbp), %rax
	movq	-4336(%rbp), %rsi
	movq	%rax, %rdi
	movq	%rax, -4376(%rbp)
.LEHB115:
	call	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELin1ELi0EEC1ERKS3_
.LEHE115:
	leaq	-1728(%rbp), %rax
	movq	-4376(%rbp), %rsi
	movq	%r14, %rcx
	movq	%r12, %rdx
	movq	%rax, %rdi
	movq	%rax, -4312(%rbp)
.LEHB116:
	call	_Z4cphfN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEES3_NS0_IS2_Lin1ELi1ELi0ELin1ELi1EEE@PLT
.LEHE116:
	movq	-1504(%rbp), %rdi
	call	free@PLT
	movq	-1376(%rbp), %rdi
	call	free@PLT
	movq	-1616(%rbp), %rdi
	call	free@PLT
	leaq	.LC112(%rip), %rdi
.LEHB117:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	leaq	.LC114(%rip), %rdi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
.LEHE117:
	movl	%r13d, %esi
	movq	%rbx, %rdi
	leaq	-1152(%rbp), %r14
	call	_ZNSt7__cxx119to_stringEi
	movl	%r15d, %esi
	movq	%r14, %rdi
	call	_ZNSt7__cxx119to_stringEi
	leaq	.LC115(%rip), %rdx
	xorl	%esi, %esi
	movq	%r14, %rdi
.LEHB118:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@PLT
.LEHE118:
	leaq	-1120(%rbp), %r14
	movq	%rax, %rsi
	movq	%r14, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	leaq	.LC116(%rip), %rsi
	movq	%r14, %rdi
.LEHB119:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendEPKc@PLT
.LEHE119:
	leaq	-1088(%rbp), %r14
	movq	%rax, %rsi
	movq	%r14, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	leaq	-576(%rbp), %rax
	movq	%r14, %rsi
	movq	%rbx, %rdx
	movq	%rax, %rdi
	movq	%rax, %r14
.LEHB120:
	call	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_
.LEHE120:
	movq	-4312(%rbp), %rsi
	movq	%r12, %rdi
.LEHB121:
	call	_ZN5Eigen12DenseStorageISt7complexIdELin1ELin1ELi1ELi0EEC1ERKS3_
.LEHE121:
	movq	%r14, %rsi
	movq	%r12, %rdi
.LEHB122:
	call	_Z10saveVectorN5Eigen6MatrixISt7complexIdELin1ELi1ELi0ELin1ELi1EEENSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE122:
	movq	-1376(%rbp), %rdi
	call	free@PLT
	movq	-576(%rbp), %rdi
	leaq	-560(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2041
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2041:
	movq	-1088(%rbp), %rdi
	leaq	-1072(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2042
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2042:
	movq	-1120(%rbp), %rdi
	leaq	-1104(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2043
	movq	-1104(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2043:
	movq	-1152(%rbp), %rdi
	leaq	-1136(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2044
	movq	-1136(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2044:
	movq	-960(%rbp), %rdi
	leaq	-944(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2045
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2045:
	leaq	.LC112(%rip), %rdi
.LEHB123:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
.LEHE123:
	movq	-1728(%rbp), %rdi
	addl	$1, %r13d
	call	free@PLT
	movq	-1840(%rbp), %rdi
	call	free@PLT
	cmpl	$3, %r13d
	jne	.L2046
	addl	$1, %r15d
	cmpl	%r15d, -4272(%rbp)
	jg	.L2039
.L2040:
	leaq	.LC117(%rip), %rdi
	leaq	-960(%rbp), %rbx
.LEHB124:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt5flushIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	movq	%rax, -4752(%rbp)
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	movl	$28, %edx
	leaq	.LC118(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	movq	%rax, -4760(%rbp)
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt5flushIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0
	movl	-4272(%rbp), %edi
	call	_Z11splitBraKeti@PLT
	movl	$6, %edx
	leaq	.LC112(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt5flushIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0
	leaq	-1376(%rbp), %rbx
	movl	-4272(%rbp), %eax
	leaq	-1504(%rbp), %r14
	movq	%rbx, %r15
	leaq	-3968(%rbp), %rdi
	movq	%r14, %rsi
	movq	%rbx, -4368(%rbp)
	leal	(%rax,%rax,2), %eax
	movq	%r15, %rdx
	movq	%r14, -4376(%rbp)
	leaq	-960(%rbp), %rbx
	movl	%eax, -1376(%rbp)
	movl	%eax, -1504(%rbp)
	movq	%rdi, -4640(%rbp)
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE124:
	movl	-4272(%rbp), %eax
	leaq	-3936(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r14, %rsi
	movq	%rdi, -4488(%rbp)
	leal	(%rax,%rax,2), %eax
	movl	%eax, -1376(%rbp)
	movl	%eax, -1504(%rbp)
.LEHB125:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE125:
	movl	-4272(%rbp), %eax
	leaq	-3904(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r14, %rsi
	movq	%rdi, -4496(%rbp)
	leal	(%rax,%rax,2), %eax
	movl	%eax, -1376(%rbp)
	movl	%eax, -1504(%rbp)
.LEHB126:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE126:
	movl	-4272(%rbp), %eax
	leaq	-3872(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r14, %rsi
	movq	%rdi, -4504(%rbp)
	leal	(%rax,%rax,2), %eax
	movl	%eax, -1376(%rbp)
	movl	%eax, -1504(%rbp)
.LEHB127:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE127:
	movl	-4272(%rbp), %eax
	leaq	-3840(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r14, %rsi
	movq	%rdi, -4544(%rbp)
	leal	(%rax,%rax,2), %eax
	movl	%eax, -1376(%rbp)
	movl	%eax, -1504(%rbp)
.LEHB128:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE128:
	movl	-4272(%rbp), %eax
	leaq	-3808(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r14, %rsi
	leal	(%rax,%rax,2), %eax
	movl	%eax, -1376(%rbp)
	movl	%eax, -1504(%rbp)
.LEHB129:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE129:
	movl	-4272(%rbp), %eax
	leaq	-3776(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r14, %rsi
	leal	(%rax,%rax,2), %eax
	movl	%eax, -1376(%rbp)
	movl	%eax, -1504(%rbp)
.LEHB130:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE130:
	movl	-4272(%rbp), %eax
	movq	%r15, %rdx
	movq	%r14, %rsi
	leal	(%rax,%rax,2), %eax
	movl	%eax, -1376(%rbp)
	movl	%eax, -1504(%rbp)
	leaq	-3744(%rbp), %rax
	movq	%rax, %rdi
	movq	%rax, -4664(%rbp)
.LEHB131:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE131:
	movl	-4272(%rbp), %r14d
	testl	%r14d, %r14d
	jle	.L2048
	movq	$0, -4608(%rbp)
	leaq	-960(%rbp), %rbx
.L2047:
	movl	-4608(%rbp), %eax
	movq	$0, -4616(%rbp)
	movl	%eax, -4628(%rbp)
	movq	-4600(%rbp), %rax
	movq	%rax, -4648(%rbp)
.L2235:
	movl	-4628(%rbp), %r15d
	movl	-4616(%rbp), %r14d
	leaq	_ZSt4cout(%rip), %rdi
	movl	%r15d, %esi
	movl	%r14d, -4632(%rbp)
.LEHB132:
	call	_ZNSolsEi@PLT
	movq	%rax, %rdi
	movl	%r14d, %esi
	call	_ZNSolsEi@PLT
	movq	%rax, %rdi
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0
.LEHE132:
	movl	%r14d, %esi
	movq	%rbx, %rdi
	call	_ZNSt7__cxx119to_stringEi
	leaq	-1152(%rbp), %rax
	movl	%r15d, %esi
	movq	%rax, %rdi
	movq	%rax, %r15
	movq	%rax, -4672(%rbp)
	call	_ZNSt7__cxx119to_stringEi
	leaq	.LC115(%rip), %rdx
	xorl	%esi, %esi
	movq	%r15, %rdi
.LEHB133:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@PLT
.LEHE133:
	movq	%rax, %rsi
	leaq	-1120(%rbp), %rax
	movq	%rax, %r15
	movq	%rax, %rdi
	movq	%rax, -4720(%rbp)
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	leaq	.LC116(%rip), %rsi
	movq	%r15, %rdi
.LEHB134:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendEPKc@PLT
.LEHE134:
	leaq	-1088(%rbp), %r14
	movq	%rax, %rsi
	movq	%r14, %rdi
	movq	%r14, -4552(%rbp)
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	leaq	-576(%rbp), %rax
	movq	%rbx, %rdx
	movq	%r14, %rsi
	movq	%rax, %rdi
	movq	%rax, -4464(%rbp)
	movq	%rax, %r15
.LEHB135:
	call	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_
.LEHE135:
	leaq	-4208(%rbp), %rdi
	movq	%r15, %rsi
.LEHB136:
	call	_Z10readVectorNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE136:
	movq	-576(%rbp), %rdi
	leaq	-560(%rbp), %rax
	movq	%rax, -4456(%rbp)
	cmpq	%rax, %rdi
	je	.L2049
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2049:
	movq	-1088(%rbp), %rdi
	leaq	-1072(%rbp), %rax
	movq	%rax, -4400(%rbp)
	cmpq	%rax, %rdi
	je	.L2050
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2050:
	movq	-1120(%rbp), %rdi
	leaq	-1104(%rbp), %rax
	movq	%rax, -4568(%rbp)
	cmpq	%rax, %rdi
	je	.L2051
	movq	-1104(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2051:
	movq	-1152(%rbp), %rdi
	leaq	-1136(%rbp), %rax
	movq	%rax, -4560(%rbp)
	cmpq	%rax, %rdi
	je	.L2052
	movq	-1136(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2052:
	movq	-960(%rbp), %rdi
	leaq	-944(%rbp), %rax
	movq	%rax, -4384(%rbp)
	cmpq	%rax, %rdi
	je	.L2053
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2053:
	movq	-4552(%rbp), %r15
	movl	-4628(%rbp), %esi
	movq	%r15, %rdi
	call	_ZNSt7__cxx119to_stringEi
	leaq	.LC119(%rip), %rdx
	xorl	%esi, %esi
	movq	%r15, %rdi
.LEHB137:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@PLT
.LEHE137:
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	-4648(%rbp), %rax
	movq	%rbx, %rdi
	movq	8(%rax), %rdx
	movq	(%rax), %rsi
.LEHB138:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendEPKcm@PLT
.LEHE138:
	movq	-4464(%rbp), %r15
	movq	%rax, %rsi
	leaq	-3712(%rbp), %r13
	movq	%r15, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	%r15, %rsi
	movq	%r13, %rdi
.LEHB139:
	call	_Z19readMatrixTransformNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE139:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2054
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2054:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2055
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2055:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2056
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2056:
	movq	-4552(%rbp), %r15
	movl	-4628(%rbp), %esi
	movq	%r15, %rdi
	call	_ZNSt7__cxx119to_stringEi
	leaq	.LC120(%rip), %rdx
	xorl	%esi, %esi
	movq	%r15, %rdi
.LEHB140:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@PLT
.LEHE140:
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	-4648(%rbp), %rax
	movq	%rbx, %rdi
	movq	8(%rax), %rdx
	movq	(%rax), %rsi
.LEHB141:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendEPKcm@PLT
.LEHE141:
	movq	-4464(%rbp), %r14
	movq	%rax, %rsi
	leaq	-3680(%rbp), %r15
	movq	%r14, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	%r14, %rsi
	movq	%r15, %rdi
.LEHB142:
	call	_Z19readMatrixTransformNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE142:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2057
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2057:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2058
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2058:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2059
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2059:
	leaq	-4260(%rbp), %rax
	leaq	-3648(%rbp), %rdi
	movq	%rax, %rsi
	movq	%rax, %rdx
	movq	%rax, -4704(%rbp)
	movq	%rdi, %r12
	movq	%rdi, -4512(%rbp)
.LEHB143:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE143:
	leaq	-1616(%rbp), %r14
	movq	%r12, %rsi
	movq	%r13, %rdx
	movq	%r14, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC1IS4_EERS4_RKNS_9DenseBaseIT_EE
	movl	-4532(%rbp), %esi
	movl	%esi, %eax
	shrl	$31, %eax
	addl	%esi, %eax
	sarl	%eax
	cltq
	movq	%rax, %xmm7
	punpcklqdq	%xmm7, %xmm7
	movaps	%xmm7, -4528(%rbp)
	movaps	%xmm7, -1504(%rbp)
	movsd	.LC121(%rip), %xmm7
	movsd	%xmm7, -1488(%rbp)
	movsd	8+.LC121(%rip), %xmm7
	movsd	%xmm7, -1480(%rbp)
	cmpl	$-1, %esi
	jl	.L2615
	movq	-4512(%rbp), %xmm7
	movq	-4376(%rbp), %rsi
	movq	%r14, %rdi
	movhps	-4392(%rbp), %xmm7
	movaps	%xmm7, -4688(%rbp)
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	movdqa	-4528(%rbp), %xmm7
	movq	-4368(%rbp), %rsi
	movq	%rax, %rdi
	movaps	%xmm7, -1376(%rbp)
	movsd	.LC121(%rip), %xmm7
	movsd	%xmm7, -1360(%rbp)
	movsd	8+.LC121(%rip), %xmm7
	movsd	%xmm7, -1352(%rbp)
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	movq	%r13, %rsi
	movq	%rax, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmIS4_EERS5_RKNS_9DenseBaseIT_EE.isra.0
	movq	%r14, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8finishedEv.isra.0
	movq	-4056(%rbp), %rax
	cmpq	-3640(%rbp), %rax
	jne	.L2063
	movq	-4392(%rbp), %rsi
	movdqa	-4688(%rbp), %xmm7
	movq	%rsi, -3328(%rbp)
	movaps	%xmm7, -3312(%rbp)
	cmpq	%rax, -3632(%rbp)
	jne	.L2064
	movq	-4704(%rbp), %rsi
	leaq	-3616(%rbp), %r13
	movq	%r13, %rdi
	movq	%rsi, %rdx
.LEHB144:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE144:
	movq	%r15, %rdx
	movq	%r13, %rsi
	movq	%r14, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC1IS4_EERS4_RKNS_9DenseBaseIT_EE
	movsd	.LC121(%rip), %xmm4
	movdqa	-4528(%rbp), %xmm7
	movq	%r14, %rdi
	movq	-4376(%rbp), %rsi
	movsd	%xmm4, -1488(%rbp)
	movsd	8+.LC121(%rip), %xmm4
	movaps	%xmm7, -1504(%rbp)
	movsd	%xmm4, -1480(%rbp)
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	movdqa	-4528(%rbp), %xmm7
	movq	-4368(%rbp), %rsi
	movq	%rax, %rdi
	movaps	%xmm7, -1376(%rbp)
	movsd	.LC121(%rip), %xmm7
	movsd	%xmm7, -1360(%rbp)
	movsd	8+.LC121(%rip), %xmm7
	movsd	%xmm7, -1352(%rbp)
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	movq	%r15, %rsi
	movq	%rax, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmIS4_EERS5_RKNS_9DenseBaseIT_EE.isra.0
	movq	%r14, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8finishedEv.isra.0
	movq	-4056(%rbp), %rax
	cmpq	-3608(%rbp), %rax
	jne	.L2063
	cmpq	%rax, -3600(%rbp)
	jne	.L2064
	movl	-4272(%rbp), %r13d
	testl	%r13d, %r13d
	jle	.L2066
	movq	$0, -4712(%rbp)
	movl	$0, -4536(%rbp)
.L2065:
	movq	-4600(%rbp), %rax
	movq	$0, -4624(%rbp)
	movq	%rax, -4656(%rbp)
.L2085:
	movl	$3, %edx
	leaq	.LC124(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	movl	-4624(%rbp), %r15d
.LEHB145:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movl	-4628(%rbp), %esi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZNSolsEi@PLT
	movl	$1, %edx
	leaq	.LC116(%rip), %rsi
	movq	%rax, %rdi
	movq	%rax, %r13
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movl	-4632(%rbp), %esi
	movq	%r13, %rdi
	call	_ZNSolsEi@PLT
	movl	$6, %edx
	leaq	.LC125(%rip), %rsi
	movq	%rax, %rdi
	movq	%rax, %r13
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movl	-4536(%rbp), %r12d
	movq	%r13, %rdi
	movl	%r12d, %esi
	call	_ZNSolsEi@PLT
	movl	$1, %edx
	leaq	.LC116(%rip), %rsi
	movq	%rax, %rdi
	movq	%rax, %r13
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movl	%r15d, %esi
	movq	%r13, %rdi
	call	_ZNSolsEi@PLT
	movq	%rax, %rdi
	leaq	.LC104(%rip), %rsi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
.LEHE145:
	movl	%r15d, %esi
	movq	%rbx, %rdi
	call	_ZNSt7__cxx119to_stringEi
	movq	-4672(%rbp), %r15
	movl	%r12d, %esi
	movq	%r15, %rdi
	call	_ZNSt7__cxx119to_stringEi
	movl	$1, %r8d
	xorl	%edx, %edx
	xorl	%esi, %esi
	leaq	.LC115(%rip), %rcx
	movq	%r15, %rdi
.LEHB146:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE7replaceEmmPKcm@PLT
.LEHE146:
	movq	-4720(%rbp), %r15
	movq	%rax, %rsi
	movq	%r15, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	leaq	.LC116(%rip), %rsi
	movq	%r15, %rdi
.LEHB147:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendEPKc@PLT
.LEHE147:
	movq	-4552(%rbp), %r15
	movq	%rax, %rsi
	movq	%r15, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	%r15, %rsi
	movq	-4464(%rbp), %r15
	movq	%rbx, %rdx
	movq	%r15, %rdi
.LEHB148:
	call	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_
.LEHE148:
	leaq	-4192(%rbp), %rdi
	movq	%r15, %rsi
.LEHB149:
	call	_Z10readVectorNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE149:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2067
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2067:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2068
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2068:
	movq	-1120(%rbp), %rdi
	movq	-4568(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2069
	movq	-1104(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2069:
	movq	-1152(%rbp), %rdi
	movq	-4560(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2070
	movq	-1136(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2070:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2071
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2071:
	movq	-4552(%rbp), %r15
	movl	-4536(%rbp), %esi
	movq	%r15, %rdi
	call	_ZNSt7__cxx119to_stringEi
	leaq	.LC119(%rip), %rdx
	xorl	%esi, %esi
	movq	%r15, %rdi
.LEHB150:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@PLT
.LEHE150:
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	-4656(%rbp), %rax
	movq	%rbx, %rdi
	movq	8(%rax), %rdx
	movq	(%rax), %rsi
.LEHB151:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendEPKcm@PLT
.LEHE151:
	movq	-4464(%rbp), %r12
	movq	%rax, %rsi
	leaq	-3584(%rbp), %r15
	movq	%r12, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	%r12, %rsi
	movq	%r15, %rdi
.LEHB152:
	call	_Z19readMatrixTransformNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE152:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2072
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2072:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2073
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2073:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2074
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2074:
	movq	-4552(%rbp), %r12
	movl	-4536(%rbp), %esi
	movq	%r12, %rdi
	call	_ZNSt7__cxx119to_stringEi
	leaq	.LC120(%rip), %rdx
	xorl	%esi, %esi
	movq	%r12, %rdi
.LEHB153:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@PLT
.LEHE153:
	movq	%rax, %rsi
	movq	%rbx, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	-4656(%rbp), %rax
	movq	%rbx, %rdi
	movq	8(%rax), %rdx
	movq	(%rax), %rsi
.LEHB154:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendEPKcm@PLT
.LEHE154:
	movq	-4464(%rbp), %r12
	movq	%rax, %rsi
	leaq	-3552(%rbp), %r13
	movq	%r12, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	%r12, %rsi
	movq	%r13, %rdi
.LEHB155:
	call	_Z19readMatrixTransformNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE155:
	leaq	-3520(%rbp), %rax
	movq	-576(%rbp), %rdi
	movq	%rax, %xmm7
	movq	%rax, -4312(%rbp)
	movq	-4456(%rbp), %rax
	movhps	-4392(%rbp), %xmm7
	movaps	%xmm7, -4336(%rbp)
	cmpq	%rax, %rdi
	je	.L2075
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2075:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2076
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2076:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2077
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2077:
	movq	-4704(%rbp), %r12
	movq	-4312(%rbp), %rdi
	movq	%r12, %rdx
	movq	%r12, %rsi
.LEHB156:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE156:
	movq	-4312(%rbp), %rsi
	movq	%r13, %rdx
	movq	%r14, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC1IS4_EERS4_RKNS_9DenseBaseIT_EE
	movsd	.LC121(%rip), %xmm4
	movdqa	-4528(%rbp), %xmm7
	movq	%r14, %rdi
	movq	-4376(%rbp), %rsi
	movsd	%xmm4, -1488(%rbp)
	movsd	8+.LC121(%rip), %xmm4
	movaps	%xmm7, -1504(%rbp)
	movsd	%xmm4, -1480(%rbp)
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	movsd	.LC121(%rip), %xmm4
	movdqa	-4528(%rbp), %xmm7
	movq	-4368(%rbp), %rsi
	movq	%rax, %rdi
	movsd	%xmm4, -1360(%rbp)
	movsd	8+.LC121(%rip), %xmm4
	movaps	%xmm7, -1376(%rbp)
	movsd	%xmm4, -1352(%rbp)
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	movq	%r13, %rsi
	movq	%rax, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmIS4_EERS5_RKNS_9DenseBaseIT_EE.isra.0
	movq	%r14, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8finishedEv.isra.0
	movq	-4056(%rbp), %rax
	cmpq	-3512(%rbp), %rax
	jne	.L2063
	movq	-4392(%rbp), %rdi
	movdqa	-4336(%rbp), %xmm4
	movq	%rdi, -3296(%rbp)
	movaps	%xmm4, -3280(%rbp)
	cmpq	%rax, -3504(%rbp)
	jne	.L2064
	leaq	-3488(%rbp), %r13
	movq	%r12, %rdx
	movq	%r12, %rsi
	movq	%r13, %rdi
.LEHB157:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE157:
	movq	%r15, %rdx
	movq	%r13, %rsi
	movq	%r14, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC1IS4_EERS4_RKNS_9DenseBaseIT_EE
	movsd	.LC121(%rip), %xmm4
	movdqa	-4528(%rbp), %xmm7
	movq	%r14, %rdi
	movq	-4376(%rbp), %rsi
	movsd	%xmm4, -1488(%rbp)
	movsd	8+.LC121(%rip), %xmm4
	movaps	%xmm7, -1504(%rbp)
	movsd	%xmm4, -1480(%rbp)
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	movdqa	-4528(%rbp), %xmm7
	movq	-4368(%rbp), %rsi
	movq	%rax, %rdi
	movaps	%xmm7, -1376(%rbp)
	movsd	.LC121(%rip), %xmm7
	movsd	%xmm7, -1360(%rbp)
	movsd	8+.LC121(%rip), %xmm7
	movsd	%xmm7, -1352(%rbp)
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	movq	%r15, %rsi
	movq	%rax, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmIS4_EERS5_RKNS_9DenseBaseIT_EE.isra.0
	movq	%r14, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8finishedEv.isra.0
	movq	-4056(%rbp), %rax
	cmpq	-3480(%rbp), %rax
	jne	.L2063
	cmpq	%rax, -3472(%rbp)
	jne	.L2064
	movq	-4552(%rbp), %r12
	movl	-4536(%rbp), %esi
	leaq	-1184(%rbp), %r13
	movq	%r12, %rdi
	call	_ZNSt7__cxx119to_stringEi
	movl	-4628(%rbp), %esi
	movq	%r13, %rdi
	call	_ZNSt7__cxx119to_stringEi
	leaq	.LC126(%rip), %rdx
	xorl	%esi, %esi
	movq	%r13, %rdi
.LEHB158:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6insertEmPKc@PLT
.LEHE158:
	movq	-4672(%rbp), %r15
	movq	%rax, %rsi
	movq	%r15, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	-4648(%rbp), %rax
	movq	%r15, %rdi
	movq	8(%rax), %rdx
	movq	(%rax), %rsi
.LEHB159:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendEPKcm@PLT
.LEHE159:
	movq	-4720(%rbp), %r15
	movq	%rax, %rsi
	movq	%r15, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	%r12, %rdx
	movq	%r15, %rsi
	movq	%rbx, %rdi
.LEHB160:
	call	_ZStplIcSt11char_traitsIcESaIcEENSt7__cxx1112basic_stringIT_T0_T1_EEOS8_S9_
.LEHE160:
	movq	-4656(%rbp), %rax
	movq	%rbx, %rdi
	movq	8(%rax), %rdx
	movq	(%rax), %rsi
.LEHB161:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE6appendEPKcm@PLT
.LEHE161:
	movq	-4464(%rbp), %r15
	movq	%rax, %rsi
	leaq	-3456(%rbp), %r13
	movq	%r15, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC1EOS4_@PLT
	movq	%r15, %rsi
	movq	%r13, %rdi
.LEHB162:
	call	_Z19readMatrixTransformNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE162:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2078
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2078:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2079
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2079:
	movq	-1120(%rbp), %rdi
	movq	-4568(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2080
	movq	-1104(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2080:
	movq	-1152(%rbp), %rdi
	movq	-4560(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2081
	movq	-1136(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2081:
	movq	-1184(%rbp), %rdi
	leaq	-1168(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2082
	movq	-1168(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2082:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2083
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2083:
	movq	-4704(%rbp), %rsi
	leaq	-3424(%rbp), %r15
	movq	%r15, %rdi
	movq	%rsi, %rdx
.LEHB163:
	call	_ZN5Eigen6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEC1IiiEERKT_RKT0_
.LEHE163:
	movq	%r13, %rdx
	movq	%r15, %rsi
	movq	%r14, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEC1IS4_EERS4_RKNS_9DenseBaseIT_EE
	movsd	.LC121(%rip), %xmm4
	movdqa	-4528(%rbp), %xmm7
	movq	%r14, %rdi
	movq	-4376(%rbp), %rsi
	movsd	%xmm4, -1488(%rbp)
	movsd	8+.LC121(%rip), %xmm4
	movaps	%xmm7, -1504(%rbp)
	movsd	%xmm4, -1480(%rbp)
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	movdqa	-4528(%rbp), %xmm7
	movq	-4368(%rbp), %rsi
	movq	%rax, %rdi
	movaps	%xmm7, -1376(%rbp)
	movsd	.LC121(%rip), %xmm7
	movsd	%xmm7, -1360(%rbp)
	movsd	8+.LC121(%rip), %xmm7
	movsd	%xmm7, -1352(%rbp)
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmINS_14CwiseNullaryOpINS_8internal18scalar_constant_opIS3_EES4_EEEERS5_RKNS_9DenseBaseIT_EE
	movq	%r13, %rsi
	movq	%rax, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEcmIS4_EERS5_RKNS_9DenseBaseIT_EE.isra.0
	movq	%r14, %rdi
	call	_ZN5Eigen16CommaInitializerINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8finishedEv.isra.0
	movq	-4056(%rbp), %rax
	cmpq	-3416(%rbp), %rax
	jne	.L2063
	movq	-4392(%rbp), %rsi
	movq	%r15, %xmm0
	movhps	-4392(%rbp), %xmm0
	movq	%rsi, -3264(%rbp)
	movaps	%xmm0, -3248(%rbp)
	cmpq	%rax, -3408(%rbp)
	jne	.L2064
	movq	-4712(%rbp), %rax
	movq	-4624(%rbp), %rsi
	movq	-4616(%rbp), %rdi
	addq	%rax, %rsi
	movq	-4608(%rbp), %rax
	movq	%rsi, %r15
	movq	%rsi, %rdx
	movq	%rsi, -4312(%rbp)
	leaq	(%rax,%rax,2), %rax
	leaq	(%rax,%rdi), %r12
	movq	-4640(%rbp), %rdi
	movq	%r12, %rsi
	movq	%r12, -4320(%rbp)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	pxor	%xmm7, %xmm7
	movq	%r12, %rsi
	movq	%r15, %rdx
	movq	-4488(%rbp), %rdi
	movups	%xmm7, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	pxor	%xmm7, %xmm7
	movq	%r12, %rsi
	movq	%r15, %rdx
	movq	-4496(%rbp), %rdi
	movups	%xmm7, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	pxor	%xmm7, %xmm7
	movq	%r12, %rsi
	movq	%r15, %rdx
	movq	-4504(%rbp), %rdi
	movups	%xmm7, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	pxor	%xmm7, %xmm7
	movq	%r12, %rsi
	movq	%r15, %rdx
	movq	-4544(%rbp), %rdi
	movups	%xmm7, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	pxor	%xmm7, %xmm7
	movq	%r12, %rsi
	movq	%r15, %rdx
	movq	-4664(%rbp), %rdi
	movups	%xmm7, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	pxor	%xmm7, %xmm7
	movl	-4268(%rbp), %r12d
	movq	$0, -4336(%rbp)
	movups	%xmm7, (%rax)
	leaq	-3264(%rbp), %rax
	movq	%rax, -4696(%rbp)
	testl	%r12d, %r12d
	jle	.L2234
.L2084:
	movq	-4336(%rbp), %rsi
	movq	-4696(%rbp), %rdi
	movq	%rsi, %rdx
	movl	%esi, -4448(%rbp)
	movl	%esi, %r15d
.LEHB164:
	call	_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll
.LEHE164:
	movq	-4312(%rbp), %rdx
	movq	-4320(%rbp), %rsi
	movsd	%xmm0, -4304(%rbp)
	movq	-4640(%rbp), %rdi
	movsd	%xmm1, -4296(%rbp)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movl	-4268(%rbp), %r13d
	movupd	(%rax), %xmm0
	addpd	-4304(%rbp), %xmm0
	movups	%xmm0, (%rax)
	cmpl	%r13d, -4532(%rbp)
	jle	.L2383
	movl	-4264(%rbp), %eax
	imull	%r15d, %eax
	cltq
	testq	%rax, %rax
	js	.L2090
	movslq	%r13d, %rsi
	movq	-4200(%rbp), %r8
	movq	%rsi, -4352(%rbp)
	leaq	-3296(%rbp), %rsi
	movq	%rsi, -4592(%rbp)
	jmp	.L2089
.L2098:
	movq	%rax, %rsi
	orq	%rcx, %rsi
	js	.L2132
	imulq	%rcx, %rax
	testq	%rax, %rax
	je	.L2119
	salq	$4, %rax
	je	.L2119
	movq	%rax, %rdx
	xorl	%esi, %esi
	call	memset@PLT
.L2119:
	movq	.LC86(%rip), %xmm7
	leaq	-3328(%rbp), %rax
	movq	%r15, %rdi
	movq	-4392(%rbp), %rdx
	leaq	-4176(%rbp), %rcx
	movq	%rax, %rsi
	movaps	%xmm7, -4176(%rbp)
.LEHB165:
	call	_ZN5Eigen8internal20generic_product_implINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_NS_10DenseShapeESG_Li8EE13scaleAndAddToISA_EEvRT_RKSF_RSB_RKS6_
.LEHE165:
.L2114:
	movq	-4352(%rbp), %rax
	imulq	-3160(%rbp), %rax
	movq	-4336(%rbp), %rsi
	movq	-3152(%rbp), %rdi
	addq	%rsi, %rax
	salq	$4, %rax
	addq	-3168(%rbp), %rax
	movupd	(%rax), %xmm1
	movaps	%xmm1, -4416(%rbp)
	call	free@PLT
	movq	-4432(%rbp), %rax
	movapd	-4416(%rbp), %xmm1
	movupd	(%rax), %xmm2
	movapd	%xmm1, %xmm3
	shufpd	$1, %xmm1, %xmm3
	movapd	%xmm2, %xmm4
	movapd	%xmm2, %xmm0
	unpckhpd	%xmm2, %xmm4
	unpcklpd	%xmm2, %xmm0
	mulpd	%xmm3, %xmm4
	mulpd	%xmm1, %xmm0
	movapd	%xmm4, %xmm3
	addpd	%xmm0, %xmm3
	subpd	%xmm4, %xmm0
	movapd	%xmm3, %xmm4
	unpckhpd	%xmm3, %xmm3
	ucomisd	%xmm0, %xmm3
	movsd	%xmm0, %xmm4
	jp	.L2616
.L2120:
	movq	-4312(%rbp), %rdx
	movq	-4320(%rbp), %rsi
	movaps	%xmm4, -4416(%rbp)
	movq	-4496(%rbp), %rdi
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movapd	-4416(%rbp), %xmm4
	movl	-4448(%rbp), %r15d
	movupd	(%rax), %xmm0
	movl	-4268(%rbp), %ecx
	addpd	%xmm4, %xmm0
	movups	%xmm0, (%rax)
	movl	-4264(%rbp), %eax
	imull	%eax, %r15d
	leal	(%r15,%r13), %esi
	subl	%ecx, %esi
	movslq	%esi, %rdx
	testq	%rdx, %rdx
	js	.L2090
	cmpq	-4184(%rbp), %rdx
	jge	.L2090
	imull	%ecx, %eax
	salq	$4, %rdx
	addq	-4192(%rbp), %rdx
	addl	%esi, %eax
	cltq
	testq	%rax, %rax
	js	.L2090
	movq	-4200(%rbp), %r8
	cmpq	%r8, %rax
	jge	.L2090
	movupd	(%rdx), %xmm2
	salq	$4, %rax
	addq	-4208(%rbp), %rax
	movupd	(%rax), %xmm1
	movapd	%xmm2, %xmm4
	movapd	%xmm2, %xmm0
	unpckhpd	%xmm2, %xmm4
	movapd	%xmm1, %xmm3
	unpcklpd	%xmm2, %xmm0
	shufpd	$1, %xmm1, %xmm3
	mulpd	%xmm3, %xmm4
	mulpd	%xmm1, %xmm0
	movapd	%xmm4, %xmm3
	addpd	%xmm0, %xmm3
	subpd	%xmm4, %xmm0
	movapd	%xmm3, %xmm4
	unpckhpd	%xmm3, %xmm3
	ucomisd	%xmm0, %xmm3
	movsd	%xmm0, %xmm4
	jp	.L2617
.L2122:
	movq	-4312(%rbp), %rdx
	movq	-4320(%rbp), %rsi
	movq	%r8, -4432(%rbp)
	addl	$1, %r13d
	movq	-4504(%rbp), %rdi
	movl	%ecx, -4416(%rbp)
	movaps	%xmm4, -4480(%rbp)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movapd	-4480(%rbp), %xmm4
	cmpl	%r13d, -4532(%rbp)
	movupd	(%rax), %xmm0
	movl	-4416(%rbp), %ecx
	movq	-4432(%rbp), %r8
	addpd	%xmm4, %xmm0
	movups	%xmm0, (%rax)
	je	.L2087
	leal	(%r15,%r13), %eax
	addq	$1, -4352(%rbp)
	subl	%ecx, %eax
	cltq
	testq	%rax, %rax
	js	.L2090
.L2089:
	cmpq	%rax, %r8
	jle	.L2090
	salq	$4, %rax
	addq	-4208(%rbp), %rax
	movq	-4336(%rbp), %rdx
	movsd	8(%rax), %xmm0
	movsd	(%rax), %xmm7
	xorpd	.LC8(%rip), %xmm0
	movq	-4352(%rbp), %rsi
	movq	-4592(%rbp), %rdi
	movsd	%xmm7, -4432(%rbp)
	movsd	%xmm0, -4416(%rbp)
.LEHB166:
	call	_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll
.LEHE166:
	movsd	-4416(%rbp), %xmm2
	movsd	%xmm0, -4304(%rbp)
	movsd	-4432(%rbp), %xmm0
	movsd	%xmm1, -4296(%rbp)
	movapd	-4304(%rbp), %xmm1
	unpcklpd	%xmm2, %xmm2
	unpcklpd	%xmm0, %xmm0
	mulpd	%xmm1, %xmm2
	mulpd	%xmm1, %xmm0
	shufpd	$1, %xmm2, %xmm2
	movapd	%xmm2, %xmm3
	addpd	%xmm0, %xmm3
	subpd	%xmm2, %xmm0
	movapd	%xmm3, %xmm2
	unpckhpd	%xmm3, %xmm3
	ucomisd	%xmm0, %xmm3
	movsd	%xmm0, %xmm2
	jp	.L2618
.L2091:
	movq	-4312(%rbp), %rdx
	movq	-4320(%rbp), %rsi
	movaps	%xmm2, -4416(%rbp)
	movq	-4488(%rbp), %rdi
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movapd	-4416(%rbp), %xmm2
	movupd	(%rax), %xmm0
	addpd	%xmm2, %xmm0
	movups	%xmm0, (%rax)
	movl	-4448(%rbp), %eax
	imull	-4264(%rbp), %eax
	addl	%r13d, %eax
	subl	-4268(%rbp), %eax
	cltq
	testq	%rax, %rax
	js	.L2090
	cmpq	-4184(%rbp), %rax
	jge	.L2090
	salq	$4, %rax
	addq	-4192(%rbp), %rax
	cmpq	$0, -4352(%rbp)
	movq	%rax, -4432(%rbp)
	movq	-4048(%rbp), %rsi
	js	.L2094
	cmpq	%rsi, -4336(%rbp)
	jge	.L2094
	cmpq	%rsi, -4352(%rbp)
	jge	.L2094
	leaq	-3152(%rbp), %r15
	pxor	%xmm7, %xmm7
	movq	%rsi, %rdx
	movq	$0, -3168(%rbp)
	movq	%r15, %rdi
	movups	%xmm7, -3144(%rbp)
	movq	$-1, -3160(%rbp)
	movq	$0, -3152(%rbp)
.LEHB167:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
.LEHE167:
	movq	-3144(%rbp), %rax
	movq	-4056(%rbp), %rsi
	movq	-3136(%rbp), %rcx
	movq	-3152(%rbp), %rdi
	leaq	(%rax,%rsi), %rdx
	movq	%rax, -3160(%rbp)
	addq	%rcx, %rdx
	movq	%rdi, -3168(%rbp)
	cmpq	$19, %rdx
	jg	.L2098
	testq	%rsi, %rsi
	jle	.L2098
	movq	-3328(%rbp), %rax
	movdqa	-4688(%rbp), %xmm7
	movq	%rax, -3232(%rbp)
	movaps	%xmm7, -3216(%rbp)
	cmpq	-3632(%rbp), %rsi
	jne	.L2128
	movq	16(%rax), %rax
	movq	-4368(%rbp), %r12
	pxor	%xmm7, %xmm7
	movq	%rsi, %rdx
	movq	$0, -1376(%rbp)
	movq	%rax, %rsi
	movq	%r12, %rdi
	movups	%xmm7, -1368(%rbp)
.LEHB168:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
	movq	-4512(%rbp), %rdx
	leaq	-3232(%rbp), %rsi
	movq	%r12, %rdi
	call	_ZN5Eigen8internal20generic_product_implINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEES9_NS_10DenseShapeESE_Li8EE6evalToIS9_EEvRT_RKSD_RSA_
.LEHE168:
	movq	-4392(%rbp), %rax
	movq	-3144(%rbp), %r10
	movq	-4056(%rbp), %xmm0
	movq	%rax, -1352(%rbp)
	movq	-1376(%rbp), %rax
	movhps	-3632(%rbp), %xmm0
	movq	%rax, -1344(%rbp)
	movq	-1368(%rbp), %rax
	movups	%xmm0, -1320(%rbp)
	movq	%rax, -1336(%rbp)
	movq	-4064(%rbp), %rax
	movq	%rax, -1328(%rbp)
	movq	-4048(%rbp), %rax
	cmpq	%r10, %rax
	jne	.L2100
	movq	-3136(%rbp), %rcx
	cmpq	%rcx, %rax
	je	.L2104
.L2100:
	movq	%rax, %rdx
	movq	%rax, %rsi
	movq	%r15, %rdi
	movq	%rax, %r12
.LEHB169:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
.LEHE169:
	movq	-3144(%rbp), %r10
	cmpq	%r10, %r12
	jne	.L2106
	movq	-3136(%rbp), %rcx
	cmpq	%r10, %rcx
	jne	.L2106
.L2104:
	movq	-3152(%rbp), %r11
	testq	%rcx, %rcx
	jle	.L2107
	movl	%r13d, -4480(%rbp)
	movq	%r10, %rsi
	xorl	%r9d, %r9d
	xorl	%r15d, %r15d
	movq	%rcx, -4416(%rbp)
	movq	-4368(%rbp), %rdx
	salq	$4, %rsi
	movq	%r14, %rdi
	.p2align 4,,10
	.p2align 3
.L2108:
	movq	%r9, %r14
	movq	%rdx, %r12
	salq	$4, %r14
	cmpq	%r9, %r10
	jle	.L2113
	.p2align 4,,10
	.p2align 3
.L2111:
	movq	-1312(%rbp), %r13
	testq	%r13, %r13
	jle	.L2384
	movq	-1320(%rbp), %rax
	movq	-1336(%rbp), %r8
	pxor	%xmm2, %xmm2
	xorl	%edx, %edx
	movq	-1344(%rbp), %rcx
	imulq	%r15, %rax
	salq	$4, %r8
	addq	%r14, %rcx
	salq	$4, %rax
	addq	-1328(%rbp), %rax
	.p2align 4,,10
	.p2align 3
.L2110:
	movdqu	(%rcx), %xmm6
	movdqu	(%rax), %xmm5
	addq	$1, %rdx
	addq	$16, %rax
	movupd	-16(%rax), %xmm4
	addq	%r8, %rcx
	pshufd	$238, %xmm6, %xmm0
	pshufd	$78, %xmm5, %xmm1
	mulpd	%xmm0, %xmm1
	pshufd	$68, %xmm6, %xmm0
	mulpd	%xmm4, %xmm0
	xorpd	.LC8(%rip), %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm0, %xmm2
	cmpq	%rdx, %r13
	jne	.L2110
.L2109:
	addq	$1, %r9
	movaps	%xmm2, (%r11,%r14)
	addq	$16, %r14
	cmpq	%r9, %r10
	jne	.L2111
	movq	%r12, %rdx
.L2113:
	xorl	%r9d, %r9d
	testq	%r10, %r10
	cmovle	%r10, %r9
	addq	$1, %r15
	addq	%rsi, %r11
	cmpq	%r15, -4416(%rbp)
	jne	.L2108
	movq	%rdx, -4368(%rbp)
	movl	-4480(%rbp), %r13d
	movq	%rdi, %r14
.L2107:
	movq	-1376(%rbp), %rdi
	call	free@PLT
	jmp	.L2114
	.p2align 4,,10
	.p2align 3
.L2384:
	pxor	%xmm2, %xmm2
	jmp	.L2109
.L2383:
	movl	%r13d, %ecx
.L2087:
	testl	%ecx, %ecx
	jle	.L2124
	leaq	-3296(%rbp), %rax
	xorl	%r15d, %r15d
	movq	%rax, -4576(%rbp)
	jmp	.L2233
.L2127:
	movq	%rcx, %rax
	orq	%r8, %rax
	js	.L2132
	movq	%rcx, %rdx
	imulq	%r8, %rdx
	testq	%rdx, %rdx
	je	.L2177
	salq	$4, %rdx
	je	.L2177
	xorl	%esi, %esi
	movq	%r13, %rdi
	movq	%r8, -4448(%rbp)
	movq	%rcx, -4432(%rbp)
	call	memset@PLT
	movq	-4432(%rbp), %rcx
	movq	-4448(%rbp), %r8
.L2177:
	movq	-4048(%rbp), %rsi
	testq	%rsi, %rsi
	movq	%rsi, %r12
	sete	%al
	cmpq	%rsi, %rcx
	jne	.L2176
	cmpq	%rsi, %r8
	jne	.L2176
	movq	-3632(%rbp), %r9
	testq	%r9, %r9
	je	.L2388
	testb	%al, %al
	je	.L2619
.L2388:
	movq	%r13, %rdi
.L2172:
	movq	%rcx, %rax
	movq	-4336(%rbp), %rsi
	movsd	-4352(%rbp), %xmm0
	imulq	%r15, %rax
	addq	%rsi, %rax
	salq	$4, %rax
	movupd	0(%r13,%rax), %xmm1
	movq	-4416(%rbp), %rax
	movapd	%xmm1, %xmm3
	movq	%rax, %xmm7
	movq	%rax, %xmm2
	movaps	%xmm1, -4592(%rbp)
	unpckhpd	%xmm1, %xmm3
	unpcklpd	%xmm7, %xmm0
	movhpd	-4352(%rbp), %xmm2
	mulpd	%xmm0, %xmm3
	movapd	%xmm1, %xmm0
	unpcklpd	%xmm1, %xmm0
	mulpd	%xmm2, %xmm0
	movapd	%xmm3, %xmm2
	addpd	%xmm0, %xmm2
	subpd	%xmm3, %xmm0
	movapd	%xmm2, %xmm7
	movaps	%xmm2, -4480(%rbp)
	movsd	%xmm0, %xmm7
	movaps	%xmm0, -4448(%rbp)
	movaps	%xmm7, -4432(%rbp)
	call	free@PLT
	movapd	-4480(%rbp), %xmm2
	movapd	-4448(%rbp), %xmm0
	movapd	-4592(%rbp), %xmm1
	unpckhpd	%xmm2, %xmm2
	ucomisd	%xmm0, %xmm2
	jp	.L2620
	movq	-4312(%rbp), %rdx
	movq	-4320(%rbp), %rsi
	movq	-4544(%rbp), %rdi
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movupd	(%rax), %xmm0
	subpd	-4432(%rbp), %xmm0
.L2600:
	movl	-4268(%rbp), %ecx
	addq	$1, %r15
	movups	%xmm0, (%rax)
	cmpl	%r15d, %ecx
	jle	.L2124
.L2233:
	movq	-4336(%rbp), %r12
	movq	-4576(%rbp), %rdi
	movq	%r15, %rsi
	movq	%r12, %rdx
.LEHB170:
	call	_ZNK5Eigen15DenseCoeffsBaseINS_7ProductINS1_INS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEESA_Li0EEELi0EEclEll
.LEHE170:
	movq	-4048(%rbp), %rsi
	movsd	%xmm0, -4416(%rbp)
	movsd	%xmm1, -4352(%rbp)
	cmpq	%rsi, %r12
	jge	.L2094
	cmpq	%rsi, %r15
	jge	.L2094
	leaq	-3104(%rbp), %rax
	pxor	%xmm0, %xmm0
	movq	%rsi, %rdx
	movq	$0, -3120(%rbp)
	movq	%rax, %rdi
	movq	%rax, -4432(%rbp)
	movq	$-1, -3112(%rbp)
	movq	$0, -3104(%rbp)
	movups	%xmm0, -3096(%rbp)
.LEHB171:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
.LEHE171:
	movq	-3096(%rbp), %rcx
	movq	-4056(%rbp), %rdx
	movq	-3088(%rbp), %r8
	movq	-3104(%rbp), %r13
	leaq	(%rcx,%rdx), %rax
	movq	%rcx, -3112(%rbp)
	addq	%r8, %rax
	movq	%r13, -3120(%rbp)
	cmpq	$19, %rax
	jg	.L2127
	testq	%rdx, %rdx
	jle	.L2127
	movq	-3328(%rbp), %rax
	movdqa	-4688(%rbp), %xmm7
	movq	%rax, -3200(%rbp)
	movaps	%xmm7, -3184(%rbp)
	cmpq	-3632(%rbp), %rdx
	jne	.L2128
	movq	16(%rax), %rsi
	pxor	%xmm0, %xmm0
	movq	-4376(%rbp), %rdi
	movq	$0, -1504(%rbp)
	movups	%xmm0, -1496(%rbp)
.LEHB172:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
.LEHE172:
	movq	-3640(%rbp), %rax
	movq	-1496(%rbp), %rcx
	movq	-1488(%rbp), %r8
	leaq	(%rax,%rcx), %rdx
	addq	%r8, %rdx
	cmpq	$19, %rdx
	jg	.L2129
	testq	%rax, %rax
	jg	.L2621
.L2129:
	movq	%rcx, %rax
	orq	%r8, %rax
	js	.L2132
	movq	%rcx, %r9
	movq	-1504(%rbp), %r13
	imulq	%r8, %r9
	testq	%r9, %r9
	je	.L2136
	movq	%r9, %rdx
	salq	$4, %rdx
	je	.L2136
	xorl	%esi, %esi
	movq	%r13, %rdi
	movq	%r9, -4592(%rbp)
	movq	%r8, -4480(%rbp)
	movq	%rcx, -4448(%rbp)
	call	memset@PLT
	movq	-4448(%rbp), %rcx
	movq	-4480(%rbp), %r8
	movq	-4592(%rbp), %r9
.L2136:
	movq	.LC86(%rip), %xmm7
	movq	-4048(%rbp), %rax
	movaps	%xmm7, -4160(%rbp)
	movq	%rax, %r12
	cmpq	%rax, %rcx
	jne	.L2135
	movq	-3632(%rbp), %rdx
	cmpq	%rdx, %r8
	jne	.L2135
	movq	-4056(%rbp), %rax
	testq	%r12, %r12
	sete	%sil
	testq	%rax, %rax
	sete	%dil
	orb	%dil, %sil
	jne	.L2131
	testq	%r8, %r8
	jne	.L2622
.L2131:
	movq	%rcx, -1464(%rbp)
	movq	-4392(%rbp), %rsi
	movq	-4064(%rbp), %rcx
	movq	%r13, -1472(%rbp)
	movq	%rsi, -1480(%rbp)
	movq	%rcx, -1456(%rbp)
	movq	%rax, -1448(%rbp)
	movq	%rdx, -1440(%rbp)
	cmpq	-3096(%rbp), %r12
	jne	.L2161
	movq	-3088(%rbp), %r10
	cmpq	%r10, %r12
	je	.L2386
.L2161:
	movq	-4432(%rbp), %rdi
	movq	%r12, %rdx
	movq	%r12, %rsi
.LEHB173:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
.LEHE173:
	movq	-3096(%rbp), %r10
	cmpq	%r10, %r12
	jne	.L2106
	movq	-3088(%rbp), %rcx
	cmpq	%rcx, %r10
	jne	.L2106
.L2164:
	movq	-3104(%rbp), %r11
	testq	%rcx, %rcx
	jle	.L2165
	movq	%rcx, -4432(%rbp)
	movq	%r10, %rsi
	xorl	%r9d, %r9d
	xorl	%r13d, %r13d
	movq	%r15, -4448(%rbp)
	movq	-4368(%rbp), %rdi
	salq	$4, %rsi
.L2166:
	movq	%r9, %r15
	movq	%rdi, %r12
	salq	$4, %r15
	cmpq	%r9, %r10
	jle	.L2171
	.p2align 4,,10
	.p2align 3
.L2169:
	movq	-1440(%rbp), %r8
	testq	%r8, %r8
	jle	.L2387
	movq	-1448(%rbp), %rax
	movq	-1464(%rbp), %rdi
	pxor	%xmm2, %xmm2
	xorl	%edx, %edx
	movq	-1472(%rbp), %rcx
	imulq	%r13, %rax
	salq	$4, %rdi
	addq	%r15, %rcx
	salq	$4, %rax
	addq	-1456(%rbp), %rax
	.p2align 4,,10
	.p2align 3
.L2168:
	movdqu	(%rax), %xmm3
	movupd	(%rax), %xmm6
	addq	$1, %rdx
	addq	$16, %rax
	pshufd	$78, %xmm3, %xmm1
	movdqu	(%rcx), %xmm3
	addq	%rdi, %rcx
	pshufd	$238, %xmm3, %xmm0
	mulpd	%xmm0, %xmm1
	pshufd	$68, %xmm3, %xmm0
	mulpd	%xmm6, %xmm0
	xorpd	.LC8(%rip), %xmm1
	addpd	%xmm1, %xmm0
	addpd	%xmm0, %xmm2
	cmpq	%rdx, %r8
	jne	.L2168
.L2167:
	addq	$1, %r9
	movaps	%xmm2, (%r11,%r15)
	addq	$16, %r15
	cmpq	%r9, %r10
	jne	.L2169
	movq	%r12, %rdi
.L2171:
	xorl	%r9d, %r9d
	movq	-4432(%rbp), %rax
	testq	%r10, %r10
	cmovle	%r10, %r9
	addq	$1, %r13
	addq	%rsi, %r11
	cmpq	%rax, %r13
	jne	.L2166
	movq	%rdi, -4368(%rbp)
	movq	-4448(%rbp), %r15
.L2165:
	movq	-1504(%rbp), %rdi
	call	free@PLT
	movq	-3120(%rbp), %r13
	movq	-3112(%rbp), %rcx
	movq	-3104(%rbp), %rdi
	jmp	.L2172
.L2387:
	pxor	%xmm2, %xmm2
	jmp	.L2167
.L2619:
	cmpq	$1, %rsi
	je	.L2623
.L2466:
.L2468:
.L2467:
	movabsq	$9223372036854775807, %rax
	pxor	%xmm0, %xmm0
	movq	$0, -2336(%rbp)
	cqto
	movups	%xmm0, -2328(%rbp)
	idivq	%r9
	cmpq	%rax, %r8
	jg	.L2570
	leaq	-2336(%rbp), %r13
	movq	%r9, %rdx
	movq	%r8, %rsi
	movq	%r13, %rdi
.LEHB174:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
	movq	-4048(%rbp), %rsi
	movq	-3632(%rbp), %rdx
	cmpq	-2328(%rbp), %rsi
	jne	.L2206
	cmpq	-2320(%rbp), %rdx
	je	.L2207
.L2206:
	movq	%r13, %rdi
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE6resizeEll
.LEHE174:
.L2207:
	movdqu	-2328(%rbp), %xmm0
	movq	-3640(%rbp), %rax
	movq	-2328(%rbp), %r8
	movq	-2320(%rbp), %rsi
	movhlps	%xmm0, %xmm7
	movq	%xmm7, %rcx
	leaq	(%rax,%r8), %rdx
	addq	%rcx, %rdx
	cmpq	$19, %rdx
	jg	.L2208
	testq	%rax, %rax
	jg	.L2624
.L2208:
	orq	%r8, %rsi
	js	.L2132
	movq	%r8, %rax
	movq	-2336(%rbp), %rdi
	imulq	%rcx, %rax
	testq	%rax, %rax
	je	.L2213
	salq	$4, %rax
	movq	%rax, %rdx
	je	.L2213
	xorl	%esi, %esi
	movq	%r8, -4448(%rbp)
	movq	%rcx, -4432(%rbp)
	movaps	%xmm0, -4480(%rbp)
	call	memset@PLT
	movq	-4432(%rbp), %rcx
	movq	-4448(%rbp), %r8
	movdqa	-4480(%rbp), %xmm0
	movq	%rax, %rdi
.L2213:
	movq	.LC86(%rip), %xmm7
	movaps	%xmm7, -4128(%rbp)
	cmpq	-4048(%rbp), %r8
	jne	.L2135
	cmpq	%rcx, -3632(%rbp)
	jne	.L2135
	movq	-4056(%rbp), %rax
	testq	%rax, %rax
	sete	%dl
	testq	%r8, %r8
	sete	%sil
	orb	%sil, %dl
	jne	.L2605
	testq	%rcx, %rcx
	jne	.L2592
.L2605:
	movsd	.LC131(%rip), %xmm7
	leaq	-1352(%rbp), %rax
	leaq	-1344(%rbp), %r13
	movq	%rax, -4448(%rbp)
	leaq	-1360(%rbp), %rax
	movsd	%xmm7, -4480(%rbp)
	movsd	8+.LC131(%rip), %xmm7
	movq	%rax, %rsi
	movq	%rax, -4432(%rbp)
	movq	%xmm7, %r12
.L2210:
	movdqu	-3096(%rbp), %xmm4
	pxor	%xmm0, %xmm0
	movq	%r13, %rdi
	movq	-4448(%rbp), %rdx
	movq	%rcx, -1344(%rbp)
	movl	$1, %ecx
	movaps	%xmm0, -1376(%rbp)
	movaps	%xmm4, -1360(%rbp)
	call	_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_
	movq	-1344(%rbp), %rax
	movq	-1360(%rbp), %rdx
	subq	$8, %rsp
	movsd	-4480(%rbp), %xmm0
	movq	-4048(%rbp), %rdi
	movq	%r12, %xmm7
	movq	%r12, -1472(%rbp)
	imulq	%rax, %rdx
	movapd	%xmm7, %xmm1
	imulq	-1352(%rbp), %rax
	movsd	%xmm0, -1480(%rbp)
	movq	%rdi, %rsi
	movq	%rdx, -1336(%rbp)
	movq	%rax, -1328(%rbp)
	pushq	-4368(%rbp)
	pushq	-3096(%rbp)
	pushq	$1
	pushq	-3104(%rbp)
	movq	-4064(%rbp), %r9
	movq	-2336(%rbp), %rcx
	pushq	-4056(%rbp)
	movq	-2328(%rbp), %r8
	movq	-2320(%rbp), %rdx
.LEHB175:
	.cfi_escape 0x2e,0x30
	call	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi0ELb0ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0
.LEHE175:
	movq	-1376(%rbp), %rdi
	addq	$48, %rsp
	call	free@PLT
	movq	-1368(%rbp), %rdi
	call	free@PLT
	movq	-2336(%rbp), %rdi
	call	free@PLT
	movq	-3120(%rbp), %r13
	movq	-3112(%rbp), %rcx
	movq	-3104(%rbp), %rdi
	jmp	.L2172
.L2124:
	addq	$1, -4336(%rbp)
	movq	-4336(%rbp), %rax
	cmpl	%eax, %ecx
	jg	.L2084
.L2234:
	movq	-4320(%rbp), %r12
	movq	-4312(%rbp), %r15
	movq	-4640(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r12, %rsi
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	-4488(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r12, %rsi
	movupd	(%rax), %xmm0
	mulpd	.LC127(%rip), %xmm0
	movups	%xmm0, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	-4496(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r12, %rsi
	movupd	(%rax), %xmm0
	mulpd	.LC127(%rip), %xmm0
	movups	%xmm0, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	-4504(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r12, %rsi
	movupd	(%rax), %xmm0
	mulpd	.LC127(%rip), %xmm0
	movups	%xmm0, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	-4544(%rbp), %r13
	movq	%r15, %rdx
	movq	%r12, %rsi
	movupd	(%rax), %xmm0
	mulpd	.LC127(%rip), %xmm0
	movq	%r13, %rdi
	movups	%xmm0, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	%r15, %rdx
	movq	%r13, %rdi
	movq	%r12, %rsi
	movupd	(%rax), %xmm0
	mulpd	.LC127(%rip), %xmm0
	movups	%xmm0, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	-4504(%rbp), %rdi
	movq	%r15, %rdx
	movq	%r12, %rsi
	movq	%r15, -4312(%rbp)
	movq	%rax, %r13
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	-4312(%rbp), %rdx
	movq	-4496(%rbp), %rdi
	movq	%r12, %rsi
	movq	%rax, %r15
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	-4488(%rbp), %rdi
	movq	-4312(%rbp), %rdx
	movq	%r12, %rsi
	movq	%rax, -4352(%rbp)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	%r12, %rsi
	movq	%rax, -4336(%rbp)
	movq	-4312(%rbp), %rdx
	movq	-4640(%rbp), %rdi
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	-4336(%rbp), %r8
	movq	%r12, %rsi
	movq	-4352(%rbp), %rcx
	movupd	(%rax), %xmm1
	movq	-4312(%rbp), %rdx
	movupd	(%r8), %xmm0
	movq	-4664(%rbp), %rdi
	addpd	%xmm1, %xmm0
	movupd	(%rcx), %xmm1
	addpd	%xmm1, %xmm0
	movupd	(%r15), %xmm1
	addpd	%xmm1, %xmm0
	movupd	0(%r13), %xmm1
	addpd	%xmm1, %xmm0
	movaps	%xmm0, -4336(%rbp)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movapd	-4336(%rbp), %xmm0
	movq	-3424(%rbp), %rdi
	movups	%xmm0, (%rax)
	call	free@PLT
	movq	-3456(%rbp), %rdi
	call	free@PLT
	movq	-3488(%rbp), %rdi
	call	free@PLT
	movq	-3520(%rbp), %rdi
	call	free@PLT
	movq	-3552(%rbp), %rdi
	call	free@PLT
	movq	-3584(%rbp), %rdi
	call	free@PLT
	movq	-4192(%rbp), %rdi
	call	free@PLT
	addq	$1, -4624(%rbp)
	movq	-4624(%rbp), %rax
	addq	$32, -4656(%rbp)
	cmpq	$3, %rax
	jne	.L2085
	addl	$1, -4536(%rbp)
	movl	-4536(%rbp), %eax
	addq	$3, -4712(%rbp)
	cmpl	%eax, -4272(%rbp)
	jg	.L2065
.L2066:
	movq	-3616(%rbp), %rdi
	call	free@PLT
	movq	-3648(%rbp), %rdi
	call	free@PLT
	movq	-3680(%rbp), %rdi
	call	free@PLT
	movq	-3712(%rbp), %rdi
	call	free@PLT
	movq	-4208(%rbp), %rdi
	call	free@PLT
	addq	$1, -4616(%rbp)
	movq	-4616(%rbp), %rax
	addq	$32, -4648(%rbp)
	cmpq	$3, %rax
	jne	.L2235
	addq	$1, -4608(%rbp)
	movq	-4608(%rbp), %rax
	cmpl	%eax, -4272(%rbp)
	jg	.L2047
.L2048:
	call	_ZNSt6chrono3_V212system_clock3nowEv@PLT
	leaq	_ZSt4cout(%rip), %rdx
	leaq	.LC132(%rip), %rsi
	movq	%rax, -4416(%rbp)
	movq	_ZSt4cout(%rip), %rax
	leaq	_ZSt4cout(%rip), %rdi
	leaq	-960(%rbp), %rbx
	addq	-24(%rax), %rdx
	movl	24(%rdx), %eax
	movq	$10, 8(%rdx)
	andl	$-261, %eax
	orl	$4, %eax
	movl	%eax, 24(%rdx)
	movl	$27, %edx
.LEHB176:
	.cfi_escape 0x2e,0
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	-4368(%rbp), %r14
	movq	-4640(%rbp), %rax
	leaq	_ZSt4cout(%rip), %rdi
	movq	%r14, %rsi
	movq	%rax, -1376(%rbp)
	call	_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE
	movq	%rax, %rdi
	leaq	.LC133(%rip), %rsi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
	movl	$27, %edx
	leaq	.LC134(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	-4488(%rbp), %rax
	movq	%r14, %rsi
	leaq	_ZSt4cout(%rip), %rdi
	movq	%rax, -1376(%rbp)
	call	_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE
	movq	%rax, %rdi
	leaq	.LC133(%rip), %rsi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
	movl	$27, %edx
	leaq	.LC135(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	-4496(%rbp), %rax
	movq	%r14, %rsi
	leaq	_ZSt4cout(%rip), %rdi
	movq	%rax, -1376(%rbp)
	call	_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE
	movq	%rax, %rdi
	leaq	.LC133(%rip), %rsi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
	movl	$27, %edx
	leaq	.LC136(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	-4504(%rbp), %rax
	movq	%r14, %rsi
	leaq	_ZSt4cout(%rip), %rdi
	movq	%rax, -1376(%rbp)
	call	_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE
	movq	%rax, %rdi
	leaq	.LC133(%rip), %rsi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
	movl	$27, %edx
	leaq	.LC137(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	-4544(%rbp), %rax
	movq	%r14, %rsi
	leaq	_ZSt4cout(%rip), %rdi
	movq	%rax, -1376(%rbp)
	call	_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE
	movq	%rax, %rdi
	leaq	.LC133(%rip), %rsi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
	movl	$26, %edx
	leaq	.LC138(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	-4664(%rbp), %rax
	movq	%r14, %rsi
	leaq	_ZSt4cout(%rip), %rdi
	movq	%rax, -1376(%rbp)
	call	_ZN5EigenlsINS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEERSoSA_RKNS_9DenseBaseIT_EE
	movq	%rax, %rdi
	leaq	.LC133(%rip), %rsi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
	leaq	.LC139(%rip), %rdi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	movslq	-4272(%rbp), %rax
	movq	$0x000000000, -1360(%rbp)
	movq	%rax, %xmm0
	punpcklqdq	%xmm0, %xmm0
	movaps	%xmm0, -1376(%rbp)
	testq	%rax, %rax
	js	.L2625
	leaq	-3424(%rbp), %r15
	movq	-4368(%rbp), %rsi
	leaq	-960(%rbp), %rbx
	movq	%r15, %rdi
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseNullaryOpINS_8internal18scalar_constant_opIdEES2_EEEERKNS_9DenseBaseIT_EE
.LEHE176:
	xorl	%r13d, %r13d
	xorl	%r12d, %r12d
.L2237:
	cmpl	%r12d, -4272(%rbp)
	jle	.L2626
	leaq	2(%r13), %rax
	xorl	%ebx, %ebx
	movq	%rax, -4336(%rbp)
.L2242:
	movq	-3736(%rbp), %rdx
	leaq	(%rbx,%rbx,2), %rcx
	movq	-4336(%rbp), %rsi
	movq	%rcx, %rax
	imulq	%rdx, %rax
	addq	%r13, %rax
	salq	$4, %rax
	addq	-3744(%rbp), %rax
	cmpq	%rsi, %rdx
	jle	.L2238
	leaq	2(%rcx), %rsi
	cmpq	-3728(%rbp), %rsi
	jge	.L2238
	movq	%rax, -1376(%rbp)
	leaq	-1616(%rbp), %r14
	movq	-4664(%rbp), %rax
	movq	-4368(%rbp), %rsi
	movq	%r14, %rdi
	movq	%r13, -1352(%rbp)
	movq	%rax, -1360(%rbp)
	movq	%rcx, -1344(%rbp)
	movq	%rdx, -1336(%rbp)
.LEHB177:
	call	_ZN5Eigen15PlainObjectBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEC2INS_14CwiseUnaryViewINS_8internal18scalar_imag_ref_opISt7complexIdEEENS_5BlockINS1_IS9_Lin1ELin1ELi0ELin1ELin1EEELi3ELi3ELb0EEEEEEERKNS_9DenseBaseIT_EE
.LEHE177:
	movq	%r14, -1496(%rbp)
	movq	-1600(%rbp), %rax
	cmpq	%rax, -1608(%rbp)
	jne	.L2627
	testq	%rax, %rax
	js	.L2628
	movq	%r14, -1344(%rbp)
	movl	$1, %edx
	movl	$2, %esi
	movq	.LC145(%rip), %rax
	movq	%r14, -1336(%rbp)
	movq	-4368(%rbp), %r14
	movq	%rax, -1352(%rbp)
	movq	%r14, %rdi
	call	_ZNK5Eigen15DenseCoeffsBaseINS_13CwiseBinaryOpINS_8internal17scalar_product_opIddEEKNS_14CwiseNullaryOpINS2_18scalar_constant_opIdEEKNS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEEKNS1_INS2_20scalar_difference_opIddEESA_KNS_9TransposeISA_EEEEEELi0EEclEll
	movl	$2, %edx
	xorl	%esi, %esi
	movq	%r14, %rdi
	mulsd	-4256(%rbp), %xmm0
	movsd	%xmm0, -4312(%rbp)
	call	_ZNK5Eigen15DenseCoeffsBaseINS_13CwiseBinaryOpINS_8internal17scalar_product_opIddEEKNS_14CwiseNullaryOpINS2_18scalar_constant_opIdEEKNS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEEKNS1_INS2_20scalar_difference_opIddEESA_KNS_9TransposeISA_EEEEEELi0EEclEll
	xorl	%edx, %edx
	movl	$1, %esi
	movq	%r14, %rdi
	mulsd	-4248(%rbp), %xmm0
	movsd	-4312(%rbp), %xmm4
	addsd	%xmm0, %xmm4
	movsd	%xmm4, -4312(%rbp)
	call	_ZNK5Eigen15DenseCoeffsBaseINS_13CwiseBinaryOpINS_8internal17scalar_product_opIddEEKNS_14CwiseNullaryOpINS2_18scalar_constant_opIdEEKNS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEEKNS1_INS2_20scalar_difference_opIddEESA_KNS_9TransposeISA_EEEEEELi0EEclEll
	movq	%rbx, %rdx
	movq	%r12, %rsi
	movq	%r15, %rdi
	mulsd	-4240(%rbp), %xmm0
	addsd	-4312(%rbp), %xmm0
	movsd	%xmm0, -4312(%rbp)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movq	%rbx, %rdx
	movq	%r12, %rsi
	movq	%r15, %rdi
	movsd	-4312(%rbp), %xmm0
	addq	$1, %rbx
	movsd	%xmm0, (%rax)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movsd	-4256(%rbp), %xmm0
	movsd	-4248(%rbp), %xmm2
	movsd	-4240(%rbp), %xmm1
	movq	-1616(%rbp), %rdi
	mulsd	%xmm2, %xmm2
	mulsd	%xmm0, %xmm0
	mulsd	%xmm1, %xmm1
	addsd	%xmm2, %xmm0
	addsd	%xmm1, %xmm0
	movsd	(%rax), %xmm1
	divsd	%xmm0, %xmm1
	movsd	%xmm1, (%rax)
	call	free@PLT
	cmpl	%ebx, -4272(%rbp)
	jg	.L2242
	movq	%r12, %rax
	addq	$3, %r13
	addq	$1, %rax
	movq	%rax, %r12
	jmp	.L2237
.L2621:
	movq	-4392(%rbp), %rsi
	movq	%rsi, -3392(%rbp)
	movq	-4512(%rbp), %rsi
	movq	%rsi, -3376(%rbp)
	cmpq	-4056(%rbp), %rax
	jne	.L2209
	movq	-4376(%rbp), %rdi
	leaq	-4274(%rbp), %rdx
	leaq	-3392(%rbp), %rsi
.LEHB178:
	call	_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_
.LEHE178:
.L2602:
	movq	-1496(%rbp), %rcx
	movq	-1504(%rbp), %r13
	movq	-4056(%rbp), %rax
	movq	-3632(%rbp), %rdx
	movq	-4048(%rbp), %r12
	jmp	.L2131
.L2622:
	cmpq	$1, %r8
	je	.L2629
	cmpq	$1, %rcx
	je	.L2630
	pxor	%xmm0, %xmm0
	movq	%r8, %xmm7
	movq	%rax, -1344(%rbp)
	leaq	-1344(%rbp), %r13
	leaq	-1360(%rbp), %rax
	movaps	%xmm0, -1376(%rbp)
	movq	%rcx, %xmm0
	movq	%r13, %rdi
	punpcklqdq	%xmm7, %xmm0
	movq	%rax, %rsi
	movl	$1, %ecx
	leaq	-1352(%rbp), %rdx
	movaps	%xmm0, -1360(%rbp)
	call	_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_
	movq	-1344(%rbp), %rax
	movq	-1360(%rbp), %rdx
	subq	$8, %rsp
	movsd	.LC131(%rip), %xmm0
	movsd	8+.LC131(%rip), %xmm1
	imulq	%rax, %rdx
	imulq	-1352(%rbp), %rax
	movsd	%xmm0, -1592(%rbp)
	movsd	%xmm1, -1584(%rbp)
	movq	%rdx, -1336(%rbp)
	movq	-4056(%rbp), %rdx
	movq	%rax, -1328(%rbp)
	pushq	-4368(%rbp)
	movq	%rdx, %r8
	pushq	-1496(%rbp)
	pushq	$1
	pushq	-1504(%rbp)
	movq	-3648(%rbp), %r9
	movq	-4064(%rbp), %rcx
	pushq	-3640(%rbp)
	movq	-3632(%rbp), %rsi
	movq	-4048(%rbp), %rdi
.LEHB179:
	.cfi_escape 0x2e,0x30
	call	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0
.LEHE179:
	movq	-1376(%rbp), %rdi
	addq	$48, %rsp
	call	free@PLT
	movq	-1368(%rbp), %rdi
	call	free@PLT
	jmp	.L2602
.L2386:
	movq	%r10, %rcx
	jmp	.L2164
.L2624:
	movq	-4392(%rbp), %rsi
	movq	%rsi, -3360(%rbp)
	movq	-4512(%rbp), %rsi
	movq	%rsi, -3344(%rbp)
	cmpq	-4056(%rbp), %rax
	jne	.L2209
	leaq	-4273(%rbp), %rdx
	leaq	-3360(%rbp), %rsi
	movq	%r13, %rdi
.LEHB180:
	.cfi_escape 0x2e,0
	call	_ZN5Eigen8internal42call_restricted_packet_assignment_no_aliasINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEENS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS4_EEKNS_9TransposeIKS5_EEEES5_Li1EEENS0_9assign_opIS4_S4_EEEEvRT_RKT0_RKT1_
.LEHE180:
	movq	-2320(%rbp), %rcx
	jmp	.L2605
.L2592:
	cmpq	$1, %rcx
	je	.L2631
	subq	$1, %r8
	je	.L2632
	movq	%rax, -1344(%rbp)
	pxor	%xmm1, %xmm1
	leaq	-1360(%rbp), %rax
	leaq	-1344(%rbp), %r13
	leaq	-1352(%rbp), %rdx
	movq	%rax, %rsi
	movl	$1, %ecx
	movq	%r13, %rdi
	movq	%rdx, -4448(%rbp)
	movq	%rax, -4432(%rbp)
	movaps	%xmm1, -1376(%rbp)
	movaps	%xmm0, -1360(%rbp)
	call	_ZN5Eigen8internal37evaluateProductBlockingSizesHeuristicISt7complexIdES3_Li1ElEEvRT2_S5_S5_S4_
	movq	-1344(%rbp), %rax
	movq	-1360(%rbp), %rdx
	movsd	.LC131(%rip), %xmm0
	movsd	8+.LC131(%rip), %xmm1
	imulq	%rax, %rdx
	imulq	-1352(%rbp), %rax
	movsd	%xmm0, -4480(%rbp)
	movq	%xmm1, %r12
	movsd	%xmm0, -1480(%rbp)
	movq	%rdx, -1336(%rbp)
	movq	-4056(%rbp), %rdx
	movq	%rax, -1328(%rbp)
	movsd	%xmm1, -1472(%rbp)
	pushq	%r8
	movq	%rdx, %r8
	pushq	-4368(%rbp)
	pushq	-2328(%rbp)
	pushq	$1
	pushq	-2336(%rbp)
	movq	-3648(%rbp), %r9
	movq	-4064(%rbp), %rcx
	pushq	-3640(%rbp)
	movq	-3632(%rbp), %rsi
	movq	-4048(%rbp), %rdi
.LEHB181:
	.cfi_escape 0x2e,0x30
	call	_ZN5Eigen8internal29general_matrix_matrix_productIlSt7complexIdELi1ELb1ES3_Li0ELb0ELi0ELi1EE3runElllPKS3_lS6_lPS3_llS3_RNS0_15level3_blockingIS3_S3_EEPNS0_16GemmParallelInfoIlEE.isra.0
.LEHE181:
	movq	-1376(%rbp), %rdi
	addq	$48, %rsp
	call	free@PLT
	movq	-1368(%rbp), %rdi
	call	free@PLT
	movq	-2320(%rbp), %rcx
	movq	-4432(%rbp), %rsi
	jmp	.L2210
.L2623:
	leaq	-2528(%rbp), %rax
	movq	-4392(%rbp), %rsi
	xorl	%edx, %edx
	movq	%r9, -4784(%rbp)
	movq	%rax, %rdi
	movq	%rax, -4592(%rbp)
	call	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS5_l
	movq	-2528(%rbp), %rax
	movq	-3328(%rbp), %rdx
	movq	-2504(%rbp), %rdi
	movq	-2520(%rbp), %rsi
	testq	%rax, %rax
	movq	%rdx, -2400(%rbp)
	movq	-4784(%rbp), %r9
	movq	%rdx, -2008(%rbp)
	movq	8(%rdi), %rdx
	movq	%rax, -4432(%rbp)
	movq	-4592(%rbp), %rax
	movq	%rsi, -4480(%rbp)
	movq	%rdi, -4448(%rbp)
	je	.L2494
	testq	%rsi, %rsi
	js	.L2182
.L2494:
	movq	%rax, %rsi
	movq	-4480(%rbp), %rax
	leaq	-1816(%rbp), %rdi
	movl	$14, %ecx
	rep movsl
	cmpq	%rax, %r9
	jne	.L2633
	movq	-2008(%rbp), %rax
	pxor	%xmm0, %xmm0
	xorl	%r10d, %r10d
	movq	%r9, -1312(%rbp)
	movq	%r9, -1832(%rbp)
	leaq	-1296(%rbp), %rdi
	leaq	-1840(%rbp), %rsi
	movl	$26, %ecx
	movq	%rax, -1360(%rbp)
	movq	%rax, -1496(%rbp)
	movq	-4512(%rbp), %rax
	movq	%r10, -1760(%rbp)
	movq	%rax, -1344(%rbp)
	movq	-4432(%rbp), %rax
	movq	%rdx, -1744(%rbp)
	movq	%rax, -1840(%rbp)
	movups	%xmm0, -1336(%rbp)
	rep movsl
	testq	%r9, %r9
	jle	.L2634
	leaq	-1360(%rbp), %rax
	leaq	-1600(%rbp), %rdi
	movq	%rax, %rsi
.LEHB182:
	.cfi_escape 0x2e,0
	call	_ZN5Eigen8internal17product_evaluatorINS_7ProductINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS6_Lin1ELin1ELi0ELin1ELin1EEEEEEESA_Li0EEELi8ENS_10DenseShapeESG_S6_S6_EC2ERKSF_
.LEHE182:
	movq	-4432(%rbp), %rax
	pxor	%xmm0, %xmm0
	xorl	%edx, %edx
	xorl	%esi, %esi
	movups	%xmm0, -1560(%rbp)
	movq	%r14, %rdi
	movq	%rax, -1536(%rbp)
	movq	-4448(%rbp), %rax
	movq	8(%rax), %rax
	movq	%rax, -1520(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll
	movsd	%xmm0, -4448(%rbp)
	movsd	%xmm1, -4432(%rbp)
	jmp	.L2186
.L2632:
	movq	-3328(%rbp), %rdx
	movq	%rdi, -1376(%rbp)
	pxor	%xmm0, %xmm0
	xorl	%r9d, %r9d
	movq	-4512(%rbp), %rdi
	movq	%rcx, -1360(%rbp)
	leaq	-2832(%rbp), %rsi
	leaq	-4128(%rbp), %rcx
	movq	%rdx, -2928(%rbp)
	movq	%rdx, -2880(%rbp)
	movq	%rdx, -2832(%rbp)
	movq	-4368(%rbp), %rdx
	movq	%rax, -2840(%rbp)
	movq	%r13, -1352(%rbp)
	movq	$1, -1328(%rbp)
	movq	%r9, -2816(%rbp)
	movq	%r9, -2808(%rbp)
	movq	%rax, -2792(%rbp)
	movaps	%xmm0, -2864(%rbp)
	movaps	%xmm0, -1344(%rbp)
.LEHB183:
	call	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0
.L2603:
	leaq	-1352(%rbp), %rax
	movsd	.LC131(%rip), %xmm4
	movsd	8+.LC131(%rip), %xmm7
	leaq	-1344(%rbp), %r13
	movq	%rax, -4448(%rbp)
	leaq	-1360(%rbp), %rax
	movq	-2320(%rbp), %rcx
	movq	%rax, -4432(%rbp)
	movq	%xmm7, %r12
	movq	%rax, %rsi
	movsd	%xmm4, -4480(%rbp)
	jmp	.L2210
.L2090:
	leaq	.LC129(%rip), %rcx
	movl	$427, %edx
	leaq	.LC4(%rip), %rsi
	leaq	.LC130(%rip), %rdi
	call	__assert_fail@PLT
.L2094:
	leaq	.LC93(%rip), %rcx
	movl	$118, %edx
	leaq	.LC4(%rip), %rsi
	leaq	.LC5(%rip), %rdi
	call	__assert_fail@PLT
.L2631:
	leaq	-2208(%rbp), %rax
	xorl	%edx, %edx
	movq	%r13, %rsi
	movq	%r8, -4448(%rbp)
	movq	%rax, %rdi
	movq	%rax, -4432(%rbp)
	leaq	-2144(%rbp), %r13
	call	_ZN5Eigen5BlockINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS4_l
	movq	-4512(%rbp), %rsi
	xorl	%edx, %edx
	movq	%r13, %rdi
	call	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS5_l
	movq	-4448(%rbp), %r8
	movq	-4432(%rbp), %rax
	subq	$1, %r8
	je	.L2635
	leaq	-2080(%rbp), %rdi
	movl	$14, %ecx
	movq	%r13, %rsi
	movq	-3328(%rbp), %rdx
	rep movsl
	movsd	.LC131(%rip), %xmm0
	movsd	8+.LC131(%rip), %xmm1
	leaq	-2080(%rbp), %rsi
	movq	%rdx, -4112(%rbp)
	leaq	-4112(%rbp), %rdi
	movq	%rax, %rdx
	call	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0
.LEHE183:
	jmp	.L2603
.L2630:
	movq	-3200(%rbp), %rdx
	movq	-4376(%rbp), %rsi
	pxor	%xmm0, %xmm0
	xorl	%r11d, %r11d
	movq	-4512(%rbp), %rdi
	leaq	-4160(%rbp), %rcx
	movq	%rax, -2984(%rbp)
	movq	%rdx, -3072(%rbp)
	movq	%rdx, -3024(%rbp)
	movq	%rsi, -2568(%rbp)
	leaq	-2976(%rbp), %rsi
	movq	%rdx, -2976(%rbp)
	leaq	-2592(%rbp), %rdx
	movq	%r13, -2592(%rbp)
	movq	%r9, -2576(%rbp)
	movq	$1, -2544(%rbp)
	movq	%r11, -2960(%rbp)
	movq	%r11, -2952(%rbp)
	movq	%rax, -2936(%rbp)
	movaps	%xmm0, -3008(%rbp)
	movaps	%xmm0, -2560(%rbp)
.LEHB184:
	call	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_9TransposeIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEENS4_IKNS_5BlockIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS7_EEKSA_EELi1ELin1ELb1EEEEENS4_INSB_IS8_Li1ELin1ELb0EEEEEEEvRKT_RKT0_RT1_RKNST_6ScalarE.isra.0
	jmp	.L2602
.L2187:
	movq	%r12, %rsi
	xorl	%edx, %edx
	movq	%r14, %rdi
	addq	$1, %r12
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS_7ProductINS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEESE_Li0EEELi1ELin1ELb0EEEEEEEKNSB_IKNSB_ISF_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESW_S5_S5_E5coeffEll
	addsd	-4448(%rbp), %xmm0
	addsd	-4432(%rbp), %xmm1
	movsd	%xmm0, -4448(%rbp)
	movsd	%xmm1, -4432(%rbp)
.L2186:
	cmpq	%r12, -4480(%rbp)
	jne	.L2187
	movq	-1584(%rbp), %rdi
	call	free@PLT
	movsd	-4432(%rbp), %xmm1
	pxor	%xmm3, %xmm3
	movsd	-4448(%rbp), %xmm0
	movsd	.LC24(%rip), %xmm2
	call	__muldc3@PLT
	movupd	0(%r13), %xmm7
	movq	-3112(%rbp), %rcx
	unpcklpd	%xmm1, %xmm0
	movq	-3104(%rbp), %rdi
	addpd	%xmm7, %xmm0
	movups	%xmm0, 0(%r13)
	movq	-3120(%rbp), %r13
	jmp	.L2172
.L2629:
	movq	-4376(%rbp), %rax
	movq	-4512(%rbp), %rsi
	pxor	%xmm0, %xmm0
	xorl	%edx, %edx
	movq	%r13, -2784(%rbp)
	movq	%rax, -2760(%rbp)
	leaq	-2720(%rbp), %rax
	movq	%rax, %rdi
	movq	%rax, -4448(%rbp)
	movq	%r9, -2776(%rbp)
	movq	%r9, -2736(%rbp)
	movaps	%xmm0, -2752(%rbp)
	call	_ZN5Eigen5BlockIKNS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELin1ELi1ELb1EEC1ERS5_l
	movq	%r12, %rax
	subq	$1, %rax
	movq	-4448(%rbp), %rax
	je	.L2636
	movq	-3200(%rbp), %rdx
	movl	$14, %ecx
	movq	%rax, %rsi
	movsd	.LC131(%rip), %xmm0
	movsd	8+.LC131(%rip), %xmm1
	leaq	-2656(%rbp), %rdi
	rep movsl
	movq	%rdx, -4144(%rbp)
	leaq	-2656(%rbp), %rsi
	leaq	-2784(%rbp), %rdx
	leaq	-4144(%rbp), %rdi
	call	_ZN5Eigen8internal19gemv_dense_selectorILi2ELi1ELb1EE3runINS_12CwiseUnaryOpINS0_19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS7_Lin1ELin1ELi0ELin1ELin1EEEEEEENS_5BlockISC_Lin1ELi1ELb1EEENSG_ISB_Lin1ELi1ELb1EEEEEvRKT_RKT0_RT1_RKNSP_6ScalarE.isra.0
.LEHE184:
	jmp	.L2602
.L2182:
	leaq	.LC75(%rip), %rcx
	movl	$176, %edx
	leaq	.LC17(%rip), %rsi
	leaq	.LC18(%rip), %rdi
	call	__assert_fail@PLT
.L2636:
	xorl	%edx, %edx
	leaq	-3200(%rbp), %rsi
	leaq	-2272(%rbp), %rdi
	call	_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC1ERSE_l
	movq	-2272(%rbp), %rax
	movq	-2696(%rbp), %rcx
	movq	-2720(%rbp), %r8
	movq	-2232(%rbp), %rdx
	movq	%rax, -1928(%rbp)
	movq	-2256(%rbp), %rax
	movq	8(%rcx), %r9
	movq	%rax, -1912(%rbp)
	movq	-2248(%rbp), %rax
	movq	%rax, -1904(%rbp)
	movq	-2712(%rbp), %rax
	testq	%r8, %r8
	je	.L2142
	testq	%rax, %rax
	js	.L2182
.L2142:
	cmpq	%rdx, %rax
	jne	.L2218
	movq	-1928(%rbp), %rdi
	movq	-1912(%rbp), %rsi
	movq	-1904(%rbp), %rcx
	movq	%rdi, -1832(%rbp)
	movq	%rsi, -1816(%rbp)
	movq	%rcx, -1808(%rbp)
	testq	%rax, %rax
	jne	.L2637
	pxor	%xmm1, %xmm1
	movapd	%xmm1, %xmm0
.L2144:
	movsd	.LC24(%rip), %xmm2
	pxor	%xmm3, %xmm3
	call	__muldc3@PLT
	movupd	0(%r13), %xmm7
	unpcklpd	%xmm1, %xmm0
	addpd	%xmm7, %xmm0
	movups	%xmm0, 0(%r13)
	jmp	.L2602
.L2209:
	leaq	.LC83(%rip), %rcx
	movl	$96, %edx
	leaq	.LC84(%rip), %rsi
	leaq	.LC85(%rip), %rdi
	call	__assert_fail@PLT
.L2637:
	jle	.L2220
	movq	8(%rdi), %rdx
	movq	%rsi, %xmm7
	movq	(%rdi), %rdi
	movq	%rax, -4448(%rbp)
	movq	-4368(%rbp), %r12
	movq	%r8, -1312(%rbp)
	movq	%rdx, %xmm0
	imulq	%rsi, %rdx
	xorl	%esi, %esi
	movq	%rdi, -1352(%rbp)
	punpcklqdq	%xmm7, %xmm0
	movq	%r12, %rdi
	movq	%r9, -1296(%rbp)
	movaps	%xmm0, -1344(%rbp)
	movq	%rcx, %xmm0
	addq	%rcx, %rdx
	movq	%rdx, %xmm7
	punpcklqdq	%xmm7, %xmm0
	movaps	%xmm0, -1328(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movq	-4448(%rbp), %rax
	movapd	%xmm0, %xmm2
	cmpq	$1, %rax
	je	.L2146
	movq	%rax, %rsi
	movq	%r12, %rdi
	movq	%rax, -4592(%rbp)
	andq	$-2, %rsi
	movaps	%xmm0, -4368(%rbp)
	movq	%rsi, -4448(%rbp)
	movl	$1, %esi
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movl	$2, %edx
	movq	%r14, -4480(%rbp)
	movq	-4592(%rbp), %r14
	movapd	%xmm0, %xmm1
	movq	%rbx, -4592(%rbp)
	movq	%rdx, %rbx
	jmp	.L2147
.L2148:
	movq	%rbx, %rsi
	movq	%r12, %rdi
	movaps	%xmm1, -4784(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	addpd	-4368(%rbp), %xmm0
	leaq	1(%rbx), %rsi
	movq	%r12, %rdi
	addq	$2, %rbx
	movaps	%xmm0, -4368(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	-4784(%rbp), %xmm1
	addpd	%xmm0, %xmm1
.L2147:
	cmpq	%rbx, -4448(%rbp)
	jg	.L2148
	movapd	-4368(%rbp), %xmm2
	movq	%r14, %rax
	movq	%r12, %rdi
	movq	-4448(%rbp), %rsi
	movq	%r12, -4368(%rbp)
	movq	-4480(%rbp), %r14
	addpd	%xmm1, %xmm2
	movq	-4592(%rbp), %rbx
	cmpq	%rsi, %rax
	jle	.L2146
	movaps	%xmm2, -4448(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	-4448(%rbp), %xmm2
	addpd	%xmm0, %xmm2
.L2146:
	movapd	%xmm2, %xmm7
	movapd	%xmm2, %xmm0
	unpckhpd	%xmm7, %xmm7
	movapd	%xmm7, %xmm1
	jmp	.L2144
.L2218:
	leaq	.LC76(%rip), %rcx
	movl	$82, %edx
	leaq	.LC77(%rip), %rsi
	leaq	.LC78(%rip), %rdi
	call	__assert_fail@PLT
.L2220:
	leaq	.LC79(%rip), %rcx
	movl	$411, %edx
	leaq	.LC67(%rip), %rsi
	leaq	.LC68(%rip), %rdi
	call	__assert_fail@PLT
.L2620:
	movsd	-4352(%rbp), %xmm3
	movsd	-4416(%rbp), %xmm2
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm1
	call	__muldc3@PLT
	movq	-4312(%rbp), %rdx
	movq	-4320(%rbp), %rsi
	movapd	%xmm0, %xmm7
	movq	-4544(%rbp), %rdi
	unpcklpd	%xmm1, %xmm7
	movaps	%xmm7, -4352(%rbp)
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movupd	(%rax), %xmm0
	subpd	-4352(%rbp), %xmm0
	jmp	.L2600
.L2176:
	leaq	.LC87(%rip), %rcx
	movl	$470, %edx
	leaq	.LC73(%rip), %rsi
	leaq	.LC74(%rip), %rdi
	call	__assert_fail@PLT
.L2135:
	leaq	.LC72(%rip), %rcx
	movl	$470, %edx
	leaq	.LC73(%rip), %rsi
	leaq	.LC74(%rip), %rdi
	call	__assert_fail@PLT
.L2634:
	leaq	.LC89(%rip), %rcx
	movl	$411, %edx
	leaq	.LC67(%rip), %rsi
	leaq	.LC68(%rip), %rdi
	call	__assert_fail@PLT
.L2633:
	leaq	.LC88(%rip), %rcx
	movl	$82, %edx
	leaq	.LC77(%rip), %rsi
	leaq	.LC78(%rip), %rdi
	call	__assert_fail@PLT
.L2635:
	leaq	-3328(%rbp), %rax
	xorl	%edx, %edx
	leaq	-2080(%rbp), %rdi
	movq	%rax, %rsi
	call	_ZN5Eigen5BlockIKNS_12CwiseUnaryOpINS_8internal19scalar_conjugate_opISt7complexIdEEEKNS_9TransposeIKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEC1ERSE_l
	movq	-2080(%rbp), %rax
	movq	-2120(%rbp), %rcx
	movq	-2144(%rbp), %r8
	movq	-2040(%rbp), %rdx
	movq	%rax, -1720(%rbp)
	movq	-2064(%rbp), %rax
	movq	-2136(%rbp), %r13
	movq	8(%rcx), %r9
	movq	%rax, -1704(%rbp)
	movq	-2056(%rbp), %rax
	movq	%rax, -1696(%rbp)
	testq	%r8, %r8
	je	.L2495
	testq	%r13, %r13
	js	.L2182
.L2495:
	cmpq	%rdx, %r13
	jne	.L2218
	movq	-1720(%rbp), %rdi
	movq	-1704(%rbp), %rsi
	movq	-1696(%rbp), %rcx
	movq	%rdi, -1608(%rbp)
	movq	%rsi, -1592(%rbp)
	movq	%rcx, -1584(%rbp)
	testq	%r13, %r13
	jne	.L2638
	pxor	%xmm1, %xmm1
	movapd	%xmm1, %xmm0
.L2219:
	movsd	.LC24(%rip), %xmm2
	pxor	%xmm3, %xmm3
	call	__muldc3@PLT
	movq	-2208(%rbp), %rax
	unpcklpd	%xmm1, %xmm0
	movupd	(%rax), %xmm7
	addpd	%xmm7, %xmm0
	movups	%xmm0, (%rax)
	jmp	.L2603
.L2618:
	movsd	-4432(%rbp), %xmm2
	movsd	-4416(%rbp), %xmm3
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm1
	call	__muldc3@PLT
	movapd	%xmm0, %xmm2
	unpcklpd	%xmm1, %xmm2
	jmp	.L2091
.L2106:
	leaq	.LC95(%rip), %rcx
	movl	$765, %edx
	leaq	.LC96(%rip), %rsi
	leaq	.LC97(%rip), %rdi
	call	__assert_fail@PLT
.L2638:
	jle	.L2220
	movq	8(%rdi), %rdx
	movq	%rsi, %xmm7
	movq	(%rdi), %rdi
	movq	%r8, -1312(%rbp)
	movq	-4368(%rbp), %r12
	movq	%r9, -1296(%rbp)
	movq	%rdx, %xmm0
	imulq	%rsi, %rdx
	movq	%rdi, -1352(%rbp)
	xorl	%esi, %esi
	punpcklqdq	%xmm7, %xmm0
	movq	%r12, %rdi
	movaps	%xmm0, -1344(%rbp)
	movq	%rcx, %xmm0
	addq	%rcx, %rdx
	movq	%rdx, %xmm7
	punpcklqdq	%xmm7, %xmm0
	movaps	%xmm0, -1328(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	%xmm0, %xmm2
	cmpq	$1, %r13
	je	.L2221
	movq	%r13, %rax
	movl	$1, %esi
	movq	%r12, %rdi
	movaps	%xmm0, -4368(%rbp)
	andq	$-2, %rax
	movq	%rax, -4432(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movl	$2, %edx
	movq	%rbx, -4448(%rbp)
	movapd	%xmm0, %xmm1
	movq	%rdx, %rbx
	jmp	.L2222
.L2223:
	movq	%rbx, %rsi
	movq	%r12, %rdi
	movaps	%xmm1, -4480(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	addpd	-4368(%rbp), %xmm0
	leaq	1(%rbx), %rsi
	movq	%r12, %rdi
	addq	$2, %rbx
	movaps	%xmm0, -4368(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	-4480(%rbp), %xmm1
	addpd	%xmm0, %xmm1
.L2222:
	cmpq	%rbx, -4432(%rbp)
	jg	.L2223
	movapd	-4368(%rbp), %xmm2
	movq	-4432(%rbp), %rsi
	movq	%r12, %rdi
	movq	%r12, -4368(%rbp)
	movq	-4448(%rbp), %rbx
	addpd	%xmm1, %xmm2
	cmpq	%rsi, %r13
	jle	.L2221
	movaps	%xmm2, -4432(%rbp)
	call	_ZNK5Eigen8internal16binary_evaluatorINS_13CwiseBinaryOpINS0_22scalar_conj_product_opISt7complexIdES5_EEKNS_9TransposeIKNS_12CwiseUnaryOpINS0_19scalar_conjugate_opIS5_EEKNS_5BlockIKNS8_ISA_KNS7_IKNS_6MatrixIS5_Lin1ELin1ELi0ELin1ELin1EEEEEEELi1ELin1ELb1EEEEEEEKNSB_IKNSB_ISE_Lin1ELi1ELb1EEELin1ELi1ELb1EEEEENS0_10IndexBasedESU_S5_S5_E6packetILi0ENS0_9Packet1cdEEET0_l
	movapd	-4432(%rbp), %xmm2
	addpd	%xmm0, %xmm2
.L2221:
	movapd	%xmm2, %xmm7
	movapd	%xmm2, %xmm0
	unpckhpd	%xmm7, %xmm7
	movapd	%xmm7, %xmm1
	jmp	.L2219
.L2128:
	leaq	.LC94(%rip), %rcx
	movl	$96, %edx
	leaq	.LC84(%rip), %rsi
	leaq	.LC85(%rip), %rdi
	call	__assert_fail@PLT
.L2617:
	movapd	%xmm2, %xmm7
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm1
	movq	%r8, -4432(%rbp)
	unpckhpd	%xmm7, %xmm7
	movl	%ecx, -4416(%rbp)
	movapd	%xmm7, %xmm3
	call	__muldc3@PLT
	movq	-4432(%rbp), %r8
	movl	-4416(%rbp), %ecx
	movapd	%xmm0, %xmm4
	unpcklpd	%xmm1, %xmm4
	jmp	.L2122
.L2064:
	leaq	.LC123(%rip), %rcx
	movl	$96, %edx
	leaq	.LC84(%rip), %rsi
	leaq	.LC85(%rip), %rdi
	call	__assert_fail@PLT
.L2063:
	leaq	.LC122(%rip), %rcx
	movl	$96, %edx
	leaq	.LC84(%rip), %rsi
	leaq	.LC85(%rip), %rdi
	call	__assert_fail@PLT
.L2615:
	call	_ZN5Eigen9DenseBaseINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEE8ConstantEllRKS3_.part.0
.L2616:
	movapd	%xmm2, %xmm7
	movapd	%xmm1, %xmm0
	unpckhpd	%xmm1, %xmm1
	unpckhpd	%xmm7, %xmm7
	movapd	%xmm7, %xmm3
	call	__muldc3@PLT
	movapd	%xmm0, %xmm4
	unpcklpd	%xmm1, %xmm4
	jmp	.L2120
.L2132:
	leaq	.LC0(%rip), %rcx
	movl	$71, %edx
	leaq	.LC1(%rip), %rsi
	leaq	.LC2(%rip), %rdi
	call	__assert_fail@PLT
.L2628:
	leaq	.LC144(%rip), %rcx
	movl	$71, %edx
	leaq	.LC1(%rip), %rsi
	leaq	.LC2(%rip), %rdi
	call	__assert_fail@PLT
.L2627:
	leaq	.LC143(%rip), %rcx
	movl	$116, %edx
	leaq	.LC64(%rip), %rsi
	leaq	.LC65(%rip), %rdi
	call	__assert_fail@PLT
.L2238:
	leaq	.LC141(%rip), %rcx
	movl	$132, %edx
	leaq	.LC20(%rip), %rsi
	leaq	.LC142(%rip), %rdi
	call	__assert_fail@PLT
.L2626:
	movl	$22, %edx
	leaq	.LC146(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
.LEHB185:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	%r15, %rsi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZN5EigenlsINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEEEERSoS3_RKNS_9DenseBaseIT_EE
	movq	%rax, %rdi
	leaq	.LC133(%rip), %rsi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
	leaq	-576(%rbp), %r14
	movl	$8, %edx
	leaq	.LC147(%rip), %rsi
	movq	%r14, %rdi
	movq	%r14, -4464(%rbp)
	call	_ZNSt14basic_ifstreamIcSt11char_traitsIcEEC1EPKcSt13_Ios_Openmode@PLT
.LEHE185:
	leaq	-1136(%rbp), %rax
	xorl	%edi, %edi
	xorl	%edx, %edx
	xorl	%ecx, %ecx
	movq	%rax, -4560(%rbp)
	leaq	-1152(%rbp), %rsi
	pxor	%xmm0, %xmm0
	leaq	-960(%rbp), %rbx
	movq	%rax, -1152(%rbp)
	leaq	-1104(%rbp), %rax
	movq	%rax, -4568(%rbp)
	movq	%rax, -1120(%rbp)
	leaq	-3328(%rbp), %rax
	movq	%rdi, -3312(%rbp)
	movq	%r14, %rdi
	movq	%rsi, -4672(%rbp)
	movq	%rdx, -1144(%rbp)
	movb	$0, -1136(%rbp)
	movq	%rcx, -1112(%rbp)
	movb	$0, -1104(%rbp)
	movq	%rax, -4336(%rbp)
	movaps	%xmm0, -3328(%rbp)
.LEHB186:
	call	_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE.isra.0
	leaq	-3328(%rbp), %rax
	xorl	%r13d, %r13d
	leaq	-960(%rbp), %rbx
	movq	%rax, -4336(%rbp)
	leaq	-1120(%rbp), %r14
	jmp	.L2244
.L2245:
	movq	-4672(%rbp), %r12
	movq	-4464(%rbp), %rdi
	movq	%r12, %rsi
	call	_ZSt7getlineIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE.isra.0
	movl	$8, %edx
	movq	%r12, %rsi
	movq	%rbx, %rdi
	call	_ZNSt7__cxx1119basic_istringstreamIcSt11char_traitsIcESaIcEEC1ERKNS_12basic_stringIcS2_S3_EESt13_Ios_Openmode@PLT
.LEHE186:
	movq	%r14, %rsi
	movq	%rbx, %rdi
.LEHB187:
	call	_ZStrsIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@PLT
	movq	%r14, %rsi
	movq	%rbx, %rdi
	call	_ZStrsIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@PLT
	movq	%r14, %rsi
	movq	%rbx, %rdi
	call	_ZStrsIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@PLT
	movq	%r14, %rsi
	movq	%rbx, %rdi
	call	_ZStrsIcSt11char_traitsIcESaIcEERSt13basic_istreamIT_T0_ES7_RNSt7__cxx1112basic_stringIS4_S5_T1_EE@PLT
	movq	-4336(%rbp), %rdi
	movq	%r14, %rsi
	call	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE9push_backERKS5_
.LEHE187:
	movq	%rbx, %rdi
	addl	$1, %r13d
	call	_ZNSt7__cxx1119basic_istringstreamIcSt11char_traitsIcESaIcEED1Ev@PLT
.L2244:
	cmpl	%r13d, -4272(%rbp)
	jg	.L2245
	movq	-4464(%rbp), %rdi
.LEHB188:
	call	_ZNSt14basic_ifstreamIcSt11char_traitsIcEE5closeEv@PLT
	leaq	.LC148(%rip), %rdi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	leaq	.LC149(%rip), %rdi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
.LEHE188:
	xorl	%eax, %eax
	pxor	%xmm0, %xmm0
	xorl	%r13d, %r13d
	movq	$0x000000000, -4320(%rbp)
	movq	%rax, -1360(%rbp)
	movq	$0x000000000, -4352(%rbp)
	movaps	%xmm0, -1376(%rbp)
.L2246:
	movslq	-4272(%rbp), %rax
	cmpl	%r13d, %eax
	jle	.L2639
	movq	$0x000000000, -4312(%rbp)
	movq	%rax, %r12
	xorl	%r14d, %r14d
.L2247:
	movq	%r14, %rdx
	movq	%r13, %rsi
	movq	%r15, %rdi
	addq	$1, %r14
	call	_ZN5Eigen15DenseCoeffsBaseINS_6MatrixIdLin1ELin1ELi0ELin1ELin1EEELi1EEclEll
	movsd	-4312(%rbp), %xmm7
	addsd	(%rax), %xmm7
	movsd	%xmm7, -4312(%rbp)
	movsd	%xmm7, -1504(%rbp)
	cmpq	%r14, %r12
	jne	.L2247
	movq	-1368(%rbp), %rsi
	cmpq	-1360(%rbp), %rsi
	je	.L2248
	movsd	%xmm7, (%rsi)
	addq	$8, %rsi
	movq	%rsi, -1368(%rbp)
.L2249:
	movq	%r13, %rax
	movq	%rbx, %rdi
	salq	$5, %rax
	addq	-3328(%rbp), %rax
	movq	%rax, %r14
	leaq	-944(%rbp), %rax
	movq	%rax, -960(%rbp)
	movq	8(%r14), %rdx
	movq	(%r14), %rsi
	movq	%rax, -4384(%rbp)
	addq	%rsi, %rdx
.LEHB189:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
.LEHE189:
	movq	%rbx, %rdi
.LEHB190:
	call	_Z8chargeOfNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
	pxor	%xmm0, %xmm0
	movsd	-4312(%rbp), %xmm7
	leaq	-1072(%rbp), %r12
	movq	8(%r14), %rdx
	cvtsi2sdl	%eax, %xmm0
	movq	%r12, -1088(%rbp)
	movq	(%r14), %rsi
	leaq	-1088(%rbp), %rax
	movq	%rax, %rdi
	movq	%r12, -4400(%rbp)
	addq	%rsi, %rdx
	movq	%rax, -4552(%rbp)
	addsd	%xmm0, %xmm7
	movsd	%xmm7, -4392(%rbp)
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
.LEHE190:
	movq	-4552(%rbp), %rdi
.LEHB191:
	call	_Z8chargeOfNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
	pxor	%xmm1, %xmm1
	movsd	-4392(%rbp), %xmm2
	movq	(%r14), %rdx
	leal	1(%r13), %esi
	cvtsi2sdl	%eax, %xmm1
	movsd	-4312(%rbp), %xmm0
	leaq	.LC150(%rip), %rdi
	movl	$3, %eax
	call	printf@PLT
.LEHE191:
	movq	-1088(%rbp), %rdi
	cmpq	%r12, %rdi
	je	.L2250
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2250:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2251
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2251:
	movq	-4384(%rbp), %r12
	movsd	-4352(%rbp), %xmm7
	movq	%rbx, %rdi
	addsd	-4312(%rbp), %xmm7
	movq	8(%r14), %rdx
	movq	%r12, -960(%rbp)
	movq	(%r14), %rsi
	addq	%rsi, %rdx
	movsd	%xmm7, -4352(%rbp)
.LEHB192:
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_M_constructIPcEEvT_S7_St20forward_iterator_tag.isra.0
.LEHE192:
	movq	%rbx, %rdi
.LEHB193:
	call	_Z8chargeOfNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE193:
	pxor	%xmm0, %xmm0
	movq	-960(%rbp), %rdi
	cvtsi2sdl	%eax, %xmm0
	addsd	-4320(%rbp), %xmm0
	movsd	%xmm0, -4320(%rbp)
	cmpq	%r12, %rdi
	je	.L2252
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2252:
	addq	$1, %r13
	jmp	.L2246
.L2625:
	leaq	.LC140(%rip), %rcx
	movl	$71, %edx
	leaq	.LC1(%rip), %rsi
	leaq	.LC2(%rip), %rdi
	call	__assert_fail@PLT
.L2639:
	leaq	.LC151(%rip), %rdi
.LEHB194:
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	movsd	-4352(%rbp), %xmm0
	movsd	-4320(%rbp), %xmm7
	movl	$3, %eax
	leaq	.LC152(%rip), %rdi
	movapd	%xmm0, %xmm2
	movapd	%xmm7, %xmm1
	addsd	%xmm7, %xmm2
	call	printf@PLT
	movq	-4728(%rbp), %rsi
	movq	-4736(%rbp), %rax
	movl	$1000000, %ecx
	leaq	.LC153(%rip), %rdi
	subq	%rsi, %rax
	movq	-4744(%rbp), %rsi
	cqto
	idivq	%rcx
	movq	%rax, %r14
	movq	-4752(%rbp), %rax
	subq	%rsi, %rax
	movq	-4760(%rbp), %rsi
	cqto
	idivq	%rcx
	movq	%rax, -4312(%rbp)
	movq	-4416(%rbp), %rax
	subq	%rsi, %rax
	cqto
	idivq	%rcx
	movq	%rax, %r13
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.constprop.0.isra.0
	pxor	%xmm0, %xmm0
	leaq	.LC155(%rip), %rdi
	movl	$1, %eax
	cvtsi2sdq	%r14, %xmm0
	mulsd	.LC154(%rip), %xmm0
	call	printf@PLT
	movq	-4312(%rbp), %r15
	pxor	%xmm0, %xmm0
	movl	$1, %eax
	leaq	.LC156(%rip), %rdi
	cvtsi2sdq	%r15, %xmm0
	mulsd	.LC154(%rip), %xmm0
	call	printf@PLT
	pxor	%xmm0, %xmm0
	leaq	.LC157(%rip), %rdi
	movl	$1, %eax
	cvtsi2sdq	%r13, %xmm0
	mulsd	.LC154(%rip), %xmm0
	call	printf@PLT
	leaq	.LC158(%rip), %rdi
	call	puts@PLT
	movq	%r15, %rax
	pxor	%xmm0, %xmm0
	leaq	.LC159(%rip), %rdi
	addq	%r14, %rax
	addq	%r13, %rax
	cvtsi2sdq	%rax, %xmm0
	mulsd	.LC154(%rip), %xmm0
	movl	$1, %eax
	call	printf@PLT
	movl	$10, %edx
	leaq	.LC160(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0
	leaq	.LC161(%rip), %rsi
	movq	%rbx, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC2IS3_EEPKcRKS3_.constprop.3
.LEHE194:
	leaq	-3296(%rbp), %rax
	movq	%rbx, %rsi
	movq	%rax, %rdi
	movq	%rax, -4576(%rbp)
.LEHB195:
	call	_Z8readHermNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE195:
	movq	-960(%rbp), %rdi
	leaq	-944(%rbp), %rax
	movq	%rax, -4384(%rbp)
	cmpq	%rax, %rdi
	je	.L2254
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2254:
	movq	_ZSt4cout(%rip), %rax
	leaq	_ZSt4cout(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	movq	-24(%rax), %rdx
	movq	$6, 8(%rsi,%rdx)
	movq	-24(%rax), %rdx
	addq	%rsi, %rdx
	leaq	.LC162(%rip), %rsi
	movl	24(%rdx), %eax
	andl	$-261, %eax
	orl	$4, %eax
	movl	%eax, 24(%rdx)
	movl	$9, %edx
.LEHB196:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	-4576(%rbp), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE
	movq	%rax, %rdi
	leaq	.LC104(%rip), %rsi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
	leaq	.LC163(%rip), %rsi
	movq	%rbx, %rdi
	call	_ZNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEC2IS3_EEPKcRKS3_.constprop.3
.LEHE196:
	leaq	-3264(%rbp), %rax
	movq	%rbx, %rsi
	movq	%rax, %rdi
	movq	%rax, -4696(%rbp)
.LEHB197:
	call	_Z10readFDEBUGNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE@PLT
.LEHE197:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2255
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2255:
	movq	_ZSt4cout(%rip), %rax
	leaq	_ZSt4cout(%rip), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	movq	-24(%rax), %rdx
	movq	$6, 8(%rsi,%rdx)
	movq	-24(%rax), %rdx
	addq	%rsi, %rdx
	leaq	.LC164(%rip), %rsi
	movl	24(%rdx), %eax
	andl	$-261, %eax
	orl	$4, %eax
	movl	%eax, 24(%rdx)
	movl	$6, %edx
.LEHB198:
	call	_ZSt16__ostream_insertIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_PKS3_l@PLT
	movq	-4696(%rbp), %rsi
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZN5EigenlsINS_6MatrixISt7complexIdELin1ELin1ELi0ELin1ELin1EEEEERSoS5_RKNS_9DenseBaseIT_EE
	movq	%rax, %rdi
	leaq	.LC104(%rip), %rsi
	call	_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc.isra.0
	leaq	_ZSt4cout(%rip), %rdi
	call	_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_.isra.0
.LEHE198:
	movq	-3264(%rbp), %rdi
	call	free@PLT
	movq	-3296(%rbp), %rdi
	call	free@PLT
	movq	-1376(%rbp), %rdi
	testq	%rdi, %rdi
	je	.L2256
	movq	-1360(%rbp), %rsi
	subq	%rdi, %rsi
	call	_ZdlPvm@PLT
.L2256:
	movq	-4336(%rbp), %rdi
	call	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED1Ev
	movq	-1120(%rbp), %rdi
	movq	-4568(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2257
	movq	-1104(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2257:
	movq	-1152(%rbp), %rdi
	movq	-4560(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2258
	movq	-1136(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2258:
	movq	-4464(%rbp), %rdi
	call	_ZNSt14basic_ifstreamIcSt11char_traitsIcEED1Ev@PLT
	movq	-3424(%rbp), %rdi
	call	free@PLT
	movq	-3744(%rbp), %rdi
	call	free@PLT
	movq	-3776(%rbp), %rdi
	call	free@PLT
	movq	-3808(%rbp), %rdi
	call	free@PLT
	movq	-3840(%rbp), %rdi
	call	free@PLT
	movq	-3872(%rbp), %rdi
	call	free@PLT
	movq	-3904(%rbp), %rdi
	call	free@PLT
	movq	-3936(%rbp), %rdi
	call	free@PLT
	movq	-3968(%rbp), %rdi
	call	free@PLT
	movq	-4224(%rbp), %rdi
	call	free@PLT
	movq	-4000(%rbp), %rdi
	call	free@PLT
	movq	-4032(%rbp), %rdi
	call	free@PLT
	movq	-4064(%rbp), %rdi
	call	free@PLT
	movq	-4096(%rbp), %rdi
	testq	%rdi, %rdi
	je	.L2261
	movq	-4080(%rbp), %rsi
	subq	%rdi, %rsi
	call	_ZdlPvm@PLT
.L2261:
	subq	$32, %rbx
	movq	(%rbx), %rdi
	leaq	16(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L2260
	movq	16(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2260:
	movq	-4600(%rbp), %rax
	cmpq	%rax, %rbx
	jne	.L2261
	movq	-56(%rbp), %rax
	subq	%fs:40, %rax
	jne	.L2640
	leaq	-40(%rbp), %rsp
	xorl	%eax, %eax
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
.L2248:
	.cfi_restore_state
	movq	-4376(%rbp), %rdx
	movq	-4368(%rbp), %rdi
.LEHB199:
	call	_ZNSt6vectorIdSaIdEE17_M_realloc_insertIJRKdEEEvN9__gnu_cxx17__normal_iteratorIPdS1_EEDpOT_
.LEHE199:
	jmp	.L2249
.L2640:
	call	__stack_chk_fail@PLT
.L2446:
	movq	%rax, %r14
	jmp	.L2348
.L2450:
	movq	%rax, %r15
	jmp	.L2354
.L2454:
	movq	%rax, %r15
	jmp	.L2358
.L2455:
	movq	%rax, %r14
	jmp	.L2361
.L2453:
	movq	%rax, %r15
	jmp	.L2360
.L2449:
	movq	%rax, %r14
	jmp	.L2349
.L2438:
	movq	%rax, %r15
	jmp	.L2336
.L2476:
	jmp	.L2038
.L2434:
	movq	%rax, %r15
	jmp	.L2323
.L2460:
.L2612:
	movq	%rax, %r15
	jmp	.L2163
.L2440:
	movq	%rax, %r15
	jmp	.L2332
.L2483:
	jmp	.L2038
.L2430:
	movq	%rax, %r15
	jmp	.L2318
.L2429:
	movq	%rax, %r14
	jmp	.L2320
.L2452:
	movq	%rax, %r14
	jmp	.L2356
.L2451:
	movq	%rax, %r14
	jmp	.L2353
.L2448:
	movq	%rax, %r15
	jmp	.L2354
.L2445:
	movq	%rax, %r14
	leaq	-960(%rbp), %rbx
	jmp	.L2365
.L2470:
	movq	%rax, %r14
	jmp	.L2227
.L2447:
	movq	%rax, %r15
	jmp	.L2347
.L2490:
	jmp	.L2038
.L2484:
	jmp	.L2038
.L2481:
	jmp	.L2038
.L2480:
	jmp	.L2038
.L2491:
	jmp	.L2038
.L2474:
	jmp	.L2038
.L2436:
	movq	%rax, %r15
	jmp	.L2343
.L2435:
	movq	%rax, %r14
	jmp	.L2344
.L2437:
	movq	%rax, %r14
	jmp	.L2338
.L2492:
	jmp	.L2038
.L2478:
	jmp	.L2038
.L2439:
	movq	%rax, %r14
	jmp	.L2334
.L2432:
	movq	%rax, %r15
	jmp	.L2327
.L2431:
	movq	%rax, %r14
	jmp	.L2316
.L2428:
	movq	%rax, %r15
	jmp	.L2306
.L2433:
	movq	%rax, %r14
	jmp	.L2325
.L2427:
	movq	%rax, %r14
	jmp	.L2308
.L2426:
	movq	%rax, %r15
	jmp	.L2310
.L2422:
	movq	%rax, %r14
	jmp	.L2346
.L2414:
	movq	%rax, %r15
	jmp	.L2284
.L2399:
	movq	%rax, %r14
	jmp	.L2280
.L2400:
	movq	%rax, %r15
	jmp	.L2277
.L2487:
	jmp	.L2038
.L2408:
	movq	%rax, %r14
	leaq	-960(%rbp), %rbx
	jmp	.L2369
.L2407:
	movq	%rax, %r15
	leaq	-960(%rbp), %rbx
	jmp	.L2370
.L2406:
	movq	%rax, %r14
	leaq	-960(%rbp), %rbx
	jmp	.L2371
.L2405:
	movq	%rax, %r15
	leaq	-960(%rbp), %rbx
	jmp	.L2372
.L2404:
	movq	%rax, %r15
	jmp	.L2270
.L2403:
	movq	%rax, %r14
	jmp	.L2271
.L2402:
	movq	%rax, %r15
	jmp	.L2273
.L2401:
	movq	%rax, %r14
	jmp	.L2275
.L2412:
	movq	%rax, %r15
	jmp	.L2288
.L2411:
	movq	%rax, %r15
	leaq	-960(%rbp), %rbx
	jmp	.L2366
.L2410:
	movq	%rax, %r14
	leaq	-960(%rbp), %rbx
	jmp	.L2367
.L2409:
	movq	%rax, %r15
	leaq	-960(%rbp), %rbx
	jmp	.L2368
.L2444:
	movq	%rax, %r15
	jmp	.L2291
.L2413:
	movq	%rax, %r14
	jmp	.L2286
.L2424:
	movq	%rax, %r15
	jmp	.L2315
.L2421:
	movq	%rax, %r14
	jmp	.L2299
.L2420:
	movq	%rax, %r15
	jmp	.L2301
.L2419:
	movq	%rax, %r14
	jmp	.L2303
.L2457:
	movq	%rax, %r14
	jmp	.L2103
.L2456:
	movq	%rax, %r14
	jmp	.L2103
.L2423:
	movq	%rax, %r15
	jmp	.L2345
.L2425:
	movq	%rax, %r14
	jmp	.L2312
.L2459:
	movq	%rax, %r15
	jmp	.L2102
.L2442:
	movq	%rax, %r14
	jmp	.L2342
.L2458:
	movq	%rax, %r14
	jmp	.L2115
.L2418:
	movq	%rax, %r15
	jmp	.L2292
.L2417:
	movq	%rax, %r14
	jmp	.L2294
.L2416:
	movq	%rax, %r15
	jmp	.L2296
.L2415:
	movq	%rax, %r14
	jmp	.L2282
.L2392:
	movq	%rax, %r15
	leaq	-960(%rbp), %rbx
	jmp	.L2374
.L2471:
	jmp	.L2038
.L2475:
	jmp	.L2038
.L2479:
	jmp	.L2038
.L2482:
	jmp	.L2038
.L2485:
	jmp	.L2038
.L2486:
	jmp	.L2038
.L2488:
	jmp	.L2038
.L2396:
	movq	%rax, %r15
	jmp	.L2268
.L2395:
	movq	%rax, %r15
	jmp	.L2269
.L2394:
	movq	%rax, %r14
	jmp	.L2281
.L2393:
	movq	%rax, %r14
	leaq	-960(%rbp), %rbx
	jmp	.L2373
.L2398:
	movq	%rax, %r14
	jmp	.L2266
.L2397:
	movq	%rax, %r14
	jmp	.L2267
.L2464:
	movq	%rax, %r14
	jmp	.L2159
.L2441:
	movq	%rax, %r14
	jmp	.L2330
.L2465:
	movq	%rax, %r15
	jmp	.L2230
.L2472:
	movq	%rax, %rdi
	jmp	.L2262
.L2473:
	movq	%rax, %rdi
	jmp	.L1995
.L2477:
	jmp	.L2038
.L2391:
	movq	%rax, %rdi
	jmp	.L1992
.L2493:
	jmp	.L2038
.L2463:
	movq	%rax, %r14
	jmp	.L2160
.L2461:
	jmp	.L2612
.L2462:
	movq	%rax, %r15
	jmp	.L2173
.L2568:
	jmp	.L2569
.L2489:
	jmp	.L2038
.L2443:
	movq	%rax, %r14
	jmp	.L2097
	.section	.gcc_except_table
.LLSDA9961:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSE9961-.LLSDACSB9961
.LLSDACSB9961:
	.uleb128 .LEHB86-.LFB9961
	.uleb128 .LEHE86-.LEHB86
	.uleb128 .L2472-.LFB9961
	.uleb128 0
	.uleb128 .LEHB87-.LFB9961
	.uleb128 .LEHE87-.LEHB87
	.uleb128 .L2391-.LFB9961
	.uleb128 0
	.uleb128 .LEHB88-.LFB9961
	.uleb128 .LEHE88-.LEHB88
	.uleb128 .L2473-.LFB9961
	.uleb128 0
	.uleb128 .LEHB89-.LFB9961
	.uleb128 .LEHE89-.LEHB89
	.uleb128 .L2474-.LFB9961
	.uleb128 0
	.uleb128 .LEHB90-.LFB9961
	.uleb128 .LEHE90-.LEHB90
	.uleb128 .L2476-.LFB9961
	.uleb128 0
	.uleb128 .LEHB91-.LFB9961
	.uleb128 .LEHE91-.LEHB91
	.uleb128 .L2493-.LFB9961
	.uleb128 0
	.uleb128 .LEHB92-.LFB9961
	.uleb128 .LEHE92-.LEHB92
	.uleb128 .L2477-.LFB9961
	.uleb128 0
	.uleb128 .LEHB93-.LFB9961
	.uleb128 .LEHE93-.LEHB93
	.uleb128 .L2478-.LFB9961
	.uleb128 0
	.uleb128 .LEHB94-.LFB9961
	.uleb128 .LEHE94-.LEHB94
	.uleb128 .L2492-.LFB9961
	.uleb128 0
	.uleb128 .LEHB95-.LFB9961
	.uleb128 .LEHE95-.LEHB95
	.uleb128 .L2480-.LFB9961
	.uleb128 0
	.uleb128 .LEHB96-.LFB9961
	.uleb128 .LEHE96-.LEHB96
	.uleb128 .L2481-.LFB9961
	.uleb128 0
	.uleb128 .LEHB97-.LFB9961
	.uleb128 .LEHE97-.LEHB97
	.uleb128 .L2491-.LFB9961
	.uleb128 0
	.uleb128 .LEHB98-.LFB9961
	.uleb128 .LEHE98-.LEHB98
	.uleb128 .L2483-.LFB9961
	.uleb128 0
	.uleb128 .LEHB99-.LFB9961
	.uleb128 .LEHE99-.LEHB99
	.uleb128 .L2484-.LFB9961
	.uleb128 0
	.uleb128 .LEHB100-.LFB9961
	.uleb128 .LEHE100-.LEHB100
	.uleb128 .L2490-.LFB9961
	.uleb128 0
	.uleb128 .LEHB101-.LFB9961
	.uleb128 .LEHE101-.LEHB101
	.uleb128 .L2489-.LFB9961
	.uleb128 0
	.uleb128 .LEHB102-.LFB9961
	.uleb128 .LEHE102-.LEHB102
	.uleb128 .L2487-.LFB9961
	.uleb128 0
	.uleb128 .LEHB103-.LFB9961
	.uleb128 .LEHE103-.LEHB103
	.uleb128 .L2488-.LFB9961
	.uleb128 0
	.uleb128 .LEHB104-.LFB9961
	.uleb128 .LEHE104-.LEHB104
	.uleb128 .L2486-.LFB9961
	.uleb128 0
	.uleb128 .LEHB105-.LFB9961
	.uleb128 .LEHE105-.LEHB105
	.uleb128 .L2485-.LFB9961
	.uleb128 0
	.uleb128 .LEHB106-.LFB9961
	.uleb128 .LEHE106-.LEHB106
	.uleb128 .L2482-.LFB9961
	.uleb128 0
	.uleb128 .LEHB107-.LFB9961
	.uleb128 .LEHE107-.LEHB107
	.uleb128 .L2479-.LFB9961
	.uleb128 0
	.uleb128 .LEHB108-.LFB9961
	.uleb128 .LEHE108-.LEHB108
	.uleb128 .L2475-.LFB9961
	.uleb128 0
	.uleb128 .LEHB109-.LFB9961
	.uleb128 .LEHE109-.LEHB109
	.uleb128 .L2471-.LFB9961
	.uleb128 0
	.uleb128 .LEHB110-.LFB9961
	.uleb128 .LEHE110-.LEHB110
	.uleb128 .L2392-.LFB9961
	.uleb128 0
	.uleb128 .LEHB111-.LFB9961
	.uleb128 .LEHE111-.LEHB111
	.uleb128 .L2393-.LFB9961
	.uleb128 0
	.uleb128 .LEHB112-.LFB9961
	.uleb128 .LEHE112-.LEHB112
	.uleb128 .L2394-.LFB9961
	.uleb128 0
	.uleb128 .LEHB113-.LFB9961
	.uleb128 .LEHE113-.LEHB113
	.uleb128 .L2395-.LFB9961
	.uleb128 0
	.uleb128 .LEHB114-.LFB9961
	.uleb128 .LEHE114-.LEHB114
	.uleb128 .L2396-.LFB9961
	.uleb128 0
	.uleb128 .LEHB115-.LFB9961
	.uleb128 .LEHE115-.LEHB115
	.uleb128 .L2397-.LFB9961
	.uleb128 0
	.uleb128 .LEHB116-.LFB9961
	.uleb128 .LEHE116-.LEHB116
	.uleb128 .L2398-.LFB9961
	.uleb128 0
	.uleb128 .LEHB117-.LFB9961
	.uleb128 .LEHE117-.LEHB117
	.uleb128 .L2399-.LFB9961
	.uleb128 0
	.uleb128 .LEHB118-.LFB9961
	.uleb128 .LEHE118-.LEHB118
	.uleb128 .L2400-.LFB9961
	.uleb128 0
	.uleb128 .LEHB119-.LFB9961
	.uleb128 .LEHE119-.LEHB119
	.uleb128 .L2401-.LFB9961
	.uleb128 0
	.uleb128 .LEHB120-.LFB9961
	.uleb128 .LEHE120-.LEHB120
	.uleb128 .L2402-.LFB9961
	.uleb128 0
	.uleb128 .LEHB121-.LFB9961
	.uleb128 .LEHE121-.LEHB121
	.uleb128 .L2403-.LFB9961
	.uleb128 0
	.uleb128 .LEHB122-.LFB9961
	.uleb128 .LEHE122-.LEHB122
	.uleb128 .L2404-.LFB9961
	.uleb128 0
	.uleb128 .LEHB123-.LFB9961
	.uleb128 .LEHE123-.LEHB123
	.uleb128 .L2399-.LFB9961
	.uleb128 0
	.uleb128 .LEHB124-.LFB9961
	.uleb128 .LEHE124-.LEHB124
	.uleb128 .L2394-.LFB9961
	.uleb128 0
	.uleb128 .LEHB125-.LFB9961
	.uleb128 .LEHE125-.LEHB125
	.uleb128 .L2405-.LFB9961
	.uleb128 0
	.uleb128 .LEHB126-.LFB9961
	.uleb128 .LEHE126-.LEHB126
	.uleb128 .L2406-.LFB9961
	.uleb128 0
	.uleb128 .LEHB127-.LFB9961
	.uleb128 .LEHE127-.LEHB127
	.uleb128 .L2407-.LFB9961
	.uleb128 0
	.uleb128 .LEHB128-.LFB9961
	.uleb128 .LEHE128-.LEHB128
	.uleb128 .L2408-.LFB9961
	.uleb128 0
	.uleb128 .LEHB129-.LFB9961
	.uleb128 .LEHE129-.LEHB129
	.uleb128 .L2409-.LFB9961
	.uleb128 0
	.uleb128 .LEHB130-.LFB9961
	.uleb128 .LEHE130-.LEHB130
	.uleb128 .L2410-.LFB9961
	.uleb128 0
	.uleb128 .LEHB131-.LFB9961
	.uleb128 .LEHE131-.LEHB131
	.uleb128 .L2411-.LFB9961
	.uleb128 0
	.uleb128 .LEHB132-.LFB9961
	.uleb128 .LEHE132-.LEHB132
	.uleb128 .L2444-.LFB9961
	.uleb128 0
	.uleb128 .LEHB133-.LFB9961
	.uleb128 .LEHE133-.LEHB133
	.uleb128 .L2412-.LFB9961
	.uleb128 0
	.uleb128 .LEHB134-.LFB9961
	.uleb128 .LEHE134-.LEHB134
	.uleb128 .L2413-.LFB9961
	.uleb128 0
	.uleb128 .LEHB135-.LFB9961
	.uleb128 .LEHE135-.LEHB135
	.uleb128 .L2414-.LFB9961
	.uleb128 0
	.uleb128 .LEHB136-.LFB9961
	.uleb128 .LEHE136-.LEHB136
	.uleb128 .L2415-.LFB9961
	.uleb128 0
	.uleb128 .LEHB137-.LFB9961
	.uleb128 .LEHE137-.LEHB137
	.uleb128 .L2416-.LFB9961
	.uleb128 0
	.uleb128 .LEHB138-.LFB9961
	.uleb128 .LEHE138-.LEHB138
	.uleb128 .L2417-.LFB9961
	.uleb128 0
	.uleb128 .LEHB139-.LFB9961
	.uleb128 .LEHE139-.LEHB139
	.uleb128 .L2418-.LFB9961
	.uleb128 0
	.uleb128 .LEHB140-.LFB9961
	.uleb128 .LEHE140-.LEHB140
	.uleb128 .L2419-.LFB9961
	.uleb128 0
	.uleb128 .LEHB141-.LFB9961
	.uleb128 .LEHE141-.LEHB141
	.uleb128 .L2420-.LFB9961
	.uleb128 0
	.uleb128 .LEHB142-.LFB9961
	.uleb128 .LEHE142-.LEHB142
	.uleb128 .L2421-.LFB9961
	.uleb128 0
	.uleb128 .LEHB143-.LFB9961
	.uleb128 .LEHE143-.LEHB143
	.uleb128 .L2422-.LFB9961
	.uleb128 0
	.uleb128 .LEHB144-.LFB9961
	.uleb128 .LEHE144-.LEHB144
	.uleb128 .L2423-.LFB9961
	.uleb128 0
	.uleb128 .LEHB145-.LFB9961
	.uleb128 .LEHE145-.LEHB145
	.uleb128 .L2424-.LFB9961
	.uleb128 0
	.uleb128 .LEHB146-.LFB9961
	.uleb128 .LEHE146-.LEHB146
	.uleb128 .L2425-.LFB9961
	.uleb128 0
	.uleb128 .LEHB147-.LFB9961
	.uleb128 .LEHE147-.LEHB147
	.uleb128 .L2426-.LFB9961
	.uleb128 0
	.uleb128 .LEHB148-.LFB9961
	.uleb128 .LEHE148-.LEHB148
	.uleb128 .L2427-.LFB9961
	.uleb128 0
	.uleb128 .LEHB149-.LFB9961
	.uleb128 .LEHE149-.LEHB149
	.uleb128 .L2428-.LFB9961
	.uleb128 0
	.uleb128 .LEHB150-.LFB9961
	.uleb128 .LEHE150-.LEHB150
	.uleb128 .L2429-.LFB9961
	.uleb128 0
	.uleb128 .LEHB151-.LFB9961
	.uleb128 .LEHE151-.LEHB151
	.uleb128 .L2430-.LFB9961
	.uleb128 0
	.uleb128 .LEHB152-.LFB9961
	.uleb128 .LEHE152-.LEHB152
	.uleb128 .L2431-.LFB9961
	.uleb128 0
	.uleb128 .LEHB153-.LFB9961
	.uleb128 .LEHE153-.LEHB153
	.uleb128 .L2432-.LFB9961
	.uleb128 0
	.uleb128 .LEHB154-.LFB9961
	.uleb128 .LEHE154-.LEHB154
	.uleb128 .L2433-.LFB9961
	.uleb128 0
	.uleb128 .LEHB155-.LFB9961
	.uleb128 .LEHE155-.LEHB155
	.uleb128 .L2434-.LFB9961
	.uleb128 0
	.uleb128 .LEHB156-.LFB9961
	.uleb128 .LEHE156-.LEHB156
	.uleb128 .L2435-.LFB9961
	.uleb128 0
	.uleb128 .LEHB157-.LFB9961
	.uleb128 .LEHE157-.LEHB157
	.uleb128 .L2436-.LFB9961
	.uleb128 0
	.uleb128 .LEHB158-.LFB9961
	.uleb128 .LEHE158-.LEHB158
	.uleb128 .L2437-.LFB9961
	.uleb128 0
	.uleb128 .LEHB159-.LFB9961
	.uleb128 .LEHE159-.LEHB159
	.uleb128 .L2438-.LFB9961
	.uleb128 0
	.uleb128 .LEHB160-.LFB9961
	.uleb128 .LEHE160-.LEHB160
	.uleb128 .L2439-.LFB9961
	.uleb128 0
	.uleb128 .LEHB161-.LFB9961
	.uleb128 .LEHE161-.LEHB161
	.uleb128 .L2440-.LFB9961
	.uleb128 0
	.uleb128 .LEHB162-.LFB9961
	.uleb128 .LEHE162-.LEHB162
	.uleb128 .L2441-.LFB9961
	.uleb128 0
	.uleb128 .LEHB163-.LFB9961
	.uleb128 .LEHE163-.LEHB163
	.uleb128 .L2442-.LFB9961
	.uleb128 0
	.uleb128 .LEHB164-.LFB9961
	.uleb128 .LEHE164-.LEHB164
	.uleb128 .L2443-.LFB9961
	.uleb128 0
	.uleb128 .LEHB165-.LFB9961
	.uleb128 .LEHE165-.LEHB165
	.uleb128 .L2456-.LFB9961
	.uleb128 0
	.uleb128 .LEHB166-.LFB9961
	.uleb128 .LEHE166-.LEHB166
	.uleb128 .L2443-.LFB9961
	.uleb128 0
	.uleb128 .LEHB167-.LFB9961
	.uleb128 .LEHE167-.LEHB167
	.uleb128 .L2457-.LFB9961
	.uleb128 0
	.uleb128 .LEHB168-.LFB9961
	.uleb128 .LEHE168-.LEHB168
	.uleb128 .L2459-.LFB9961
	.uleb128 0
	.uleb128 .LEHB169-.LFB9961
	.uleb128 .LEHE169-.LEHB169
	.uleb128 .L2458-.LFB9961
	.uleb128 0
	.uleb128 .LEHB170-.LFB9961
	.uleb128 .LEHE170-.LEHB170
	.uleb128 .L2443-.LFB9961
	.uleb128 0
	.uleb128 .LEHB171-.LFB9961
	.uleb128 .LEHE171-.LEHB171
	.uleb128 .L2461-.LFB9961
	.uleb128 0
	.uleb128 .LEHB172-.LFB9961
	.uleb128 .LEHE172-.LEHB172
	.uleb128 .L2463-.LFB9961
	.uleb128 0
	.uleb128 .LEHB173-.LFB9961
	.uleb128 .LEHE173-.LEHB173
	.uleb128 .L2462-.LFB9961
	.uleb128 0
	.uleb128 .LEHB174-.LFB9961
	.uleb128 .LEHE174-.LEHB174
	.uleb128 .L2568-.LFB9961
	.uleb128 0
	.uleb128 .LEHB175-.LFB9961
	.uleb128 .LEHE175-.LEHB175
	.uleb128 .L2465-.LFB9961
	.uleb128 0
	.uleb128 .LEHB176-.LFB9961
	.uleb128 .LEHE176-.LEHB176
	.uleb128 .L2444-.LFB9961
	.uleb128 0
	.uleb128 .LEHB177-.LFB9961
	.uleb128 .LEHE177-.LEHB177
	.uleb128 .L2445-.LFB9961
	.uleb128 0
	.uleb128 .LEHB178-.LFB9961
	.uleb128 .LEHE178-.LEHB178
	.uleb128 .L2463-.LFB9961
	.uleb128 0
	.uleb128 .LEHB179-.LFB9961
	.uleb128 .LEHE179-.LEHB179
	.uleb128 .L2464-.LFB9961
	.uleb128 0
	.uleb128 .LEHB180-.LFB9961
	.uleb128 .LEHE180-.LEHB180
	.uleb128 .L2568-.LFB9961
	.uleb128 0
	.uleb128 .LEHB181-.LFB9961
	.uleb128 .LEHE181-.LEHB181
	.uleb128 .L2470-.LFB9961
	.uleb128 0
	.uleb128 .LEHB182-.LFB9961
	.uleb128 .LEHE182-.LEHB182
	.uleb128 .L2460-.LFB9961
	.uleb128 0
	.uleb128 .LEHB183-.LFB9961
	.uleb128 .LEHE183-.LEHB183
	.uleb128 .L2568-.LFB9961
	.uleb128 0
	.uleb128 .LEHB184-.LFB9961
	.uleb128 .LEHE184-.LEHB184
	.uleb128 .L2463-.LFB9961
	.uleb128 0
	.uleb128 .LEHB185-.LFB9961
	.uleb128 .LEHE185-.LEHB185
	.uleb128 .L2445-.LFB9961
	.uleb128 0
	.uleb128 .LEHB186-.LFB9961
	.uleb128 .LEHE186-.LEHB186
	.uleb128 .L2446-.LFB9961
	.uleb128 0
	.uleb128 .LEHB187-.LFB9961
	.uleb128 .LEHE187-.LEHB187
	.uleb128 .L2447-.LFB9961
	.uleb128 0
	.uleb128 .LEHB188-.LFB9961
	.uleb128 .LEHE188-.LEHB188
	.uleb128 .L2446-.LFB9961
	.uleb128 0
	.uleb128 .LEHB189-.LFB9961
	.uleb128 .LEHE189-.LEHB189
	.uleb128 .L2451-.LFB9961
	.uleb128 0
	.uleb128 .LEHB190-.LFB9961
	.uleb128 .LEHE190-.LEHB190
	.uleb128 .L2448-.LFB9961
	.uleb128 0
	.uleb128 .LEHB191-.LFB9961
	.uleb128 .LEHE191-.LEHB191
	.uleb128 .L2449-.LFB9961
	.uleb128 0
	.uleb128 .LEHB192-.LFB9961
	.uleb128 .LEHE192-.LEHB192
	.uleb128 .L2451-.LFB9961
	.uleb128 0
	.uleb128 .LEHB193-.LFB9961
	.uleb128 .LEHE193-.LEHB193
	.uleb128 .L2450-.LFB9961
	.uleb128 0
	.uleb128 .LEHB194-.LFB9961
	.uleb128 .LEHE194-.LEHB194
	.uleb128 .L2451-.LFB9961
	.uleb128 0
	.uleb128 .LEHB195-.LFB9961
	.uleb128 .LEHE195-.LEHB195
	.uleb128 .L2452-.LFB9961
	.uleb128 0
	.uleb128 .LEHB196-.LFB9961
	.uleb128 .LEHE196-.LEHB196
	.uleb128 .L2453-.LFB9961
	.uleb128 0
	.uleb128 .LEHB197-.LFB9961
	.uleb128 .LEHE197-.LEHB197
	.uleb128 .L2454-.LFB9961
	.uleb128 0
	.uleb128 .LEHB198-.LFB9961
	.uleb128 .LEHE198-.LEHB198
	.uleb128 .L2455-.LFB9961
	.uleb128 0
	.uleb128 .LEHB199-.LFB9961
	.uleb128 .LEHE199-.LEHB199
	.uleb128 .L2451-.LFB9961
	.uleb128 0
.LLSDACSE9961:
	.section	.text.startup
	.cfi_endproc
	.section	.text.unlikely
	.cfi_startproc
	.cfi_personality 0x9b,DW.ref.__gxx_personality_v0
	.cfi_lsda 0x1b,.LLSDAC9961
	.type	main.cold, @function
main.cold:
.LFSB9961:
.L2347:
	.cfi_def_cfa 6, 16
	.cfi_offset 3, -56
	.cfi_offset 6, -16
	.cfi_offset 12, -48
	.cfi_offset 13, -40
	.cfi_offset 14, -32
	.cfi_offset 15, -24
	movq	%rbx, %rdi
	movq	%r15, %r14
	call	_ZNSt7__cxx1119basic_istringstreamIcSt11char_traitsIcESaIcEED1Ev@PLT
.L2348:
	movq	-4336(%rbp), %rdi
	call	_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED1Ev
	movq	-1120(%rbp), %rdi
	movq	-4568(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2363
	movq	-1104(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2363:
	movq	-1152(%rbp), %rdi
	movq	-4560(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2641
.L2364:
	movq	-4464(%rbp), %rdi
	call	_ZNSt14basic_ifstreamIcSt11char_traitsIcEED1Ev@PLT
.L2365:
	movq	-3424(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
.L2291:
	movq	-3744(%rbp), %rdi
	call	free@PLT
.L2366:
	movq	-3776(%rbp), %rdi
	movq	%r15, %r14
	call	free@PLT
.L2367:
	movq	-3808(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
.L2368:
	movq	-3840(%rbp), %rdi
	movq	%r15, %r14
	call	free@PLT
.L2369:
	movq	-3872(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
.L2370:
	movq	-3904(%rbp), %rdi
	movq	%r15, %r14
	call	free@PLT
.L2371:
	movq	-3936(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
.L2372:
	movq	-3968(%rbp), %rdi
	movq	%r15, %r14
	call	free@PLT
.L2281:
	movq	-4224(%rbp), %rdi
	call	free@PLT
.L2373:
	movq	-4000(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
	movq	-4032(%rbp), %rdi
	call	free@PLT
	movq	-4064(%rbp), %rdi
	call	free@PLT
.L2374:
	movq	-4096(%rbp), %rdi
	movq	-4080(%rbp), %rsi
	subq	%rdi, %rsi
	testq	%rdi, %rdi
	je	.L2375
	call	_ZdlPvm@PLT
.L2375:
	movq	%r15, %r14
.L2377:
	subq	$32, %rbx
	movq	(%rbx), %rdi
	leaq	16(%rbx), %rdx
	cmpq	%rdx, %rdi
	je	.L2376
	movq	16(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2376:
	movq	-4600(%rbp), %rax
	cmpq	%rax, %rbx
	jne	.L2377
	movq	%r14, %rdi
.LEHB200:
	call	_Unwind_Resume@PLT
.L2349:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2642
.L2350:
	movq	%r14, %r15
.L2354:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2611
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2611:
	movq	%r15, %r14
.L2353:
	movq	-1376(%rbp), %rdi
	movq	-1360(%rbp), %rsi
	subq	%rdi, %rsi
	testq	%rdi, %rdi
	je	.L2348
	call	_ZdlPvm@PLT
	jmp	.L2348
.L2361:
	movq	-3264(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
.L2360:
	movq	-3296(%rbp), %rdi
	call	free@PLT
	jmp	.L2611
.L2358:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2360
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2360
.L2641:
	movq	-1136(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2364
.L2316:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2643
.L2317:
	movq	%r14, %r15
.L2318:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2319
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2319:
	movq	%r15, %r14
.L2320:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2608
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2608
.L2356:
	movq	-960(%rbp), %rdi
	leaq	-944(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2353
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2353
.L2323:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2324
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2324:
	movq	%r15, %r14
.L2325:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2326
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2326:
	movq	%r14, %r15
.L2327:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2328
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2328:
	movq	%r15, %r14
	jmp	.L2329
.L2227:
	movq	-1376(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
	movq	-1368(%rbp), %rdi
	call	free@PLT
.L2606:
	movq	-2336(%rbp), %rdi
	call	free@PLT
.L2163:
	movq	-3104(%rbp), %rdi
	movq	%r15, %r14
	call	free@PLT
.L2097:
	movq	-3424(%rbp), %rdi
	call	free@PLT
.L2342:
	movq	-3456(%rbp), %rdi
	call	free@PLT
.L2607:
	movq	-3488(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
.L2343:
	movq	-3520(%rbp), %rdi
	movq	%r15, %r14
	call	free@PLT
.L2344:
	movq	-3552(%rbp), %rdi
	call	free@PLT
.L2329:
	movq	-3584(%rbp), %rdi
	call	free@PLT
.L2608:
	movq	-4192(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
.L2315:
	movq	-3616(%rbp), %rdi
	call	free@PLT
.L2345:
	movq	-3648(%rbp), %rdi
	movq	%r15, %r14
	call	free@PLT
.L2346:
	movq	-3680(%rbp), %rdi
	call	free@PLT
.L2609:
	movq	-3712(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
.L2610:
	movq	-4208(%rbp), %rdi
	call	free@PLT
	jmp	.L2291
.L2330:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2644
.L2331:
	movq	%r14, %r15
.L2332:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2333
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2333:
	movq	%r15, %r14
.L2334:
	movq	-1120(%rbp), %rdi
	movq	-4568(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2335
	movq	-1104(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2335:
	movq	%r14, %r15
.L2336:
	movq	-1152(%rbp), %rdi
	movq	-4560(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2337
	movq	-1136(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2337:
	movq	%r15, %r14
.L2338:
	movq	-1184(%rbp), %rdi
	leaq	-1168(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2339
	movq	-1168(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2339:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2607
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2607
.L2642:
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2350
.L2306:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2645
.L2307:
	movq	%r15, %r14
.L2308:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2309
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2309:
	movq	%r14, %r15
.L2310:
	movq	-1120(%rbp), %rdi
	movq	-4568(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2311
	movq	-1104(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2311:
	movq	%r15, %r14
.L2312:
	movq	-1152(%rbp), %rdi
	movq	-4560(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2313
	movq	-1136(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2313:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2314
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2314:
	movq	%r14, %r15
	jmp	.L2315
.L2643:
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2317
.L2038:
	leaq	-960(%rbp), %rbx
	movq	%rax, %r14
	jmp	.L2377
.L2645:
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2307
.L2282:
	movq	-576(%rbp), %rdi
	leaq	-560(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2646
.L2283:
	movq	%r14, %r15
.L2284:
	movq	-1088(%rbp), %rdi
	leaq	-1072(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2285
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2285:
	movq	%r15, %r14
.L2286:
	movq	-1120(%rbp), %rdi
	leaq	-1104(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2287
	movq	-1104(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2287:
	movq	%r14, %r15
.L2288:
	movq	-1152(%rbp), %rdi
	leaq	-1136(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2289
	movq	-1136(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2289:
	movq	-960(%rbp), %rdi
	leaq	-944(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2291
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2291
.L2270:
	movq	-1376(%rbp), %rdi
	movq	%r15, %r14
	call	free@PLT
.L2271:
	movq	-576(%rbp), %rdi
	leaq	-560(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2647
.L2272:
	movq	%r14, %r15
.L2273:
	movq	-1088(%rbp), %rdi
	leaq	-1072(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2648
.L2274:
	movq	%r15, %r14
.L2275:
	movq	-1120(%rbp), %rdi
	leaq	-1104(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2649
.L2276:
	movq	%r14, %r15
.L2277:
	movq	-1152(%rbp), %rdi
	leaq	-1136(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2650
.L2278:
	movq	-960(%rbp), %rdi
	leaq	-944(%rbp), %rax
	cmpq	%rax, %rdi
	jne	.L2651
.L2279:
	movq	%r15, %r14
.L2280:
	movq	-1728(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
	jmp	.L2269
.L2648:
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2274
.L2647:
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2272
.L2649:
	movq	-1104(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2276
.L2299:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2300
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2300:
	movq	%r14, %r15
.L2301:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2302
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2302:
	movq	%r15, %r14
.L2303:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2609
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2609
.L2102:
	movq	-1376(%rbp), %rdi
	movq	%r15, %r14
	call	free@PLT
.L2103:
	movq	-3152(%rbp), %rdi
	call	free@PLT
	jmp	.L2097
.L2115:
	movq	-1376(%rbp), %rdi
	call	free@PLT
	jmp	.L2103
.L2292:
	movq	-576(%rbp), %rdi
	movq	-4456(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2293
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2293:
	movq	%r15, %r14
.L2294:
	movq	-960(%rbp), %rdi
	movq	-4384(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2295
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
.L2295:
	movq	%r14, %r15
.L2296:
	movq	-1088(%rbp), %rdi
	movq	-4400(%rbp), %rax
	cmpq	%rax, %rdi
	je	.L2610
	movq	-1072(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2610
.L2646:
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2283
.L2266:
	movq	-1504(%rbp), %rdi
	call	free@PLT
.L2267:
	movq	-1376(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
.L2268:
	movq	-1616(%rbp), %rdi
	call	free@PLT
.L2269:
	movq	-1840(%rbp), %rdi
	movq	%r15, %r14
	call	free@PLT
	jmp	.L2281
.L2159:
	movq	-1376(%rbp), %rdi
	call	free@PLT
	movq	-1368(%rbp), %rdi
	call	free@PLT
.L2160:
	movq	-1504(%rbp), %rdi
	movq	%r14, %r15
	call	free@PLT
	jmp	.L2163
.L2230:
	movq	-1376(%rbp), %rdi
	call	free@PLT
	movq	-1368(%rbp), %rdi
	call	free@PLT
	jmp	.L2606
.L2262:
	movl	$2, %eax
.L1993:
	movl	$2, %ebx
	movq	%rdi, %r15
	subq	%rax, %rbx
	movq	-4600(%rbp), %rax
	salq	$5, %rbx
	addq	%rax, %rbx
.L2265:
	movq	-4600(%rbp), %rax
	cmpq	%rax, %rbx
	jne	.L2652
	movq	%r15, %rdi
	call	_Unwind_Resume@PLT
.LEHE200:
.L2644:
	movq	-560(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2331
.L1992:
	movl	$1, %eax
	jmp	.L1993
.L2652:
	subq	$32, %rbx
	movq	(%rbx), %rdi
	leaq	16(%rbx), %rax
	cmpq	%rax, %rdi
	je	.L2265
	movq	16(%rbx), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2265
.L1995:
	xorl	%eax, %eax
	jmp	.L1993
.L2570:
.LEHB201:
	call	_ZN5Eigen8internal19throw_std_bad_allocEv
.LEHE201:
.L2173:
	movq	-1504(%rbp), %rdi
	call	free@PLT
	jmp	.L2163
.L2469:
.L2569:
	movq	%rax, %r15
	jmp	.L2606
.L2650:
	movq	-1136(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2278
.L2651:
	movq	-944(%rbp), %rax
	leaq	1(%rax), %rsi
	call	_ZdlPvm@PLT
	jmp	.L2279
	.cfi_endproc
.LFE9961:
	.section	.gcc_except_table
.LLSDAC9961:
	.byte	0xff
	.byte	0xff
	.byte	0x1
	.uleb128 .LLSDACSEC9961-.LLSDACSBC9961
.LLSDACSBC9961:
	.uleb128 .LEHB200-.LCOLDB166
	.uleb128 .LEHE200-.LEHB200
	.uleb128 0
	.uleb128 0
	.uleb128 .LEHB201-.LCOLDB166
	.uleb128 .LEHE201-.LEHB201
	.uleb128 .L2469-.LCOLDB166
	.uleb128 0
.LLSDACSEC9961:
	.section	.text.unlikely
	.section	.text.startup
	.size	main, .-main
	.section	.text.unlikely
	.size	main.cold, .-main.cold
.LCOLDE166:
	.section	.text.startup
.LHOTE166:
	.p2align 4
	.type	_GLOBAL__sub_I_main, @function
_GLOBAL__sub_I_main:
.LFB14814:
	.cfi_startproc
	pushq	%rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	leaq	_ZStL8__ioinit(%rip), %rbx
	movq	%rbx, %rdi
	call	_ZNSt8ios_base4InitC1Ev@PLT
	movq	_ZNSt8ios_base4InitD1Ev@GOTPCREL(%rip), %rdi
	movq	%rbx, %rsi
	popq	%rbx
	.cfi_def_cfa_offset 8
	leaq	__dso_handle(%rip), %rdx
	jmp	__cxa_atexit@PLT
	.cfi_endproc
.LFE14814:
	.size	_GLOBAL__sub_I_main, .-_GLOBAL__sub_I_main
	.section	.init_array,"aw"
	.align 8
	.quad	_GLOBAL__sub_I_main
	.weak	_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits
	.section	.rodata._ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits,"aG",@progbits,_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits,comdat
	.align 32
	.type	_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits, @gnu_unique_object
	.size	_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits, 201
_ZZNSt8__detail18__to_chars_10_implIjEEvPcjT_E8__digits:
	.string	"00010203040506070809101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899"
	.weak	_ZGVZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes
	.section	.bss._ZGVZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes,"awG",@nobits,_ZGVZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes,comdat
	.align 8
	.type	_ZGVZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes, @gnu_unique_object
	.size	_ZGVZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes, 8
_ZGVZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes:
	.zero	8
	.weak	_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes
	.section	.bss._ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes,"awG",@nobits,_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes,comdat
	.align 16
	.type	_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes, @gnu_unique_object
	.size	_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes, 24
_ZZN5Eigen8internal20manage_caching_sizesENS_6ActionEPlS2_S2_E12m_cacheSizes:
	.zero	24
	.local	_ZStL8__ioinit
	.comm	_ZStL8__ioinit,1,1
	.section	.rodata.cst16,"aM",@progbits,16
	.align 16
.LC7:
	.long	0
	.long	0
	.long	0
	.long	-2147483648
	.align 16
.LC8:
	.long	0
	.long	-2147483648
	.long	0
	.long	0
	.set	.LC23,.LC7+8
	.set	.LC24,.LC86
	.align 16
.LC45:
	.quad	3
	.quad	3
	.section	.data.rel.ro,"aw"
	.align 8
.LC48:
	.quad	_ZTVNSt7__cxx1115basic_stringbufIcSt11char_traitsIcESaIcEEE+16
	.align 8
.LC49:
	.quad	_ZTVSt15basic_streambufIcSt11char_traitsIcEE+16
	.section	.rodata.cst2,"aM",@progbits,2
	.align 2
.LC50:
	.byte	10
	.byte	0
	.align 2
.LC51:
	.byte	32
	.byte	0
	.section	.rodata.cst16
	.align 16
.LC86:
	.long	0
	.long	1072693248
	.long	0
	.long	0
	.align 16
.LC127:
	.long	0
	.long	-1073741824
	.long	0
	.long	-1073741824
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC145:
	.long	0
	.long	1071644672
	.align 8
.LC154:
	.long	-755914244
	.long	1062232653
	.hidden	DW.ref.__gxx_personality_v0
	.weak	DW.ref.__gxx_personality_v0
	.section	.data.rel.local.DW.ref.__gxx_personality_v0,"awG",@progbits,DW.ref.__gxx_personality_v0,comdat
	.align 8
	.type	DW.ref.__gxx_personality_v0, @object
	.size	DW.ref.__gxx_personality_v0, 8
DW.ref.__gxx_personality_v0:
	.quad	__gxx_personality_v0
	.hidden	__dso_handle
	.ident	"GCC: (GNU) 12.2.1 20230111"
	.section	.note.GNU-stack,"",@progbits
