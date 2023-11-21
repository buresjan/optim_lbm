// SpMatrix
// (c) Radek Fučík
// 2013, 2014, 2015, 2016
// Main class for the sparse matrix storage format
// spMatrix: symmetric structure CSC vs. CSR format problems

#ifndef __SpRectMatrix_H
#define __SpRectMatrix_H

#include "defs.h"

#ifdef __CUDACC__

    // CUDA Runtime
    #include <cuda_runtime.h>

    // Using updated (v2) interfaces for CUBLAS and CUSPARSE
    #include <cusparse_v2.h>
    #include <cublas_v2.h>
    
    // debug
//    #include "helper_cuda.h"
    

#endif

#ifdef USE_CUDA

//#define CSC_FORMAT 1
//#define CSR_FORMAT 2

template < typename DREAL > 
struct SpRectMatrix
{
	enum { CSC_FORMAT, CSR_FORMAT };

#ifdef __CUDACC__
    // cublas and cusparse wrapper
	cusparseStatus_t 
	cusparseTcsrmv(cusparseHandle_t handle, cusparseOperation_t transA, 
	int m, int n, int nnz, const float           *alpha, 
	const cusparseMatDescr_t descrA, 
	const float           *csrValA, 
	const int *csrRowPtrA, const int *csrColIndA,
	const float           *x, const float           *beta, 
	float           *y) { 
		return cusparseScsrmv(handle, transA, m, n, nnz, alpha, descrA,  csrValA,  csrRowPtrA,  csrColIndA, x,  beta,  y);
	}
	cusparseStatus_t 
	cusparseTcsrmv(cusparseHandle_t handle, cusparseOperation_t transA, 
	int m, int n, int nnz, const double          *alpha, 
	const cusparseMatDescr_t descrA, 
	const double          *csrValA, 
	const int *csrRowPtrA, const int *csrColIndA,
	const double          *x, const double          *beta, 
	double          *y) { 
		return cusparseDcsrmv(handle, transA, m, n, nnz, alpha, descrA,  csrValA,  csrRowPtrA,  csrColIndA, x,  beta,  y);
	}
	cublasStatus_t cublasTaxpy(cublasHandle_t handle, int n,
	const float           *alpha,
	const float           *x, int incx,
	float                 *y, int incy) {
		return cublasSaxpy( handle, n, alpha, x,incx, y, incy);
	}
	cublasStatus_t cublasTaxpy(cublasHandle_t handle, int n,
	const double          *alpha,
	const double          *x, int incx,
	double                *y, int incy)	{
		return cublasDaxpy( handle, n, alpha, x,incx, y, incy);
	}
	cublasStatus_t cublasTdot (cublasHandle_t handle, int n,
	const float           *x, int incx,
	const float           *y, int incy,
	float           *result) {
		return cublasSdot( handle, n, x, incx, y, incy, result);
	}
	cublasStatus_t cublasTdot (cublasHandle_t handle, int n,
	const double          *x, int incx,
	const double          *y, int incy,
	double          *result)    {
		return cublasDdot( handle, n, x, incx, y, incy, result);
	}
	cublasStatus_t  cublasTscal(cublasHandle_t handle, int n,
	const float           *alpha,
	float           *x, int incx)	{
		return cublasSscal( handle, n, alpha, x, incx);
	}
	cublasStatus_t  cublasTscal(cublasHandle_t handle, int n,
	const double          *alpha,
	double          *x, int incx) 	{
		return cublasDscal( handle, n, alpha, x, incx);
	}
	cublasStatus_t cublasTcopy(cublasHandle_t handle, int n,
	const float           *x, int incx,
	float                 *y, int incy)	{
		return cublasScopy( handle, n, x, incx, y, incy);
	}
	cublasStatus_t cublasTcopy(cublasHandle_t handle, int n,
	const double          *x, int incx,
	double                *y, int incy)	{
		return cublasDcopy( handle, n, x, incx, y, incy);
	}
#endif

	int m_nr; // number of rows
	int m_nc; // number of cols
	int m_nz; // number of all non-zeros entries of the matrix
	DREAL null; // dummy value to return
	
	////	// CSR or CSC format based on FLAG_CSC_format=true: CSC format or FLAG_CSC_format=false: CSR format
	int *Ai;
	int *Ap;
	DREAL *Ax;
	
	
#ifdef __CUDACC__
	// CUDA
	DREAL *dAx;
	int *dAi;
	int *dAp;
	DREAL *dy, *dp;

	int *dunitAi; // unit matrix in CSR format: cublas cannot multiply 2 vectors per components ....ugh
	int *dunitAp; // unit matrix in CSR format 
	DREAL *dunitV; // auxiliary vector
	
	cusparseMatDescr_t descr;
	cusparseHandle_t cusparseHandle;
	cublasHandle_t cublasHandle;
	
	void cusparseInit();
	void unitMatrixInit(int N);
#endif
	int xpos(int i, int j);
	// solver setup

	DREAL solverTolerance;
	int solverMaxIter;

	inline int Nr() { return m_nr; }
	inline int Nc() { return m_nc; }
	inline int Nz() { return m_nz; }

	inline DREAL &get(int i, int j) { return get(i,j,1); }
	DREAL &get(int i, int j, int verbose);		/// returns exact position of a matrix entry
	
	inline DREAL& operator()(int i, int j) { return get(i,j,1); }
	inline DREAL& operator()(int i, int j, int verbose) { return get(i,j,verbose); }
	
	// matrix * input = output
//	void multiply(DREAL *input, DREAL*output) { return multiply(input, output, 1.0); }
//	void multiply(DREAL *input, DREAL*output, DREAL multiplicator);

//	void multiply(DREAL *input, real*output) { return multiply(input, output, 1.0); }
//	void multiply(DREAL *input, real*output, real multiplicator);

	// CUDA CUSPARSE
	void dmultiply(DREAL *input, DREAL*output) { return dmultiply(input, output, (DREAL)1.0); }
	void dmultiply(DREAL *input, DREAL*output, DREAL multiplicator);
	
	void dVectorMultiply(int N, DREAL*a, DREAL *b, DREAL alpha); // b[i] *= a[i]
	void dVectorAddition(int N, DREAL*a, DREAL* b, DREAL alpha); // a[i] += b[i]
	void dVectorMultiplyAndAdd(int N, DREAL*a, DREAL*b, DREAL alpha, DREAL*output); // output[i] += a[i]*b[i]*alpha

	int empty();
	int print_nonzeros();
	int print();
	int print2();
	int print2d();

	int dsolve(DREAL *b, DREAL *x);
	int dsolveCG(DREAL *b, DREAL *x);

	int matrixFormat();

	SpRectMatrix();
	virtual ~SpRectMatrix();
};

#include "sprectmatrix.hpp"

#endif // USE_CUDA
#endif