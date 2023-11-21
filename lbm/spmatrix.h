#ifndef __SpMatrix_H
#define __SpMatrix_H

#include "defs.h"

#ifdef USE_UMFPACK
#include "suitesparse/umfpack.h"
#endif

#ifdef USE_PETSC
#include <petscksp.h>
#endif

//#define CSC_FORMAT 1
//#define CSR_FORMAT 2

//typedef real T;

template < typename REAL >
struct SpMatrix
{
	enum { CSC_FORMAT, CSR_FORMAT };

	int m_n;    // size of the matrix Nb*NbEq
	int m_nz;   // number of all non-zeros entries of the matrix
	REAL null;	// dummy value to return
	
	bool doFactorizeOnce; // factorize the matrix once only and store the information for future use
	bool factorized; // is the matrix factorized?
	
	int Solver; // SOLVER_UMFPACK

	#ifdef USE_PETSC
	Vec            PETSC_x, PETSC_b;      /* approx solution, RHS, exact solution */
	Mat            PETSC_A;            /* linear system matrix */
	KSP            PETSC_ksp;         /* linear solver context */
	bool PETSC_INITIALIZED;
	#endif
	
	#ifdef USE_UMFPACK
	bool UMFPACK_INITIALIZED;
	void *Symbolic, *Numeric;
	double umf_Info [UMFPACK_INFO], umf_Control [UMFPACK_CONTROL];
	#endif

////	// CSR or CSC format based on FLAG_CSC_format=true: CSC format or FLAG_CSC_format=false: CSR format
	int *Ai;
	int *Ap;
	REAL  *Ax;

	int xpos(int i, int j);

	inline int N() { return m_n; }
	inline int Nz() { return m_nz; }

	inline REAL  &get(int i, int j) { return get(i,j,1); }
	REAL  &get(int i, int j, int verbose);		/// returns exact position of a matrix entry
	
	inline REAL & operator()(int i, int j) { return get(i,j,1); }
	inline REAL & operator()(int i, int j, int verbose) { return get(i,j,verbose); }

	int empty();
	int print_nonzeros();
	int print();
	int print2();
	int print2d();

	int solve(REAL  *b, REAL  *x);
	
	int initializePETSC();
	int initializeUMFPACK();
	
	int solveUMFPACK(REAL  *b, REAL  *x);
	int solvePETSC(REAL  *b, REAL  *x);

	bool isEntryInArray(int *array, int entry, int pos);
	int arrayAddEntry(int *array, int entry, int &pos);
	int arrayAddEntry(int *array, int entry, int &pos, int max);

	int transpose();

	int matrixFormat();

	SpMatrix();
	~SpMatrix();
};

#include "spmatrix.hpp"

#endif
