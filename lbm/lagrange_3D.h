#pragma once

#include <vector>
#include <algorithm>
#include <math.h>
#include "defs.h"
#include "lbm.h"
#include "spmatrix.h"


#ifdef USE_TNL
	#include <TNL/Matrices/SparseMatrix.h>
	#include <TNL/Algorithms/Segments/SlicedEllpack.h>
	#include <TNL/Pointers/DevicePointer.h>
	#include <TNL/Pointers/SharedPointer.h>
	#include <TNL/Algorithms/ParallelFor.h>
	#include <TNL/Solvers/Linear/CG.h>
	#include <TNL/Solvers/Linear/Preconditioners/Diagonal.h>
	#include <TNL/Solvers/Linear/Preconditioners/ILU0.h>
	#include <TNL/Allocators/CudaHost.h>

	template< typename Device, typename Index, typename IndexAlocator >
	using SlicedEllpackSegments = TNL::Algorithms::Segments::SlicedEllpack< Device, Index, IndexAlocator >;
	template< typename Real, typename Device, typename Index >
	using SlicedEllpack = TNL::Matrices::SparseMatrix< Real,
	                                                   Device,
	                                                   Index,
	                                                   TNL::Matrices::GeneralMatrix,
	                                                   SlicedEllpackSegments
	                                                 >;
#endif

#ifdef USE_CUSPARSE
	#include "sprectmatrix.h"
#endif

enum { 
	ws_computeCPU,		// use UMFPACK or PETSC ... TODO, FIXME
	ws_computeGPU_CUSPARSE,
	ws_computeHybrid_CUSPARSE,
	ws_computeCPU_TNL,
	ws_computeGPU_TNL,
	ws_computeHybrid_TNL,
	ws_computeHybrid_TNL_zerocopy
};

// cyclic vector
template < typename T >
class CyclicVector : public std::vector<T>
{
public:
    T& operator[](int n)
    {
                if (n<0)
                        return operator[](n+std::vector<T>::size());
                if (n>=std::vector<T>::size())
                        return operator[](n-std::vector<T>::size());
                return std::vector<T>::operator[](n);
//        int n=in;
//        if (in<0) n=std::vector<T>::size()+in; else
//        if (in>=std::vector<T>::size()) n=in-std::vector<T>::size();
//        return std::vector<T>::operator[](n);
    }
};

template < typename REAL>
struct LagrangePoint3D
{
	REAL x=0,y=0,z=0;
	REAL x_ref=0, y_ref=0, z_ref=0;
	REAL vx=0, vy=0, vz=0;
	REAL fx=0, fy=0, fz=0;
	// Lagrangian coordinates of the surface (as a grid)
	int lag_x, lag_y;
};


template< typename LBM >
struct Lagrange3D
{
	using T_TRAITS = typename LBM::T_TRAITS;
	using T_LBM = LBM;
	
	using idx = typename T_TRAITS::idx;
	using dreal = typename T_TRAITS::dreal;
	using real = typename T_TRAITS::real;

	T_LBM &lbm;

#ifdef USE_TNL
//	using hVector = TNL::Containers::Vector< real, TNL::Devices::Host, idx, TNL::Allocators::CudaHost<real> >;
	using hVector = TNL::Containers::Vector< dreal, TNL::Devices::Host, idx, TNL::Allocators::CudaHost<dreal> >;
	using dVector = TNL::Containers::Vector< dreal, TNL::Devices::Cuda, idx >;
//	using hEllpack = SlicedEllpack< real, TNL::Devices::Host, idx >;
	using hEllpack = SlicedEllpack< dreal, TNL::Devices::Host, idx >;
	using dEllpack = SlicedEllpack< dreal, TNL::Devices::Cuda, idx >;
	using hEllpackPtr = TNL::Pointers::SharedPointer< hEllpack >;
	using dEllpackPtr = TNL::Pointers::SharedPointer< dEllpack >;
	
	// ws_ using sparse matrices
	hEllpackPtr ws_tnl_hA;
	hEllpack ws_tnl_hM; // matrix realizing projection of u* to lagrange desc.
	hEllpack ws_tnl_hMT; // matrix realizing projection of uB from lagrange desc. to Euler desc. .... basially transpose of M
	hVector ws_tnl_hx[3], ws_tnl_hb[3];

	using hSolver = TNL::Solvers::Linear::CG< hEllpack >;
	using hPreconditioner = TNL::Solvers::Linear::Preconditioners::Diagonal< hEllpack >;
//	using hPreconditioner = TNL::Solvers::Linear::Preconditioners::ILU0< hEllpack >;
	hSolver ws_tnl_hsolver;
	typename std::shared_ptr< hPreconditioner > ws_tnl_hprecond;

	#ifdef USE_CUDA
	dEllpackPtr ws_tnl_dA; // square matrix A
	dEllpack ws_tnl_dM; // matrix realizing projection of u* to lagrange desc.
	dEllpack ws_tnl_dMT; // matrix realizing projection of uB from lagrange desc. to Euler desc. .... basially transpose of M
	dVector ws_tnl_dx[3], ws_tnl_db[3];
	// for ws_computeHybrid_TNL_zerocopy
	using hVectorPinned = TNL::Containers::Vector< dreal, TNL::Devices::Host, idx, TNL::Allocators::CudaHost<dreal> >;
	hVectorPinned ws_tnl_hxz[3], ws_tnl_hbz[3];

	using dSolver = TNL::Solvers::Linear::CG< dEllpack >;
	using dPreconditioner = TNL::Solvers::Linear::Preconditioners::Diagonal< dEllpack >;
//	using dPreconditioner = TNL::Solvers::Linear::Preconditioners::ILU0< dEllpack >;
	dSolver ws_tnl_dsolver;
	typename std::shared_ptr< dPreconditioner > ws_tnl_dprecond;
	#endif
	bool ws_tnl_constructed=false;
#endif // USE_TNL

	// WuShu using sparse matrices
	SpMatrix<real> *ws_A;
	#ifdef USE_CUSPARSE
	SpRectMatrix<dreal> *ws_dA; // square matrix ws_A
	SpRectMatrix<dreal> *ws_M; // matrix realizing projection of u* to lagrange desc.
	SpRectMatrix<dreal> *ws_MT; // matrix realizing projection of uB from lagrange desc. to Euler desc. .... basially transpose of M
	#endif
	real *ws_x[3], *ws_b[3];
	dreal *ws_hx[3], *ws_hb[3];
	dreal *ws_dx[3], *ws_db[3], *ws_du;
	int Solver = SOLVER_UMFPACK;

	// sparse vectors
	std::vector<idx> *d_i=0;
	std::vector<real> *d_x=0;
	
	bool indexed=false; // is i,j lagrangian index created?
	int **index_array=0;
	
	// Lagrange surface dimensions
	int lag_X=-1; // size
	int lag_Y=-1; // size

	bool ws_regularDirac=true;		// use continuous ws_ trick with 2 dirac functions
	int ws_compute=ws_computeCPU;		// ws_computeCPU, ws_computeGPU, ws_computeHybrid
	bool ws_speedUpAllocation=false;	// choose neighbors based on lag_x and lag_y proximity !!! experimental
	int ws_speedUpAllocationSupport=1000000; // very big

	void constructWuShuMatricesSparse();
	void constructWuShuMatricesSparse_TNL();
	void computeWuShuForcesSparse(real time);

	bool ws_constructed=false;	// Wu Shu matrices constructed?

	bool implicitWuShuForcing=true; // use implicit WuShu Forcing (= run extra kernels etc...)

	real maxDist;			// maximal distance between points
	real minDist;			// minimal distance between points
	int diracDeltaType=2;

	inline int N() { return LL.size(); }	// return the amount of points in the filament

	CyclicVector<LagrangePoint3D<real>> LL;		// by this

	real diracDelta(int i, real r);
	real diracDelta(real r) { return diracDelta(diracDeltaType, r); }

	real dist(LagrangePoint3D<real> &A, LagrangePoint3D<real> &B) { return NORM( A.x - B.x, A.y - B.y, A.z - B.z ); }
	real dist_ref(LagrangePoint3D<real> &A, LagrangePoint3D<real> &B) { return NORM( A.x_ref - B.x_ref, A.y_ref - B.y_ref, A.z_ref - B.z_ref ); }
	real distx(LagrangePoint3D<real> &A, LagrangePoint3D<real> &B) { return fabs(A.x - B.x); }
	real disty(LagrangePoint3D<real> &A, LagrangePoint3D<real> &B) { return fabs(A.y - B.y); }
	real distz(LagrangePoint3D<real> &A, LagrangePoint3D<real> &B) { return fabs(A.z - B.z); }
	real distx_ref(LagrangePoint3D<real> &A, LagrangePoint3D<real> &B) { return fabs(A.x_ref - B.x_ref); }
	real disty_ref(LagrangePoint3D<real> &A, LagrangePoint3D<real> &B) { return fabs(A.y_ref - B.y_ref); }
	real distz_ref(LagrangePoint3D<real> &A, LagrangePoint3D<real> &B) { return fabs(A.z_ref - B.z_ref); }
	
	int findIndex(int i, int j);
	int createIndexArray();
	int findIndexOfNearestX(real x);

	void computeMaxMinDist();	// computes min and max distance between neinghboring nodes
	real computeMinDist();		// computes min and max distance between neinghboring nodes
	real computeMaxDistFromMinDist(real mindist);		// computes min and max distance between neighboring nodes
//	void integrateForce(real &Fx, real &Fy, real &Fz, real surface_element_size) { printf("integrateForce not implemented yet."); }

	// special log file for the linear system solvers
	char logfile[FILENAME_CHARS];

	template < typename... ARGS >
	void log(const char* fmt, ARGS... args);

	// constructors
	Lagrange3D(T_LBM &inputLBM, const char* resultsDir);
	~Lagrange3D();

	// disable copy-constructor and copy-assignment, leave only move-constructor and move-assignment
	// (because this class has a "T_LBM &lbm;" member)
	Lagrange3D(const Lagrange3D&) = delete;
	Lagrange3D(Lagrange3D&&) = default;
	Lagrange3D& operator=(const Lagrange3D&) = delete;
	Lagrange3D& operator=(Lagrange3D&&) = default;
};

#include "lagrange_3D.hpp"
