#ifndef __DEFS_H
#define __DEFS_H

#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <iostream>
#include <png.h>
#include "ciselnik.h"

#ifdef __CUDACC__
	#define CUDA_HOSTDEV __host__ __device__
	#define CUDA_HOSTDEV_NOINLINE CUDA_HOSTDEV __noinline__
#else
	#define CUDA_HOSTDEV
	#define CUDA_HOSTDEV_NOINLINE
#endif

template < 
	typename _dreal = float,	// real number representation on GPU
	typename _real = double,	// real number representation on CPU
	typename _idx = long int,	// array index on CPU and GPU (can be very large)
	typename _map_t = short int 
>
struct Traits
{
	using real = _real;
	using dreal = _dreal;
	using idx = _idx;
	using map_t = _map_t;
};

using TraitsSP = Traits<float>; //_dreal is float only
using TraitsDP = Traits<double>;

template < typename REAL >
struct KernelStruct 
{
	REAL fz=0,fx=0,fy=0;
	REAL f[27];
	REAL vz=0,vx=0, vy=0, rho=1.0, lbmViscosity=1.0;
};


//#define USE_HIGH_PRECISION_RHO // use num value ordering to compute rho inlbm_common.h .. slow!!!
//#define USE_GALILEAN_CORRECTION // Geier 2015: use Gal correction in BKG and CUM?
//#define USE_GEIER_CUM_2017 // use Geier 2017 Cummulant improvement A,B terms
//#define USE_GEIER_CUM_ANTIALIAS // use antialiasing Dxu, Dyv, Dzw from Geier 2015/2017

#define MAX( a , b) (((a)>(b))?(a):(b))
#define MIN( a , b) (((a)<(b))?(a):(b))

#define FILENAME_CHARS 500

#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code
#define NORM(x, y, z) sqrt(SQ(x) + SQ(y) + SQ(z))

#define POS(x,y,z,X,Y,Z) (y+(z)*(Y)+(x)*(Y)*(Z))
#define Fxyz(q,x,y,z,X,Y,Z) (y+(z)*(Y)+(q)*(Y)*(Z)+(x)*(Y)*(Z)*27)

// i=0,.., number of macro variables
#define mPOS(i, x, y, z, X, Y, Z) (POS(x,y,z,X,Y,Z)+(i)*(X)*(Y)*(Z))

enum { SOLVER_UMFPACK, SOLVER_PETSC };

// number of dist. functions, default=2 
// quick fix, use templates to define DFMAX ... through TRAITS maybe ?
#ifdef USE_DFMAX3
enum { df_cur, df_out, df_prev, DFMAX }; // special 3 dfs 
#else
enum { df_cur, df_out, DFMAX }; // default 2 dfs
#endif

// opposite directions have IDs different by 13
// (first half)
enum 
{
pzz=0,
zpz=1,
zzp=2,
ppz=3,
pzp=4,
zpp=5,
ppp=6,
ppm=7,
pmp=8,
mpp=9,
zpm=10,
pzm=11,
pmz=12,
// (second half)
mzz=13,
zmz=14,
zzm=15,
mmz=16,
mzm=17,
zmm=18,
mmm=19,
mmp=20,
mpm=21,
pmm=22,
zmp=23,
mzp=24,
mpz=25,
// (central)
zzz=26
};

#if defined(USE_CUDA) && defined(__CUDACC__)
	#define max_cuda_streams 27
	static cudaStream_t cuda_streams[max_cuda_streams];
#endif

#define Main main // TNL fix when LBM is included into TNL
	
	
#ifdef USE_CUDA
	#define checkCudaDevice checkDevice( __FILE__, __LINE__, cudaGetLastError() )

	bool checkDevice( const char* file_name, int line, cudaError error ) 
	{
		if( error == cudaSuccess ) return true;
		std::cerr << "CUDA ERROR(" << error << ") at line " << line << " in " << file_name << ":" << std::endl;
		std::cerr << cudaGetErrorString( error )  << std::endl;
		abort();
		return true;
	}
	#include <cuda_profiler_api.h>
#endif // USE_CUDA

#endif
