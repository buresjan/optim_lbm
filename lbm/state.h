#ifndef __STATE_H__
#define __STATE_H__ 

#include <deque>
#include <limits>
#include <chrono>

//#include <stdarg.h>
#include <sys/stat.h>

// spawn process vtkZip
#include <unistd.h>
#include <spawn.h>


#include "defs.h"
#include "lbm.h"
#include "vtk_writer.h"
#include "png_tool.h"

// ibm: lagrangian filament/surface
#include "lagrange_3D.h"

#include <sys/wait.h>

// sparse box origin+length where to plot - can be just a part of the domain
template < typename IDX >
struct probe3Dcut
{
	IDX ox,oy,oz; // lower left front point
	IDX lx,ly,lz; // length
	IDX step; // 1: every voxel 2: every 3 voxels etc.
	char name[FILENAME_CHARS];
	int cycle;
};


template < typename IDX >
struct probe2Dcut
{
	int type; // 0=X, 1=Y, 2=Z
	char name[FILENAME_CHARS];
	IDX position; // x/y/z ... LBM units ... int
	int cycle;
};

template < typename IDX >
struct probe1Dcut
{
	int type; // 0=X, 1=Y, 2=Z
	char name[FILENAME_CHARS];
	IDX pos1; // x/y/z
	IDX pos2; // y/z
	int cycle;
};

template < typename REAL >
struct probe1Dlinecut
{
	char name[FILENAME_CHARS];
	REAL from[3]; // physical units
	REAL to[3];   // physical units
	int cycle;
};

// for print/stat/write/reset counters
template < typename REAL >
struct counter
{
	int count=0;
	REAL period=-1.0;
//	bool copyMemoryInSimupdate=true;
	bool action(REAL time) { return (period>0 && time >= count * period) ? true : false; }
};

enum { STAT_RESET, STAT2_RESET, PRINT, VTK1D, VTK2D, VTK3D, PROBE1, PROBE2, PROBE3, SAVESTATE, VTK3DCUT, MAX_COUNTER };
enum { MemoryToFile, FileToMemory };
enum { vtk3DsingleFile, vtk3DmanyFiles, vtk3DmanyFilesExtraHeader };


template<
	typename LBM_TYPE,
	typename MACRO,
	typename CPU_MACRO,
	typename LBM_DATA,
	typename LBM_BC
>
struct State 
{
	using T_MACRO = MACRO;
	using T_CPU_MACRO = CPU_MACRO;
	using T_LBM_TYPE = LBM_TYPE;
	using T_LBM_DATA = LBM_DATA;
	using T_LBM = LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>;
	using T_LBM_BC = LBM_BC;
	using T_LBM_EQ = typename LBM_TYPE::T_LBM_EQ;
	using T_TRAITS = typename LBM_TYPE::T_TRAITS;
	using T_Lagrange3D = Lagrange3D< T_LBM >;

	using dreal = typename  T_TRAITS::dreal;
	using idx = typename T_TRAITS::idx;
	using real = typename T_TRAITS::real;

	using T_PROBE3DCUT = probe3Dcut<idx>;
	using T_PROBE2DCUT = probe2Dcut<idx>;
	using T_PROBE1DCUT = probe1Dcut<idx>;
	using T_PROBE1DLINECUT = probe1Dlinecut<real>;
	using T_COUNTER = counter<real>;

	T_LBM lbm;
//	VTKWriter vtk;
	
	std::vector< T_PROBE3DCUT > probe3Dvec;
	std::vector< T_PROBE2DCUT > probe2Dvec;
	std::vector< T_PROBE1DCUT > probe1Dvec;
	std::vector< T_PROBE1DLINECUT > probe1Dlinevec;
	
	// Lagrange
	std::deque<T_Lagrange3D> FF;			// array of filaments (std::deque instead of std::vector because T_Lagrange3D is not copy-constructible)
	int addLagrange3D();				// add a filament into the array and returns its index 
	void computeAllLagrangeForces();
	
	// vtk surface rotational
	void writeVTK_Surface(const char* name, real time, int cycle, T_Lagrange3D &fil);
	void writeVTK_Points(const char* name, real time, int cycle, T_Lagrange3D &fil);

//	template < typename... ARGS >
//	int run_cmd(const char *fmt, ARGS... args);

	// how often to probe/print/write/stat
	T_COUNTER cnt[MAX_COUNTER];
	virtual void probe1() {  }
	virtual void probe2() {  }
	virtual void probe3() {  }
	virtual void statReset() { }
	virtual void stat2Reset() { }

	// vtk export
	template< typename real1, typename real2 >
	bool vtk_helper(const char*iid, real1 ivalue, int idofs, char*id, real2 &value, int &dofs) /// simplifies data output routine
	{
		sprintf(id,"%s",iid);
		dofs=idofs;
		value=ivalue;
		return true;
	}
//	void probe2D() { writeVTKs_2D(); }; // replaced by writeVTKs_2D
	void writeVTKs_2D();
	void writeVTK_2DcutX(const char* name, real time, int cycle, idx XPOS);
	void writeVTK_2DcutY(const char* name, real time, int cycle, idx YPOS);
	void writeVTK_2DcutZ(const char* name, real time, int cycle, idx ZPOS);
	template < typename... ARGS >
	void add2Dcut_X(idx x, const char* fmt, ARGS... args);
	template < typename... ARGS >
	void add2Dcut_Y(idx y, const char* fmt, ARGS... args);
	template < typename... ARGS >
	void add2Dcut_Z(idx z, const char* fmt, ARGS... args);

	//bool VTK3D_write_separate_header=true; // true = split 3D vtks to ".prevtk" files header + data ... false = single standalone vtks
	int vtk3Dstyle = vtk3DsingleFile; // enum { vtk3DsingleFile, vtk3DmanyFiles, vtk3DmanyFilesExtraHeader };

	void writeVTKs_3D();
//	void writeVTK(const char* name, real time, int cycle) { writeVTK_3D(name, time, cycle); } // replaced by writeVTK_3D
	void writeVTK_3D(const char* name, real time, int cycle); 
	void writeVTK_3D_singlefile(const char* name, real time, int cycle); 

	// 3D cuts
	void writeVTKs_3Dcut();
	template < typename... ARGS >
	void add3Dcut(idx ox, idx oy, idx oz, idx lx, idx ly, idx lz, idx step, const char* fmt, ARGS... args);
//	void writeVTK_3Dcut(const char* name, real time, int cycle); 
	void writeVTK_3Dcut(const char* name, real time, int cycle, idx ox, idx oy, idx oz, idx lx, idx ly, idx lz, idx step); 

//	void probe1D() { writeVTKs_1D(); } // replaced by writeVTKs_1D
	void writeVTKs_1D();

	template < typename... ARGS >
	void add1Dcut(real fromx, real fromy, real fromz, real tox, real toy, real toz, const char* fmt, ARGS... args);
	template < typename... ARGS >
	void add1Dcut_X(real y, real z, const char* fmt, ARGS... args);
	template < typename... ARGS >
	void add1Dcut_Y(real x, real z, const char* fmt, ARGS... args);
	template < typename... ARGS >
	void add1Dcut_Z(real x, real y, const char* fmt, ARGS... args);
	void write1Dcut(real fromx, real fromy, real fromz, real tox, real toy, real toz, const char * desc);
	void write1Dcut_X(idx y, idx z, const char * desc);
	void write1Dcut_Y(idx y, idx z, const char * desc);
	void write1Dcut_Z(idx x, idx y, const char * desc);

	int verbosity=1;
	char id[FILENAME_CHARS] = {'d','e','f','a','u','l','t','\0' };		    // simulation id = default this is for 4.9 gcc
//	char id[FILENAME_CHARS] = "default"; // gcc 5.0 + allows that

	virtual bool outputData(int index, int dof, char *desc, idx x, idx y, idx z, real &value, int &dofs) { return false; }
	
	bool fileExists(const char*filename);
	bool getPNGdimensions(const char * filename, int &w, int &h);

	bool projectPNG_X(const char * filename, idx x0, int rotate=0, bool mirror=false, bool flip=false);
	bool projectPNG_Y(const char * filename, idx y0, int rotate=0, bool mirror=false, bool flip=false);
	bool projectPNG_Z(const char * filename, idx z0, int rotate=0, bool mirror=false, bool flip=false);

	virtual void setupBoundaries() { }
	virtual void updateKernelVelocities(T_LBM &local_lbm) { } // setup current velocity profile for the Kernel

	template < typename... ARGS >
	void mark(const char* fmt, ARGS... args);
	bool isMark();
	
	void flagCreate(const char*flagname);
	void flagDelete(const char*flagname);
	bool flagExists(const char*flagname);


	template < typename... ARGS >
	void setid(const char* fmt, ARGS... args);

	template < typename... ARGS >
	void log(const char* fmt, ARGS... args);

	// sim lbm
	void reset();
	virtual void resetLattice(real irho, real ivx, real ivy, real ivz);
	void setEqLat(idx x, idx y, idx z, real irho, real ivx, real ivy, real ivz); // prescribe rho,vx,vy,vz at a given point into "lat" array
	void setEqLat(dreal*lat, idx x, idx y, idx z, real irho, real ivx, real ivy, real ivz); // prescribe rho,vx,vy,vz at a given point into "lat" array

	// save & load state
	template < typename VARTYPE >
	int saveloadBinaryData(int direction, const char*subdirname, const char*filename, VARTYPE*data, idx length);

	template< typename... ARGS >
	int saveLoadTextData(int direction, const char*subdirname, const char*filename, ARGS&... args);
	
	// old version
	//template < typename... ARGS >
	//int saveLoadTextData(int direction, const char*subdirname, const char*filename, const char*fmt, ARGS&... args);

	virtual void saveAndLoadState(int direction, const char*subdirname);
//	void loadInitState(int direction, const char*subdirname);

	// called periodically thru cnt[SAVESTATE]
	bool check_savestate_flag=true; 	// false = output savestate every cnt[SAVESTATE].period, true = output savestate only if "savestate" file exists
	bool delete_savestate_flag=true;	// true = delete "savestate" flag file after savestate is completed
	virtual void saveState(bool forced=false);	
	virtual void loadState(bool forced=false);
	
	// JK magic starts from here
	std::string getFmt(short)               { return "%hd"; }
	std::string getFmt(int)                 { return "%d"; }
	std::string getFmt(long)                { return "%ld"; }
	std::string getFmt(long long)           { return "%lld"; }
	std::string getFmt(unsigned short)      { return "%hu"; }
	std::string getFmt(unsigned int)        { return "%u"; }
	std::string getFmt(unsigned long)       { return "%lu"; }
	std::string getFmt(unsigned long long)  { return "%llu"; }
	std::string getFmt(float)               { return "%e"; }
	std::string getFmt(double)              { return "%le"; }

	template< typename ARG0 >
	std::string getSaveLoadFmt(ARG0 arg0)	{ return getFmt(arg0) + "\n"; }

	template< typename ARG0, typename... ARGS >
	std::string getSaveLoadFmt(ARG0 arg0, ARGS... args) {	return getFmt(arg0) + "\n" + getSaveLoadFmt(args...); }
	// JK magic ends here

	timespec t_init;
	long wallTime=-1; //wallTime in seconds, use negative value to disable wall time check
	bool wallTimeReached();
    
	// vopicarna
	int directionIndex(int i, int j, int k);
    
	// constructors
//	void StateCommonConstructor();
	State(idx iX, idx iY, idx iZ, real iphysViscosity, real iphysDl, real iphysDt) : lbm(iX, iY, iZ, iphysViscosity, iphysDl, iphysDt) {};
};

#include "state.hpp"

#endif
