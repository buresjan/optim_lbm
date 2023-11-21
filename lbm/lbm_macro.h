// empty Macro containing required forcing quantities for IBM (see lbm.h -> hfx() etc.)
template < typename TRAITS >
struct MacroBase
{
	using dreal = typename TRAITS::dreal;
	using idx = typename TRAITS::idx;
	
	static const bool use_kernelWorker = false;
	static const bool use_syncMacro = false;

	// compulsory method, void here
	CUDA_HOSTDEV static void zeroForces(dreal*macro, idx size, idx pos) {}
//	CUDA_HOSTDEV static void getMacro(dreal*macro, idx size, KernelStruct<dreal> &KS, idx pos) {}

	template < 
		typename LBM_TYPE, 
		typename STREAMING, 
		typename LBM_DATA,
		typename LBM_BC
	>
	CUDA_HOSTDEV static void kernelWorker(LBM_DATA &SD, idx x, idx y, idx z) {}

	CUDA_HOSTDEV static dreal* fx(dreal*macro, idx size) { return macro; } // empty
	CUDA_HOSTDEV static dreal* fy(dreal*macro, idx size) { return macro; } // empty
	CUDA_HOSTDEV static dreal* fz(dreal*macro, idx size) { return macro; } // empty
};

// empty Macro containing required empty rho, vx, vy, vz quantities for IBM (see lbm.h -> hfx() etc.)
template < typename TRAITS >
struct MacroBaseAll : MacroBase < TRAITS >
{
	using dreal = typename TRAITS::dreal;
	using idx = typename TRAITS::idx;

	CUDA_HOSTDEV static dreal* rho(dreal*macro, idx size) { return macro; } // empty
	CUDA_HOSTDEV static dreal* vx(dreal*macro, idx size) { return macro; } // empty
	CUDA_HOSTDEV static dreal* vy(dreal*macro, idx size) { return macro; } // empty
	CUDA_HOSTDEV static dreal* vz(dreal*macro, idx size) { return macro; } // empty
};


template < typename TRAITS >
struct MacroDefault : MacroBase< TRAITS >
{
	using dreal = typename TRAITS::dreal;
	using idx = typename TRAITS::idx;

	enum { e_rho, e_vx, e_vy, e_vz, N };
	CUDA_HOSTDEV static void outputMacro(dreal*macro, idx size, KernelStruct<dreal> &KS, idx pos)
	{
		rho(macro,size)[pos] = KS.rho;
		vx(macro,size)[pos] = KS.vx;
		vy(macro,size)[pos] = KS.vy;
		vz(macro,size)[pos] = KS.vz;
	}
	
	template < typename LBM_DATA >
	CUDA_HOSTDEV static void copyQuantities(LBM_DATA &SD, KernelStruct<dreal> &KS, idx i)
	{
		KS.lbmViscosity = SD.lbmViscosity;
		KS.fx=SD.fx;
		KS.fy=SD.fy;
		KS.fz=SD.fz;
	}

	CUDA_HOSTDEV static dreal* rho(dreal*macro, idx size) { return macro+e_rho*size; }
	CUDA_HOSTDEV static dreal* vx(dreal*macro, idx size) { return macro+e_vx*size; }
	CUDA_HOSTDEV static dreal* vy(dreal*macro, idx size) { return macro+e_vy*size; }
	CUDA_HOSTDEV static dreal* vz(dreal*macro, idx size) { return macro+e_vz*size; }
};


template < typename TRAITS >
struct MacroVoid : MacroBase < TRAITS >
{
	using dreal = typename TRAITS::dreal;
	using idx = typename TRAITS::idx;

	static const int N=0;
	CUDA_HOSTDEV static void outputMacro(dreal*macro, idx size, KernelStruct<dreal> &KS, idx pos)
	{
	}

	template < typename LBM_DATA >
	CUDA_HOSTDEV static void copyQuantities(LBM_DATA &SD, KernelStruct<dreal> &KS, idx i)
	{
	}
};

