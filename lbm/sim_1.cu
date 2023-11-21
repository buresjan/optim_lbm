#include "core.h"

// 3D test domain
template < 
	typename LBM_TYPE,
	typename MACRO = MacroDefault< typename LBM_TYPE::T_TRAITS >,
	typename CPU_MACRO = MacroVoid< typename LBM_TYPE::T_TRAITS >,
	typename LBM_DATA = LBM_Data_ConstInflow< typename LBM_TYPE::T_TRAITS >,
	typename LBM_BC = LBM_BC_All< typename LBM_TYPE::T_TRAITS >
>
struct StateLocal : State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>
{
	using State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::lbm;
	using State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::vtk_helper;
	
	using idx = typename LBM_TYPE::T_TRAITS::idx;
	using real = typename LBM_TYPE::T_TRAITS::real;
	using dreal = typename LBM_TYPE::T_TRAITS::dreal;

	void generateObject()
    /*
     generates entire object from file fname
    */
    {
        std::ifstream infile("mesh/output");
        int x, y, z;
        real type;
        while (infile >> x >> y >> z >> type)
        {
            if (type == no1)
            {
                lbm.map(x, y) = LBM_BC::GEO_WALL;
            }
        }
    }

	virtual void setupBoundaries()
	{
	    generateObject();

        for (idx y = 0; y <= lbm.Y; y++)
            for (idx x = 0; x <= lbm.X; x++)
            {
                if (lbm.map(x, y, 0) == LBM_BC::GEO_FLUID)
                    lbm.map(x, y, 0) = LBM_BC::GEO_INFLOW;
                if (lbm.map(x, y, lbm.Z) == LBM_BC::GEO_FLUID)
                    lbm.map(x, y, 0) = LBM_BC::GEO_OUTFLOW_TOP;
        }
	}

	virtual bool outputData(int index, int dof, char *desc, idx x, idx y, idx z, real &value, int &dofs)
	{
		int k=0;
		if (index==k++) return vtk_helper("lbm_density", lbm.macro(MACRO::e_rho,x,y,z), 1, desc, value, dofs);
		if (index==k++)
		{
			switch (dof)
			{
				case 0: return vtk_helper("velocity", lbm.macro(MACRO::e_vx,x,y,z), 3, desc, value, dofs);
				case 1: return vtk_helper("velocity", lbm.macro(MACRO::e_vy,x,y,z), 3, desc, value, dofs);
				case 2: return vtk_helper("velocity", lbm.macro(MACRO::e_vz,x,y,z), 3, desc, value, dofs);
			}
		}
		return false;
	}
	
	StateLocal(int iX, int iY, int iZ, real iphysViscosity, real iphysVelocity, real iphysDl, real iphysDt) 
		: State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>(iX, iY, iZ, iphysViscosity, iphysDl, iphysDt)
	{
		lbm.data.inflow_rho = no1;
		lbm.data.inflow_vx = 0;
		lbm.data.inflow_vy = 0;
		lbm.data.inflow_vz = lbm.phys2lbmVelocity(iphysVelocity);
	}
};

template < typename LBM_TYPE >
int sim01_test(int RESOLUTION = 2)
{
	using real = typename LBM_TYPE::T_TRAITS::real;

	int block_size=32;
	int X = 128*RESOLUTION;// width in pixels --- product of 128.
	int Y = 82;// height in pixels --- top and bottom walls 1px
	int Z = 82;// height in pixels --- top and bottom walls 1px
	real LBM_VISCOSITY = 0.001;//1.0/6.0; /// GIVEN: optimal is 1/6
	real PHYS_HEIGHT = 0.1; // [m] domain height (physical)
	real PHYS_VISCOSITY = 0.001;// [m^2/s] fluid viscosity .... blood?
	real PHYS_VELOCITY = 1.0; // this is only average velocity .... will be multiplied by 9/4 to get the correct value Um from Schafer Turek
	real PHYS_DL = PHYS_HEIGHT/((real)Y-2);
	real PHYS_DT = LBM_VISCOSITY / PHYS_VISCOSITY*PHYS_DL*PHYS_DL;//PHYS_HEIGHT/(real)LBM_HEIGHT;
	
	StateLocal< LBM_TYPE > state(X, Y, Z, PHYS_VISCOSITY, PHYS_VELOCITY, PHYS_DL, PHYS_DT);
	state.setid("sim01");
	state.lbm.block_size = 32;
	state.lbm.physCharLength = 0.1; // [m]
	state.lbm.use_multiple_gpus = false;
	state.lbm.physFinalTime = 10;//2520.0;
//	state.printIter = 100;
//	state.printIter = 100;
//	state.vtk.writePeriod = 0.01;
	state.cnt[SAVESTATE].period = 0.01; // test
	state.cnt[PRINT].period = 0.01;
	state.cnt[VTK2D].period = 0.1;
	state.cnt[VTK3D].period = 0.1;
	state.cnt[VTK3DCUT].period = 0.01;
	state.wallTime=-10;
	
	state.add2Dcut_Y(Y/2,"cut_Y");


	state.add3Dcut(0,0,0,X,Y,Z,10,"box10");

	execute(state);
	
	return 0;
}

template < typename TRAITS=TraitsSP >
void run()
{
	sim01_test<LBM_CUM< TRAITS > >();
//	sim01_test<LBM_CUM< TRAITS, LBM_EQ_INV_CUM<TRAITS> > >();
}

int Main(int argc, char **argv)
{
	run();
	return 0;
}
