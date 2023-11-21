#include "lbm_common.h"

template <
	typename TRAITS,
	typename LBM_EQ=LBM_EQ_DEFAULT<TRAITS>
>
struct LBM_MRT : LBM_COMMON< TRAITS, LBM_EQ >
{
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	static constexpr const char* id = "MRT_LES";
	CUDA_HOSTDEV static void collision(KernelStruct<dreal> &KS)
	{
		dreal Pxx = KS.f[mmm]+KS.f[mmz]+KS.f[mmp]+KS.f[mzm]+KS.f[mzz]+KS.f[mzp]+KS.f[mpm]+KS.f[mpz]+KS.f[mpp]+KS.f[pmm]+KS.f[pmz]+KS.f[pmp]+KS.f[pzm]+KS.f[pzz]+KS.f[pzp]+KS.f[ppm]+KS.f[ppz]+KS.f[ppp]; //Second order moment Pi
		dreal Pyy = KS.f[mmm]+KS.f[mmz]+KS.f[mmp]+KS.f[mpm]+KS.f[mpz]+KS.f[mpp]+KS.f[zmm]+KS.f[zmz]+KS.f[zmp]+KS.f[zpm]+KS.f[zpz]+KS.f[zpp]+KS.f[pmm]+KS.f[pmz]+KS.f[pmp]+KS.f[ppm]+KS.f[ppz]+KS.f[ppp];
		dreal Pzz = KS.f[mmm]+KS.f[mmp]+KS.f[mzm]+KS.f[mzp]+KS.f[mpm]+KS.f[mpp]+KS.f[zmm]+KS.f[zmp]+KS.f[zzm]+KS.f[zzp]+KS.f[zpm]+KS.f[zpp]+KS.f[pmm]+KS.f[pmp]+KS.f[pzm]+KS.f[pzp]+KS.f[ppm]+KS.f[ppp];
		dreal Pxy = KS.f[mmm]+KS.f[mmz]+KS.f[mmp]-KS.f[mpm]-KS.f[mpz]-KS.f[mpp]-KS.f[pmm]-KS.f[pmz]-KS.f[pmp]+KS.f[ppm]+KS.f[ppz]+KS.f[ppp];
		dreal Pxz = KS.f[mmm]-KS.f[mmp]+KS.f[mzm]-KS.f[mzp]+KS.f[mpm]-KS.f[mpp]-KS.f[pmm]+KS.f[pmp]-KS.f[pzm]+KS.f[pzp]-KS.f[ppm]+KS.f[ppp];
		dreal Pyz = KS.f[mmm]-KS.f[mmp]-KS.f[mpm]+KS.f[mpp]+KS.f[zmm]-KS.f[zmp]-KS.f[zpm]+KS.f[zpp]+KS.f[pmm]-KS.f[pmp]-KS.f[ppm]+KS.f[ppp];
		const dreal Pnxx = Pxx - KS.rho*(n1o3 + KS.vx*KS.vx); //Non-equilibrium part of Pi
		const dreal Pnyy = Pyy - KS.rho*(n1o3 + KS.vy*KS.vy);
		const dreal Pnzz = Pzz - KS.rho*(n1o3 + KS.vz*KS.vz);
		const dreal Pnxz = Pxz - KS.rho*KS.vx*KS.vz;
		const dreal Pnxy = Pxy - KS.rho*KS.vx*KS.vy;
		const dreal Pnyz = Pyz - KS.rho*KS.vy*KS.vz;

		//FIXME check all LES sometime ];)
		const dreal Q = no2*(Pnxx*Pnxx+Pnyy*Pnyy+Pnzz*Pnzz+no2*(Pnxy*Pnxy+Pnxz*Pnxz+Pnyz*Pnyz)); //Strain rate tensor magnitude
		const dreal tau = no3*KS.lbmViscosity+n1o2;
		const dreal C = (dreal)0.0342; //.18^2 Smagorinsky constant .18-.19
		const dreal omega = no2/(sqrt(tau*tau+no2*C*no3*no3*sqrt(Q)/KS.rho)+tau); //Smagorinsky LES
//const dreal omega = no1/tau;

		/*
		T	S = (sqrtf(visclb*visclb+(T)18.*C*C*sqrtf(Q))-visclb)/((T)6.*C);
		omega = (T)1./(iCsq*(visclb+C*S)+(T).5);
		tau = (T)1./omega;
		*/
		/*
		tau = iCsq*visclb+(T).5;
		tau = (T).5*(sqrtf(tau*tau+C*sqrtf((T)648.*Q)/(KS.rho+(T)1.))+tau); //Smagorinsky LES
		omega = (T)1./tau;
		*/
		Pxx -= omega*Pnxx; //Relaxation, only Pi is relaxed, rho&u are conserved, rest is equilibrated (P.Dellar 2013(talk) -> a'la Ladd 1994)
		Pyy -= omega*Pnyy;
		Pzz -= omega*Pnzz;
		Pxy -= omega*Pnxy;
		Pxz -= omega*Pnxz;
		Pyz -= omega*Pnyz;

		KS.f[mmm]=n1o216*(KS.rho*(n5o2-n3o2*3+no3*(KS.vx*-no1+KS.vy*-no1+KS.vz*(-no1)))+n9o2*(Pxx*no1+Pyy*no1+Pzz*no1+no2*(Pxy*no1+Pxz*no1+Pyz*no1))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[mmz]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*-no1+KS.vy*-no1+KS.vz*0))+n9o2*(Pxx*no1+Pyy*no1+Pzz*0+no2*(Pxy*no1+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[mmp]=n1o216*(KS.rho*(n5o2-n3o2*3+no3*(KS.vx*-no1+KS.vy*-no1+KS.vz*no1))+n9o2*(Pxx*no1+Pyy*no1+Pzz*no1+no2*(Pxy*no1+Pxz*-no1+Pyz*(-no1)))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[mzm]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*-no1+KS.vy*0+KS.vz*(-no1)))+n9o2*(Pxx*no1+Pyy*0+Pzz*no1+no2*(Pxy*0+Pxz*no1+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[mzz]=n2o27*(KS.rho*(n5o2-n3o2*no1+no3*(KS.vx*-no1+KS.vy*0+KS.vz*0))+n9o2*(Pxx*no1+Pyy*0+Pzz*0+no2*(Pxy*0+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[mzp]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*-no1+KS.vy*0+KS.vz*no1))+n9o2*(Pxx*no1+Pyy*0+Pzz*no1+no2*(Pxy*0+Pxz*-no1+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[mpm]=n1o216*(KS.rho*(n5o2-n3o2*3+no3*(KS.vx*-no1+KS.vy*no1+KS.vz*(-no1)))+n9o2*(Pxx*no1+Pyy*no1+Pzz*no1+no2*(Pxy*-no1+Pxz*no1+Pyz*(-no1)))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[mpz]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*-no1+KS.vy*no1+KS.vz*0))+n9o2*(Pxx*no1+Pyy*no1+Pzz*0+no2*(Pxy*-no1+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[mpp]=n1o216*(KS.rho*(n5o2-n3o2*3+no3*(KS.vx*-no1+KS.vy*no1+KS.vz*no1))+n9o2*(Pxx*no1+Pyy*no1+Pzz*no1+no2*(Pxy*-no1+Pxz*-no1+Pyz*no1))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[zmm]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*0+KS.vy*-no1+KS.vz*(-no1)))+n9o2*(Pxx*0+Pyy*no1+Pzz*no1+no2*(Pxy*0+Pxz*0+Pyz*no1))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[zmz]=n2o27*(KS.rho*(n5o2-n3o2*no1+no3*(KS.vx*0+KS.vy*-no1+KS.vz*0))+n9o2*(Pxx*0+Pyy*no1+Pzz*0+no2*(Pxy*0+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[zmp]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*0+KS.vy*-no1+KS.vz*no1))+n9o2*(Pxx*0+Pyy*no1+Pzz*no1+no2*(Pxy*0+Pxz*0+Pyz*(-no1)))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[zzm]=n2o27*(KS.rho*(n5o2-n3o2*no1+no3*(KS.vx*0+KS.vy*0+KS.vz*(-no1)))+n9o2*(Pxx*0+Pyy*0+Pzz*no1+no2*(Pxy*0+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[zzz]=n8o27*(KS.rho*(n5o2-n3o2*0+no3*(KS.vx*0+KS.vy*0+KS.vz*0))+n9o2*(Pxx*0+Pyy*0+Pzz*0+no2*(Pxy*0+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[zzp]=n2o27*(KS.rho*(n5o2-n3o2*no1+no3*(KS.vx*0+KS.vy*0+KS.vz*no1))+n9o2*(Pxx*0+Pyy*0+Pzz*no1+no2*(Pxy*0+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[zpm]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*0+KS.vy*no1+KS.vz*(-no1)))+n9o2*(Pxx*0+Pyy*no1+Pzz*no1+no2*(Pxy*0+Pxz*0+Pyz*(-no1)))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[zpz]=n2o27*(KS.rho*(n5o2-n3o2*no1+no3*(KS.vx*0+KS.vy*no1+KS.vz*0))+n9o2*(Pxx*0+Pyy*no1+Pzz*0+no2*(Pxy*0+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[zpp]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*0+KS.vy*no1+KS.vz*no1))+n9o2*(Pxx*0+Pyy*no1+Pzz*no1+no2*(Pxy*0+Pxz*0+Pyz*no1))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[pmm]=n1o216*(KS.rho*(n5o2-n3o2*3+no3*(KS.vx*no1+KS.vy*-no1+KS.vz*(-no1)))+n9o2*(Pxx*no1+Pyy*no1+Pzz*no1+no2*(Pxy*-no1+Pxz*-no1+Pyz*no1))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[pmz]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*no1+KS.vy*-no1+KS.vz*0))+n9o2*(Pxx*no1+Pyy*no1+Pzz*0+no2*(Pxy*-no1+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[pmp]=n1o216*(KS.rho*(n5o2-n3o2*3+no3*(KS.vx*no1+KS.vy*-no1+KS.vz*no1))+n9o2*(Pxx*no1+Pyy*no1+Pzz*no1+no2*(Pxy*-no1+Pxz*no1+Pyz*(-no1)))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[pzm]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*no1+KS.vy*0+KS.vz*(-no1)))+n9o2*(Pxx*no1+Pyy*0+Pzz*no1+no2*(Pxy*0+Pxz*-no1+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[pzz]=n2o27*(KS.rho*(n5o2-n3o2*no1+no3*(KS.vx*no1+KS.vy*0+KS.vz*0))+n9o2*(Pxx*no1+Pyy*0+Pzz*0+no2*(Pxy*0+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[pzp]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*no1+KS.vy*0+KS.vz*no1))+n9o2*(Pxx*no1+Pyy*0+Pzz*no1+no2*(Pxy*0+Pxz*no1+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[ppm]=n1o216*(KS.rho*(n5o2-n3o2*3+no3*(KS.vx*no1+KS.vy*no1+KS.vz*(-no1)))+n9o2*(Pxx*no1+Pyy*no1+Pzz*no1+no2*(Pxy*no1+Pxz*-no1+Pyz*(-no1)))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[ppz]=n1o54*(KS.rho*(n5o2-n3o2*2+no3*(KS.vx*no1+KS.vy*no1+KS.vz*0))+n9o2*(Pxx*no1+Pyy*no1+Pzz*0+no2*(Pxy*no1+Pxz*0+Pyz*0))-n3o2*(Pxx+Pyy+Pzz));
		KS.f[ppp]=n1o216*(KS.rho*(n5o2-n3o2*3+no3*(KS.vx*no1+KS.vy*no1+KS.vz*no1))+n9o2*(Pxx*no1+Pyy*no1+Pzz*no1+no2*(Pxy*no1+Pxz*no1+Pyz*no1))-n3o2*(Pxx+Pyy+Pzz));
	}
};
