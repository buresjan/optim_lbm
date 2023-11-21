#ifndef _LBM_COMMON_H_
#define _LBM_COMMON_H_

template <
	typename TRAITS,
	typename LBM_EQ
>
struct LBM_COMMON
{
	using T_TRAITS = TRAITS;
	using T_LBM_EQ = LBM_EQ;

	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	CUDA_HOSTDEV static void computeDensityAndVelocity(KernelStruct<dreal> &KS)
	{ 
		#ifdef USE_HIGH_PRECISION_RHO
		// src: https://en.wikipedia.org/wiki/Kahan_summation_algorithm
		KS.rho = 0;
		dreal c = 0;                 // A running compensation for lost low-order bits.
		for (int i=0;i<27;i++)
		{
			dreal y = KS.f[i] - c;
			dreal t = KS.rho + y;
			c = (t-KS.rho)-y;
			KS.rho = t;
		}
		#else
		// based on Geier 2015 Appendix J
		KS.rho= ((((KS.f[ppp]+KS.f[mmm]) + (KS.f[pmp]+KS.f[mpm])) + ((KS.f[ppm]+KS.f[mmp])+(KS.f[mpp]+KS.f[pmm])))
			+(((KS.f[zpp]+KS.f[zmm]) + (KS.f[zpm]+KS.f[zmp])) + ((KS.f[pzp]+KS.f[mzm])+(KS.f[pzm]+KS.f[mzp])) + ((KS.f[ppz]+KS.f[mmz]) + (KS.f[pmz]+KS.f[mpz])))
			+((KS.f[pzz]+KS.f[mzz]) + (KS.f[zpz]+KS.f[zmz]) + (KS.f[zzp]+KS.f[zzm]))) + KS.f[zzz];
		#endif
		
		KS.vz=((((KS.f[ppp]-KS.f[mmm])+(KS.f[mpp]-KS.f[pmm]))+((KS.f[pmp]-KS.f[mpm])+(KS.f[mmp]-KS.f[ppm])))+(((KS.f[zpp]-KS.f[zmm])+(KS.f[zmp]-KS.f[zpm]))+((KS.f[pzp]-KS.f[mzm])+(KS.f[mzp]-KS.f[pzm])))+(KS.f[zzp]-KS.f[zzm])+KS.fz*n1o2)/KS.rho;
		KS.vx=((((KS.f[ppp]-KS.f[mmm])+(KS.f[pmp]-KS.f[mpm]))+((KS.f[ppm]-KS.f[mmp])+(KS.f[pmm]-KS.f[mpp])))+(((KS.f[pzp]-KS.f[mzm])+(KS.f[pzm]-KS.f[mzp]))+((KS.f[ppz]-KS.f[mmz])+(KS.f[pmz]-KS.f[mpz])))+(KS.f[pzz]-KS.f[mzz])+KS.fx*n1o2)/KS.rho;
		KS.vy=((((KS.f[ppp]-KS.f[mmm])+(KS.f[ppm]-KS.f[mmp]))+((KS.f[mpp]-KS.f[pmm])+(KS.f[mpm]-KS.f[pmp])))+(((KS.f[ppz]-KS.f[mmz])+(KS.f[mpz]-KS.f[pmz]))+((KS.f[zpp]-KS.f[zmm])+(KS.f[zpm]-KS.f[zmp])))+(KS.f[zpz]-KS.f[zmz])+KS.fy*n1o2)/KS.rho;
	}

	CUDA_HOSTDEV static void setEquilibrium(KernelStruct<dreal> &KS)
	{
		KS.f[mmm] = LBM_EQ::feq_mmm(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[mmz] = LBM_EQ::feq_mmz(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[mmp] = LBM_EQ::feq_mmp(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[mzm] = LBM_EQ::feq_mzm(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[mzz] = LBM_EQ::feq_mzz(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[mzp] = LBM_EQ::feq_mzp(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[mpm] = LBM_EQ::feq_mpm(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[mpz] = LBM_EQ::feq_mpz(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[mpp] = LBM_EQ::feq_mpp(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[zmm] = LBM_EQ::feq_zmm(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[zmz] = LBM_EQ::feq_zmz(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[zmp] = LBM_EQ::feq_zmp(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[zzm] = LBM_EQ::feq_zzm(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[zzz] = LBM_EQ::feq_zzz(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[zzp] = LBM_EQ::feq_zzp(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[zpm] = LBM_EQ::feq_zpm(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[zpz] = LBM_EQ::feq_zpz(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[zpp] = LBM_EQ::feq_zpp(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[pmm] = LBM_EQ::feq_pmm(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[pmz] = LBM_EQ::feq_pmz(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[pmp] = LBM_EQ::feq_pmp(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[pzm] = LBM_EQ::feq_pzm(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[pzz] = LBM_EQ::feq_pzz(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[pzp] = LBM_EQ::feq_pzp(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[ppm] = LBM_EQ::feq_ppm(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[ppz] = LBM_EQ::feq_ppz(KS.rho,KS.vx,KS.vy,KS.vz);
		KS.f[ppp] = LBM_EQ::feq_ppp(KS.rho,KS.vx,KS.vy,KS.vz);
	}
	
	template < typename LBM_DATA > 
	CUDA_HOSTDEV static void copyDFcur2KS(LBM_DATA &SD, KernelStruct<dreal> &KS, idx x, idx y, idx z)
	{
		for (int i=0;i<27;i++) KS.f[i] = SD.df(df_cur,i,x,y,z);
//		for (int i=0;i<27;i++) KS.f[i] = SD.cdf(i,x,y,z);
//		KS.f[i] = SD.cdf[Fxyz(i,x,y,z,SD.X,SD.Y,SD.Z)];
	}

	template < typename LBM_DATA > 
	CUDA_HOSTDEV static void copyKS2DFout(LBM_DATA &SD, KernelStruct<dreal> &KS, idx x, idx y, idx z)
	{
		for (int i=0;i<27;i++) SD.df(df_out,i,x,y,z) = KS.f[i]; 
//		for (int i=0;i<27;i++) KS.f[i] = SD.cdf(i,x,y,z);
//		KS.f[i] = SD.cdf[Fxyz(i,x,y,z,SD.X,SD.Y,SD.Z)];
	}


};

#endif
