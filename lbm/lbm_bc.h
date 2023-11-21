#ifndef _LBM_BC_H_
#define _LBM_BC_H_

template< typename TRAITS >
struct LBM_BC_All
{
	using map_t = typename TRAITS::map_t;
	using idx = typename TRAITS::idx;
	using dreal = typename TRAITS::dreal;

	enum GEO : map_t {
		GEO_FLUID, 		// compulsory
		GEO_WALL, 		// compulsory
		GEO_INFLOW, 
		GEO_INFLOW_FREE_RHO, 
		GEO_OUTFLOW_RIGHT, 
		GEO_OUTFLOW_RIGHT_INTERP, 
		GEO_OUTFLOW_LEFT, 
		GEO_OUTFLOW_TOP, 
		GEO_OUTFLOW_EQ, 
		GEO_OUTFLOW_RIGHT_EQ, 
		GEO_PERIODIC, 
		GEO_NOTHING,
		GEO_SYM_TOP,
		GEO_SYM_BOTTOM,
		GEO_SYM_BACK,
		GEO_SYM_FRONT,
		GEO_SYM_LEFT,
		GEO_SYM_RIGHT,
		GEO_INFLOW_TOP
	};

	CUDA_HOSTDEV static bool isPeriodic(map_t mapgi) 
	{ 
		return (mapgi==GEO_PERIODIC); 
	}

	CUDA_HOSTDEV static bool isFluid(map_t mapgi) 
	{ 
		return (mapgi==GEO_FLUID); 
	}

	CUDA_HOSTDEV static bool isStreaming(map_t mapgi) 
	{ 
		return (
			mapgi != GEO_WALL &&
			mapgi != GEO_OUTFLOW_RIGHT &&
			mapgi != GEO_OUTFLOW_RIGHT_INTERP &&
			mapgi != GEO_OUTFLOW_LEFT &&
			mapgi != GEO_INFLOW_TOP
		);
	}
	
	CUDA_HOSTDEV static bool isComputeDensityAndVelocity(map_t mapgi) 
	{ 
		return (
			mapgi != GEO_INFLOW && 
			mapgi != GEO_INFLOW_FREE_RHO && 
			mapgi != GEO_WALL &&
			mapgi != GEO_SYM_TOP &&
			mapgi != GEO_SYM_FRONT &&
			mapgi != GEO_SYM_BOTTOM &&
			mapgi != GEO_SYM_BACK && 
			mapgi != GEO_OUTFLOW_RIGHT &&
			mapgi != GEO_OUTFLOW_RIGHT_INTERP &&
			mapgi != GEO_OUTFLOW_LEFT &&
			mapgi != GEO_INFLOW_TOP
		);
	}
	
	template< typename L, typename S, typename LBM_DATA >
	CUDA_HOSTDEV static bool BC(LBM_DATA &SD, KernelStruct<dreal> &KS, map_t mapgi, idx xm, idx x, idx xp, idx ym, idx y, idx yp, idx zm, idx z, idx zp)
	{
		switch (mapgi)
		{
		case GEO_OUTFLOW_RIGHT:
			S::streaming(SD,KS,xm,x,xm,ym,y,yp,zm,z,zp);
			L::computeDensityAndVelocity(KS);
			KS.rho=no1;
			L::collision(KS);
			break;
		case GEO_OUTFLOW_RIGHT_INTERP:
			// streaming: interpolation from Geier - CuLBM (2015)
			// NOTE: velocity is neglected (for the case velocity << speed of sound)
			// TODO: move to lbm_streaming_*.h, implement it for esotwist as well
			#define SpeedOfSound 0.5773502691896257
			KS.f[mmm] = SpeedOfSound*SD.df(df_cur,mmm,xm,yp,zp) + (1.0-SpeedOfSound)*SD.df(df_cur,mmm,x,yp,zp);
			KS.f[mmz] = SpeedOfSound*SD.df(df_cur,mmz,xm,yp, z) + (1.0-SpeedOfSound)*SD.df(df_cur,mmz,x,yp, z);
			KS.f[mmp] = SpeedOfSound*SD.df(df_cur,mmp,xm,yp,zm) + (1.0-SpeedOfSound)*SD.df(df_cur,mmp,x,yp,zm);
			KS.f[mzm] = SpeedOfSound*SD.df(df_cur,mzm,xm, y,zp) + (1.0-SpeedOfSound)*SD.df(df_cur,mzm,x, y,zp);
			KS.f[mzz] = SpeedOfSound*SD.df(df_cur,mzz,xm, y, z) + (1.0-SpeedOfSound)*SD.df(df_cur,mzz,x, y, z);
			KS.f[mzp] = SpeedOfSound*SD.df(df_cur,mzp,xm, y,zm) + (1.0-SpeedOfSound)*SD.df(df_cur,mzp,x, y,zm);
			KS.f[mpm] = SpeedOfSound*SD.df(df_cur,mpm,xm,ym,zp) + (1.0-SpeedOfSound)*SD.df(df_cur,mpm,x,ym,zp);
			KS.f[mpz] = SpeedOfSound*SD.df(df_cur,mpz,xm,ym, z) + (1.0-SpeedOfSound)*SD.df(df_cur,mpz,x,ym, z);
			KS.f[mpp] = SpeedOfSound*SD.df(df_cur,mpp,xm,ym,zm) + (1.0-SpeedOfSound)*SD.df(df_cur,mpp,x,ym,zm);
			#undef SpeedOfSound
			KS.f[zmm] = SD.df(df_cur,zmm, x,yp,zp);
			KS.f[zmz] = SD.df(df_cur,zmz, x,yp, z);
			KS.f[zmp] = SD.df(df_cur,zmp, x,yp,zm);
			KS.f[zzm] = SD.df(df_cur,zzm, x, y,zp);
			KS.f[zzz] = SD.df(df_cur,zzz, x, y, z);
			KS.f[zzp] = SD.df(df_cur,zzp, x, y,zm);
			KS.f[zpm] = SD.df(df_cur,zpm, x,ym,zp);
			KS.f[zpz] = SD.df(df_cur,zpz, x,ym, z);
			KS.f[zpp] = SD.df(df_cur,zpp, x,ym,zm);
			KS.f[pmm] = SD.df(df_cur,pmm,xm,yp,zp);
			KS.f[pmz] = SD.df(df_cur,pmz,xm,yp, z);
			KS.f[pmp] = SD.df(df_cur,pmp,xm,yp,zm);
			KS.f[pzm] = SD.df(df_cur,pzm,xm, y,zp);
			KS.f[pzz] = SD.df(df_cur,pzz,xm, y, z);
			KS.f[pzp] = SD.df(df_cur,pzp,xm, y,zm);
			KS.f[ppm] = SD.df(df_cur,ppm,xm,ym,zp);
			KS.f[ppz] = SD.df(df_cur,ppz,xm,ym, z);
			KS.f[ppp] = SD.df(df_cur,ppp,xm,ym,zm);

			L::computeDensityAndVelocity(KS);
			KS.rho=no1;
			L::collision(KS);
			break;
		case GEO_OUTFLOW_LEFT:
			S::streaming(SD,KS,xp,x,xp,ym,y,yp,zm,z,zp);
			L::computeDensityAndVelocity(KS);
			KS.rho=no1;
			L::collision(KS);
			break;
		case GEO_INFLOW:
			SD.inflow(KS,x,y,z);
			L::setEquilibrium(KS);
			break;
		case GEO_INFLOW_FREE_RHO:
			S::streaming(SD,KS,xm,x,xp,ym,y,yp,zm,z,zp);
			L::computeDensityAndVelocity(KS);
			SD.inflow(KS,x,y,z);
			L::setEquilibrium(KS);
			break;
//		case GEO_INFLOW_CONST:
			// use initial distributions == do nothing
//			L::copyDFtoKS(SD,KS,x,y,z);
//			L::computeDensityAndVelocity(KS);
//			break;
		case GEO_OUTFLOW_EQ: // FIXME
//			S::streaming(SD,KS,xm,x,xp,ym,y,yp,zm,z,zp);
			L::computeDensityAndVelocity(KS);
			KS.rho=no1;
			L::setEquilibrium(KS);
			break;
		case GEO_OUTFLOW_RIGHT_EQ:
			S::streaming(SD,KS,xm-1,xm,x,ym,y,yp,zm,z,zp);
			L::computeDensityAndVelocity(KS);
			KS.rho=no1;
			L::setEquilibrium(KS);
			break;
		case GEO_WALL:
			S::streamingBounceBack(SD,KS,xm,x,xp,ym,y,yp,zm,z,zp);
			// nema zadny vliv na vypocet, jen pro output
			KS.rho=no1; 
			KS.vx=(dreal)0;
			KS.vy=(dreal)0;
			KS.vz=(dreal)0;
			break;
		case GEO_INFLOW_TOP: // PE GEO_WALL_W_VEL
			S::streaming(SD,KS,xm,x,xp,ym,y,y,zm,z,zp);
			SD.inflow(KS,x,y,z);
			L::collision(KS);

			dreal t;
			t = KS.f[ppp];
			KS.f[ppp] = KS.f[mmm] - no6*KS.rho*n1o216*(-KS.vx - KS.vy - KS.vz);
			KS.f[mmm] = t - no6*KS.rho*n1o216*(KS.vx + KS.vy + KS.vz);

			t = KS.f[ppz];
			KS.f[ppz] = KS.f[mmz] - no6*KS.rho*n1o54*(-KS.vx - KS.vy);
			KS.f[mmz] = t - no6*KS.rho*n1o54*(KS.vx + KS.vy);

			t = KS.f[ppm];
			KS.f[ppm] = KS.f[mmp] - no6*KS.rho*n1o216*(-KS.vx - KS.vy + KS.vz);
			KS.f[mmp] = t - no6*KS.rho*n1o216*(KS.vx + KS.vy - KS.vz);

			t = KS.f[pzp];
			KS.f[pzp] = KS.f[mzm] - no6*KS.rho*n1o54*(-KS.vx - KS.vz);
			KS.f[mzm] = t - no6*KS.rho*n1o54*(KS.vx + KS.vz);

			t = KS.f[pzz];
			KS.f[pzz] = KS.f[mzz] - no6*KS.rho*n2o27*(-KS.vx);
			KS.f[mzz] = t - no6*KS.rho*n2o27*(KS.vx);

			t = KS.f[pzm];
			KS.f[pzm] = KS.f[mzp] - no6*KS.rho*n1o54*(-KS.vx + KS.vz);
			KS.f[mzp] = t - no6*KS.rho*n1o54*(KS.vx - KS.vz);

			t = KS.f[pmp];
			KS.f[pmp] = KS.f[mpm] - no6*KS.rho*n1o216*(-KS.vx + KS.vy - KS.vz);
			KS.f[mpm] = t - no6*KS.rho*n1o216*(KS.vx - KS.vy + KS.vz);

			t = KS.f[pmz];
			KS.f[pmz] = KS.f[mpz] - no6*KS.rho*n1o54*(-KS.vx + KS.vy);
			KS.f[mpz] = t - no6*KS.rho*n1o54*(KS.vx - KS.vy);

			t = KS.f[pmm];
			KS.f[pmm] = KS.f[mpp] - no6*KS.rho*n1o216*(-KS.vx + KS.vy + KS.vz);
			KS.f[mpp] = t - no6*KS.rho*n1o216*(KS.vx - KS.vy - KS.vz);

			t = KS.f[zpp];
			KS.f[zpp] = KS.f[zmm] - no6*KS.rho*n1o54*(- KS.vy - KS.vz);
			KS.f[zmm] = t - no6*KS.rho*n1o54*(KS.vy + KS.vz);

			t = KS.f[zpz];
			KS.f[zpz] = KS.f[zmz] - no6*KS.rho*n2o27*(- KS.vy);
			KS.f[zmz] = t - no6*KS.rho*n2o27*(KS.vy);

			t = KS.f[zpm];
			KS.f[zpm] = KS.f[zmp] - no6*KS.rho*n1o54*(- KS.vy + KS.vz);
			KS.f[zmp] = t - no6*KS.rho*n1o54*(KS.vy - KS.vz);

			t = KS.f[zzp];
			KS.f[zzp] = KS.f[zzm] - no6*KS.rho*n2o27*(- KS.vz);
			KS.f[zzm] = t - no6*KS.rho*n2o27*(KS.vz);
			break;
		case GEO_SYM_BACK:
			KS.f[mpm] = KS.f[mmm];
			KS.f[mpz] = KS.f[mmz];
			KS.f[mpp] = KS.f[mmp];
			KS.f[zpm] = KS.f[zmm];
			KS.f[zpz] = KS.f[zmz];
			KS.f[zpp] = KS.f[zmp];
			KS.f[ppm] = KS.f[pmm];
			KS.f[ppz] = KS.f[pmz];
			KS.f[ppp] = KS.f[pmp];
			L::computeDensityAndVelocity(KS);
			break;
		case GEO_SYM_FRONT:
			KS.f[mmm] = KS.f[mpm];
			KS.f[mmz] = KS.f[mpz];
			KS.f[mmp] = KS.f[mpp];
			KS.f[zmm] = KS.f[zpm];
			KS.f[zmz] = KS.f[zpz];
			KS.f[zmp] = KS.f[zpp];
			KS.f[pmm] = KS.f[ppm];
			KS.f[pmz] = KS.f[ppz];
			KS.f[pmp] = KS.f[ppp];
			L::computeDensityAndVelocity(KS);
			break;
		case GEO_SYM_BOTTOM:
			KS.f[mmp] = KS.f[mmm];
			KS.f[mzp] = KS.f[mzm];
			KS.f[mpp] = KS.f[mpm];
			KS.f[zmp] = KS.f[zmm];
			KS.f[zzp] = KS.f[zzm];
			KS.f[zpp] = KS.f[zpm];
			KS.f[pmp] = KS.f[pmm];
			KS.f[pzp] = KS.f[pzm];
			KS.f[ppp] = KS.f[ppm];
			L::computeDensityAndVelocity(KS);
			break;
		case GEO_SYM_TOP:
			KS.f[mmm] = KS.f[mmp];
			KS.f[mzm] = KS.f[mzp];
			KS.f[mpm] = KS.f[mpp];
			KS.f[zmm] = KS.f[zmp];
			KS.f[zzm] = KS.f[zzp];
			KS.f[zpm] = KS.f[zpp];
			KS.f[pmm] = KS.f[pmp];
			KS.f[pzm] = KS.f[pzp];
			KS.f[ppm] = KS.f[ppp];
			L::computeDensityAndVelocity(KS);
			break;
		case GEO_SYM_LEFT:
			KS.f[pmm] = KS.f[mmm];
			KS.f[pmz] = KS.f[mmz];
			KS.f[pmp] = KS.f[mmp];
			KS.f[pzm] = KS.f[mzm];
			KS.f[pzz] = KS.f[mzz];
			KS.f[pzp] = KS.f[mzp];
			KS.f[ppm] = KS.f[mpm];
			KS.f[ppz] = KS.f[mpz];
			KS.f[ppp] = KS.f[mpp];
			L::computeDensityAndVelocity(KS);
			break;
		case GEO_SYM_RIGHT:
			KS.f[mmm] = KS.f[pmm];
			KS.f[mmz] = KS.f[pmz];
			KS.f[mmp] = KS.f[pmp];
			KS.f[mzm] = KS.f[pzm];
			KS.f[mzz] = KS.f[pzz];
			KS.f[mzp] = KS.f[pzp];
			KS.f[mpm] = KS.f[ppm];
			KS.f[mpz] = KS.f[ppz];
			KS.f[mpp] = KS.f[ppp];
			L::computeDensityAndVelocity(KS);
			break;
		default:
			return false;
		}
		return true;
	}
};

#endif
