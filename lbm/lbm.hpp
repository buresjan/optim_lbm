template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC >
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::resetForces(real ifx, real ify, real ifz)
{
        /// Reset forces - This is necessary since '+=' is used afterwards.
        for (idx i=0;i<X*Y*Z;i++)
        {
                hfx(i) = ifx;
                hfy(i) = ify;
                hfz(i) = ifz;
        }
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC >
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::copyForcesToDevice()
{
        #ifdef USE_CUDA
        cudaMemcpy(dfx(),  hfx(), X*Y*Z*sizeof(dreal), cudaMemcpyHostToDevice);
        cudaMemcpy(dfy(),  hfy(), X*Y*Z*sizeof(dreal), cudaMemcpyHostToDevice);
        cudaMemcpy(dfz(),  hfz(), X*Y*Z*sizeof(dreal), cudaMemcpyHostToDevice);
        #endif
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::resetMap(map_t geo_type)
{
	for (idx i=0; i<X*Y*Z; i++) map(i)=geo_type;
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::defineWall(idx x, idx y, idx z, bool value)
{
//	if (x>0 && x<X-1 && y > 0 && y<Y-1 && z > 0 && z<Z-1) wall[pos(x,y,z)] = value;
	if (x>=0 && x<=X-1 && y >= 0 && y<=Y-1 && z>=0 && z<=Z-1) wall[pos(x,y,z)] = value;
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void  LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::copyMapToDevice()
{
	#ifdef USE_CUDA
	if (!use_multiple_gpus)
		cudaMemcpy(dmap, hmap, X*Y*Z*sizeof(map_t), cudaMemcpyHostToDevice);
	else
		cudaMemcpy(dmap, hmap + offset_X*Y*Z, local_X*Y*Z*sizeof(map_t), cudaMemcpyHostToDevice);
	#endif
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void  LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::copyMapToHost()
{
	#ifdef USE_CUDA
	if (!use_multiple_gpus)
		cudaMemcpy(hmap, dmap, X*Y*Z*sizeof(map_t), cudaMemcpyDeviceToHost);
	else
		cudaMemcpy(hmap + offset_X*Y*Z, dmap, local_X*Y*Z*sizeof(map_t), cudaMemcpyDeviceToHost);
	#endif
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::copyDFsToDevice(dreal*source, dreal*target)
{
	#ifdef USE_CUDA
	if (!use_multiple_gpus)
	{
		cudaMemcpy(target, source, 27*X*Y*Z*sizeof(dreal), cudaMemcpyHostToDevice);
	}
	else
	{
		cudaMemcpyAsync(
			target + overlap_left*27*Y*Z,
			source + 27*offset_X*Y*Z,
			27*local_X*Y*Z*sizeof(dreal),
			cudaMemcpyHostToDevice,
			cuda_streams[0]
		);
	}
	#endif
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::copyDFsToDevice()
{
	for (int dfty=0;dfty<DFMAX;dfty++)
		copyDFsToDevice(hfs[dfty], dfs[dfty]);
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::copyDFsToHost(dreal*source, dreal*target)
{
	#ifdef USE_CUDA
	if (!use_multiple_gpus)
		cudaMemcpy(target, source, 27*X*Y*Z*sizeof(dreal), cudaMemcpyDeviceToHost);
	else
	{
		cudaMemcpy(
			target + 27*offset_X*Y*Z,
			source + overlap_left*27*Y*Z, // JK hotfix email 2019.05.01
			27*local_X*Y*Z*sizeof(dreal),
			cudaMemcpyDeviceToHost);
	}
	#endif
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::copyDFsToHost()
{
	for (int dfty=0;dfty<DFMAX;dfty++) 
		copyDFsToHost(dfs[dfty], hfs[dfty]);
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::copyMacroToHost()
{
	#ifdef USE_CUDA
	if (MACRO::N>0)
	{
		if (!use_multiple_gpus)
			cudaMemcpy(hmacro, dmacro, MACRO::N*X*Y*Z*sizeof(dreal), cudaMemcpyDeviceToHost);
		else
		{
			idx macro_overlap_left = (MACRO::use_syncMacro) ? overlap_left : 0;
			idx macro_overlap_right = (MACRO::use_syncMacro) ? overlap_right : 0;
			for (int i=0; i<MACRO::N; i++)
				cudaMemcpyAsync(
					hmacro + i*X*Y*Z + offset_X*Y*Z,
					dmacro + macro_overlap_left*Y*Z + i*(macro_overlap_left + local_X + macro_overlap_right)*Y*Z, // FIXME
					local_X*Y*Z*sizeof(dreal),
					cudaMemcpyDeviceToHost,
					cuda_streams[i % max_cuda_streams]
				);
		}
	}
	#endif
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::copyMacroToDevice()
{
	#ifdef USE_CUDA
	if (MACRO::N>0)
	{
		if (!use_multiple_gpus)
			cudaMemcpy(dmacro, hmacro, MACRO::N*X*Y*Z*sizeof(dreal), cudaMemcpyHostToDevice);
		else
		{
			idx macro_overlap_left = (MACRO::use_syncMacro) ? overlap_left : 0;
			idx macro_overlap_right = (MACRO::use_syncMacro) ? overlap_right : 0;
			for (int i=0; i<MACRO::N; i++)
				cudaMemcpyAsync(
					dmacro + macro_overlap_left*Y*Z + i*(macro_overlap_left + local_X + macro_overlap_right)*Y*Z, // FIXME
					hmacro + i*X*Y*Z + offset_X*Y*Z,
					local_X*Y*Z*sizeof(dreal),
					cudaMemcpyHostToDevice,
					cuda_streams[i % max_cuda_streams]
				);
//			cudaDeviceSynchronize();
		}
	}
	#endif
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::computeCPUMacroFromLat()
{
	// take Lat, compute KS and then CPU_MACRO
	if (CPU_MACRO::N > 0)
	{
		LBM_DATA SD;
		for (int dfty=0;dfty<DFMAX;dfty++) SD.dfs[dfty] = hfs[df_out]; //todo check if correct dfs are used for cpu_macro RF FIXME
		SD.X = X;
		SD.Y = Y;
		SD.Z = Z;
//		SD.offset_X = 0;
		SD.overlap_left=0;
		SD.overlap_right=0;
//		SD.X_plus_overlaps = X;
		
		#pragma omp parallel for schedule(static) collapse(2)
		for (idx x=0;x<X;x++)
		for (idx z=0;z<Z;z++)
		for (idx y=0;y<Y;y++)
		{
			idx gi = POS(x, y, z, X, Y, Z);
			KernelStruct<dreal> KS;
			KS.fx=0;
			KS.fy=0;
			KS.fz=0;
			LBM_TYPE::copyDFcur2KS(SD, KS, x, y, z);
			LBM_TYPE::computeDensityAndVelocity(KS);
			CPU_MACRO::outputMacro(cpumacro, X*Y*Z, KS, gi);
//			if (x==128 && y==23 && z==103)
//			printf("KS: %e %e %e %e vs. cpumacro %e %e %e %e [at %d %d %d]\n", KS.vx, KS.vy, KS.vz, KS.rho, cpumacro[mpos(CPU_MACRO::e_vx,x,y,z)], cpumacro[mpos(CPU_MACRO::e_vy,x,y,z)], cpumacro[mpos(CPU_MACRO::e_vz,x,y,z)],cpumacro[mpos(CPU_MACRO::e_rho,x,y,z)],x,y,z);
                }
//                printf("computeCPUMAcroFromLat done.\n");
        }
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::synchronizeDFsToRight(int gpu_id, int gpu_id_right, LBM& lbm_right)
{
	#ifdef USE_CUDA
	idx offset_src = 27*(overlap_left + local_X - 1)*Y*Z;
	if (overlap_left >=2 || overlap_left>=2) printf("warning: overlap >= 2 unexpected & unhandled in synchronizeDFsToRight!!!\n");
	for (int dfty=0;dfty<DFMAX;dfty++)
		cudaMemcpyPeerAsync(
		    lbm_right.dfs[dfty], 
		    gpu_id_right, 
		    dfs[dfty] + offset_src, 
		    gpu_id, 
		    27*Y*Z*sizeof(dreal),  		// FIXME: * overlap  ---> posunout lbm_right.overlap_left
		    cuda_streams[dfty % max_cuda_streams]
		);
	#endif
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::synchronizeDFsToLeft(int gpu_id, int gpu_id_left, LBM& lbm_left)
{
	#ifdef USE_CUDA
	idx offset_src = 27*overlap_left*Y*Z;
	idx offset_dst = 27*(lbm_left.overlap_left + lbm_left.local_X)*Y*Z;
	if (overlap_left >=2 || overlap_right >=2) printf("warning: overlap >= 2 unexpected & unhandled in synchronizeDFsToLeft!!!\n");
	for (int dfty=0;dfty<DFMAX;dfty++)
		cudaMemcpyPeerAsync(
		    lbm_left.dfs[dfty] + offset_dst, 
		    gpu_id_left, 
		    dfs[dfty] + offset_src, 
		    gpu_id, 
		    27*Y*Z*sizeof(dreal), 		// FIXME: * overlap
		    cuda_streams[(DFMAX+dfty) % max_cuda_streams]
		); //FIXME cuda_streams
	#endif
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::synchronizeMacroToRight(int gpu_id, int gpu_id_right, LBM& lbm_right)
{
	#ifdef USE_CUDA
	idx offset_src = (overlap_left + local_X - 1)*Y*Z;
	if (overlap_left >=2 || overlap_left>=2) printf("warning: overlap >= 2 unexpected & unhandled in synchronizeMacroToRight!!!\n");
	for (int i=0; i<MACRO::N; i++)
		cudaMemcpyPeerAsync(lbm_right.dmacro + i*(lbm_right.overlap_left + lbm_right.local_X + lbm_right.overlap_right)*Y*Z, 
		gpu_id_right, dmacro + i*(overlap_left + local_X + overlap_right)*Y*Z + offset_src, gpu_id, Y*Z*sizeof(dreal), cuda_streams[i%max_cuda_streams]);
	#endif
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::synchronizeMacroToLeft(int gpu_id, int gpu_id_left, LBM& lbm_left)
{
	#ifdef USE_CUDA
	idx offset_src = overlap_left*Y*Z;
	idx offset_dst = (lbm_left.overlap_left + lbm_left.local_X)*Y*Z;
	if (overlap_left >=2 || overlap_left>=2) printf("warning: overlap >= 2 unexpected & unhandled in synchronizeMacroToLeft!!!\n");
	for (int i=0; i<MACRO::N; i++)
		cudaMemcpyPeerAsync(lbm_left.dmacro + i*(lbm_left.overlap_left + lbm_left.local_X + lbm_left.overlap_right)*Y*Z + offset_dst ,
		gpu_id_left, dmacro + i*(overlap_left + local_X + overlap_right)*Y*Z + offset_src, gpu_id, Y*Z*sizeof(dreal), cuda_streams[(MACRO::N + i)%max_cuda_streams]); //FIXME cuda_streams
	#endif
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::clearHostMacro()
{
	if (MACRO::N>0)
	for (idx i = 0; i < X*Y*Z*MACRO::N; i++) hmacro[i]=0;
//	copyMacroToDevice();
//	printf("LBM clear macro done\n");
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::allocateDeviceData()
{
	if (local_X == 0) local_X = X;
//	local_X*Y*Z*sizeof(dreal) = local_X*Y*Z*sizeof(dreal);

	#ifdef USE_CUDA
	cudaMalloc((void**)&dmap, local_X*Y*Z*sizeof(map_t));
	dmacro=0;

	idx macro_overlap_left = (MACRO::use_syncMacro) ? overlap_left : 0;
	idx macro_overlap_right = (MACRO::use_syncMacro) ? overlap_right : 0;
	if (MACRO::N>0) cudaMalloc((void**)&dmacro, MACRO::N*(macro_overlap_left + macro_overlap_right + local_X)*Y*Z*sizeof(dreal)); 
	for (int dfty=0;dfty<DFMAX;dfty++)
		cudaMalloc((void**)&dfs[dfty], 27 * (overlap_left + overlap_right + local_X)*Y*Z*sizeof(dreal));
	#else
	dmap=hmap;
	dmacro=hmacro;
	for (int dfty=0;dfty<DFMAX;dfty++) dfs[dfty] = hfs[dfty];
	#endif
	
//	copyMacroToDevice(); // in core.h

	// initialize data pointers
	for (int dfty=0;dfty<DFMAX;dfty++) data.dfs[dfty] = dfs[dfty];
	data.X = local_X;
	data.overlap_left = overlap_left;
	data.overlap_right = overlap_right;
	data.macro_overlap_left = (MACRO::use_syncMacro) ? overlap_left : 0;
	data.macro_overlap_right = (MACRO::use_syncMacro) ? overlap_right : 0;
	data.Y = Y;
	data.Z = Z;
	data.dmacro = dmacro;
	data.dmap = dmap;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::updateKernelData()
{
	data.lbmViscosity = (dreal)lbmViscosity();

	// rotation
	int i = iterations % DFMAX; 			// i = 0, 1, 2, ... DMAX-1
	
	for (int k=0;k<DFMAX;k++)
	{
		int knew = (k-i)<=0 ? (k-i+DFMAX) % DFMAX : k-i;
//		data.dfs[k] = dfs[knew]; // FIXME !!!
		data.dfs[knew] = dfs[k]; // FIXME !!!
//		printf("updateKernelData:: assigning data.dfs[%d] to dfs[%d]\n",k, knew);
	}
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::LBM(idx iX, idx iY, idx iZ, real iphysViscosity, real iphysDl, real iphysDt)
{
	X = iX;
	Y = iY;
	Z = iZ;

	physDl = iphysDl;
	physDt = iphysDt;

	physCharLength = physDl * (real)Y;

	physViscosity = iphysViscosity;
	physFluidDensity = 1000.0;		// override this to your fluid

	iterations = 0;

	physFinalTime = 1e10;
	terminate=false;

	hmap = (map_t*)calloc(X*Y*Z, sizeof(map_t));
	wall = (bool*)calloc(X*Y*Z, sizeof(bool));
	hmacro = (MACRO::N>0) ? (dreal*)malloc(MACRO::N*X*Y*Z*sizeof(dreal)) : 0;
	clearHostMacro();
	cpumacro = (CPU_MACRO::N>0) ? (dreal*)malloc(CPU_MACRO::N*X*Y*Z*sizeof(dreal)) : 0;
//	lat = (dreal*)malloc(27*size_dreal);
	for (int dfty=0;dfty<DFMAX;dfty++) hfs[dfty] = (dreal*)malloc(27*X*Y*Z*sizeof(dreal));
//	lat = hsf[0];
	for (idx i = 0; i < X*Y*Z; i++) wall[i] = false;

	// CUDA
	dmap = NULL;
	dmacro = NULL;
	for (int dfty=0;dfty<DFMAX;dfty++) dfs[dfty] = NULL;

	block_size = 128;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
LBM<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::~LBM()
{
	if (hmap)	free(hmap);
	if (wall)	free(wall);
	if (hmacro)	free(hmacro);
	if (cpumacro)	free(cpumacro);
	for (int dfty=0;dfty<DFMAX;dfty++) 
		if (hfs[dfty]) free(hfs[dfty]);

	#ifdef USE_CUDA
	if (dmap) cudaFree(dmap);
	if (dmacro) cudaFree(dmacro);
	for (int dfty=0;dfty<DFMAX;dfty++)
		cudaFree(dfs[dfty]);
	#endif
}

