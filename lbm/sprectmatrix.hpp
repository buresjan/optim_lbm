template < typename DREAL >
int SpRectMatrix<DREAL>::matrixFormat()
{
// TODO
//	if (Solver == SOLVER_UMFPACK) return CSC_FORMAT;
	return CSR_FORMAT;
}

template < typename DREAL >
int SpRectMatrix<DREAL>::xpos(int i, int j) // i,j, in [0, n-1] !!
{
	switch (matrixFormat())
	{
		case CSC_FORMAT:
			// CSC format:
			// based on Ai and Ap structures get x
			// i ... row index - they are in Ai[ ] 
			// j ... column index
			if (j>=m_nc || j<0) { printf("SpRectMatrix<DREAL>::xpos xpos(%d,%d) : column index out of range! max %d \n",i,j,m_nc); return 0;}
			for (int k=Ap[j];k<Ap[j+1];k++) if (Ai[k]==i) return k; // this is the xpos
		break;
		case CSR_FORMAT:
			// CSR format:
			// based on Ai and Ap structures get x
			// i ... row index
			// j ... column index - they are in Ai[ ] 
			if (i>=m_nr || i<0) { printf("SpRectMatrix<DREAL>::xpos xpos(%d,%d) : row index out of range! max %d \n",i,j,m_nr); return 0; }
			for (int k=Ap[i];k<Ap[i+1];k++) if (Ai[k]==j) return k; // this is the xpos
		break;
	}
	return -1;//CTRL->error("SpRectMatrixxpos","unhandled case i %d j %d  n %d nz %d ",i ,j, m_n, m_nz);
}

template < typename DREAL >
DREAL& SpRectMatrix<DREAL>::get(int i, int j, int verbose)
{
	int xp = xpos(i,j);
	if (xp < 0) return null; 
	return Ax[xp];
}

template < typename DREAL >
int SpRectMatrix<DREAL>::empty()
{
	for (int i=0;i<m_nz;i++) Ax[i] = 0;
	return 1;
}

template < typename DREAL >
int SpRectMatrix<DREAL>::print_nonzeros()
{
	printf("Matrix print_nonzeros\n");
	for (int i=0;i<m_nz;i++) if (Ax[i]!=0) printf(" i %d Ax %e \n", i, Ax[i]);
	return 1;
}

template < typename DREAL >
int SpRectMatrix<DREAL>::print()
{
	DREAL val;
	printf("Matrix print \n");
	for (int i=0;i<m_nr;i++)
	{
//			        printf("row %d : ",i);
//			        printf("\t %d",i);
		for (int j=0;j<m_nc;j++)
		{
			val = get(i,j);
//			if (val!=0) printf("\t %d \t %d \t %1.20e\n",i,j,val);
			if (xpos(i,j)!=-1) printf("\t %d \t %d \t %1.20e\n",i,j,val);
//						if (val==0) printf(" "); else
//						if (val==1) printf(""); else
//						if (val>0) printf("+"); else
//						printf("x");
//						printf("%s ",(val>=0) ? "+" : "", val);
		}
		//		printf(" \n");
	}
	return 1;
}

template < typename DREAL >
int SpRectMatrix<DREAL>::print2()
{
	DREAL val;
	printf("Matrix print \n");
	for (int i=0;i<m_nr;i++)
	{
	        printf("Row: %d ->  ",i);
//			        printf("\t %d",i);
		for (int j=0;j<m_nc;j++)
		{
			val = get(i,j);
//			if (val!=0) printf("\t %d \t %d \t %1.20e\n",i,j,val);
			if (val!=0)
			if (xpos(i,j)!=-1) printf("Col:%d->%1.5e\t ",j,val);
//						if (val==0) printf(" "); else
//						if (val==1) printf(""); else
//						if (val>0) printf("+"); else
//						printf("x");
//						printf("%s ",(val>=0) ? "+" : "", val);
		}
		printf(" \n");
	}
	return 1;
}

template < typename DREAL >
int SpRectMatrix<DREAL>::print2d()
{
	DREAL val;
	printf("Matrix print \n");
	for (int i=0;i<m_nr;i++)
	{
			        printf("row %d : ",i);
			        printf("\t %d",i);
		for (int j=0;j<m_nc;j++)
		{
			val = get(i,j);
			//if (val!=0) printf("\t %d \t %d \t %g\n",i,j,val);
			if (val==0) printf("  "); else
//			if (val==1) printf(" "); else
			if (val>0) printf("++"); else
			printf("--");
//			printf("%s ",(val>=0) ? "+" : " ");
		}
		printf(" \n");
	}
	return 1;
}

template < typename DREAL >
void SpRectMatrix<DREAL>::cusparseInit()
{
    // Get handle to the CUBLAS context 
//    cublasHandle_t cublasHandle = 0;
//    cublasStatus_t cublasStatus;
//    cublasStatus = 
    cublasCreate(&cublasHandle);

//    checkCudaErrors(cublasStatus);

    // Get handle to the CUSPARSE context 
//    cusparseHandle_t cusparseHandle = 0;
//    cusparseStatus_t cusparseStatus;
//    cusparseStatus = 
    cusparseCreate(&cusparseHandle);

//    checkCudaErrors(cusparseStatus);

//    cusparseMatDescr_t descr = 0;
//    cusparseStatus = 
    cusparseCreateMatDescr(&descr);
//    checkCudaErrors(cusparseStatus);

    cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

    cudaMalloc((void **)&dAi, m_nz*sizeof(int));
    cudaMalloc((void **)&dAp, (m_nr+1)*sizeof(int));
    cudaMalloc((void **)&dAx, m_nz*sizeof(DREAL));
    cudaMalloc((void **)&dp, m_nr*sizeof(DREAL));
    cudaMalloc((void **)&dy, m_nr*sizeof(DREAL));

    cudaMemcpy(dAi, Ai, m_nz*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dAp, Ap, (m_nr+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dAx, Ax, m_nz*sizeof(DREAL), cudaMemcpyHostToDevice);
}

/*
void SpRectMatrix<DREAL>::multiply(DREAL *input, DREAL*output, DREAL multiplicator)
{
	for (int i=0;i<m_nr;i++)
	{
		output[i]=0;
		for (int j=Ap[i];j<Ap[i+1];j++) output[i] +=Ax[j] * input[Ai[j]] * multiplicator;
	}
}

void SpRectMatrix<DREAL>::multiply(DREAL *input, real*output, real multiplicator)
{
	for (int i=0;i<m_nr;i++)
	{
		output[i]=0;
		for (int j=Ap[i];j<Ap[i+1];j++) output[i] += (real)Ax[j] * (real)input[Ai[j]] * multiplicator;
	}
}
*/

template < typename DREAL >
void SpRectMatrix<DREAL>::dmultiply(DREAL *input, DREAL*output, DREAL multiplicator)
{
	if (!dAi || !dAp || !dAx) cusparseInit();
	DREAL beta=0;
	cusparseTcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, m_nr, m_nc, m_nz, &multiplicator, descr, dAx, dAp, dAi, input, &beta, output);
}

// d_r = rhs
// d_x = output vector
template < typename DREAL >
int SpRectMatrix<DREAL>::dsolveCG(DREAL *d_r, DREAL *d_x)
{
	if (m_nr!=m_nc)
	{
		printf("SpRectMatrix error: nr %d != nc %d\n", m_nr, m_nc);
		return 0;
	}
	int N=m_nr;
	DREAL alpha = 1.0;
	DREAL alpham1 = -1.0;
	DREAL beta = 0.0;
	DREAL r0 = 0.;
	DREAL a,b,na, r1, dot;
	int k;

	// CG stopping criterion: ||r|| / ||b|| < solverTolerance
	// (note the ||b|| to be consistent with TNL)
	DREAL normb;
	cublasTdot(cublasHandle, N, d_r, 1, d_r, 1, &normb);
	normb = sqrt(normb);

	cusparseTcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, m_nz, &alpha, descr, dAx, dAp, dAi, d_x, &beta, dy);
/*	
	// print dy
	DREAL *w,*q;
	w = new DREAL[N];
	q = new DREAL[N];
	cudaMemcpy(w, dy, N*sizeof(DREAL), cudaMemcpyDeviceToHost);
	cudaMemcpy(q, d_x, N*sizeof(DREAL), cudaMemcpyDeviceToHost);
	for (int i=0;i<N;i++) printf("dy[%d]=%e rhs[%d]=%e\n",i,w[i],i,q[i]);
	delete [] w;
	delete [] q;
	w = new DREAL[N];
	cudaMemcpy(w, dAx, N*sizeof(DREAL), cudaMemcpyDeviceToHost);
	for (int i=0;i<N;i++) printf("dAx[%d]=%e | Ax[%d]=%e\n ",i,w[i],i,Ax[i] );
	delete [] w;
*/	
	
    
	cublasTaxpy(cublasHandle, N, &alpham1, dy, 1, d_r, 1);
	cublasTdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
	
//	printf("r1 = %e\n",r1);

	k = 1;

        while (sqrt(r1) / normb > solverTolerance && k <= solverMaxIter)
//        while (k <= 10)
        {
            if (k > 1)
            {
                b = r1 / r0;
                cublasTscal(cublasHandle, N, &b, dp, 1);
                cublasTaxpy(cublasHandle, N, &alpha, d_r, 1, dp, 1);
            }
            else
            {
                cublasTcopy(cublasHandle, N, d_r, 1, dp, 1);
            }
            
            cusparseTcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, m_nz, &alpha, descr, dAx, dAp, dAi, dp, &beta, dy);
            
            // analyze
            /*
		w = new DREAL[N];
		q = new DREAL[N];
		cudaMemcpy(w, dy, N*sizeof(DREAL), cudaMemcpyDeviceToHost);
		cudaMemcpy(q, dp, N*sizeof(DREAL), cudaMemcpyDeviceToHost);
		for (int i=0;i<N;i++) printf("dy[%d]=%e dp[%d]=%e\n",i,w[i],i,q[i]);
		delete [] w;
		delete [] q;
	int *v;
	v = new int[m_nr+1];
	cudaMemcpy(v, dAp, (m_nr+1)*sizeof(int), cudaMemcpyDeviceToHost);
	for (int i=0;i<=m_nr;i++) printf("dAp[%d]=%d vs Ap[%d]=%d\n",i,v[i], i,Ap[i]);
	delete [] v;
	v = new int[N];
	cudaMemcpy(v, dAi, (N)*sizeof(int), cudaMemcpyDeviceToHost);
	for (int i=0;i<N;i++) printf("dAi[%d]=%d vs Ai[%d]=%d\n",i,v[i], i,Ai[i]);
	delete [] v;
            */
            
            
            cublasTdot(cublasHandle, N, dp, 1, dy, 1, &dot);
            a = r1 / dot;
//            printf("a %e r1 %e dot %e\n",a,r1,dot);

            cublasTaxpy(cublasHandle, N, &a, dp, 1, d_x, 1);
            na = -a;
            cublasTaxpy(cublasHandle, N, &na, dy, 1, d_r, 1);

            r0 = r1;
            cublasTdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
            cudaDeviceSynchronize();
//            printf("iteration = %3d, residual = %e\n", k, sqrt(r1));
            k++;
        }
//        printf("CG solver: end: iteration = %3d, residual = %e\n", k, sqrt(r1));
        return 1;
}

template < typename DREAL >
int SpRectMatrix<DREAL>::dsolve(DREAL *b, DREAL *x)
{
	if (!dAi || !dAp || !dAx) cusparseInit();
	if (m_nr!=m_nc)
	{
		printf("SpRectMatrix error: nr %d != nc %d\n", m_nr, m_nc);
		return 0;
	}
	return dsolveCG(b, x);
}

template < typename DREAL >
void SpRectMatrix<DREAL>::unitMatrixInit(int N)
{
	int *unitAi = new int[N];
	int *unitAp = new int[N+1];
	for (int i=0;i<=N;i++) unitAp[i]=i;
	for (int i=0;i<N;i++) unitAi[i]=i;
	
	cudaMalloc((void **)&dunitAi, N*sizeof(int));
	cudaMalloc((void **)&dunitAp, (N+1)*sizeof(int));
//	cudaMalloc((void **)&dunitV, (N)*sizeof(DREAL));
	
	cudaMemcpy(dunitAi, unitAi, N*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dunitAp, unitAp, (N+1)*sizeof(int), cudaMemcpyHostToDevice);

	delete [] unitAi;
	delete [] unitAp;
}

template < typename DREAL >
void SpRectMatrix<DREAL>::dVectorMultiply(int N, DREAL*a, DREAL *b, DREAL alpha) // b[i] *= a[i]
{
	if (!dAi || !dAp || !dAx) cusparseInit();
	if (!dunitAi || !dunitAp) unitMatrixInit(N);
	DREAL beta=0;
//	cublasTcopy(cublasHandle, N, b, 1, dunitV, 1);
	cusparseTcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, N, &alpha, descr, a, dunitAp, dunitAi, b, &beta, b);
}

template < typename DREAL >
void SpRectMatrix<DREAL>::dVectorMultiplyAndAdd(int N, DREAL*a, DREAL*b, DREAL alpha, DREAL*output)
{
	if (!dAi || !dAp || !dAx) cusparseInit();
	if (!dunitAi || !dunitAp) unitMatrixInit(N);
	DREAL beta=(DREAL)1.0;
	//cublasTcopy(cublasHandle, N, b, 1, dunitV, 1);
	cusparseTcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, N, &alpha, descr, a, dunitAp, dunitAi, b, &beta, output);
}

template < typename DREAL >
void SpRectMatrix<DREAL>::dVectorAddition(int N, DREAL*a, DREAL* b, DREAL alpha) // a[i] += b[i]
{
	cublasTaxpy(cublasHandle, N, &alpha, a, 1, b, 1);
}

template < typename DREAL >
SpRectMatrix<DREAL>::SpRectMatrix()
{
	m_nr=0;
	m_nc=0;
	m_nz=0;
	Ai=0;
	Ap=0;
	Ax=0;
	null=0; // dummy variable that always returns zero (zero entry in the matrix)

	dAi=0;
	dAp=0;
	dAx=0;
	dunitAi=0;
	dunitAp=0;
	dunitV=0;
	dy=0; // dummy unneeded
	dp=0;
	descr=0;
	
	cusparseHandle=0;
	cublasHandle=0;
	
	solverTolerance = (DREAL)(3e-4);
	solverMaxIter = 10000;
}


template < typename DREAL >
SpRectMatrix<DREAL>::~SpRectMatrix()
{
	if (Ai) delete [] Ai;
	if (Ap) delete [] Ap;
	if (Ax) delete [] Ax;

	if (dAi) cudaFree(dAi);
	if (dAp) cudaFree(dAp);
	if (dAx) cudaFree(dAx);

	if (dunitAi) cudaFree(dunitAi);
	if (dunitAp) cudaFree(dunitAp);
	if (dunitV) cudaFree(dunitV);

	if (dy) cudaFree(dy);
	if (dp) cudaFree(dp);
	
	if (cusparseHandle) cusparseDestroy(cusparseHandle);
	if (cublasHandle) cublasDestroy(cublasHandle);
}
