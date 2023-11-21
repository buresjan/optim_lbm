template < typename REAL >
int SpMatrix< REAL >::matrixFormat()
{
	if (Solver == SOLVER_UMFPACK) return CSC_FORMAT;
	return CSR_FORMAT;
}

template < typename REAL >
int SpMatrix< REAL >::xpos(int i, int j) // i,j, in [0, n-1] !!
{
	switch (matrixFormat())
	{
		case CSC_FORMAT:
			// CSC format:
			// based on Ai and Ap structures get x
			// i ... row index - they are in Ai[ ] 
			// j ... column index
			if (j>=m_n || j<0) { printf("SpMatrix< REAL >::xpos xpos(%d,%d) : column index out of range! max %d \n",i,j,m_n); return 0;}
			for (int k=Ap[j];k<Ap[j+1];k++) if (Ai[k]==i) return k; // this is the xpos
		break;
		case CSR_FORMAT:
			// CSR format:
			// based on Ai and Ap structures get x
			// i ... row index
			// j ... column index - they are in Ai[ ] 
			if (i>=m_n || i<0) { printf("SpMatrix< REAL >::xpos xpos(%d,%d) : row index out of range! max %d \n",i,j,m_n); return 0; }
			for (int k=Ap[i];k<Ap[i+1];k++) if (Ai[k]==j) return k; // this is the xpos
		break;
	}
	return -1;//CTRL->error("SpMatrix< REAL >xpos","unhandled case i %d j %d  n %d nz %d ",i ,j, m_n, m_nz);
}


template < typename REAL >
REAL& SpMatrix< REAL >::get(int i, int j, int verbose)
{
	int xp = xpos(i,j);
	if (xp < 0) return null; 
	return Ax[xp];
}


template < typename REAL >
int SpMatrix< REAL >::solve(REAL *b, REAL *x)
{
	switch (Solver)
	{
		case SOLVER_UMFPACK:  		return solveUMFPACK(b, x); 
		case SOLVER_PETSC: 		return solvePETSC(b, x); 
		default: 
			printf("SpMatrix< REAL >::solve Solver %d undefined.\n", Solver);
			break;
	}
	return 1;
}


template < typename REAL >
int SpMatrix< REAL >::initializeUMFPACK()
{
	#ifdef USE_UMFPACK
	if (UMFPACK_INITIALIZED) return 1;
	umfpack_di_defaults (umf_Control) ;
	umf_Control[UMFPACK_PRL] = 0;//6;
		// new RFFIXME
		umf_Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_NONE;
//		umf_Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_AMD;
//	umf_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
	umf_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;
	umf_Control[UMFPACK_SINGLETONS] = 0;
	umf_Control[UMFPACK_AGGRESSIVE] = 0;
//		umf_Control[UMFPACK_SCALE] = UMFPACK_SCALE_SUM;
	umf_Control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE;
//		umf_Control[UMFPACK_SCALE] = UMFPACK_SCALE_MAX;
	umfpack_di_report_control (umf_Control) ;
	umf_Control[UMFPACK_PRL] = 0;
	UMFPACK_INITIALIZED=true;
	return 1;
	#else
	printf("SpMatrix< REAL >::initializeUMFPACK UMFPACK not included in the current build.");
	return 0;
	#endif
}



template < typename REAL >
int SpMatrix< REAL >::solveUMFPACK(REAL *b, REAL *x)
{
	#ifdef USE_UMFPACK
	initializeUMFPACK();
	int status;
	
//	if (doFactorizeOnce && !)
	if (!factorized)
	{
	
		status =  umfpack_di_symbolic (m_n, m_n, Ap, Ai, Ax, &Symbolic, umf_Control, umf_Info) ;
		umfpack_di_report_status (umf_Control, status) ;
		if (status < 0) 
		{ 
			umfpack_di_report_status (umf_Control, status); 
			printf("SpMatrix< REAL >::solveUMFPACK umfpack_di_symbolic error status %d",status);
			return 0; 
		}
		if (status > 0) printf("SpMatrix< REAL >::solveUMFPACK umfpack_di_symbolic warning status %d",status);
	
		status = umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, umf_Control, umf_Info) ;
		umfpack_di_report_status (umf_Control, status) ;
	
		if (status < 0) 
		{ 
			umfpack_di_report_status (umf_Control, status); 
			printf("SpMatrix< REAL >::solveUMFPACK umfpack_di_numeric error status %d",status);
			return 0; 
		}
		if (status > 0) printf("SpMatrix< REAL >::solveUMFPACK umfpack_di_numeric warning status %d",status);
		umfpack_di_free_symbolic (&Symbolic) ;
	}
	
	if (doFactorizeOnce) factorized=true;

	status = umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, umf_Control, umf_Info) ;
	umfpack_di_report_status (umf_Control, status);
	if (status < 0) 
	{
		umfpack_di_report_status (umf_Control, status); 
		printf("SpMatrix< REAL >::solveUMFPACK umfpack_di_solve error status %d",status);
		return 0;
	}
	if (status > 0) printf("SpMatrix< REAL >::solveUMFPACK umfpack_di_solve warning status %d",status);

	if (!factorized) umfpack_di_free_numeric (&Numeric);
	return 1;
	#else
	printf("SpMatrix< REAL >::solveUMFPACK UMFPACK not included in the current build.");
	return 0;
	#endif
}


template < typename REAL >
int SpMatrix< REAL >::initializePETSC()
{
	#ifdef USE_PETSC
	PetscScalar ierr;
	if (PETSC_INITIALIZED) return 1;
	ierr = VecCreate(PETSC_COMM_WORLD,&PETSC_x);CHKERRQ(ierr);
	ierr = VecSetSizes(PETSC_x,PETSC_DECIDE,m_n);CHKERRQ(ierr);
	ierr = VecSetFromOptions(PETSC_x);CHKERRQ(ierr);
	ierr = VecDuplicate(PETSC_x,&PETSC_b);CHKERRQ(ierr);

// 	Create matrix.  When using MatCreate(), the matrix format can be specified at runtime.
// 	Performance tuning note:  For problems of substantial size,
// 	pTlocation of matrix memory is crucial for attaining good
// 	performance. See the matrix chapter of the users manual for details.
	PetscInt *nnz = new PetscInt[m_n];
	for (int i=0;i<m_n;i++) nnz[i] = Ap[i+1]-Ap[i];
	ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD,m_n,m_n,0,nnz,&PETSC_A);CHKERRQ(ierr);
	delete [] nnz;
	
// 	Create linear solver context
	ierr = KSPCreate(PETSC_COMM_WORLD,&PETSC_ksp);CHKERRQ(ierr);
	
	PETSC_INITIALIZED=true;
	return 1;
	#else
	printf("SpMatrix< REAL >::initializePETSC PETSc not included in the current build.");
	return 0;
	#endif
}

template < typename REAL >
int SpMatrix< REAL >::solvePETSC(REAL *i_b, REAL *i_x)
{
	#ifdef USE_PETSC
	PC             pc;
	PetscErrorCode ierr;
	PetscInt       i, its;
	PetscMPIInt    size;
	PetscScalar    *bufx;

	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	if (size != 1)
	{
		printf("SpMatrix< REAL >::solvePETSC MPI size %d. Uniprocessor computation is supported only.", size);
		SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor code only!");
	}

	initializePETSC();

	// copy i_b --> b
	for (i=0;i<m_n;i++) VecSetValue(PETSC_b, i, (PetscScalar)i_b[i], INSERT_VALUES);
	// copy i_x --> x as non-zero guess
	for (i=0;i<m_n;i++) VecSetValue(PETSC_x, i, (PetscScalar)i_x[i], INSERT_VALUES);

	// copy matrix data from this to PETSc
	// TODO; optimiza memory access
	for (int i1=0;i1<m_n;i1++)
	for (int j1=Ap[i1];j1<Ap[i1+1];j1++)
	{
		MatSetValue(PETSC_A, i1, Ai[j1], (PetscScalar)Ax[j1], INSERT_VALUES);CHKERRQ(ierr);
//		MatSetValue(PETSC_A, Ai[j1], i1, (PetscScalar)get(Ai[j1],i1), INSERT_VALUES);CHKERRQ(ierr);
	}
	ierr = MatAssemblyBegin(PETSC_A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	ierr = MatAssemblyEnd(PETSC_A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	
//	if (CTRL->verbosity(V_DEBUG_SOLVER_MATVIEW)) MatView(PETSC_A,PETSC_VIEWER_DRAW_WORLD);

	ierr = KSPSetOperators(PETSC_ksp,PETSC_A,PETSC_A);CHKERRQ(ierr);

	// configure linear solver -- only default options, override using command line parameters
	ierr = KSPGetPC(PETSC_ksp,&pc);CHKERRQ(ierr);
//	ierr = PCSetType(pc,PCILU);CHKERRQ(ierr);
	ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr); 
	ierr = KSPSetType(PETSC_ksp, KSPBCGS);CHKERRQ(ierr); 
	
//	ierr = KSPSetTolerances(PETSC_ksp,CTRL->SolverTolerance,PETSC_DEFAULT,PETSC_DEFAULT,CTRL->SolverMaxIter);CHKERRQ(ierr);

// 	Set runtime options, e.g.,
// 	-ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
// 	These options will override those specified above as long as
// 	KSPSetFromOptions() is called _after_ any other customization
// 	routines.
	ierr = KSPSetFromOptions(PETSC_ksp);CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(PETSC_ksp,PETSC_TRUE);CHKERRQ(ierr);

// 	Solve linear system
	ierr = KSPSolve(PETSC_ksp,PETSC_b,PETSC_x);CHKERRQ(ierr);

// 	View solver info; we could instead use the option -ksp_view to
// 	print this info to the screen at the conclusion of KSPSolve().
//	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

//	Check solution and clean up
	ierr = KSPGetIterationNumber(PETSC_ksp,&its);CHKERRQ(ierr);
//	if (CTRL->verbosity(V_DEBUG_SOLVER)) CTRL->log("PETSC iteration solver: number of iterations %d", its);

	// copy x --> i_x
	VecGetArray(PETSC_x, &bufx);
	for (i=0;i<m_n;i++) i_x[i] = bufx[i];

	return 1;
	#else
	printf("SpMatrix< REAL >::solvePETSC PETSc not included in the current build.");
	return 0;
	#endif
}


template < typename REAL >
int SpMatrix< REAL >::transpose()
{
	// transpose the matrix --- note,, a stupid nxn algorithm is implemented, the matrix MUST be symmetrically allocated
	REAL temp;
	for (int i=0;i<m_n;i++)
		for (int j=i+1;j<m_n;j++)
		{
			temp = get(i,j,0);
			get(i,j,0) = get(j,i,0);
			get(j,i,0) = temp;
		}
	return 1;
}

template < typename REAL >
int SpMatrix< REAL >::empty()
{
	for (int i=0;i<m_nz;i++) Ax[i] = 0;
	return 1;
}

template < typename REAL >
int SpMatrix< REAL >::print_nonzeros()
{
	printf("Matrix print_nonzeros\n");
	for (int i=0;i<m_nz;i++) if (Ax[i]!=0) printf(" i %d Ax %e \n", i, Ax[i]);
	return 1;
}

template < typename REAL >
int SpMatrix< REAL >::print()
{
	REAL val;
	printf("Matrix print \n");
	for (int i=0;i<m_n;i++)
	{
//			        printf("row %d : ",i);
//			        printf("\t %d",i);
		for (int j=0;j<m_n;j++)
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


template < typename REAL >
int SpMatrix< REAL >::print2()
{
	REAL val;
	printf("Matrix print \n");
	for (int i=0;i<m_n;i++)
	{
	        printf("Row: %d ->  ",i);
//			        printf("\t %d",i);
		for (int j=0;j<m_n;j++)
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



template < typename REAL >
int SpMatrix< REAL >::print2d()
{
	REAL val;
	printf("Matrix print \n");
	for (int i=0;i<m_n;i++)
	{
			        printf("row %d : ",i);
			        printf("\t %d",i);
		for (int j=0;j<m_n;j++)
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

template < typename REAL >
bool SpMatrix< REAL >::isEntryInArray(int *array, int entry, int pos)
{
	for (int i=0;i<pos;i++) if (array[i] == entry) return true;
	return false;
}

template < typename REAL >
int SpMatrix< REAL >::arrayAddEntry(int *array, int entry, int &pos)
// ordered array - add entry
{
	if (isEntryInArray(array,entry,pos)) return 1;
	array[pos] = entry;
	// reorder
	for (int i=pos;i>0;i--)
	{
		if (array[i] < array[i-1]) 
		{
			array[i] = array[i-1];
			array[i-1]=entry;
		}  //else return 1;
	}
	pos++;
	return 1;
}

template < typename REAL >
int SpMatrix< REAL >::arrayAddEntry(int *array, int entry, int &pos, int max)
// ordered array - add entry
{
	if (isEntryInArray(array,entry,pos)) return 1;
	if (pos >= max) { printf("SpMatrix< REAL >::arrayAddEntry out of range pos %d max %d",pos, max); return 0; }
	return arrayAddEntry(array, entry, pos);
}

template < typename REAL >
SpMatrix<REAL>::SpMatrix()
{
//	m_n = i_n;
	m_n=0;
	Ai=0;
	Ap=0;
	Ax=0;

	Solver = SOLVER_UMFPACK;

	null=0; // dummy variable that always returns zero (zero entry in the matrix)
	
	#ifdef USE_UMFPACK
	UMFPACK_INITIALIZED=false;
	#endif
	doFactorizeOnce=true;
	factorized=false;
	
	
	#ifdef USE_PETSC
	PETSC_INITIALIZED=false;
	#endif
}

template < typename REAL >
SpMatrix< REAL >::~SpMatrix()
{
	#ifdef USE_PETSC
	if (PETSC_INITIALIZED)
	{
		VecDestroy(&PETSC_x);//CHKERRQ(ierr);
		VecDestroy(&PETSC_b);//CHKERRQ(ierr);
		MatDestroy(&PETSC_A);//CHKERRQ(ierr);
		KSPDestroy(&PETSC_ksp);//CHKERRQ(ierr);
	}
	#endif
	
	#ifdef USE_UMFPACK
		if (factorized) umfpack_di_free_numeric (&Numeric);
	#endif

	if (Ai) delete [] Ai;
	if (Ap) delete [] Ap;
	if (Ax) delete [] Ax;
}
