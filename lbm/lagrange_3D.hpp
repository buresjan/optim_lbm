#include "timeutils.h"

template< typename LBM >
typename LBM::T_TRAITS::real Lagrange3D<LBM>::computeMinDist()
{
	minDist=1e10;
	maxDist=-1e10;
	for (int i=0;i<LL.size()-1;i++)
	{
	for (int j=i+1;j<LL.size();j++)
	{
		real d = dist(LL[j],LL[i]);
		if (d>maxDist) maxDist=d;
		if (d<minDist) minDist=d;
	}
	if (i%1000==0)
	printf("computeMinDist %d of %d\n", i, LL.size());
	}
	return minDist;
}

template< typename LBM >
typename LBM::T_TRAITS::real Lagrange3D<LBM>::computeMaxDistFromMinDist(typename LBM::T_TRAITS::real mindist)
{
	maxDist=-1e10;
	for (int i=0;i<LL.size()-1;i++)
	{
		real search_dist = 2.0*mindist;
		for (int j=i+1;j<LL.size();j++)
		{
			real d = dist(LL[j],LL[i]);
			if (d < search_dist)
				if (d>maxDist) maxDist=d;
		}
	}
	return maxDist;
}


template< typename LBM >
void Lagrange3D<LBM>::computeMaxMinDist()
{
	maxDist=-1e10;
	minDist=1e10;
	if (lag_X<=0 || lag_Y<=0) return;
	for (int i=0;i<lag_X-1;i++)
	for (int j=0;j<lag_Y-1;j++)
	{
		int index = findIndex(i,j);
		for (int i1=0;i1<=1;i1++)
		for (int j1=0;j1<=1;j1++)
		if (j1!=0 || i1!=0)
		{
			int index1 = findIndex(i+i1,j+j1);
//			int index1 = findIndex(i+i1,j);
			real d = dist(LL[index1],LL[index]);
			if (d>maxDist) maxDist=d;
			if (d<minDist) minDist=d;
		}
	}
}


template< typename LBM >
int Lagrange3D<LBM>::createIndexArray()
{
	if (lag_X<=0 || lag_Y<=0) return 0;
	index_array = new int*[lag_X];
	for (int i=0;i<lag_X;i++) index_array[i] = new int[lag_Y];
	for (int k=0;k<LL.size();k++) index_array[LL[k].lag_x][LL[k].lag_y] = k;
	indexed=true;
	return 1;
}


template< typename LBM >
int Lagrange3D<LBM>::findIndex(int i, int j)
{
	if (!indexed) createIndexArray();
	if (!indexed) return 0;
	return index_array[i][j];
	// brute force
//	for (int k=0;k<LL.size();k++) if (LL[k].lag_x == i && LL[k].lag_y == j) return k;
//	printf("findIndex(%d,%d): not found\n",i,j);
//	return 0;
}

template< typename LBM >
int Lagrange3D<LBM>::findIndexOfNearestX(typename LBM::T_TRAITS::real x)
{
	int imin=0;
	real xmindist=fabs(LL[imin].x-x);
	// brute force
	for (int k=1;k<LL.size();k++) 
	if (fabs(LL[k].x - x) < xmindist)
	{
		imin=k;
		xmindist = fabs(LL[imin].x-x);
	}
	return imin;
}

template< typename LBM >
typename LBM::T_TRAITS::real Lagrange3D<LBM>::diracDelta(int i, typename LBM::T_TRAITS::real r)
{
	switch (i)
	{
		case 1: // VU: phi3
			if(fabs(r) < 1.0)
				return (1.0 - fabs(r));
			else
				return 0;
		case 2: // VU: phi2
			if(fabs(r) < 2.0)
				return 0.25*((1.0+cos((PI*r)/2.0)));
			else
				return 0;
		case 3: // VU: phi1
			if (fabs(r)>2.0) 
				return 0;
			else if (fabs(r)>1.0) 
				return (5.0 - 2.0*fabs(r) - sqrt(-7.0 + 12.0*fabs(r) - 4.0*r*r))/8.0;
			else 
				return (3.0 - 2.0*fabs(r) + sqrt(1.0 + 4.0*fabs(r) - 4.0*r*r))/8.0;
/*	
			if(r > -2.0 && r <= -1.0 )
				result = (5.0 + 2.0*r - sqrt(-7.0 - 12.0*r - 4.0*r*r))/8.0;
			else if(r > -1.0 && r<= 0)
				result = (3.0 + 2.0*r + sqrt(1.0 - 4.0*r - 4.0*r*r))/8.0;
			else if(r > 0 && r <= 1.0)
				result = (3.0 - 2.0*r + sqrt(1.0 + 4.0*r - 4.0*r*r))/8.0;
			else if(r > 1.0 && r <2.0)
				result = (5.0 - 2.0*r - sqrt(-7.0 + 12.0*r - 4.0*r*r))/8.0;
			else result = 0;
*/			
		case 4: // VU: phi4
			if (fabs(r)>1.5)
				return 0;
			else if (fabs(r)>0.5)
				return (5.0 - 3.0*fabs(r) - sqrt(-2.0+6.0*fabs(r)-3.0*r*r))/6.0;
			else
				return (1.0 + sqrt(1.0 - 3.0*r*r))/3.0;
	}
	printf("warning: zero Dirac delta: type=%d\n",i);
	return 0;
}


template< typename LBM >
void Lagrange3D<LBM>::constructWuShuMatricesSparse()
{

	if (ws_constructed) return;
	ws_constructed=true;

	
	int rDirac=1;
	// count non zero elements in matrix ws_A
	int m=LL.size();	// number of lagrangian nodes
	int n=lbm.X*lbm.Y*lbm.Z;	// number of eulerian nodes
	// projdi veskery filament a najdi sousedni body (ve vzdaelenosti mensi nez je prekryv delta funkci)
	struct DD_struct
	{
		real DD;
		int ka;
	};


	typedef std::vector<DD_struct> VECDD;
	typedef std::vector<int> VEC;
	typedef std::vector<real> VECR;
	VEC *v = new VEC[m];
	VECDD *vr = new VECDD[m];


	printf("wushu construct loop 1: start\n");
//	int resrv = 2*((ws_speedUpAllocation) ? ws_speedUpAllocationSupport : 10); // we have plenty of memory, we need performance
	real *LLx = new real[m];
	real *LLy = new real[m];
	real *LLz = new real[m];
	int *LLlagx = new int[m];
	for (int i=0;i<m;i++)
	{
		LLx[i]=LL[i].x;
		LLy[i]=LL[i].y;
		LLz[i]=LL[i].z;
		LLlagx[i]=LL[i].lag_x;
	}

	printf("wushu construct loop 1: cont\n");
	#pragma omp parallel for schedule(dynamic)
	for (int el=0;el<m;el++)
	{
		if (el%100==0)printf("progress %5.2f %%    \r",100.0*el/(real)m);
		for (int ka=0;ka<m;ka++)
		{
			bool proceed=true;
			real d1,d2,ddd;
			if (ws_speedUpAllocation)
			{
				if (abs(LLlagx[el] - LLlagx[ka]) > ws_speedUpAllocationSupport) proceed=false;
			}
			if (!proceed) continue;
			if (ws_regularDirac)
			{
				d1 = diracDelta(rDirac,(LLx[el] - LLx[ka])/lbm.physDl);
				if (d1>0)
				{
					d2 = diracDelta(rDirac, (LLy[el] - LLy[ka])/lbm.physDl);
					if (d2>0)
					{
						ddd=d1*d2*diracDelta(rDirac, (LLz[el] - LLz[ka])/lbm.physDl);
						if (ddd>0)
						{
							v[el].push_back(ka);
							DD_struct sdd;
							sdd.DD = ddd;
							sdd.ka = ka;
							vr[el].push_back(sdd);
						}
					}
				}
			} else
			{
				d1 = diracDelta((LLx[el] - LLx[ka])/lbm.physDl/2.0);
				if (d1>0)
				{
					d2 = diracDelta((LLy[el] - LLy[ka])/lbm.physDl/2.0);
					if (d2>0)
					{
						ddd=d1*d2*diracDelta((LLz[el] - LLz[ka])/lbm.physDl/2.0);
						if (ddd>0)
						{
							v[el].push_back(ka);
							DD_struct sdd;
							sdd.DD = ddd;
							sdd.ka = ka;
							vr[el].push_back(sdd);
						}
					}
				}
			}
		}
	}
	delete [] LLx;
	delete [] LLy;
	delete [] LLz;
	delete [] LLlagx;
	// count non zero
	int nz=0;
	for (int el=0;el<m;el++) nz += vr[el].size();
//	printf("non-zeros: %d\n",nz);

	printf("wushu construct loop 1: end\n");



	// create spmatrix
	ws_A = new SpMatrix<real>();
	ws_A->Solver = Solver;
	ws_A->Ap = new int[m+1];
	ws_A->Ai = new int[nz];
	ws_A->Ax = new real[nz];
	ws_A->m_nz = nz;
	ws_A->m_n = m;
	#ifdef USE_CUSPARSE
	// create sprectmatrix ... for cuda
	ws_dA = new SpRectMatrix<dreal>();
	ws_dA->Ap = new int[m+1];
	ws_dA->Ai = new int[nz];
	ws_dA->Ax = new dreal[nz];
	ws_dA->m_nz = nz;
	ws_dA->m_nr = m;
	ws_dA->m_nc = m;
	#endif 

	ws_A->Ap[0]=0;
	for (int i=0;i<nz;i++) ws_A->Ax[i] = 0; // empty 
	
	#ifdef USE_CUSPARSE
	ws_dA->Ap[0]=0;
	for (int i=0;i<nz;i++) ws_dA->Ax[i] = 0; // empty 
	#endif
	
//	printf("Ai construct\n");
	int count=0;
	for (int i=0;i<m;i++)
	{
		ws_A->Ap[i+1] = ws_A->Ap[i] + v[i].size();
		#ifdef USE_CUSPARSE
		ws_dA->Ap[i+1] = ws_dA->Ap[i] + v[i].size();
		#endif
//		printf("Ap[%d]=%d (%d)\n",i+1,ws_A->Ap[i+1],nz);
		for (int j=0;j<v[i].size();j++)
		{
			ws_A->Ai[count]=v[i][j];
			#ifdef USE_CUSPARSE
			ws_dA->Ai[count]=v[i][j];
			#endif
			count++;
		}
	}
	delete [] v;

	// fill vectors delta_el
	// sparse vector of deltas
	d_i = new std::vector<idx>[m];
	d_x = new std::vector<real>[m];
	// fill only non zero elements-relevant
/*
	 // brute force
	printf("wushu construct loop 2: start\n");
	#pragma omp parallel for schedule(static)
	for (int i=0;i<m;i++)
	{
//		if (i%20==0) printf("\r progress: %04d of %04d",i,m);
		for (int gz=0;gz<lbm.Z;gz++)
		for (int gy=0;gy<lbm.Y;gy++)
		for (int gx=0;gx<lbm.X;gx++)
		{
			real dd = diracDelta((real)(gx + 0.5) - LL[i].x/lbm.physDl) * diracDelta((real)(gy + 0.5) - LL[i].y/lbm.physDl) * diracDelta((real)(gz + 0.5) - LL[i].z/lbm.physDl);
			if (dd>0)
			{
				d_i[i].push_back(lbm.pos(gx,gy,gz));
				d_x[i].push_back(dd);
			}
		}
	}
	printf("wushu construct loop 2: end\n");
*/

	printf("wushu construct loop 2: start\n");
	idx support=5; // search in this support
	#pragma omp parallel for schedule(static)
	for (int i=0;i<m;i++)
	{
		idx fi_x = floor(LL[i].x/lbm.physDl);
		idx fi_y = floor(LL[i].y/lbm.physDl);
		idx fi_z = floor(LL[i].z/lbm.physDl);
		
		for (int gz=MAX( 0, fi_z - support);gz<MIN(lbm.Z, fi_z + support);gz++)
		for (int gy=MAX( 0, fi_y - support);gy<MIN(lbm.Y, fi_y + support);gy++)
		for (int gx=MAX( 0, fi_x - support);gx<MIN(lbm.X, fi_x + support);gx++)
		{
			real dd = diracDelta((real)(gx + 0.5) - LL[i].x/lbm.physDl) * diracDelta((real)(gy + 0.5) - LL[i].y/lbm.physDl) * diracDelta((real)(gz + 0.5) - LL[i].z/lbm.physDl);
			if (dd>0)
			{
				d_i[i].push_back(lbm.pos(gx,gy,gz));
				d_x[i].push_back(dd);
			}
		}
	}
	printf("wushu construct loop 2: end\n");


	printf("wushu construct loop 3: start\n");
	#pragma omp parallel for schedule(static)
	for (int i=0;i<m;i++)
	{
		if (i%100==0)printf("progress %5.2f %%    \r",100.0*i/(real)m);
		for (int ka=0;ka<vr[i].size();ka++)
		{
			int j=vr[i][ka].ka;
			real ddd = vr[i][ka].DD;
			if (ws_regularDirac)
			{
				ws_A->get(i,j) = ddd;
				#ifdef USE_CUSPARSE
				ws_dA->get(i,j) = ddd;
				#endif
			} else
			{
				if (ddd>0) // we have non-zero element at i,j
				{
					real val=0;
					for (int in1=0;in1<d_i[i].size();in1++)
					{
						for (int in2=0;in2<d_i[j].size();in2++)
						{
							if (d_i[i][in1]==d_i[j][in2]) 
							{
								val += d_x[i][in1]*d_x[j][in2];
								break;
							}
						}
					}
					ws_A->get(i,j) = val;
					#ifdef USE_CUSPARSE
					ws_dA->get(i,j) = (dreal)val;
					#endif
				}
			}
		}
	}
	delete [] vr; // free
	printf("wushu construct loop 3: end\n");

	printf("wushu construct loop 4: start\n");
	for (int k=0;k<3;k++)
	{
		ws_x[k] = new real[m];
		ws_b[k] = new real[m]; // right hand side
		
		ws_hx[k] = new dreal[m];
		ws_hb[k] = new dreal[m]; // right hand side
	}

//	ws_ds = new real[m]; // delta s_ell
	#ifdef USE_CUDA
	for (int k=0;k<3;k++) 
	{
		cudaMalloc((void **)&ws_dx[k], m*sizeof(dreal));
		cudaMalloc((void **)&ws_db[k], m*sizeof(dreal));
	}
	cudaMalloc((void **)&ws_du,  n*sizeof(dreal));

	// copy zero to x1, x2 (init)
	dreal* zero = new dreal[m];
	for (int i=0;i<m;i++) zero[i]=0;
	for (int k=0;k<3;k++) cudaMemcpy(ws_dx[k], zero, m*sizeof(dreal), cudaMemcpyHostToDevice); // TODO use setCudaValue ... 
	delete [] zero;
	#endif

	// create Matrix M: matrix realizing projection of u* to lagrange desc.
	nz=0;
	for (int el=0;el<m;el++)
	for (int in1=0;in1<d_i[el].size();in1++)
		nz++;
	#ifdef USE_CUSPARSE
	ws_M = new SpRectMatrix<dreal>();
	ws_M->Ap = new int[m+1];
	ws_M->Ai = new int[nz];
	ws_M->Ax = new dreal[nz];
	ws_M->m_nz = nz;
	ws_M->m_nr = m;
	ws_M->m_nc = n;

	ws_M->Ap[0]=0;
	for (int i=0;i<nz;i++) ws_M->Ax[i] = 0; // empty 
	
//	printf("Ai construct\n");
	count=0;
	for (int i=0;i<m;i++)
	{
		ws_M->Ap[i+1] = ws_M->Ap[i] + d_i[i].size();
		for (int j=0;j<d_i[i].size();j++)
		{
			ws_M->Ai[count]=d_i[i][j];
			ws_M->Ax[count]=(dreal)d_x[i][j];
			count++;
		}
	}
	
	
	// its transpose
	ws_MT = new SpRectMatrix<dreal>();
	ws_MT->Ap = new int[n+1];
	ws_MT->Ai = new int[nz];
	ws_MT->Ax = new dreal[nz];
	ws_MT->m_nz = nz;
	ws_MT->m_nr = n;
	ws_MT->m_nc = m;

	ws_MT->Ap[0]=0;
	for (int i=0;i<nz;i++) ws_MT->Ax[i] = 0; // empty 

	// for each Euler node, assign 
	VEC *vn = new VEC[n];
	VECR *vx = new VECR[n];
	for (int i=0;i<m;i++) 
	for (int j=0;j<d_i[i].size();j++) 
	{
		vn[ d_i[i][j] ].push_back( i );
		vx[ d_i[i][j] ].push_back( d_x[i][j] );
	}

	count=0;
	for (int i=0;i<n;i++)
	{
		ws_MT->Ap[i+1] = ws_MT->Ap[i] + vn[i].size();
		for (int j=0;j<vn[i].size();j++)
		{
			ws_MT->Ai[count]=vn[i][j];
			ws_MT->Ax[count]=vx[i][j];
			count++;
		}
	}
	delete [] vn;
	delete [] vx;
	printf("wushu construct loop 4: end\n");
	#endif
}


template< typename LBM >
void Lagrange3D<LBM>::constructWuShuMatricesSparse_TNL()
{
#ifdef USE_TNL
	if (ws_tnl_constructed) return;
	ws_tnl_constructed=true;
	int rDirac=1;
	// count non zero elements in matrix A
	int m=LL.size();	// number of lagrangian nodes
	int n=lbm.X*lbm.Y*lbm.Z;	// number of eulerian nodes
	// projdi veskery filament a najdi sousedni body (ve vzdaelenosti mensi nez je prekryv delta funkci)
	struct DD_struct
	{
		real DD;
		int ka;
	};
	typedef std::vector<DD_struct> VECDD;
	typedef std::vector<int> VEC;
	typedef std::vector<real> VECR;
	VEC *v = new VEC[m];
	VECDD *vr = new VECDD[m];

	printf("tnl wushu construct loop 1: start\n");
//	int resrv = 2*((ws_speedUpAllocation) ? ws_speedUpAllocationSupport : 10); // we have plenty of memory, we need performance
	real *LLx = new real[m];
	real *LLy = new real[m];
	real *LLz = new real[m];
	int *LLlagx = new int[m];
	for (int i=0;i<m;i++)
	{
		LLx[i]=LL[i].x;
		LLy[i]=LL[i].y;
		LLz[i]=LL[i].z;
		LLlagx[i]=LL[i].lag_x;
	}

	printf("tnl wushu construct loop 1: cont\n");
	#pragma omp parallel for schedule(dynamic)
	for (int el=0;el<m;el++)
	{
		if (el%100==0)printf("progress %5.2f %%    \r",100.0*el/(real)m);
		for (int ka=0;ka<m;ka++)
		{
			bool proceed=true;
			real d1,d2,ddd;
			if (ws_speedUpAllocation)
			{
				if (abs(LLlagx[el] - LLlagx[ka]) > ws_speedUpAllocationSupport) proceed=false;
			}
			if (!proceed) continue;
			if (ws_regularDirac)
			{
				d1 = diracDelta(rDirac,(LLx[el] - LLx[ka])/lbm.physDl);
				if (d1>0)
				{
					d2 = diracDelta(rDirac, (LLy[el] - LLy[ka])/lbm.physDl);
					if (d2>0)
					{
						ddd=d1*d2*diracDelta(rDirac, (LLz[el] - LLz[ka])/lbm.physDl);
						if (ddd>0)
						{
							v[el].push_back(ka);
							DD_struct sdd;
							sdd.DD = ddd;
							sdd.ka = ka;
							vr[el].push_back(sdd);
						}
					}
				}
			} else
			{
				d1 = diracDelta((LLx[el] - LLx[ka])/lbm.physDl/2.0);
				if (d1>0)
				{
					d2 = diracDelta((LLy[el] - LLy[ka])/lbm.physDl/2.0);
					if (d2>0)
					{
						ddd=d1*d2*diracDelta((LLz[el] - LLz[ka])/lbm.physDl/2.0);
						if (ddd>0)
						{
							v[el].push_back(ka);
							DD_struct sdd;
							sdd.DD = ddd;
							sdd.ka = ka;
							vr[el].push_back(sdd);
						}
					}
				}
			}
		}
	}
	delete [] LLx;
	delete [] LLy;
	delete [] LLz;
	delete [] LLlagx;

	printf("tnl wushu construct loop 1: end\n");

	// allocate matrix A
	ws_tnl_hA->setDimensions(m, m);
//	int max_nz_per_row=0;
//	for (int el=0;el<m;el++) {
//		max_nz_per_row = TNL::max(max_nz_per_row, vr[el].size());
//	}
//	ws_tnl_hA->setConstantCompressedRowLengths(max_nz_per_row);
	typename hEllpack::CompressedRowLengthsVector hA_row_lengths( m );
	for (int el=0; el<m; el++) hA_row_lengths[el] = vr[el].size();
	ws_tnl_hA->setCompressedRowLengths(hA_row_lengths);
	
	for (int i=0;i<m;i++)
	{
		auto row = ws_tnl_hA->getRow(i);
		for (int j=0;j<v[i].size();j++)
			row.setElement(j, v[i][j], 1);
	}
	delete [] v;

	// fill vectors delta_el
	// sparse vector of deltas
	d_i = new std::vector<idx>[m];
	d_x = new  std::vector<real>[m];
	// fill only non zero elements-relevant

	printf("tnl wushu construct loop 2: start\n");
	idx support=5; // search in this support
	#pragma omp parallel for schedule(static)
	for (int i=0;i<m;i++)
	{
		idx fi_x = floor(LL[i].x/lbm.physDl);
		idx fi_y = floor(LL[i].y/lbm.physDl);
		idx fi_z = floor(LL[i].z/lbm.physDl);
		
		for (int gz=MAX( 0, fi_z - support);gz<MIN(lbm.Z, fi_z + support);gz++)
		for (int gy=MAX( 0, fi_y - support);gy<MIN(lbm.Y, fi_y + support);gy++)
		for (int gx=MAX( 0, fi_x - support);gx<MIN(lbm.X, fi_x + support);gx++)
		{
			real dd = diracDelta((real)(gx + 0.5) - LL[i].x/lbm.physDl) * diracDelta((real)(gy + 0.5) - LL[i].y/lbm.physDl) * diracDelta((real)(gz + 0.5) - LL[i].z/lbm.physDl);
			if (dd>0)
			{
				d_i[i].push_back(lbm.pos(gx,gy,gz));
				d_x[i].push_back(dd);
			}
		}
	}
	printf("tnl wushu construct loop 2: end\n");

	printf("tnl wushu construct loop 3: start\n");
	#pragma omp parallel for schedule(static)
	for (int i=0;i<m;i++)
	{
		if (i%100==0)printf("progress %5.2f %%    \r",100.0*i/(real)m);
		for (int ka=0;ka<vr[i].size();ka++)
		{
			int j=vr[i][ka].ka;
			real ddd = vr[i][ka].DD;
			if (ws_regularDirac)
			{
				ws_tnl_hA->setElement(i,j, ddd);
			} else
			{
				if (ddd>0) // we have non-zero element at i,j
				{
					real val=0;
					for (int in1=0;in1<d_i[i].size();in1++)
					{
						for (int in2=0;in2<d_i[j].size();in2++)
						{
							if (d_i[i][in1]==d_i[j][in2]) 
							{
								val += d_x[i][in1]*d_x[j][in2];
								break;
							}
						}
					}
					ws_tnl_hA->setElement(i,j, val);
				}
			}
		}
	}
	delete [] vr; // free
	printf("tnl wushu construct loop 3: end\n");

	// create vectors for the solution of the linear system
	for (int k=0;k<3;k++)
	{
		ws_tnl_hx[k].setSize(m);
		ws_tnl_hb[k].setSize(m);
		#ifdef USE_CUDA
		ws_tnl_dx[k].setSize(m);
		ws_tnl_db[k].setSize(m);
		ws_tnl_hxz[k].setSize(m);
		ws_tnl_hbz[k].setSize(m);
		#endif
	}

	// zero-initialize x1, x2, x3
	for (int k=0;k<3;k++) ws_tnl_hx[k].setValue(0);
	#ifdef USE_CUDA
	for (int k=0;k<3;k++) ws_tnl_dx[k].setValue(0);
	for (int k=0;k<3;k++) ws_tnl_hxz[k].setValue(0);
	#endif

	#ifdef USE_CUDA
	// create Matrix M: matrix realizing projection of u* to lagrange desc.
	ws_tnl_hM.setDimensions(m, n);
//	max_nz_per_row = 0;
//	for (int el=0;el<m;el++)
//		max_nz_per_row = TNL::max(max_nz_per_row, d_i[el].size());
//	ws_tnl_hM.setConstantCompressedRowLengths(max_nz_per_row);
	typename hEllpack::CompressedRowLengthsVector hM_row_lengths( m );
	for (int el=0; el<m; el++) hM_row_lengths[el] = d_i[el].size();
	ws_tnl_hM.setCompressedRowLengths(hM_row_lengths);

//	printf("Ai construct\n");
	for (int i=0;i<m;i++)
	{
		auto row = ws_tnl_hM.getRow(i);
		for (int j=0;j<d_i[i].size();j++)
			row.setElement(j, d_i[i][j], (dreal)d_x[i][j]);
	}
	
	// its transpose
	ws_tnl_hMT.setDimensions(n, m);

	// for each Euler node, assign 
	VEC *vn = new VEC[n];
	VECR *vx = new VECR[n];
	for (int i=0;i<m;i++) 
	for (int j=0;j<d_i[i].size();j++) 
	{
		vn[ d_i[i][j] ].push_back( i );
		vx[ d_i[i][j] ].push_back( d_x[i][j] );
	}

//	max_nz_per_row = 0;
//	for (int el=0;el<n;el++)
//		max_nz_per_row = TNL::max(max_nz_per_row, vn[el].size());
//	ws_tnl_hMT.setConstantCompressedRowLengths(max_nz_per_row);
	typename hEllpack::CompressedRowLengthsVector hMT_row_lengths( n );
	for (int el=0; el<n; el++) hMT_row_lengths[el] = vn[el].size();
	ws_tnl_hMT.setCompressedRowLengths(hMT_row_lengths);

	for (int i=0;i<n;i++)
	{
		auto row = ws_tnl_hMT.getRow(i);
		for (int j=0;j<vn[i].size();j++)
			row.setElement(j, vn[i][j], vx[i][j]);
	}
	delete [] vn;
	delete [] vx;
	#endif

	// update the preconditioner
	ws_tnl_hprecond->update(ws_tnl_hA);

	#ifdef USE_CUDA
	// copy matrices from host to the GPU
	*ws_tnl_dA = *ws_tnl_hA;
	ws_tnl_dA.synchronize();
	ws_tnl_dM = ws_tnl_hM;
	ws_tnl_dMT = ws_tnl_hMT;

	// update the preconditioner
	ws_tnl_dprecond->update(ws_tnl_dA);
	#endif
	printf("tnl wushu lagrange_3D_end\n");
	printf("number of lagrangian points: %d\n", m);

	const char* compute_desc;
	switch (ws_compute)
	{
		case ws_computeCPU:                    compute_desc = "ws_computeCPU"; break;
		case ws_computeGPU_CUSPARSE:           compute_desc = "ws_computeGPU_CUSPARSE"; break;
		case ws_computeHybrid_CUSPARSE:        compute_desc = "ws_computeHybrid_CUSPARSE"; break;
		case ws_computeCPU_TNL:                compute_desc = "ws_computeCPU_TNL"; break;
		case ws_computeGPU_TNL:                compute_desc = "ws_computeGPU_TNL"; break;
		case ws_computeHybrid_TNL:             compute_desc = "ws_computeHybrid_TNL"; break;
		case ws_computeHybrid_TNL_zerocopy:    compute_desc = "ws_computeHybrid_TNL_zerocopy"; break;
	}
	log("constructed WuShu matrices for ws_compute=%s", compute_desc);
#endif // USE_TNL
}

#ifdef USE_TNL
template< typename Matrix, typename Vector >
__cuda_callable__
typename Matrix::RealType
rowVectorProduct( const Matrix& matrix, typename Matrix::IndexType i, const Vector& vector )
{
    typename Matrix::RealType result = 0;
    const auto row = matrix.getRow( i );

    for( typename Matrix::IndexType c = 0; c < row.getSize(); c++ ) {
        const typename Matrix::IndexType column = row.getColumnIndex( c );
        if( column != matrix.getPaddingIndex() )
            result += row.getValue( c ) * vector[ column ];
    }

    return result;
}
#endif

//require: rho, vx, vy, vz
template< typename LBM >
void Lagrange3D<LBM>::computeWuShuForcesSparse(real time)
{
	switch (ws_compute)
	{
		case ws_computeCPU_TNL:
		case ws_computeGPU_TNL:
		case ws_computeHybrid_TNL:
		case ws_computeHybrid_TNL_zerocopy:
			constructWuShuMatricesSparse_TNL();
			break;
		default:
			constructWuShuMatricesSparse();
			break;
	}

	int m=LL.size();
	idx n=lbm.X*lbm.Y*lbm.Z;

	#ifdef USE_TNL
	using VectorView = TNL::Containers::VectorView< dreal, TNL::Devices::Cuda, idx >;
	using ConstVectorView = TNL::Containers::VectorView< const dreal, TNL::Devices::Cuda, idx >;
	#endif

	switch (ws_compute)
	{
		#ifdef USE_CUDA
		case ws_computeGPU_TNL:
		{
			#ifdef USE_TNL
			// no Device--Host copy is required
			ws_tnl_dM.vectorProduct(ConstVectorView(lbm.dvx(), n), ws_tnl_db[0], -1.0);
			ws_tnl_dM.vectorProduct(ConstVectorView(lbm.dvy(), n), ws_tnl_db[1], -1.0);
			ws_tnl_dM.vectorProduct(ConstVectorView(lbm.dvz(), n), ws_tnl_db[2], -1.0);
			// solver
			for (int k=0;k<3;k++) {
				auto start = std::chrono::steady_clock::now();
				ws_tnl_dsolver.solve(ws_tnl_db[k], ws_tnl_dx[k]);
				auto end = std::chrono::steady_clock::now();
				auto int_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
				real WT = int_ms * 1e-6;
				log("t=%es k=%d TNL CG solver: WT=%e iterations=%d residual=%e", time, k, WT, ws_tnl_dsolver.getIterations(), ws_tnl_dsolver.getResidue());
			}
			ConstVectorView x1(ws_tnl_dx[0].getData(), ws_tnl_dx[0].getSize());
			ConstVectorView x2(ws_tnl_dx[1].getData(), ws_tnl_dx[1].getSize());
			ConstVectorView x3(ws_tnl_dx[2].getData(), ws_tnl_dx[2].getSize());
			ConstVectorView rho(lbm.drho(), n);
			VectorView fx(lbm.dfx(), n);
			VectorView fy(lbm.dfy(), n);
			VectorView fz(lbm.dfz(), n);
			//TNL::Pointers::DevicePointer<dEllpack> MT(ws_tnl_dMT);
			TNL::Pointers::DevicePointer<dEllpack> MT_dptr(ws_tnl_dMT);
			const dEllpack* MT = &MT_dptr.template getData<TNL::Devices::Cuda>();
			auto kernel = [=] CUDA_HOSTDEV (idx i) mutable
			{
				// skipping empty rows explicitly is much faster
				if( MT->getRowCapacity(i) > 0 ) {
					fx[i] += 2*rho[i]*rowVectorProduct(*MT, i, x1);
					fy[i] += 2*rho[i]*rowVectorProduct(*MT, i, x2);
					fz[i] += 2*rho[i]*rowVectorProduct(*MT, i, x3);
				}
			};
			TNL::Algorithms::ParallelFor< TNL::Devices::Cuda >::exec((idx) 0, n, kernel);
			#else
			printf("ws_tnl_computeGPU_TNL failed: TNL not included in the build. \n");
			#endif
			break;
		}

		case ws_computeHybrid_TNL:
		{
			#ifdef USE_TNL
			ws_tnl_dM.vectorProduct(ConstVectorView(lbm.dvx(), n), ws_tnl_db[0], -1.0);
			ws_tnl_dM.vectorProduct(ConstVectorView(lbm.dvy(), n), ws_tnl_db[1], -1.0);
			ws_tnl_dM.vectorProduct(ConstVectorView(lbm.dvz(), n), ws_tnl_db[2], -1.0);
			// copy to Host
			for (int k=0;k<3;k++) ws_tnl_hb[k] = ws_tnl_db[k];
			// solve on CPU
			for (int k=0;k<3;k++) {
				auto start = std::chrono::steady_clock::now();
				ws_tnl_hsolver.solve(ws_tnl_hb[k], ws_tnl_hx[k]);
				auto end = std::chrono::steady_clock::now();
				auto int_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
				real WT = int_ms * 1e-6;
				log("t=%es k=%d TNL CG solver: WT=%e iterations=%d residual=%e", time, k, WT, ws_tnl_hsolver.getIterations(), ws_tnl_hsolver.getResidue());
			}
			// copy to GPU
			for (int k=0;k<3;k++) ws_tnl_dx[k] = ws_tnl_hx[k];
			// continue on GPU
			ConstVectorView x1(ws_tnl_dx[0].getData(), ws_tnl_dx[0].getSize());
			ConstVectorView x2(ws_tnl_dx[1].getData(), ws_tnl_dx[1].getSize());
			ConstVectorView x3(ws_tnl_dx[2].getData(), ws_tnl_dx[2].getSize());
			ConstVectorView rho(lbm.drho(), n);
			VectorView fx(lbm.dfx(), n);
			VectorView fy(lbm.dfy(), n);
			VectorView fz(lbm.dfz(), n);
//			TNL::Pointers::DevicePointer<dEllpack> MT(ws_tnl_dMT);
			TNL::Pointers::DevicePointer<dEllpack> MT_dptr(ws_tnl_dMT);
			const dEllpack* MT = &MT_dptr.template getData<TNL::Devices::Cuda>();
			auto kernel = [=] CUDA_HOSTDEV (idx i) mutable
			{
				// skipping empty rows explicitly is much faster
				if( MT->getRowCapacity(i) > 0 ) {
					fx[i] += 2*rho[i]*rowVectorProduct(*MT, i, x1);
					fy[i] += 2*rho[i]*rowVectorProduct(*MT, i, x2);
					fz[i] += 2*rho[i]*rowVectorProduct(*MT, i, x3);
				}
			};
			TNL::Algorithms::ParallelFor< TNL::Devices::Cuda >::exec((idx) 0, n, kernel);
			#else
			printf("ws_tnl_computeHybrid_TNL failed: TNL not included in the build. \n");
			#endif
			break;
		}

		case ws_computeHybrid_TNL_zerocopy:
		{
			#ifdef USE_TNL
			ws_tnl_dM.vectorProduct(ConstVectorView(lbm.dvx(), n), ws_tnl_hbz[0], -1.0);
			ws_tnl_dM.vectorProduct(ConstVectorView(lbm.dvy(), n), ws_tnl_hbz[1], -1.0);
			ws_tnl_dM.vectorProduct(ConstVectorView(lbm.dvz(), n), ws_tnl_hbz[2], -1.0);
			// solve on CPU
			for (int k=0;k<3;k++) {
				auto start = std::chrono::steady_clock::now();
				ws_tnl_hsolver.solve(ws_tnl_hbz[k].getView(), ws_tnl_hxz[k]);
				auto end = std::chrono::steady_clock::now();
				auto int_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
				real WT = int_ms * 1e-6;
				log("t=%es k=%d TNL CG solver: WT=%e iterations=%d residual=%e", time, k, WT, ws_tnl_hsolver.getIterations(), ws_tnl_hsolver.getResidue());
			}
			// continue on GPU
			ConstVectorView x1(ws_tnl_hxz[0].getData(), ws_tnl_hxz[0].getSize());
			ConstVectorView x2(ws_tnl_hxz[1].getData(), ws_tnl_hxz[1].getSize());
			ConstVectorView x3(ws_tnl_hxz[2].getData(), ws_tnl_hxz[2].getSize());
			ConstVectorView rho(lbm.drho(), n);
			VectorView fx(lbm.dfx(), n);
			VectorView fy(lbm.dfy(), n);
			VectorView fz(lbm.dfz(), n);
//			TNL::Pointers::DevicePointer<dEllpack> MT(ws_tnl_dMT);
			TNL::Pointers::DevicePointer<dEllpack> MT_dptr(ws_tnl_dMT);
			const dEllpack* MT = &MT_dptr.template getData<TNL::Devices::Cuda>();
			auto kernel = [=] CUDA_HOSTDEV (idx i) mutable
			{
				// skipping empty rows explicitly is much faster
				if( MT->getRowCapacity(i) > 0 ) {
					fx[i] += 2*rho[i]*rowVectorProduct(*MT, i, x1);
					fy[i] += 2*rho[i]*rowVectorProduct(*MT, i, x2);
					fz[i] += 2*rho[i]*rowVectorProduct(*MT, i, x3);
				}
			};
			TNL::Algorithms::ParallelFor< TNL::Devices::Cuda >::exec((idx) 0, n, kernel);
			#else
			printf("ws_tnl_computeHybrid_TNL_zerocopy failed: TNL not included in the build. \n");
			#endif
			break;
		}

		case ws_computeGPU_CUSPARSE:
		{
			#ifdef USE_CUSPARSE
			// no Device--Host copy is required
			// HOKUS POKUS
			ws_M->dmultiply(lbm.dvx(), ws_db[0], -1.0);
			ws_M->dmultiply(lbm.dvy(), ws_db[1], -1.0);
			ws_M->dmultiply(lbm.dvz(), ws_db[2], -1.0);
			// solver
			for (int k=0;k<3;k++) ws_dA->dsolve(ws_db[k], ws_dx[k]);
			// now compute delta u
			ws_MT->dmultiply(ws_dx[0], ws_du, 2.0);
			ws_MT->dVectorMultiplyAndAdd(n, ws_du, lbm.drho(), 1.0, lbm.dfx());
			ws_MT->dmultiply(ws_dx[1], ws_du, 2.0);
			ws_MT->dVectorMultiplyAndAdd(n, ws_du, lbm.drho(), 1.0, lbm.dfy());
			ws_MT->dmultiply(ws_dx[2], ws_du, 2.0);
			ws_MT->dVectorMultiplyAndAdd(n, ws_du, lbm.drho(), 1.0, lbm.dfz());
			#else
			printf("ws_computeHybrid_CUSPARSE failed: CUSPARSE not included in the build. \n");
			#endif
			break;
		}
		case ws_computeHybrid_CUSPARSE:
		{
			#ifdef USE_CUSPARSE
			ws_M->dmultiply(lbm.dvx(), ws_db[0], -1.0);
			ws_M->dmultiply(lbm.dvy(), ws_db[1], -1.0);
			ws_M->dmultiply(lbm.dvy(), ws_db[2], -1.0);
			// copy to Host
			for (int k=0;k<3;k++) cudaMemcpy(ws_hb[k],ws_db[k],m*sizeof(dreal), cudaMemcpyDeviceToHost);
			// retype
			for (int k=0;k<3;k++)
			for (int i=0;i<m;i++)
				ws_b[k][i]=(real)ws_hb[k][i];
			// solve on CPU
			for (int k=0;k<3;k++)
				ws_A->solve(ws_b[k],ws_x[k]);
			// retype and copy to GPU
			for (int k=0;k<3;k++)
			for (int i=0;i<m;i++)
				ws_hx[k][i]=(dreal)ws_x[k][i];
			for (int k=0;k<3;k++)
				cudaMemcpy(ws_dx[k],ws_hx[k],m*sizeof(dreal), cudaMemcpyHostToDevice);
			// continue on GPU
			ws_MT->dmultiply(ws_dx[0], ws_du, 2.0);
			ws_MT->dVectorMultiplyAndAdd(n,ws_du,lbm.drho(),1.0, lbm.dfx());
			ws_MT->dmultiply(ws_dx[1], ws_du, 2.0);
			ws_MT->dVectorMultiplyAndAdd(n,ws_du,lbm.drho(),1.0, lbm.dfy());
			ws_MT->dmultiply(ws_dx[2], ws_du, 2.0);
			ws_MT->dVectorMultiplyAndAdd(n,ws_du,lbm.drho(),1.0, lbm.dfz());
			#else
			printf("ws_computeHybrid_CUSPARSE failed: CUSPARSE not included in the build. \n");
			#endif
			break;
		}
		#endif // USE_CUDA
		case ws_computeCPU_TNL:
		{
			#ifdef USE_TNL
			// vx, vy, vz, rho must be copied from the device
			ws_tnl_hM.vectorProduct(ConstVectorView(lbm.hvx(), n), ws_tnl_hb[0], -1.0);
			ws_tnl_hM.vectorProduct(ConstVectorView(lbm.hvy(), n), ws_tnl_hb[1], -1.0);
			ws_tnl_hM.vectorProduct(ConstVectorView(lbm.hvz(), n), ws_tnl_hb[2], -1.0);
			// solver
			for (int k=0;k<3;k++) {
				auto start = std::chrono::steady_clock::now();
				ws_tnl_hsolver.solve(ws_tnl_hb[k], ws_tnl_hx[k]);
				auto end = std::chrono::steady_clock::now();
				auto int_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
				real WT = int_ms * 1e-6;
				log("t=%es k=%d TNL CG solver: WT=%e iterations=%d residual=%e", time, k, WT, ws_tnl_hsolver.getIterations(), ws_tnl_hsolver.getResidue());
			}
			auto kernel = [&] (idx i) mutable
			{
				// skipping empty rows explicitly is much faster
				if( ws_tnl_hMT.getRowCapacity(i) > 0 ) {
					lbm.hfx()[i] += 2 * lbm.hrho()[i] * rowVectorProduct(ws_tnl_hMT, i, ws_tnl_hx[0]);
					lbm.hfy()[i] += 2 * lbm.hrho()[i] * rowVectorProduct(ws_tnl_hMT, i, ws_tnl_hx[1]);
					lbm.hfz()[i] += 2 * lbm.hrho()[i] * rowVectorProduct(ws_tnl_hMT, i, ws_tnl_hx[2]);
				}
			};
			TNL::Algorithms::ParallelFor< TNL::Devices::Host >::exec((idx) 0, n, kernel);
			#else
			printf("ws_tnl_computeCPU_TNL failed: TNL not included in the build. \n");
			#endif
			break;
		}
		case ws_computeCPU:
		{
			// vx, vy, rho must be copied from the device
			// fx, fy must be zero
			for (int el=0;el<m;el++)
			{
				for (int k=0;k<3;k++)
					ws_b[k][el]=0;
				for (int in1=0;in1<d_i[el].size();in1++)
				{
					int gi=d_i[el][in1];
					ws_b[0][el] -= (real)lbm.hvx()[gi] * d_x[el][in1];
					ws_b[1][el] -= (real)lbm.hvy()[gi] * d_x[el][in1];
					ws_b[2][el] -= (real)lbm.hvz()[gi] * d_x[el][in1];
				}
			}
			//solver
			for (int k=0;k<3;k++) ws_A->solve(ws_b[k],ws_x[k]);
			// transfer to fx fy
			for (int el=0;el<m;el++)
			{
				for (int in1=0;in1<d_i[el].size();in1++)
				{
					int gi=d_i[el][in1];
					lbm.hfx()[gi] += (dreal)(ws_x[0][el]*d_x[el][in1] * 2.0 * lbm.hrho()[gi]);
					lbm.hfy()[gi] += (dreal)(ws_x[1][el]*d_x[el][in1] * 2.0 * lbm.hrho()[gi]);
					lbm.hfz()[gi] += (dreal)(ws_x[2][el]*d_x[el][in1] * 2.0 * lbm.hrho()[gi]);
				}
			}
			break;
		}
		default:
			printf("lagrange_3D: Wu Shu compute flag %d unrecognized.\n", ws_compute);
			break;
	}
}

template< typename LBM >
template< typename... ARGS >
void Lagrange3D<LBM>::log(const char* fmt, ARGS... args)
{
	FILE*f = fopen(logfile,"at"); // append information
	if (f==0)
	{
		printf("unable to create/access file %s",logfile);
		return;
	}
	// insert time stamp
	char tname[FILENAME_CHARS];
	timestamp(tname);
	fprintf(f, "%s ", tname);
	fprintf(f,fmt, args...);
	fprintf(f,"\n");
	fclose(f);

//	printf(fmt, args...);
//	printf("\n");
}

template< typename LBM >
Lagrange3D<LBM>::Lagrange3D(LBM &inputLBM, const char* resultsDir) : lbm(inputLBM)
{
	sprintf(logfile, "%s/ibm_solver.log", resultsDir);

	#ifdef USE_TNL
	ws_tnl_hsolver.setMatrix(ws_tnl_hA);
	ws_tnl_hsolver.setMaxIterations(10000);
	ws_tnl_hsolver.setConvergenceResidue(3e-4);
	ws_tnl_hprecond = std::make_shared< hPreconditioner >();
//	ws_tnl_hsolver.setPreconditioner(ws_tnl_hprecond);
	#ifdef USE_CUDA
	ws_tnl_dsolver.setMatrix(ws_tnl_dA);
	ws_tnl_dsolver.setMaxIterations(10000);
	ws_tnl_dsolver.setConvergenceResidue(3e-4);
	ws_tnl_dprecond = std::make_shared< dPreconditioner >();
//	ws_tnl_dsolver.setPreconditioner(ws_tnl_dprecond);
	#endif
	#endif

	#ifdef USE_CUSPARSE
	ws_M=0;
	ws_MT=0;
	ws_dA=0;
	#endif
	ws_A=0;
	for (int k=0;k<3;k++)
	{
		ws_x[k]=0;
		ws_b[k]=0;
		ws_hx[k]=0;
		ws_hb[k]=0;
		ws_dx[k]=0;
		ws_db[k]=0;
	}
	ws_du=0;
	d_i=0;
	d_x=0;
}

template< typename LBM >
Lagrange3D<LBM>::~Lagrange3D()
{
	if (index_array)
	{
		for (int i=0;i<lag_X;i++) delete [] index_array[i];
		delete [] index_array;
	}
	// WuShu
	if (d_i) delete [] d_i;
	if (d_x) delete [] d_x;

	// WuShu
	if (ws_constructed)
	{
		for (int k=0;k<3;k++)
		{
			if (ws_x[k]) delete [] ws_x[k];
			if (ws_b[k]) delete [] ws_b[k];
			if (ws_hx[k]) delete [] ws_hx[k];
			if (ws_hb[k]) delete [] ws_hb[k];
		}
		#ifdef USE_CUSPARSE
		if (ws_du) cudaFree(ws_du);
		for (int k=0;k<3;k++)
		{
			if (ws_dx[k]) cudaFree(ws_dx[k]);
			if (ws_db[k]) cudaFree(ws_db[k]);
		}
		if (ws_dA) delete ws_dA;
		if (ws_M) delete ws_M;
		if (ws_MT) delete ws_MT;
		#endif
		if (ws_A) delete ws_A;
	}
}
