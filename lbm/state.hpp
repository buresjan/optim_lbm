#include "state.h"

#include "timeutils.h"

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC >
int State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::addLagrange3D()
{
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	FF.emplace_back(lbm, dir);
	return FF.size()-1;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::computeAllLagrangeForces()
{
	for (int i=0;i<FF.size();i++)
		if (FF[i].implicitWuShuForcing)
			FF[i].computeWuShuForcesSparse(lbm.physTime());
}

// input/output log and mark
/// returns current time and date formated for the log files
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
bool State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::fileExists(const char*filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) return false;
    fclose(fp);
    return true;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
bool State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::getPNGdimensions(const char * filename, int &w, int &h)
{
	if (!fileExists(filename)) { printf("file %s does not exist\n",filename); return false; }
	FILE *fp = fopen(filename, "rb");

	png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if(!png) { printf("file %s png read error\n",filename); return false; }

	png_infop info = png_create_info_struct(png);
	if(!png) { printf("file %s png read error\n",filename); return false; }

	if(setjmp(png_jmpbuf(png))) { printf("file %s png read error\n",filename); return false; }

	png_init_io(png, fp);

	png_read_info(png, info);

	w      = png_get_image_width(png, info);
	h     = png_get_image_height(png, info);
	//  color_type = png_get_color_type(png, info);
	//  bit_depth  = png_get_bit_depth(png, info);
	fclose(fp);
	if (w>0 && h>0) return true;
	return false;
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC> 
template< typename... ARGS >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::log(const char* fmt, ARGS... args)
{
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	char fname[FILENAME_CHARS];
	sprintf(fname,"%s/log",dir);

	FILE*f = fopen(fname,"at"); // append information
	if (f==0)
	{
		printf("unable to create/access file %s",fname);
		return;
	}
	// insert time stamp
	char tname[FILENAME_CHARS];
	timestamp(tname);
	fprintf(f, "%s ", tname);
	fprintf(f,fmt, args...);
	fprintf(f,"\n");
	fclose(f);

	printf(fmt, args...);
	printf("\n");

}


/// outputs information into log file "type"
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::setid(const char* fmt, ARGS... args)
{
	sprintf(id, fmt, args...);
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::flagCreate(const char*flagname)
{
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	char fname[FILENAME_CHARS];
	sprintf(fname,"%s/%s",dir,flagname);
	FILE*f = fopen(fname,"at"); // append information
	if (f==0)
	{
		printf("unable to create/access file %s",fname);
		return;
	}
	// insert time stamp
	char tname[FILENAME_CHARS];
	timestamp(tname);
	fprintf(f, "%s\n", tname);
	fclose(f);
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::flagDelete(const char*flagname)
{
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	char fname[FILENAME_CHARS];
	sprintf(fname,"%s/%s",dir,flagname);
	if (fileExists(fname)) remove(fname);
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
bool State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::flagExists(const char*flagname)
{
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	char fname[FILENAME_CHARS];
	sprintf(fname,"%s/%s",dir,flagname);
	return fileExists(fname);
}



template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::mark(const char* fmt, ARGS... args)
{
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	char fname[FILENAME_CHARS];
	sprintf(fname,"%s/mark",dir);
	FILE*f = fopen(fname,"at"); // append information
	if (f==0)
	{
		printf("unable to create/access file %s",fname);
		return;
	}
	// insert time stamp
	char tname[FILENAME_CHARS];
	timestamp(tname);
	fprintf(f, "%s ", tname);
	fprintf(f,fmt, args...);
	fprintf(f,"\n");
	fclose(f);
}

/// checks/creates mark and return status
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
bool State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::isMark()
{
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	char fname[FILENAME_CHARS];
	sprintf(fname,"%s/mark", dir);
	FILE*out=fopen(fname,"rt");
	if (out==0) 
	{
		printf("Mark %s does not exist. Creating new mark.\n",fname);
		out = fopen(fname,"wt");
		fclose(out);
		return false;
	}
	fclose(out);
	printf("Mark %s already exists.\n",fname);
	return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
// VTK SURFACE
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTK_Surface(const char* name, real time, int cycle, T_Lagrange3D &fil)
{
	VTKWriter vtk;
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	sprintf(dir,"results_%s/vtk3D",id);
	mkdir(dir,0777);
        char fname[FILENAME_CHARS];
	sprintf(fname,"%s/%s.vtk",dir,name);
    
        FILE* fp = fopen(fname, "w+");
        vtk.writeHeader(fp);
	
	fprintf(fp, "DATASET POLYDATA\n");
    
	fprintf(fp, "POINTS %d float\n", fil.LL.size());
	for (int i=0;i<fil.LL.size();i++)
	{
		vtk.writeFloat(fp, fil.LL[i].x);
		vtk.writeFloat(fp, fil.LL[i].y);
		vtk.writeFloat(fp, fil.LL[i].z);
	}
	vtk.writeBuffer(fp);

	fprintf(fp, "POLYGONS %d %d\n", (fil.lag_X-1)*fil.lag_Y , 5*(fil.lag_X-1)*fil.lag_Y ); // first number: number of polygons, second number: total integers describing the list
	for (int i=0;i<fil.lag_X-1;i++)
	for (int j=0;j<fil.lag_Y;j++)
	{
		int ip = i+1;
		int jp = (j==fil.lag_Y-1) ? 0 : j+1;
		vtk.writeInt(fp,4);
		vtk.writeInt(fp,fil.findIndex(i,j));
		vtk.writeInt(fp,fil.findIndex(ip,j));
		vtk.writeInt(fp,fil.findIndex(ip,jp));
		vtk.writeInt(fp,fil.findIndex(i,jp));
	}
	fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
// VTK POINTS
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTK_Points(const char* name, real time, int cycle, T_Lagrange3D &fil)
{
	VTKWriter vtk;

	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	sprintf(dir,"results_%s/vtk3D",id);
	mkdir(dir,0777);
        char fname[FILENAME_CHARS];
	sprintf(fname,"%s/%s.vtk",dir,name);
    
        FILE* fp = fopen(fname, "w+");
        vtk.writeHeader(fp);
	
	fprintf(fp, "DATASET POLYDATA\n");
    
	fprintf(fp, "POINTS %d float\n", fil.LL.size());
	for (int i=0;i<fil.LL.size();i++)
	{
		vtk.writeFloat(fp, fil.LL[i].x);
		vtk.writeFloat(fp, fil.LL[i].y);
		vtk.writeFloat(fp, fil.LL[i].z);
	}
	vtk.writeBuffer(fp);
	fclose(fp);
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
// VTK 1D CUT
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::add1Dcut(real fromx, real fromy, real fromz, real tox, real toy, real toz, const char* fmt, ARGS... args)
{
	probe1Dlinevec.push_back( T_PROBE1DLINECUT() );
	int last = probe1Dlinevec.size()-1;
	sprintf(probe1Dlinevec[last].name, fmt, args...);
	probe1Dlinevec[last].from[0] = fromx;
	probe1Dlinevec[last].from[1] = fromy;
	probe1Dlinevec[last].from[2] = fromz;

	probe1Dlinevec[last].to[0] = tox;
	probe1Dlinevec[last].to[1] = toy;
	probe1Dlinevec[last].to[2] = toz;

	probe1Dlinevec[last].cycle = 0;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::add1Dcut_X(real y, real z, const char* fmt, ARGS... args)
{
	probe1Dvec.push_back( T_PROBE1DCUT() );
	int last = probe1Dvec.size()-1;
	sprintf(probe1Dvec[last].name, fmt, args...);
	probe1Dvec[last].type = 0;
	probe1Dvec[last].pos1 = (int)(y/lbm.physDl);
	probe1Dvec[last].pos2 = (int)(z/lbm.physDl);
	probe1Dvec[last].cycle = 0;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::add1Dcut_Y(real x, real z, const char* fmt, ARGS... args)
{
	probe1Dvec.push_back( T_PROBE1DCUT() );
	int last = probe1Dvec.size()-1;
	sprintf(probe1Dvec[last].name, fmt, args...);
	probe1Dvec[last].type = 1;
	probe1Dvec[last].pos1 = (int)(x/lbm.physDl);
	probe1Dvec[last].pos2 = (int)(z/lbm.physDl);
	probe1Dvec[last].cycle = 0;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::add1Dcut_Z(real x, real y, const char* fmt, ARGS... args)
{
	probe1Dvec.push_back( T_PROBE1DCUT() );
	int last = probe1Dvec.size()-1;
	sprintf(probe1Dvec[last].name, fmt, args...);
	probe1Dvec[last].type = 2;
	probe1Dvec[last].pos1 = (int)(x/lbm.physDl);
	probe1Dvec[last].pos2 = (int)(y/lbm.physDl);
	probe1Dvec[last].cycle = 0;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTKs_1D()
{
	if (probe1Dvec.size()>0)
	{
		// check whether dir exists
		char dir[FILENAME_CHARS];
		sprintf(dir,"results_%s",id);
		mkdir(dir,0777);
		sprintf(dir,"results_%s/probes1D",id);
		mkdir(dir,0777);
		// browse all 1D probeline cuts
		for (int i=0;i<probe1Dvec.size(); i++)
		{
			char fname[FILENAME_CHARS];
//			sprintf(fname,"%s/%s_%06d_t%f", dir, probe1Dvec[i].name, probe1Dvec[i].cycle, lbm.physTime());
			sprintf(fname,"%s/%s_%06d", dir, probe1Dvec[i].name, probe1Dvec[i].cycle);
//			probeLine(probe1Dvec[i].from[0],probe1Dvec[i].from[1],probe1Dvec[i].from[2],probe1Dvec[i].to[0],probe1Dvec[i].to[1],probe1Dvec[i].to[2],fname);
			switch (probe1Dvec[i].type)
			{
				case 0: write1Dcut_X(probe1Dvec[i].pos1, probe1Dvec[i].pos2, fname);
					break;
				case 1: write1Dcut_Y(probe1Dvec[i].pos1, probe1Dvec[i].pos2, fname);
					break;
				case 2: write1Dcut_Z(probe1Dvec[i].pos1, probe1Dvec[i].pos2, fname);
					break;
			}
			probe1Dvec[i].cycle++;
		}
	}
	
	// browse all 1D probe cuts
	for (int i=0;i<probe1Dlinevec.size(); i++)
	{
		char dir[FILENAME_CHARS];
		sprintf(dir,"results_%s",id);
		mkdir(dir,0777);
		sprintf(dir,"results_%s/probes1D",id);
		mkdir(dir,0777);
		char fname[FILENAME_CHARS];
//		sprintf(fname,"%s/%s_%06d_t%f", dir, probe1Dvec[i].name, probe1Dvec[i].cycle, lbm.physTime());
		sprintf(fname,"%s/%s_%06d", dir, probe1Dlinevec[i].name, probe1Dlinevec[i].cycle);
		write1Dcut(probe1Dlinevec[i].from[0],probe1Dlinevec[i].from[1],probe1Dlinevec[i].from[2],probe1Dlinevec[i].to[0],probe1Dlinevec[i].to[1],probe1Dlinevec[i].to[2],fname);
		probe1Dlinevec[i].cycle++;
	}
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::write1Dcut_X(idx y, idx z, const char * fname)
{
	FILE*fout = fopen(fname,"wt"); // append information
	log("[probe %s]",fname);
	// probe vertical profile at x_m
	char idd[500];
	real value;
	int dofs;
	fprintf(fout,"#time %f s\n", lbm.physTime());
	fprintf(fout,"#1:rel_pos");
	int count=2, index=0;
	while (outputData(index++, 0, idd, 0,0,0, value, dofs))
	{
		if (dofs==1) fprintf(fout,"\t%d:%s",count++,idd);
		else
		for (idx i=0;i<dofs;i++) fprintf(fout,"\t%d:%s[%d]",count++,idd,(int)i);
	}
	fprintf(fout,"\n");

	for (idx i=0;i<lbm.X;i++)
	{
		if (lbm.isFluid(i,y,z))
		{
			fprintf(fout, "%e",(i-0.5)*lbm.physDl);
			int index=0;
			while (outputData(index++, 0, idd, 0,0,0, value, dofs))
			{
				for (int dof=0;dof<dofs;dof++)
				{
					outputData(index-1,dof,idd,i,y,z,value,dofs);
					fprintf(fout, "\t%e", value);
				}
			}
			fprintf(fout, "\n");
		}
	}
	fclose(fout);
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::write1Dcut_Y(idx x, idx z, const char * fname)
{
	FILE*fout = fopen(fname,"wt"); // append information
	log("[probe %s]",fname);
	// probe vertical profile at x_m
	char idd[500];
	real value;
	int dofs;
	fprintf(fout,"#time %f s\n", lbm.physTime());
	fprintf(fout,"#1:rel_pos");
	int count=2, index=0;
	while (outputData(index++, 0, idd, 0,0,0, value, dofs))
	{
		if (dofs==1) fprintf(fout,"\t%d:%s",count++,idd);
		else
		for (idx i=0;i<dofs;i++) fprintf(fout,"\t%d:%s[%d]",count++,idd,(int)i);
	}
	fprintf(fout,"\n");

	for (idx i=0;i<lbm.Y;i++)
	{
		if (lbm.isFluid(x,i,z))
		{
			fprintf(fout, "%e",(i-0.5)*lbm.physDl);
			int index=0;
			while (outputData(index++, 0, idd, 0,0,0, value, dofs))
			{
				for (int dof=0;dof<dofs;dof++)
				{
					outputData(index-1,dof,idd,x,i,z,value,dofs);
					fprintf(fout, "\t%e", value);
				}
			}
			fprintf(fout, "\n");
		}
	}
	fclose(fout);
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::write1Dcut_Z(idx x, idx y, const char * fname)
{
	FILE*fout = fopen(fname,"wt"); // append information
	log("[probe %s]",fname);
	// probe vertical profile at x_m
	char idd[500];
	real value;
	int dofs;
	fprintf(fout,"#time %f s\n", lbm.physTime());
	fprintf(fout,"#1:rel_pos");
	int count=2, index=0;
	while (outputData(index++, 0, idd, 0,0,0, value, dofs))
	{
		if (dofs==1) fprintf(fout,"\t%d:%s",count++,idd);
		else
		for (idx i=0;i<dofs;i++) fprintf(fout,"\t%d:%s[%d]",count++,idd,(int)i);
	}
	fprintf(fout,"\n");

	for (idx i=0;i<lbm.Z;i++)
	{
		if (lbm.isFluid(x,y,i))
		{
			fprintf(fout, "%e",(i-0.5)*lbm.physDl);
			index=0;
			while (outputData(index++, 0, idd, 0,0,0, value, dofs))
			{
				for (int dof=0;dof<dofs;dof++)
				{
					outputData(index-1,dof,idd,x,y,i,value,dofs);
					fprintf(fout, "\t%e", value);
				}
			}
			fprintf(fout, "\n");
		}
	}
	fclose(fout);
}


// line projection from[3] to[3] 
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::write1Dcut(real fromx, real fromy, real fromz, real tox, real toy, real toz, const char * fname)
{
//		char dir[FILENAME_CHARS];
//		sprintf(dir,"results_%s",id);
//		mkdir(dir,0777);
//		sprintf(dir,"results_%s/probes1D",id);
//		mkdir(dir,0777);
//		char fname[FILENAME_CHARS];
//		sprintf(fname,"%s/%s_it%08d_t%f",dir,desc,lbm.iterations,lbm.physTime());
	FILE*fout = fopen(fname,"wt"); // append information
	log("[probe %s]",fname);
	// probe vertical profile at x_m
	real i[3],f[3],p[3];
	i[0]=fromx/lbm.physDl;
	i[1]=fromy/lbm.physDl;
	i[2]=fromz/lbm.physDl;
	f[0]=tox/lbm.physDl;
	f[1]=toy/lbm.physDl;
	f[2]=toz/lbm.physDl;
	real dist = NORM(i[0]-f[0],i[1]-f[1],i[2]-f[2]);
	real ds=1.0/(dist*2.0); // rozliseni najit
	
	char idd[500];
	real value;
	int dofs;
	fprintf(fout,"#time %f s\n", lbm.physTime());
	fprintf(fout,"#1:rel_pos");
	
	int count=2, index=0;
	while (outputData(index++, 0, idd, 0,0,0, value, dofs))
	{
		if (dofs==1) fprintf(fout,"\t%d:%s",count++,idd);
		else
		for (int i=0;i<dofs;i++) fprintf(fout,"\t%d:%s[%d]",count++,idd,i);
	}
	fprintf(fout,"\n");
	
	for (real s=0;s<=1.0;s+=ds)
	{
		for (int k=0;k<3;k++) p[k] = i[k] + s*(f[k]-i[k]);
		if (lbm.isFluid(p[0],p[1],p[2]))
		{
			fprintf(fout, "%e",s*dist*lbm.physDl);
			index=0;
			while (outputData(index++, 0, idd, 0,0,0, value, dofs))
			{
				for (int dof=0;dof<dofs;dof++)
				{
					outputData(index-1,dof,idd,p[0],p[1],p[2],value,dofs);
					fprintf(fout, "\t%e", value);
				}
			}
			fprintf(fout, "\n");
		}
	}
	fclose(fout);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
// VTK 3D
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTKs_3D()
{
	// only one 3D vtk is written now
	char dir1[250], dir2[250], filename[250];
	sprintf(dir1,"results_%s", id);
//		sprintf(dir2,"%s/vtk_%s",dir1,State::T_LBM_TYPE::id);
	sprintf(dir2,"%s/vtk3D",dir1);
	mkdir(dir1,0755);
	mkdir(dir2,0755);
	
//	int vtk3Dstyle = vtk3DsingleFile; // enum { vtk3DsingleFile, vtk3DmanyFiles, vtk3DmanyFilesExtraHeader };
	if (vtk3Dstyle == vtk3DsingleFile)
	{
		sprintf(filename,"%s/data_%d.vtk",dir2, cnt[VTK3D].count);
		writeVTK_3D_singlefile(filename,lbm.physTime(),cnt[VTK3D].count);
	} else
	{
		sprintf(filename,"%s/",dir2);//, cnt[VTK3D].count);
		writeVTK_3D(filename,lbm.physTime(),cnt[VTK3D].count);
	}
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTK_3D(const char* name, real time, int cycle)
{
	bool VTK3D_write_separate_header=(vtk3Dstyle==vtk3DmanyFiles) ? false : true;
	// determine max objects to write
	int max_objects=0;
	{
		int index=0;
		char idd[500];
		real value;
		int dofs;
		while (outputData(index++, 0, idd, 0,0,0, value, dofs)) max_objects++;
		log("[vtk preparing to write %d objects]",max_objects);
	}

//	state.log("number of CPU cores:\t%d", omp_get_num_procs());
//	omp_set_num_threads(num_cpu_threads);
	int max_thr = MIN(max_objects, omp_get_num_procs());
	log("[vtk3d %d threads for vtk 3D write]", max_thr);
	
	if (VTK3D_write_separate_header)
	{
		VTKWriter vtk; // local vtkwriter
		char fname[500];
		sprintf(fname,"%sheader_%d.prevtk",name,cycle);
		// browse output
		FILE* fp = fopen(fname, "w+");
		vtk.writeHeader(fp);
		fprintf(fp,"DATASET RECTILINEAR_GRID\n");
		fprintf(fp,"DIMENSIONS %d %d %d\n", (int)lbm.X, (int)lbm.Y, (int)lbm.Z);
		fprintf(fp,"X_COORDINATES %d float\n", (int)lbm.X);
		for (idx x = 0; x < lbm.X; x++) vtk.writeFloat(fp, x*lbm.physDl);
		vtk.writeBuffer(fp);
	
		fprintf(fp,"Y_COORDINATES %d float\n", (int)lbm.Y);
		for (idx y = 0; y < lbm.Y; y++) vtk.writeFloat(fp, y*lbm.physDl);
		vtk.writeBuffer(fp);
	
		fprintf(fp,"Z_COORDINATES %d float\n", (int)lbm.Z);
		for (idx z = 0; z < lbm.Z; z++) vtk.writeFloat(fp, z*lbm.physDl);
		vtk.writeBuffer(fp);
	
		fprintf(fp,"FIELD FieldData %d\n",2);
		fprintf(fp,"TIME %d %d float\n",1,1);
		vtk.writeFloat(fp, time);
		vtk.writeBuffer(fp);
	
		fprintf(fp,"CYCLE %d %d float\n",1,1);
		vtk.writeFloat(fp, cycle);
		vtk.writeBuffer(fp);
	
		fprintf(fp,"POINT_DATA %d\n", (int)(lbm.X*lbm.Y*lbm.Z));

		fprintf(fp,"SCALARS wall int 1\n");
		fprintf(fp,"LOOKUP_TABLE default\n");
		for (idx z = 0; z < lbm.Z; z++)
		for (idx y = 0; y < lbm.Y; y++) 
		for (idx x = 0; x < lbm.X; x++) 
			vtk.writeInt(fp, lbm.map(x,y,z));

		fclose(fp);
	}
	
	#pragma omp parallel for schedule(static) num_threads(max_thr)
	for (int i=0;i<max_objects;i++)
	{
		VTKWriter vtk; // local vtkwriter
		char idd[500];
		real value;
		int dofs;
//		int index=0;
		char fname[500];
		outputData(i, 0, idd, 0,0,0, value, dofs); // read idd and dofs
//		{
		sprintf(fname,"%s%s_%d.%s",name,idd,cycle,(VTK3D_write_separate_header)?"prevtk":"vtk");
		log("[vtk3d: writing %s, time %f, cycle %d] ", fname, time, cycle);
		// browse output
		FILE* fp = fopen(fname, "w+");
		if (!VTK3D_write_separate_header)
		{
			vtk.writeHeader(fp);
			fprintf(fp,"DATASET RECTILINEAR_GRID\n");
			fprintf(fp,"DIMENSIONS %d %d %d\n", (int)lbm.X, (int)lbm.Y, (int)lbm.Z);
			fprintf(fp,"X_COORDINATES %d float\n", (int)lbm.X);
			for (idx x = 0; x < lbm.X; x++) vtk.writeFloat(fp, x*lbm.physDl);
			vtk.writeBuffer(fp);
	
			fprintf(fp,"Y_COORDINATES %d float\n", (int)lbm.Y);
			for (idx y = 0; y < lbm.Y; y++) vtk.writeFloat(fp, y*lbm.physDl);
			vtk.writeBuffer(fp);
	
			fprintf(fp,"Z_COORDINATES %d float\n", (int)lbm.Z);
			for (idx z = 0; z < lbm.Z; z++) vtk.writeFloat(fp, z*lbm.physDl);
			vtk.writeBuffer(fp);
	
			fprintf(fp,"FIELD FieldData %d\n",2);
			fprintf(fp,"TIME %d %d float\n",1,1);
			vtk.writeFloat(fp, time);
			vtk.writeBuffer(fp);
	
			fprintf(fp,"CYCLE %d %d float\n",1,1);
			vtk.writeFloat(fp, cycle);
			vtk.writeBuffer(fp);
	
			fprintf(fp,"POINT_DATA %d\n", (int)(lbm.X*lbm.Y*lbm.Z));

			fprintf(fp,"SCALARS wall int 1\n");
			fprintf(fp,"LOOKUP_TABLE default\n");
			for (idx z = 0; z < lbm.Z; z++)
			for (idx y = 0; y < lbm.Y; y++) 
			for (idx x = 0; x < lbm.X; x++) 
				vtk.writeInt(fp, lbm.map(x,y,z));
		}

		// insert description
		if (dofs==1)
		{
			fprintf(fp,"SCALARS %s float 1\n",idd);
			fprintf(fp,"LOOKUP_TABLE default\n");
		}
		else
		{
			fprintf(fp,"VECTORS %s float\n",idd);
		}

		for (idx z = 0; z < lbm.Z; z++)
		for (idx y = 0; y < lbm.Y; y++)
		for (idx x = 0; x < lbm.X; x++) 
		{
			// determine triangle center
			for (int dof=0;dof<dofs;dof++)
			{
				outputData(i,dof,idd,x,y,z,value,dofs);
				vtk.writeFloat(fp, value);
			}
		}
		vtk.writeBuffer(fp);
		fclose(fp);
//		log("[vtk3d %s written, time %f, cycle %d] ", fname, time, cycle);
	}
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTK_3D_singlefile(const char* name, real time, int cycle)
{
	VTKWriter vtk;

	FILE* fp = fopen(name, "w+");
	vtk.writeHeader(fp);
	fprintf(fp,"DATASET RECTILINEAR_GRID\n");
	fprintf(fp,"DIMENSIONS %d %d %d\n", (int)lbm.X, (int)lbm.Y, (int)lbm.Z);
	fprintf(fp,"X_COORDINATES %d float\n", (int)lbm.X);
	for (idx x = 0; x < lbm.X; x++) vtk.writeFloat(fp, x*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"Y_COORDINATES %d float\n", (int)lbm.Y);
	for (idx y = 0; y < lbm.Y; y++) vtk.writeFloat(fp, y*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"Z_COORDINATES %d float\n", (int)lbm.Z);
	for (idx z = 0; z < lbm.Z; z++) vtk.writeFloat(fp, z*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"FIELD FieldData %d\n",2);
	fprintf(fp,"TIME %d %d float\n",1,1);
	vtk.writeFloat(fp, time);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"CYCLE %d %d float\n",1,1);
	vtk.writeFloat(fp, cycle);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"POINT_DATA %d\n", (int)(lbm.X*lbm.Y*lbm.Z));
	
	
	fprintf(fp,"SCALARS wall int 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");
	for (idx z = 0; z < lbm.Z; z++)
	for (idx y = 0; y < lbm.Y; y++) 
	for (idx x = 0; x < lbm.X; x++) 
		vtk.writeInt(fp, lbm.map(x,y,z));
	

	char idd[500];
	real value;
	int dofs;
	int count=0, index=0;
	while (outputData(index++, 0, idd, 0,0,0, value, dofs))
	{
		// insert description
		if (dofs==1)
		{
			fprintf(fp,"SCALARS %s float 1\n",idd);
			fprintf(fp,"LOOKUP_TABLE default\n");
		}
		else
			fprintf(fp,"VECTORS %s float\n",idd);

		for (idx z = 0; z < lbm.Z; z++)
		for (idx y = 0; y < lbm.Y; y++)
		for (idx x = 0; x < lbm.X; x++) 
		{
			// determine triangle center
			for (int dof=0;dof<dofs;dof++)
			{
				outputData(index-1,dof,idd,x,y,z,value,dofs);
				vtk.writeFloat(fp, value);
			}
		}
		vtk.writeBuffer(fp);
		count++;
	}
        
	fclose(fp);
	log("[vtk %s written, time %f, cycle %d] ", name, time, cycle);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
// VTK 3D CUT
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTK_3Dcut(const char* name, real time, int cycle, idx ox, idx oy, idx oz, idx lx, idx ly, idx lz, idx step)
{
	VTKWriter vtk;
	
	// local dimensions
	idx X = floor(lx/step);
	idx Y = floor(ly/step);
	idx Z = floor(lz/step);
//	log("debug: writeVTK3Dcut X %d Y %d Z %d",(int)X,(int)Y,(int)Z);

	FILE* fp = fopen(name, "w+");
	vtk.writeHeader(fp);
	fprintf(fp,"DATASET RECTILINEAR_GRID\n");
	fprintf(fp,"DIMENSIONS %d %d %d\n", (int)X, (int)Y, (int)Z);
	fprintf(fp,"X_COORDINATES %d float\n", (int)X);
	for (idx x = 0; x < X; x++) vtk.writeFloat(fp, (ox + x*step)*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"Y_COORDINATES %d float\n", (int)Y);
	for (idx y = 0; y < Y; y++) vtk.writeFloat(fp, (oy + y*step)*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"Z_COORDINATES %d float\n", (int)Z);
	for (idx z = 0; z < Z; z++) vtk.writeFloat(fp, (oz + z*step)*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"FIELD FieldData %d\n",2);
	fprintf(fp,"TIME %d %d float\n",1,1);
	vtk.writeFloat(fp, time);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"CYCLE %d %d float\n",1,1);
	vtk.writeFloat(fp, cycle);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"POINT_DATA %d\n", (int)(X*Y*Z));
	
	
	fprintf(fp,"SCALARS wall int 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");
	for (idx z = 0; z < Z; z++)
	for (idx y = 0; y < Y; y++) 
	for (idx x = 0; x < X; x++) 
		vtk.writeInt(fp, lbm.map(ox+x*step,oy+y*step,oz+z*step));
	

	char idd[500];
	real value;
	int dofs;
	int count=0, index=0;
	while (outputData(index++, 0, idd, 0,0,0, value, dofs))
	{
		// insert description
		if (dofs==1)
		{
			fprintf(fp,"SCALARS %s float 1\n",idd);
			fprintf(fp,"LOOKUP_TABLE default\n");
		}
		else
			fprintf(fp,"VECTORS %s float\n",idd);

		for (idx z = 0; z < Z; z++)
		for (idx y = 0; y < Y; y++)
		for (idx x = 0; x < X; x++) 
		{
			// determine triangle center
			for (int dof=0;dof<dofs;dof++)
			{
				outputData(index-1,dof,idd,ox+x*step,oy+y*step,oz+z*step,value,dofs);
				vtk.writeFloat(fp, value);
			}
		}
		vtk.writeBuffer(fp);
		count++;
	}
        
	fclose(fp);
	log("[vtk %s written, time %f, cycle %d] ", name, time, cycle);
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::add3Dcut(idx ox, idx oy, idx oz, idx lx, idx ly, idx lz, idx step, const char* fmt, ARGS... args)
{
	probe3Dvec.push_back( T_PROBE3DCUT() );
	int last = probe3Dvec.size()-1;

	sprintf(probe3Dvec[last].name, fmt, args...);

	probe3Dvec[last].ox = ox;
	probe3Dvec[last].oy = oy;
	probe3Dvec[last].oz = oz;
	probe3Dvec[last].lx = lx;
	probe3Dvec[last].ly = ly;
	probe3Dvec[last].lz = lz;
	probe3Dvec[last].step = step;
	probe3Dvec[last].cycle = 0;
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTKs_3Dcut()
{
	if (probe3Dvec.size()<=0) return;
	// check whether dir exists
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0755);
	sprintf(dir,"results_%s/vtk3Dcut",id);
	mkdir(dir,0755);
	// browse all 3D vtk cuts
	for (int i=0;i<probe3Dvec.size(); i++)
	{
		char fname[FILENAME_CHARS];
		sprintf(fname,"%s/%s_%d.vtk", dir, probe3Dvec[i].name, probe3Dvec[i].cycle);
		writeVTK_3Dcut(
			fname, 
			lbm.physTime(), 
			probe3Dvec[i].cycle, 
			probe3Dvec[i].ox,
			probe3Dvec[i].oy,
			probe3Dvec[i].oz,
			probe3Dvec[i].lx,
			probe3Dvec[i].ly,
			probe3Dvec[i].lz,
			probe3Dvec[i].step
			);
		probe3Dvec[i].cycle++;
	}

// 	// only one 3D vtk is written now
// 	char dir1[250], dir2[250], filename[250];
// 	sprintf(dir1,"results_%s", id);
// //		sprintf(dir2,"%s/vtk_%s",dir1,State::T_LBM_TYPE::id);
// 	sprintf(dir2,"%s/vtk3Dcut",dir1);
// 	mkdir(dir1,0755);
// 	mkdir(dir2,0755);
// 	
// 	sprintf(filename,"%s/data_%d.vtk",dir2, cnt[VTK3DCUT].count);
// 	writeVTK_3D_singlefile(filename,lbm.physTime(),cnt[VTK3DCUT].count);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
// VTK 2D CUT
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::add2Dcut_X(idx x, const char* fmt, ARGS... args)
{
	probe2Dvec.push_back( T_PROBE2DCUT() );
	int last = probe2Dvec.size()-1;

	sprintf(probe2Dvec[last].name, fmt, args...);

	probe2Dvec[last].type = 0;
	probe2Dvec[last].cycle = 0;
	probe2Dvec[last].position = x;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::add2Dcut_Y(idx y, const char* fmt, ARGS... args)
{
	probe2Dvec.push_back( T_PROBE2DCUT() );
	int last = probe2Dvec.size()-1;

	sprintf(probe2Dvec[last].name, fmt, args...);

	probe2Dvec[last].type = 1;
	probe2Dvec[last].cycle = 0;
	probe2Dvec[last].position = y;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS >
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::add2Dcut_Z(idx z, const char* fmt, ARGS... args)
{
	probe2Dvec.push_back( T_PROBE2DCUT() );
	int last = probe2Dvec.size()-1;

	sprintf(probe2Dvec[last].name, fmt, args...);
	
	probe2Dvec[last].type = 2;
	probe2Dvec[last].cycle = 0;
	probe2Dvec[last].position = z;
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTKs_2D()
{
	if (probe2Dvec.size()<=0) return;
	// check whether dir exists
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	sprintf(dir,"results_%s/vtk2D",id);
	mkdir(dir,0777);
	// browse all 2D vtk cuts
	for (int i=0;i<probe2Dvec.size(); i++)
	{
		char fname[FILENAME_CHARS];
		sprintf(fname,"%s/%s_%d.vtk", dir, probe2Dvec[i].name, probe2Dvec[i].cycle);
		switch (probe2Dvec[i].type)
		{
			case 0: writeVTK_2DcutX(fname,lbm.physTime(),probe2Dvec[i].cycle,probe2Dvec[i].position);
				break;
			case 1: writeVTK_2DcutY(fname,lbm.physTime(),probe2Dvec[i].cycle,probe2Dvec[i].position);
				break;
			case 2: writeVTK_2DcutZ(fname,lbm.physTime(),probe2Dvec[i].cycle,probe2Dvec[i].position);
				break;
		}
		probe2Dvec[i].cycle++;
	}
}

// X-Z plane for Y=YPOS 
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTK_2DcutY(const char* name, real time, int cycle, idx YPOS)
{
	VTKWriter vtk;

	FILE* fp = fopen(name, "w+");
	vtk.writeHeader(fp);
	fprintf(fp,"DATASET RECTILINEAR_GRID\n");
	fprintf(fp,"DIMENSIONS %d %d %d\n", (int)lbm.X, 1, (int)lbm.Z);
	fprintf(fp,"X_COORDINATES %d float\n", (int)lbm.X);
	for (idx x = 0; x < lbm.X; x++) vtk.writeFloat(fp, x*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"Y_COORDINATES 1 float\n");
	vtk.writeFloat(fp, YPOS*lbm.physDl);
	vtk.writeBuffer(fp);

	fprintf(fp,"Z_COORDINATES %d float\n", (int)lbm.Z);
	for (idx z = 0; z < lbm.Z; z++) vtk.writeFloat(fp, z*lbm.physDl);
	vtk.writeBuffer(fp);
	
	
	fprintf(fp,"FIELD FieldData %d\n",2);
	fprintf(fp,"TIME %d %d float\n",1,1);
	vtk.writeFloat(fp, time);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"CYCLE %d %d float\n",1,1);
	vtk.writeFloat(fp, cycle);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"POINT_DATA %d\n", (int)(1*lbm.X*lbm.Z));
	
	fprintf(fp,"SCALARS wall int 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");
	idx y=YPOS;
	for (idx z = 0; z < lbm.Z; z++)
	for (idx x = 0; x < lbm.X; x++) 
		vtk.writeInt(fp, lbm.map(x,y,z));
	

	int count=0, index=0;
	char idd[500];
	real value;
	int dofs;
	while (outputData(index++, 0, idd, 0,0,0, value, dofs))
	{
		// insert description
		if (dofs==1)
		{
			fprintf(fp,"SCALARS %s float 1\n",idd);
			fprintf(fp,"LOOKUP_TABLE default\n");
		} else
			fprintf(fp,"VECTORS %s float\n",idd);

		for (idx z = 0; z < lbm.Z; z++)
		for (idx x = 0; x < lbm.X; x++) 
		{
			// determine triangle center
			for (int dof=0;dof<dofs;dof++)
			{
				outputData(index-1,dof,idd,x,y,z,value,dofs);
				vtk.writeFloat(fp, value);
			}
		}
		vtk.writeBuffer(fp);
		count++;
		
	}

	fclose(fp);
	log("[vtk %s written, time %f, cycle %d] ", name, time, cycle);
}

// Y-Z plane for X=XPOS 
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTK_2DcutX(const char* name, real time, int cycle, idx XPOS)
{
	VTKWriter vtk;

	FILE* fp = fopen(name, "w+");
	vtk.writeHeader(fp);
	fprintf(fp,"DATASET RECTILINEAR_GRID\n");
	fprintf(fp,"DIMENSIONS %d %d %d\n",1, (int)lbm.Y, (int)lbm.Z);

	fprintf(fp,"X_COORDINATES %d float\n", 1);
	vtk.writeFloat(fp, XPOS*lbm.physDl);
	vtk.writeBuffer(fp);

	fprintf(fp,"Y_COORDINATES %d float\n", (int)lbm.Y);
	for (idx y = 0; y < lbm.Y; y++) vtk.writeFloat(fp, y*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"Z_COORDINATES %d float\n", (int)lbm.Z);
	for (idx z = 0; z < lbm.Z; z++) vtk.writeFloat(fp, z*lbm.physDl);
	vtk.writeBuffer(fp);
	
	
	fprintf(fp,"FIELD FieldData %d\n",2);
	fprintf(fp,"TIME %d %d float\n",1,1);
	vtk.writeFloat(fp, time);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"CYCLE %d %d float\n",1,1);
	vtk.writeFloat(fp, cycle);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"POINT_DATA %d\n", (int)(1*lbm.Y*lbm.Z));
	
	fprintf(fp,"SCALARS wall int 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");
	idx x=XPOS;
	for (idx z = 0; z < lbm.Z; z++)
	for (idx y = 0; y < lbm.Y; y++) 
		vtk.writeInt(fp, lbm.map(x,y,z));
	

	int count=0, index=0;
	char idd[500];
	real value;
	int dofs;
	while (outputData(index++, 0, idd, 0,0,0, value, dofs))
	{
		// insert description
		if (dofs==1)
		{
			fprintf(fp,"SCALARS %s float 1\n",idd);
			fprintf(fp,"LOOKUP_TABLE default\n");
		} else
		{
			fprintf(fp,"VECTORS %s float\n",idd);
		}
		for (idx z = 0; z < lbm.Z; z++)
		for (idx y = 0; y < lbm.Y; y++)
		{
			// determine triangle center
			for (int dof=0;dof<dofs;dof++)
			{
				outputData(index-1,dof,idd,x,y,z,value,dofs);
				vtk.writeFloat(fp, value);
			}
		}
		vtk.writeBuffer(fp);
		count++;
	}

	fclose(fp);
	log("[vtk %s written, time %f, cycle %d] ", name, time, cycle);
}

// X-Y plane for Z=ZPOS 
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::writeVTK_2DcutZ(const char* name, real time, int cycle, idx ZPOS)
{
	VTKWriter vtk;

	FILE* fp = fopen(name, "w+");
	vtk.writeHeader(fp);
	fprintf(fp,"DATASET RECTILINEAR_GRID\n");
	fprintf(fp,"DIMENSIONS %d %d %d\n", (int)lbm.X, (int)lbm.Y, 1);
	fprintf(fp,"X_COORDINATES %d float\n", (int)lbm.X);
	for (idx x = 0; x < lbm.X; x++) vtk.writeFloat(fp, x*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"Y_COORDINATES %d float\n", (int)lbm.Y);
	for (idx y = 0; y < lbm.Y; y++) vtk.writeFloat(fp, y*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"Z_COORDINATES %d float\n", 1);
	vtk.writeFloat(fp, ZPOS*lbm.physDl);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"FIELD FieldData %d\n",2);
	fprintf(fp,"TIME %d %d float\n",1,1);
	vtk.writeFloat(fp, time);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"CYCLE %d %d float\n",1,1);
	vtk.writeFloat(fp, cycle);
	vtk.writeBuffer(fp);
	
	fprintf(fp,"POINT_DATA %d\n", (int)(1*lbm.X*lbm.Y));
	
	fprintf(fp,"SCALARS wall int 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");
	idx z=ZPOS;
	for (idx y = 0; y < lbm.Y; y++) 
	for (idx x = 0; x < lbm.X; x++) 
		vtk.writeInt(fp, lbm.map(x,y,z));
	

	int count=0, index=0;
	char idd[500];
	real value;
	int dofs;
	while (outputData(index++, 0, idd, 0,0,0, value, dofs))
	{
		// insert description
		if (dofs==1)
		{
			fprintf(fp,"SCALARS %s float 1\n",idd);
			fprintf(fp,"LOOKUP_TABLE default\n");
		} else
			fprintf(fp,"VECTORS %s float\n",idd);
		
		for (idx y = 0; y < lbm.Y; y++)
		for (idx x = 0; x < lbm.X; x++) 
		{
		// determine triangle center
			for (int dof=0;dof<dofs;dof++)
			{
				outputData(index-1,dof,idd,x,y,z,value,dofs);
				vtk.writeFloat(fp, value);
			}
		}
		vtk.writeBuffer(fp);
		count++;
	}
        
	fclose(fp);
	log("[vtk %s written, time %f, cycle %d] ", name, time, cycle);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
// PNG PROJECTION
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
bool State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::projectPNG_X(const char * filename, idx x0, int rotate, bool mirror, bool flip)
{
	if (!fileExists(filename)) { printf("file %s does not exist\n",filename); return false; }
	PNGTool *P = new PNGTool(filename);
	
	real a,b;
	
	// plane y-z
	idx x=x0;
	for (idx z=0;z<lbm.Z;z++)
	{
		a = (real)z/(real)(lbm.Z-1); // a in [0,1]
		if (mirror) a = 1.0-a;
		for (idx y=0;y<lbm.Y;y++)
		{
			b = (real)y/(real)(lbm.Y-1); // b in [0,1]
			if (flip) b=1.0-b;
			switch (rotate)
			{
				case 1:
					if (P->intensity(b,a) > 0) lbm.defineWall(x, y, z, true);
				break;
				default:
					if (P->intensity(a,b) > 0) lbm.defineWall(x, y, z, true);
			}
		}
	}
	delete P;
	return true;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
bool State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::projectPNG_Y(const char * filename, idx y0, int rotate, bool mirror, bool flip)
{
	if (!fileExists(filename)) { printf("file %s does not exist\n",filename); return false; }
	PNGTool *P = new PNGTool(filename);
	
	real a,b;
	
	// plane x-z
	idx y=y0;
	for (idx z=0;z<lbm.Z;z++)
	{
		a = (real)z/(real)(lbm.Z-1); // a in [0,1]
		if (mirror) a = 1.0-a;
		for (idx x=0;x<lbm.X;x++)
		{
			b = (real)x/(real)(lbm.X-1); // b in [0,1]
			if (flip) b=1.0-b;
			switch (rotate)
			{
				case 1:
					if (P->intensity(b,a) > 0) lbm.defineWall(x, y, z, true);
				break;
				default:
					if (P->intensity(a,b) > 0) lbm.defineWall(x, y, z, true);
			}
		}
	}
	delete P;
	return true;
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
bool State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::projectPNG_Z(const char * filename, idx z0, int rotate, bool mirror, bool flip)
{
	if (!fileExists(filename)) { printf("file %s does not exist\n",filename); return false; }
	PNGTool *P = new PNGTool(filename);
	
	real a,b;
	
	// plane x-y
	idx z=z0;
	for (idx x=0;x<lbm.X;x++)
	{
		a = (real)x/(real)(lbm.X-1); // a in [0,1]
		if (mirror) a = 1.0-a;
		for (idx y=0;y<lbm.Y;y++)
		{
			b = (real)y/(real)(lbm.Y-1); // b in [0,1]
			if (flip) b=1.0-b;
			switch (rotate)
			{
				case 1:
					if (P->intensity(b,a) > 0) lbm.defineWall(x, y, z, true);
				break;
				default:
					if (P->intensity(a,b) > 0) lbm.defineWall(x, y, z, true);
			}
		}
	}
	delete P;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
// SAVE & LOAD STATE
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// old version
// template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
// template< typename... ARGS >
// int State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::saveLoadTextData(int direction, const char*subdirname, const char*filename, const char*fmt, ARGS&... args)
// {
// 	// check if main dir exists
// 	char dir[FILENAME_CHARS];
// 	sprintf(dir,"results_%s",id);
// 	mkdir(dir,0777);
// 	char subdir[FILENAME_CHARS];
// 	sprintf(subdir,"%s/%s", dir, subdirname);
// 	mkdir(subdir,0777);
// 	char fname[FILENAME_CHARS];
// 	sprintf(fname,"%s/%s", subdir, filename);
// 
// 	if (direction==MemoryToFile)
// 	{
// 		FILE*f = fopen(fname,"wt");
// 		if (f==0)
// 		{
// 			log("unable to create file %s",fname);
// 			return 0;
// 		}
// 		fprintf(f,fmt, args...);
// 		fclose(f);
// 		log("[saveLoadTextData: saved data into %s]",fname);
// 	}
// 	if (direction==FileToMemory)
// 	{
// 		FILE*f = fopen(fname,"rt");
// 		if (f==0)
// 		{
// 			log("unable to access file %s",fname);
// 			return 0;
// 		}
// 		fscanf(f,fmt, &args...);
// 		fclose(f);
// 		log("[saveLoadTextData: read data from %s]",fname);
// 	}
// 	return 1;
// }

// new version
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename... ARGS >
int State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::saveLoadTextData(int direction, const char*subdirname, const char*filename, ARGS&... args)
{
	// check if main dir exists
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	char subdir[FILENAME_CHARS];
	sprintf(subdir,"%s/%s", dir, subdirname);
	mkdir(subdir,0777);
	char fname[FILENAME_CHARS];
	sprintf(fname,"%s/%s", subdir, filename);

	const std::string fmt = getSaveLoadFmt(args...);

	if (direction==MemoryToFile)
	{
		FILE*f = fopen(fname,"wt");
		if (f==0)
		{
			log("unable to create file %s",fname);
			return 0;
		}
		fprintf(f,fmt.c_str(), args...);
		fclose(f);
		log("[saveLoadTextData: saved data into %s]",fname);
	}
	if (direction==FileToMemory)
	{
		FILE*f = fopen(fname,"rt");
		if (f==0)
		{
			log("unable to access file %s",fname);
			return 0;
		}
		fscanf(f,fmt.c_str(), &args...);
		fclose(f);
		log("[saveLoadTextData: read data from %s]",fname);
	}
	return 1;
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
template< typename VARTYPE >
int State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::saveloadBinaryData(int direction, const char*subdirname, const char*filename, VARTYPE*data, idx length)
{
	// check if main dir exists
	char dir[FILENAME_CHARS];
	sprintf(dir,"results_%s",id);
	mkdir(dir,0777);
	char subdir[FILENAME_CHARS];
	sprintf(subdir,"%s/%s", dir, subdirname);
	mkdir(subdir,0777);
	char fname[FILENAME_CHARS];
	sprintf(fname,"%s/%s", subdir, filename);


	if (direction==MemoryToFile)
	{
		FILE*f = fopen(fname,"wb");
		if (f==0)
		{
			log("unable to create file %s",fname);
			return 0;
		}
		fwrite(data, sizeof(VARTYPE), length, f);
		fclose(f);
		log("[saveLoadBinaryData: saved data into %s]",fname);
	}
	if (direction==FileToMemory)
	{
		FILE*f = fopen(fname,"rb");
		if (f==0)
		{
			log("unable to access file %s",fname);
			return 0;
		}
		fread(data, sizeof(VARTYPE), length, f);
		fclose(f);
		log("[saveLoadTBinaryData: read data from %s]",fname);
	}
	return 1;
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::saveAndLoadState(int direction, const char*subdirname)
{
	char nid[200];
	//saveLoadTextData(direction, subdirname, "config", "%d\n%d\n%d\n%d\n%le\n", lbm.iterations, lbm.X, lbm.Y, lbm.Z, lbm.physFinalTime);
	saveLoadTextData(direction, subdirname, "config", lbm.iterations, lbm.X, lbm.Y, lbm.Z, lbm.physFinalTime);
	// save all counter states
	for (int c=0;c<MAX_COUNTER;c++)
	{
		sprintf(nid,"cnt_%d",c);
//		saveLoadTextData(direction, subdirname, nid, "%d\n%le\n", cnt[c].count, cnt[c].period);
		saveLoadTextData(direction, subdirname, nid, cnt[c].count, cnt[c].period);
	}
	// save probes
	for (int i=0;i<probe1Dvec.size();i++)
	{
		sprintf(nid,"probe1D_%d",i);
//		saveLoadTextData(direction, subdirname, nid, "%d\n", probe1Dvec[i].cycle);
		saveLoadTextData(direction, subdirname, nid, probe1Dvec[i].cycle);
	}
	for (int i=0;i<probe1Dlinevec.size();i++)
	{
		sprintf(nid,"probe1Dline_%d",i);
//		saveLoadTextData(direction, subdirname, nid, "%d\n", probe1Dlinevec[i].cycle);
		saveLoadTextData(direction, subdirname, nid, probe1Dlinevec[i].cycle);
	}
	for (int i=0;i<probe2Dvec.size();i++)
	{
		sprintf(nid,"probe2D_%d",i);
//		saveLoadTextData(direction, subdirname, nid, "%d\n", probe2Dvec[i].cycle);
		saveLoadTextData(direction, subdirname, nid, probe2Dvec[i].cycle);
	}
	// save DFs
	for (int dfty=0;dfty<DFMAX;dfty++)
	{
		sprintf(nid,"df_%d",dfty);
		saveloadBinaryData(direction, subdirname, nid, lbm.hfs[dfty], 27*lbm.X*lbm.Y*lbm.Z);
	}
	// save map
	sprintf(nid,"map");
	saveloadBinaryData(direction, subdirname, nid, lbm.hmap, lbm.X*lbm.Y*lbm.Z);
	// save macro
	if (MACRO::N>0)
	{
		sprintf(nid,"macro");
		saveloadBinaryData(direction, subdirname, nid, lbm.hmacro, MACRO::N*lbm.X*lbm.Y*lbm.Z);
	}
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::saveState(bool forced)
{
//		flagCreate("do_save_state");
    if (flagExists("savestate") || !check_savestate_flag || forced)
    {
        log("[saveState invoked]");
        saveAndLoadState(MemoryToFile, "current_state");
        if (delete_savestate_flag && !forced)
        {
            flagDelete("savestate");
//				flagRename("savestate","savestate_done");
            flagCreate("savestate_done");
        }
        if (forced) flagCreate("loadstate");
    }
    // debug
    // saveAndLoadState(FileToMemory, "current_state");
}


template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::loadState(bool forced)
{
//		flagCreate("do_save_state");
    if (flagExists("loadstate") || forced)
    {
        log("[loadState invoked]");
//			printf("Provadim cteni df\n");
        saveAndLoadState(FileToMemory, "current_state");
// 			if (delete_savestate_flag)
// 				flagDelete("savestate");
//			flagRename("savestate","savestate_saved");
    }
    // debug
    // saveAndLoadState(FileToMemory, "current_state");
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
bool State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::wallTimeReached()
{
    if (wallTime <=0) return false;
    timespec t_actual;
    clock_gettime(CLOCK_REALTIME, &t_actual);
    long actualtimediff = (t_actual.tv_sec - t_init.tv_sec);
    if(actualtimediff >= wallTime) 
    {
        log("wallTime reached: %ld / %ld [sec]", actualtimediff, wallTime);
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
// LBM RELATED
//                                                                                                                                                                                                                //                         
//                                                                                                                                                                                                                //                         
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::setEqLat(dreal*lat, idx x, idx y, idx z, real rho, real vx, real vy, real vz)
{
	lat[Fxyz(mmm,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_mmm(rho,vx,vy,vz);
	lat[Fxyz(zmm,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_zmm(rho,vx,vy,vz);
	lat[Fxyz(pmm,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_pmm(rho,vx,vy,vz);
	lat[Fxyz(mzm,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_mzm(rho,vx,vy,vz);
	lat[Fxyz(zzm,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_zzm(rho,vx,vy,vz);
	lat[Fxyz(pzm,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_pzm(rho,vx,vy,vz);
	lat[Fxyz(mpm,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_mpm(rho,vx,vy,vz);
	lat[Fxyz(zpm,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_zpm(rho,vx,vy,vz);
	lat[Fxyz(ppm,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_ppm(rho,vx,vy,vz);

	lat[Fxyz(mmz,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_mmz(rho,vx,vy,vz);
	lat[Fxyz(zmz,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_zmz(rho,vx,vy,vz);
	lat[Fxyz(pmz,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_pmz(rho,vx,vy,vz);
	lat[Fxyz(mzz,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_mzz(rho,vx,vy,vz);
	lat[Fxyz(zzz,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_zzz(rho,vx,vy,vz);
	lat[Fxyz(pzz,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_pzz(rho,vx,vy,vz);
	lat[Fxyz(mpz,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_mpz(rho,vx,vy,vz);
	lat[Fxyz(zpz,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_zpz(rho,vx,vy,vz);
	lat[Fxyz(ppz,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_ppz(rho,vx,vy,vz);

	lat[Fxyz(mmp,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_mmp(rho,vx,vy,vz);
	lat[Fxyz(zmp,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_zmp(rho,vx,vy,vz);
	lat[Fxyz(pmp,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_pmp(rho,vx,vy,vz);
	lat[Fxyz(mzp,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_mzp(rho,vx,vy,vz);
	lat[Fxyz(zzp,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_zzp(rho,vx,vy,vz);
	lat[Fxyz(pzp,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_pzp(rho,vx,vy,vz);
	lat[Fxyz(mpp,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_mpp(rho,vx,vy,vz);
	lat[Fxyz(zpp,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_zpp(rho,vx,vy,vz);
	lat[Fxyz(ppp,x,y,z,lbm.X, lbm.Y, lbm.Z)] = T_LBM_EQ::feq_ppp(rho,vx,vy,vz);
}
// clear Lattice + boundary setup
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::setEqLat(idx x, idx y, idx z, real rho, real vx, real vy, real vz)
{
	for (int dfty=0;dfty<DFMAX;dfty++) setEqLat(lbm.hfs[dfty], x, y, z, rho, vx, vy, vz);
}

template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
int State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::directionIndex(int i, int j, int k)
{
	if (i==-1 && j==-1 && k==-1) return mmm;
	if (i==-1 && j==-1 && k== 0) return mmz;
	if (i==-1 && j==-1 && k== 1) return mmp;
	if (i==-1 && j== 0 && k==-1) return mzm;
	if (i==-1 && j== 0 && k== 0) return mzz;
	if (i==-1 && j== 0 && k== 1) return mzp;
	if (i==-1 && j== 1 && k==-1) return mpm;
	if (i==-1 && j== 1 && k== 0) return mpz;
	if (i==-1 && j== 1 && k== 1) return mpp;

	if (i== 0 && j==-1 && k==-1) return zmm;
	if (i== 0 && j==-1 && k== 0) return zmz;
	if (i== 0 && j==-1 && k== 1) return zmp;
	if (i== 0 && j== 0 && k==-1) return zzm;
	if (i== 0 && j== 0 && k== 0) return zzz;
	if (i== 0 && j== 0 && k== 1) return zzp;
	if (i== 0 && j== 1 && k==-1) return zpm;
	if (i== 0 && j== 1 && k== 0) return zpz;
	if (i== 0 && j== 1 && k== 1) return zpp;

	if (i== 1 && j==-1 && k==-1) return pmm;
	if (i== 1 && j==-1 && k== 0) return pmz;
	if (i== 1 && j==-1 && k== 1) return pmp;
	if (i== 1 && j== 0 && k==-1) return pzm;
	if (i== 1 && j== 0 && k== 0) return pzz;
	if (i== 1 && j== 0 && k== 1) return pzp;
	if (i== 1 && j== 1 && k==-1) return ppm;
	if (i== 1 && j== 1 && k== 0) return ppz;
	if (i== 1 && j== 1 && k== 1) return ppp;
	log("directionindex: wrong i %d j %d k %d",i,j,k);
	return zzz;
}
    

// clear Lattice + boundary setup
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::resetLattice(real rho, real vx, real vy, real vz)
{
	#pragma omp parallel for schedule(static) collapse(2)
	for (idx x = 0; x < lbm.X;x++)
	for (idx z = 0; z < lbm.Z;z++)
	for (idx y = 0; y < lbm.Y;y++)
		setEqLat(x,y,z,rho,vx,vy,vz);
}


// clear Lattice and boundary setup
template< typename LBM_TYPE, typename MACRO, typename CPU_MACRO, typename LBM_DATA, typename LBM_BC>
void State<LBM_TYPE, MACRO, CPU_MACRO, LBM_DATA, LBM_BC>::reset()
{
    //pokud existuje soubor s ulozenymi df_0, tak ho nactu
//    if(flagExists("current_state/df_0"))
	if(flagExists("loadstate"))
	{
		loadState(); // load saved state into CPU memory
	} else
	{
		lbm.resetMap(LBM_BC::GEO_FLUID);
		resetLattice(1.0, 0, 0, 0);
		setupBoundaries();		// this can be virtualized
		lbm.projectWall();
	}
	
//	resetLattice(1.0, lbmInputVelocityX(), lbmInputVelocityY(),lbmInputVelocityZ());
//	lbm.copyMapToDevice();

	//initial time of current simulation
	clock_gettime(CLOCK_REALTIME, &t_init);
}
