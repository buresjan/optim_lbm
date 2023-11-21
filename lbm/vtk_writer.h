#ifndef __VTK_WRITER__
#define __VTK_WRITER__ 

#include "defs.h"

struct VTKWriter
{
//	bool binary=true;
	bool zip=false;
//	int numInColumn=0;
	long buffer_len=64*1024*1024;
	long buffer_pos=0;
	float *buffer=0;

	void forceBigEndian(unsigned char *bytes);
	void writeHeader(FILE*fp);
	void writeEndLine(FILE*fp);
	void writeInt(FILE*fp, int val);
	void writeFloat(FILE*fp, float val);
	
	void writeBuffer(FILE*fp);

/*
	template< typename real1, typename real2 >
	bool helper(const char*iid, real1 ivalue, int idofs, char*id, real2 &value, int &dofs) /// simplifies data output routine
	{
		sprintf(id,"%s",iid);
		dofs=idofs;
		value=ivalue;
		return true;
	}
*/
	VTKWriter()
	{
		buffer = (float*)calloc(buffer_len, sizeof(float));
		buffer_pos=0;
	}
	
	~VTKWriter()
	{
		if (buffer) free(buffer);
	}
};

#include "vtk_writer.hpp"

#endif
