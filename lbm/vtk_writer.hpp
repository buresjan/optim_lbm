void VTKWriter::forceBigEndian(unsigned char *bytes)
{
	static int doneTest = 0;
	static int shouldSwap = 0;
	if (!doneTest)
	{
		int tmp1 = 1;
		unsigned char *tmp2 = (unsigned char *) &tmp1;
		if (*tmp2 != 0)
			shouldSwap = 1;
		doneTest = 1;
	}
	
	if (shouldSwap)
	{
		unsigned char tmp = bytes[0];
		bytes[0] = bytes[3];
		bytes[3] = tmp;
		tmp = bytes[1];
		bytes[1] = bytes[2];
		bytes[2] = tmp;
	}
}

void VTKWriter::writeHeader(FILE*fp)
{
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "Written using RF writer\n");
	fprintf(fp, "BINARY\n");
}

void VTKWriter::writeInt(FILE*fp, int val)
{
	forceBigEndian((unsigned char *) &val);
	fwrite(&val, sizeof(int), 1, fp);
}

void VTKWriter::writeFloat(FILE*fp, float val)
{
	if (!fp) return;
	forceBigEndian((unsigned char *) &val);
	//fwrite(&val, sizeof(float), 1, fp);
	if (buffer_pos>=buffer_len) 
	{ 
		printf("vtk.writeFloat::unexpected pos %d vs. max %d\n",buffer_pos, buffer_len-1);
		return;
	}
	buffer[buffer_pos] = val;
	buffer_pos++;
	if (buffer_pos == buffer_len)
	{
		// we added last element to the buffer + flush it + reset buffer_pos to 0
		writeBuffer(fp);
	}
}

void VTKWriter::writeBuffer(FILE*fp)
{
	if (buffer_pos>0) fwrite(buffer, sizeof(float), buffer_pos, fp);
	buffer_pos=0;
}