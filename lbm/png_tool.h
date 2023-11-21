#pragma once
#include "defs.h"

struct PNGTool
{
	int width;
	int height;
	bool allocated;
	png_bytep *row_pointers;
	
	// coords a,b in [0,1]
	// rule:
	// a=0 ... x=0
	// a=1 ... x=width-1
	int intensity(double ia, double ib)
	{
		if (!allocated) { printf("PNGTool::intensity() png file not allocated returning 0\n"); return 0; }
		// spocitej x a y
		int y=ib*(height-1);
		int x=ia*(width-1);
		
		png_bytep row = row_pointers[y];
		png_bytep px = &(row[x * 4]);

		unsigned char r=px[0];
		unsigned char g=px[1];
		unsigned char b=px[2];

		return ((r & 0xff) << 16) + ((g & 0xff) << 8) + (b & 0xff);
	}
	
	bool readPNG(const char *filename)
	{
		FILE *fp = fopen(filename, "rb");
		if (!fp) { printf("file %s png does not exist\n",filename); return false; }


		png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		if(!png) { printf("file %s png read error\n",filename); return false; }
	
		png_infop info = png_create_info_struct(png);
		if(!png) { printf("file %s png read error\n",filename); return false; }

		if(setjmp(png_jmpbuf(png))) { printf("file %s png read error\n",filename); return false;}
    
		png_init_io(png, fp);

		png_read_info(png, info);

		width = png_get_image_width(png, info);
		height = png_get_image_height(png, info);
		png_byte color_type  = png_get_color_type(png, info);
		png_byte bit_depth = png_get_bit_depth(png, info);

		if(bit_depth == 16) png_set_strip_16(png);
		if(color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);

		// PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
		if(color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png);

		if(png_get_valid(png, info, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);

		// These color_type don't have an alpha channel then fill it with 0xff.
		if(color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_PALETTE) png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

		if(color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) png_set_gray_to_rgb(png);

		png_read_update_info(png, info);

		row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
		for(int y = 0; y < height; y++) row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png,info));
		png_read_image(png, row_pointers);
		fclose(fp);
		return true;
	}
	

	PNGTool(const char*filename)
	{
		allocated = readPNG(filename);
	}
	
	~PNGTool()
	{
		// free
		if (allocated)
		{
			for(int y = 0; y < height; y++) free(row_pointers[y]);
			free(row_pointers);
		}
	}
};
