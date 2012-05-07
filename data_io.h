#ifndef _data_io_included
#define _data_io_included

#include "interface.h"
#include <fstream>
#include <string>
//#include <fcntl.h> //Used for posix_fallocate

int load_ascii_data(const char filename[], float_2d &elevations);

#define OUTPUT_DEM	1
#define OUTPUT_OMG	2


template <class T>
int output_ascii_data(const std::string filename, const array2d<T> &output_grid, int precision=8){
	std::ofstream fout;
	std::string outputsep=" ";
	int output_type=OUTPUT_DEM;

	diagnostic_arg("Opening ASCII output file \"%s\"...",filename.c_str());
	fout.open(filename.c_str());
	if(!fout.is_open()){
		diagnostic("failed!\n");
		exit(-1);	//TODO: Need to make this safer! Don't just close after all that work!
	}
	diagnostic("succeeded.\n");

/*	diagnostic_arg("Reserving %ldMB of disk space...", output_grid.estimated_output_size()/1024/1024);
	if(posix_fallocate(fileno(fout),0,output_grid.estimated_output_size())){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");
*/

	//OmniGlyph output
	if(filename.substr(filename.length()-4)==".omg"){
		outputsep="|";
		output_type=OUTPUT_OMG;
		diagnostic("Writing OmniGlyph file header...");
		fout<<"Contents: Pixel array"<<std::endl;
		fout<<std::endl;
		fout<<"Width:    "<<output_grid.width()<<std::endl;
		fout<<"Height:   "<<output_grid.height()<<std::endl;
		fout<<std::endl;
		fout<<"Spectral bands:   1"<<std::endl;
		fout<<"Bits per band:   32"<<std::endl;
		fout<<"Range of values:   "<<output_grid.min()<<","<<output_grid.max()<<std::endl;
		fout<<"Actual range:   "<<output_grid.no_data<<","<<output_grid.max()<<std::endl;	//TODO: Assumes no_data is a small negative value
		fout<<"Gamma exponent:   0."<<std::endl;
		fout<<"Resolution:   100 pixels per inch"<<std::endl;
		fout<<std::endl;
		fout<<"|"<<std::endl;
	} else {
		diagnostic("Writing ArcGrid ASCII file header...");
		fout<<"ncols\t\t"<<output_grid.width()<<std::endl;
		fout<<"nrows\t\t"<<output_grid.height()<<std::endl;
		fout<<"xllcorner\t"<<std::fixed<<std::setprecision(precision)<<output_grid.xllcorner<<std::endl;
		fout<<"yllcorner\t"<<std::fixed<<std::setprecision(precision)<<output_grid.yllcorner<<std::endl;
		fout<<"cellsize\t"<<std::fixed<<std::setprecision(precision)<<output_grid.cellsize<<std::endl;
		fout<<"NODATA_value\t"<<std::fixed<<std::setprecision(precision)<<output_grid.no_data<<std::endl;
	}
	diagnostic("succeeded.\n");

	for(int y=0;y<output_grid.height();y++){
		if(output_type==OUTPUT_OMG)
			fout<<"|";
		for(int x=0;x<output_grid.width();x++)
			if(sizeof(T)==1)	//TODO: This is a crude way of detecting chars and bools
				fout<<(int)output_grid(x,y)<<outputsep;
			else
				fout<<std::fixed<<std::setprecision(precision)<<output_grid(x,y)<<outputsep;
		fout<<std::endl;
	}

//	diagnostic("Writing file data...");
//	fout<<std::setprecision(precision)<<output_grid;
//	diagnostic("succeeded.\n");

	fout.close();

	return 0;
}

#endif
