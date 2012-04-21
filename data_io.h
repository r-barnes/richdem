#ifndef _data_io_included
#define _data_io_included

#include "interface.h"
#include <fstream>
//#include <fcntl.h> //Used for posix_fallocate

int load_ascii_data(const char filename[], float_2d &elevations);



template <class T>
int output_ascii_data(const char filename[], const array2d<T> &output_grid, int precision=3){
	std::ofstream fout;

	diagnostic_arg("Opening ASCII output file \"%s\"...",filename);
	fout.open(filename);
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

	diagnostic("Writing file header...");
	fout<<"ncols\t\t"<<output_grid.width()<<std::endl;
	fout<<"nrows\t\t"<<output_grid.height()<<std::endl;
	fout<<"xllcorner\t"<<output_grid.xllcorner<<std::endl;
	fout<<"yllcorner\t"<<output_grid.yllcorner<<std::endl;
	fout<<"cellsize\t"<<output_grid.cellsize<<std::endl;
	fout<<"NODATA_value\t"<<std::fixed<<std::setprecision(precision)<<output_grid.no_data<<std::endl;
	diagnostic("succeeded.\n");

	diagnostic("Writing file data...");
	fout<<std::setprecision(precision)<<output_grid;
	diagnostic("succeeded.\n");

	fout.close();

	return 0;
}

#endif
