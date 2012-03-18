#ifndef _data_io_included
#define _data_io_included

#include "interface.h"
//#include <fcntl.h> //Used for posix_fallocate

int load_ascii_data(char filename[], float_2d &elevations);



template <class T>
int output_ascii_data(const char filename[], const array2d<T> &output_grid){
	FILE *fout;
//	long long file_size;

	diagnostic_arg("Opening ASCII output file \"%s\"...",filename);
	fout=fopen(filename,"w");
	if(fout==NULL){
		diagnostic("failed!\n");
		throw -1;
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
	if(fprintf(fout,"ncols\t\t%ld\nnrows\t\t%ld\nxllcorner\t%lf\nyllcorner\t%lf\ncellsize\t%d\nNODATA_value\t%f\n", output_grid.width(), output_grid.height(), output_grid.xllcorner, output_grid.yllcorner, output_grid.cellsize, output_grid.no_data)<0){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Writing file data...\n");
	progress_bar(-1);
	for(int y=0;y<output_grid.height();y++){
//		progress_bar(ftell(out)*100/file_size); //Todo: Should I check to see if ftell fails here?
		for(int x=0;x<output_grid.width();x++)
			if (!output_grid.print(fout,x,y)){
				diagnostic("failed! (Couldn't write a value!)\n");
				return -1;
			}
		fprintf(fout,"\n");
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	fclose(fout);

	return 0;
}

#endif
