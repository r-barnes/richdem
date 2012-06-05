#include <stdio.h>
#include <stdlib.h>
#include "interface.h"
#include "data_structures.h"
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include "utility.h"

int load_ascii_data(const char filename[], float_2d &elevations){
	FILE *fin;
	long long file_size;
	int rows,columns;
	Timer load_time;

	load_time.start();

	errno=0;
	diagnostic_arg("Opening input ASCII-DEM file \"%s\"...",filename);
	fin=fopen(filename,"r");
	if(fin==NULL){
		diagnostic_arg("failed with error %d: \"%s\"!\n",errno,strerror(errno));
		exit(-1);
	}
	diagnostic("succeeded.\n");

	diagnostic("Calculating file size...");
	if(fseek(fin,-1,SEEK_END)!=0){
		diagnostic("failed! (Couldn't jump to end of file.)\n");
		exit(-1);
	}
	if((file_size=ftell(fin))==-1){
		diagnostic("failed! (Couldn't determine file size.)\n");
		exit(-1);
	}
	if(fseek(fin,0,SEEK_SET)!=0){
		diagnostic("failed! (Couldn't jump back to beginning of file.)\n");
		exit(-1);
	}
	diagnostic("succeeded.\n");

//	posix_fadvise(fileno(fin),0,0,POSIX_FADV_SEQUENTIAL);

	diagnostic("Reading DEM header...");
	if(fscanf(fin,"ncols %d nrows %d xllcorner %lf yllcorner %lf cellsize %lf NODATA_value %f",&columns, &rows, &elevations.xllcorner, &elevations.yllcorner, &elevations.cellsize, &elevations.no_data)!=6){
		diagnostic("failed!\n");
		exit(-1);
	}
	diagnostic("succeeded.\n");

	diagnostic_arg("The loaded DEM will require approximately %ldMB of RAM.\n",columns*rows*((long)sizeof(float))/1024/1024);

	diagnostic("Resizing elevation matrix...");	//TODO: Consider abstracting this block
	elevations.resize(columns,rows);
	diagnostic("succeeded.\n");

	diagnostic("%%Reading elevation matrix...\n");
	progress_bar(-1);
	float temp;
	elevations.data_cells=0;
	for(int y=0;y<rows;y++){
		progress_bar(((long long)ftell(fin))*((long long)100)/file_size); //Todo: Should I check to see if ftell fails here?
		for(int x=0;x<columns;x++){
			if (fscanf(fin,"%f", &temp)!=1){
				diagnostic("\n\tFailed! (Couldn't read or convert a value!)\n");
				exit(-1);
			}
			elevations(x,y)=temp;
			if(temp!=elevations.no_data)
				elevations.data_cells++;
		}
	}
	diagnostic_arg(SUCCEEDED_IN,progress_bar(-1));

	fclose(fin);

	diagnostic_arg("Read %ld cells, of which %ld contained data (%ld%%).\n", elevations.width()*elevations.height(), elevations.data_cells, elevations.data_cells*100/elevations.width()/elevations.height());

	load_time.stop();
	diagnostic_arg("Read time was: %lf\n", load_time.get());

	return 0;
}
