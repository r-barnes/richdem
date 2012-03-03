#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include "interface.h"
#include "data_structures.h"
#include <fcntl.h>

int load_ascii_data(char filename[], float_2d &elevations){
	FILE *fin;
	long long file_size;
	int rows,columns;

	diagnostic_arg("Opening input ASCII-DEM file \"%s\"...",filename);
	fin=fopen(filename,"r");
	if(fin==NULL){
		diagnostic("failed!\n");
		throw -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Calculating file size...");
	if(fseek(fin,-1,SEEK_END)!=0){
		diagnostic("failed! (Couldn't jump to end of file.)\n");
		return -1;
	}
	if((file_size=ftell(fin))==-1){
		diagnostic("failed! (Couldn't determine file size.)\n");
		return -1;
	}
	if(fseek(fin,0,SEEK_SET)!=0){
		diagnostic("failed! (Couldn't jump back to beginning of file.)\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Reading DEM header...");
	if(fscanf(fin,"ncols %d nrows %d xllcorner %lf yllcorner %lf cellsize %d NODATA_value %f",&columns, &rows, &elevations.xllcorner, &elevations.yllcorner, &elevations.cellsize, &elevations.no_data)!=6){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic_arg("The loaded DEM will require approximately %ldMB of RAM.\n",columns*rows*sizeof(float)/1024/1024);

	diagnostic("Resizing elevation matrix...");	//TODO: Consider abstracting this block
	try{
		elevations.resize(columns,rows);
	} catch (std::exception &e) {
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Reading elevation matrix...\n");
	progress_bar(-1);
	float temp;
	elevations.data_cells=0;
	for(int y=0;y<rows;y++){
		progress_bar(ftell(fin)*100/file_size); //Todo: Should I check to see if ftell fails here?
		for(int x=0;x<columns;x++){
			if (fscanf(fin,"%f", &temp)!=1){
				diagnostic("failed! (Couldn't read or convert a value!)\n");
				return -1;
			}
			elevations(x,y)=temp;
			if(temp!=elevations.no_data)
				elevations.data_cells++;
		}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	fclose(fin);

	diagnostic_arg("Read %ld cells, of which %ld contained data (%ld%%).\n", elevations.width()*elevations.height(), elevations.data_cells, elevations.data_cells*100/elevations.width()/elevations.height());

	return 0;
}


/*
template <class T>
int load_ascii_data(char filename[], array2d<T> &output_grid){
	FILE *fout;
	long long file_size;

	diagnostic_arg("Opening ASCII output file \"%s\"...",filename);
	fout=fopen(filename,"w");
	if(fout==NULL){
		diagnostic("failed!\n");
		throw -1;
	}
	diagnostic("succeeded.\n");

	diagnostic_arg("Reserving %ldMB of disk space...", output_grid.estimated_output_size()/1024/1024);
	if(posix_fallocate(fout,0,output_grid.estimated_output_size())){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Writing file header...");
	if(fprintf(fout,"ncols\t\t%d\nnrows\t\t%d\nxllcorner\t\t%lf\nyllcorner\t\t%lf\ncellsize\t\t%d\nNODATA_value\t\t%f\n", output_grid.width(), output_grid.height(), output_grid.xllcorner, output_grid.yllcorner, output_grid.cellsize, output_grid.no_data)!=6){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Writing file data...\n");
	progress_bar(-1);
	for(int y=0;y<output_grid.height();y++){
//		progress_bar(ftell(fin)*100/file_size); //Todo: Should I check to see if ftell fails here?
		for(int x=0;x<output_grid.width();x++)
			if (!output_grid.print(x,y)){
				diagnostic("failed! (Couldn't write a value!)\n");
				return -1;
			}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	fclose(fout);

	return 0;
}*/
