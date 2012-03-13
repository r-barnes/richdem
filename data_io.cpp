#include <stdio.h>
#include <stdlib.h>
#include "interface.h"
#include "data_structures.h"

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
	elevations.resize(columns,rows);
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
