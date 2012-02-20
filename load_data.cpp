#include <stdio.h>
#include <stdlib.h>
#include <exception>
#include "interface.h"
#include "data_structures.h"

int load_ascii_data(char filename[], float_2d &elevations){
	FILE *fin;
	long long file_size;
	int rows,columns,cellsize;
	double xllcorner,yllcorner,no_data;

	diagnostic_arg("Opening input ASCII-DEM file \"%s\"...",filename);
	fin=fopen(filename,"r");
	if(fin==NULL){
		diagnostic("failed!\n");
		return -1;
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
	if(fscanf(fin,"ncols %d nrows %d xllcorner %lf yllcorner %lf cellsize %d NODATA_value %lf",&columns,&rows,&xllcorner,&yllcorner,&cellsize,&no_data)!=6){
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
	for(int y=0;y<rows;y++){
		progress_bar(100*ftell(fin)/file_size); //Todo: Should I check to see if ftell fails here?
		for(int x=0;x<columns;x++){
			if (fscanf(fin,"%f", &temp)!=1){
				diagnostic("failed! (Couldn't read or convert a value!)\n");
				return -1;
			}
			elevations(x,y)=temp;
		}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	fclose(fin);

    return 0;
}
