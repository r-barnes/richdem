#include <stdlib.h>
#include <fstream>
#include "interface.h"
#include "data_structures.h"
#include "utility.h"
#include <string>
#include <boost/lexical_cast.hpp>
#include <limits>

int load_ascii_data(char filename[], float_2d &elevations){
	std::ifstream fin;
	long file_size;
	int rows,columns;

	diagnostic_arg("Opening input ASCII-DEM file \"%s\"...",filename);
	fin.open(filename);
	if(!fin.is_open()){
		diagnostic("failed!\n");
		exit(-1);
	}
	diagnostic("succeeded.\n");

	diagnostic("Calculating file size...");
	file_size=fin.tellg();
	fin.seekg(0,std::ios::end);
	file_size=fin.tellg()-file_size;
	fin.seekg(0,std::ios::beg);
	diagnostic("succeeded.\n");

	diagnostic("Reading DEM header...");
	fin>>must_be("ncols")>>columns;
	fin>>must_be("nrows")>>rows;
	fin>>must_be("xllcorner")>>elevations.xllcorner;
	fin>>must_be("yllcorner")>>elevations.yllcorner;
	fin>>must_be("cellsize")>>elevations.cellsize;
	fin>>must_be("NODATA_value")>>elevations.no_data;
	diagnostic("succeeded.\n");

	diagnostic_arg("The loaded DEM will require approximately %ldMB of RAM.\n",columns*rows*sizeof(float)/1024/1024);

	diagnostic("Resizing elevation matrix...");	//TODO: Consider abstracting this block
	elevations.resize(columns,rows);
	diagnostic("succeeded.\n");

	diagnostic("Reading elevation matrix...\n");
	progress_bar(-1);
	elevations.data_cells=0;
	unsigned int precision_max=0;
	unsigned int number_length=0;
	float max=-std::numeric_limits<float>::infinity();
	float min=std::numeric_limits<float>::infinity();
	for(int y=0;y<rows;y++){
		progress_bar(fin.tellg()*100/file_size); //Todo: Should I check to see if ftell fails here?
		for(int x=0;x<columns;x++){
			std::string temp;
			fin>>temp;
			elevations(x,y)=boost::lexical_cast<float>(temp);

			if(elevations(x,y)==elevations.no_data)
				continue;

			if(elevations(x,y)<min)
				min=elevations(x,y);
			if(elevations(x,y)>max)
				max=elevations(x,y);

			elevations.data_cells++;
			unsigned int decimal_point_at=temp.find('.');
			if(decimal_point_at!=std::string::npos && temp.length()-decimal_point_at-1>precision_max)
				precision_max=temp.length()-decimal_point_at-1;
			if(temp.length()>number_length)
				number_length=temp.length();
		}
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));

	fin.close();

	diagnostic_arg("Read %ld cells, of which %ld contained data (%ld%%).\n", elevations.width()*elevations.height(), elevations.data_cells, elevations.data_cells*100/elevations.width()/elevations.height());
	diagnostic_arg("The number of digits (including decimal points) in an elevation was %d.\n",number_length);
	diagnostic_arg("The maximal precision in an elevation was %d decimal digits.\n",precision_max);
	diagnostic_arg("Elevations ranged from %f to %f.\n",min,max);

	return 0;
}
