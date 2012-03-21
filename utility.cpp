#include "utility.h"
#include "interface.h"
#include "data_structures.h"
#include <cstdio>

void print_dem(const float_2d &elevations, int mark_x, int mark_y, int colour){
	for(int y=0;y<elevations.height();y++){
		for(int x=0;x<elevations.width();x++){
			if(x!=-1 && x==mark_x && y==mark_y) printf("\033[%dm",colour);
			printf("%.0f ",elevations(x,y));
			if(x!=-1 && x==mark_x && y==mark_y) printf("\033[39m");
		}
		printf("\n");
	}
	printf("\n");
}

std::istream& operator>>(std::istream& inputstream, const must_be &a){
    std::string in_val;
    inputstream >> in_val;
	if(in_val!=a.val){
		diagnostic_arg("failed!\n\tCould not match \"%s\"!\n",a.val.c_str());
		throw std::string("Input file didn't match expectations!");
	}
    return inputstream;
}
