#include "utility.h"
#include "data_structures.h"
#include <cstdio>

void print_dem(float_2d &elevations, int mark_x=-1, int mark_y=-1){
	for(int y=0;y<elevations.size2();y++){
		for(int x=0;x<elevations.size1();x++){
			if(x==mark_x && y==mark_y) printf("\033[91m");
			printf("%.0f ",elevations(x,y));
			if(x==mark_x && y==mark_y) printf("\033[39m");
		}
		printf("\n");
	}
	printf("\n");
}

void print_flow(float_2d &flowdirs){
	for(int y=0;y<flowdirs.size2();y++){
		for(int x=0;x<flowdirs.size1();x++)
			printf("%5.2f ",flowdirs(x,y));
		printf("\n");
	}
	printf("\n");
}
