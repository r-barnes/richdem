#include "utility.h"
#include "data_structures.h"
#include "data_io.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"
#include "flat_resolution.h"
//#include "visualize.h"
#include "debug.h"
#include <sys/time.h>
#include <iostream>
#include <iomanip>

int main(int argc, char **argv){
	if(argc!=3){ //TODO
		printf("RichDEM was built by Richard Barnes (rbarnes@umn.edu, http://finog.org)\nIt was designed to:\n\t*Use the fastest available algorithms.\n\t*Run in parallel whenever possible.\n\t*Use straight-forward, easy-to-debug code.\n\t*Advance the state-of-the-art of DEM processing.\n\nIt is suggested you edit main.cpp to suit your needs.\nSyntax: ./richdem <INPUT FILE> <OUTPUT FILE>\n");
		return -1;
	}

	float_2d elevations;
	load_ascii_data(argv[1],elevations);

	pit_fill_wang(elevations);

	float_2d flowdirs(elevations);
	dinf_flow_directions(elevations,flowdirs);

	timeval startTime;
	gettimeofday(&startTime, NULL);
	int_2d flat_resolution_mask(elevations), groups;
	resolve_flats(elevations,flowdirs,flat_resolution_mask,groups);
	printf("\033[96mResolve time: %lf\033[39m\n",timediff(startTime));

	dinf_flow_flats(flat_resolution_mask,groups,flowdirs);

//	PRINT(flat_resolution_mask,0,2);
//	PRINT(flowdirs,0,2);

	float_2d area(elevations);
	dinf_upslope_area(flowdirs, area);

//	visualize(area,true,(float)0,"D8 Upslope Area w/ Flats Resolved");

	output_ascii_data(argv[2],area);

	return 0;
}
