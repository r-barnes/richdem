#include "utility.h"
#include "data_structures.h"
#include "data_io.h"
#include "d8_methods.h"
//#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"
#include "flat_resolution.h"
#include "visualize.h"
#include "debug.h"
#include <sys/time.h>
#include <iostream>
#include <iomanip>

int main(int argc, char **argv){
	if(argc!=2){ //TODO
		printf("RichDEM was built by Richard Barnes (rbarnes@umn.edu, http://finog.org)\nIt was designed to:\n\t*Use the fastest available algorithms.\n\t*Run in parallel whenever possible.\n\t*Use straight-forward, easy-to-debug code.\n\t*Advance the state-of-the-art of DEM processing.\n\nIt is suggested you edit main.cpp to suit your needs.\nSyntax: ./richdem <INPUT FILE>\n");
		return -1;
	}

	float_2d elevations;
	load_ascii_data(argv[1],elevations);

	char_2d flowdirs(elevations);
	d8_flow_directions(elevations,flowdirs);

	visualize(flowdirs,true,NO_FLOW);

	int_2d flat_resolution_mask(elevations,true);
	resolve_flats(elevations,flowdirs,flat_resolution_mask);

	d8_flow_directions(flat_resolution_mask,flowdirs,false);

	visualize(flowdirs,true,NO_FLOW);

//	PRINT(flowdirs,0,3)

//	PRINT(elevations,0,2);
//	visualize(diff, false, 0);
//	visualize(elevations, false, 0);

//	pit_fill_wang(elevations);
//	print2d("%2d ",elevations);

//	visualize(elevations, false, 0);
	return 0;


//	PRINT(flowdirs,0,2);

	timeval startTime;
	gettimeofday(&startTime, NULL);


	printf("\033[96mResolve time: %lf\033[39m\n",timediff(startTime));

//	PRINT(elevations,3,5);

//	print2d("%4.1f ",elevations);

//	visualize(elevations,false,0);

//	visualize(flowdirs,true,NO_FLOW);

//	print2d("%4.2f ",elevations);
//	PRINT(flowdirs,0,2);

//	d8_flow_directions(elevations,flowdirs,false);

//	float_2d flowdirs(elevations);
//	pit_fill_wang(elevations);

//	float_2d area(elevations);

//	dinf_flow_directions(elevations,flowdirs);
///	dinf_upslope_area(flowdirs, area);

//	output_ascii_data("zout",area);


	return 0;
}
