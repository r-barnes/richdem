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

	timeval startTime;

	float_2d elevations,elevND,elevbar1,elevbar2,elevbar3;
	load_ascii_data(argv[1],elevations);
	elevND=elevations;
	elevbar1=elevations;
	elevbar2=elevations;
	elevbar3=elevations;

	gettimeofday(&startTime, NULL);
	pit_fill_wang(elevations);
	printf("\033[96mWang pit fill time: %lf\033[39m\n",timediff(startTime));
/*
	gettimeofday(&startTime, NULL);
	pit_fill_wangND(elevND);
	printf("\033[96mWangND pit fill time: %lf\033[39m\n",timediff(startTime));
*/
	gettimeofday(&startTime, NULL);
	pit_fill_barnes1(elevbar1);
	printf("\033[96mBarnes v1 pit fill time: %lf\033[39m\n",timediff(startTime));
	printf("Good? %d\n",elevations==elevbar1);

	gettimeofday(&startTime, NULL);
	pit_fill_barnes2(elevbar2);
	printf("\033[96mBarnes v2 pit fill time: %lf\033[39m\n",timediff(startTime));
	printf("Good? %d\n",elevations==elevbar2);

	gettimeofday(&startTime, NULL);
	pit_fill_barnes3(elevbar3);
	printf("\033[96mBarnes v3 pit fill time: %lf\033[39m\n",timediff(startTime));
	printf("Good? %d\n",elevations==elevbar3);

//	gettimeofday(&startTime, NULL);
//	float_2d flowdirs(elevations);
//	dinf_flow_directions(elevations,flowdirs);
//	printf("\033[96Flow dirs time: %lf\033[39m\n",timediff(startTime));

//	gettimeofday(&startTime, NULL);
//	int_2d flat_resolution_mask(elevations), groups;
//	resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);
//	printf("\033[96mResolve time: %lf\033[39m\n",timediff(startTime));

//	gettimeofday(&startTime, NULL);
//	dinf_flow_flats(flat_resolution_mask,groups,flowdirs);
//	printf("\033[96Flow flats time: %lf\033[39m\n",timediff(startTime));

//	PRINT(flat_resolution_mask,0,2);
//	PRINT(flowdirs,0,2);

//	gettimeofday(&startTime, NULL);
//	float_2d area(elevations);
//	dinf_upslope_area(flowdirs, area);
//	printf("\033[96Upslope area time: %lf\033[39m\n",timediff(startTime));

//	visualize(area,true,(float)0,"D8 Upslope Area w/ Flats Resolved");

//	output_ascii_data(argv[2],area);

	return 0;
}
