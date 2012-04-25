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
#include <map>

int main(int argc, char **argv){
/*	if(argc!=4){ //TODO
		printf("RichDEM was built by Richard Barnes (rbarnes@umn.edu, http://finog.org)\nIt was designed to:\n\t*Use the fastest available algorithms.\n\t*Run in parallel whenever possible.\n\t*Use straight-forward, easy-to-debug code.\n\t*Advance the state-of-the-art of DEM processing.\n\nIt is suggested you edit main.cpp to suit your needs.\nSyntax: ./richdem <INPUT FILE> <OUTPUT FILE>\n");
		return -1;
	}*/

	float_2d elevations,angbar,angtau,angle_diff,flowdirs;
	load_ascii_data(argv[1],elevations);
	load_ascii_data(argv[2],angbar);
	load_ascii_data(argv[3],angtau);
	load_ascii_data(argv[4],angle_diff);

	barnes_flood(elevations);
	dinf_flow_directions(elevations,flowdirs);

	uint_2d pitloc(elevations);
	flat_mask(flowdirs,pitloc);

	uint_2d odd_flows(elevations,true);
	odd_flows.no_data=3;

	for(int x=0;x<angle_diff.width();x++)
	for(int y=0;y<angle_diff.height();y++)
		if(pitloc(x,y)==pitloc.no_data)
			continue;//odd_flows(x,y)=odd_flows.no_data;
		else if(angle_diff(x,y)>10.0*M_PI/180.0 && pitloc(x,y)!=1){
			diagnostic_arg("(%d,%d)\tBarnes: %f, Tau: %f, Diff: %f\n",x,y,DEG(angbar(x,y)),DEG(angtau(x,y)),DEG(angle_diff(x,y)));
			elevations.surroundings(x,y);
		}
//		else
//			odd_flows(x,y)=0;

//	output_ascii_data(argv[3],odd_flows,0);

/*	barnes_flood(elevations);
//	pit_fill_barneslehman(elevations);

	float_2d flowdirs(elevations);
	dinf_flow_directions(elevations,flowdirs);

	int_2d flat_resolution_mask(elevations), groups;
	resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);

	dinf_flow_flats(flat_resolution_mask,groups,flowdirs);

//	float_2d area(elevations);
//	dinf_upslope_area(flowdirs, area);

	output_ascii_data(argv[2],flowdirs,8);

/*
//	PRINT(flat_resolution_mask,0,2);
//	PRINT(flowdirs,0,2);

	gettimeofday(&startTime, NULL);
	float_2d area(elevations);
	dinf_upslope_area(flowdirs, area);
	printf("\033[96mUpslope area time: %lf\033[39m\n",timediff(startTime));

//	visualize(area,true,(float)0,"D8 Upslope Area w/ Flats Resolved");
*/
	return 0;
}
