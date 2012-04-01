#include "utility.h"
#include "data_structures.h"
#include "data_io.h"
#include "d8_methods.h"
#include "dinf_methods.h"
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

//	pit_fill_wang(elevations);
//	elevations.print_block(std::cerr,cx-5,cx+5,cy-5,cy+5,3,4);

	float_2d flowdirs(elevations);
	dinf_flow_directions(elevations,flowdirs);
//	visualize(flowdirs,true,NO_FLOW,"D8 Flowdirs");
//	flowdirs.print_block(std::cerr,cx-5,cx+5,cy-5,cy+5,0,4);

	timeval startTime;
	gettimeofday(&startTime, NULL);
	int_2d flat_resolution_mask(elevations), groups;
	resolve_flats(elevations,flowdirs,flat_resolution_mask,groups);
	printf("\033[96mResolve time: %lf\033[39m\n",timediff(startTime));

//	flat_resolution_mask.print_block(std::cerr,cx-5,cx+5,cy-5,cy+5,0,4);

	dinf_flow_flats(flat_resolution_mask,groups,flowdirs);

//	visualize(flowdirs,true,(float)NO_FLOW,"D8 Flow Directions With Flats Resolved");

//	flowdirs.print_block(std::cerr,cx-5,cx+5,cy-5,cy+5,0,4);

//	visualize(flowdirs,true,NO_FLOW,"D8 Flow Directions w/ Flat Resolution");

	float_2d area(elevations);
	dinf_upslope_area(flowdirs, area);

//	area.print_block(std::cerr,cx-5,cx+5,cy-5,cy+5,0,4);

	int nfs=0;
	for(int x=0;x<flowdirs.width();x++)
	for(int y=0;y<flowdirs.height();y++)
		if(flowdirs(x,y)!=flowdirs.no_data && flowdirs(x,y)==NO_FLOW){
			nfs++;
/*		if(area(x,y)==-1){
			diagnostic_arg("Area is -1 at (%d,%d). Flowdir=%d. Flows into it=",x,y,flowdirs(x,y));
			for(int n=1;n<=8;n++)
				if(flowdirs(x,y)!=NO_FLOW && flowdirs(x,y)!=flowdirs.no_data && n==inverse_flow[flowdirs(x,y)])
					diagnostic_arg("%d ",n);
			diagnostic("\n");
/*			elevations.print_block(std::cerr,x-5,x+5,y-5,y+5,3,4);
			diagnostic("---\n");
			flat_resolution_mask.print_block(std::cerr,x-5,x+5,y-5,y+5,0,4);
			diagnostic("---\n");
			groups.print_block(std::cerr,x-5,x+5,y-5,y+5,0,4);
			diagnostic("---\n");
			flowdirs.print_block(std::cerr,x-5,x+5,y-5,y+5,0,4);
			diagnostic("---\n");
			area.print_block(std::cerr,x-5,x+5,y-5,y+5,0,4);
			diagnostic("==========\n==========\n");*/
		}

	int cx=1198,cy=1541;
	elevations.print_block(std::cerr,cx-5,cx+5,cy-5,cy+5,4,8);
	diagnostic("---\n");
	flat_resolution_mask.print_block(std::cerr,cx-5,cx+5,cy-5,cy+5,0,4);
	diagnostic("---\n");
	groups.print_block(std::cerr,cx-5,cx+5,cy-5,cy+5,0,4);
	diagnostic("---\n");
	flowdirs.print_block(std::cerr,cx-5,cx+5,cy-5,cy+5,4,6);
	diagnostic("---\n");
	area.print_block(std::cerr,cx-5,cx+5,cy-5,cy+5,0,4);
	diagnostic("==========\n==========\n");

	diagnostic_arg("No flows: %d\n",nfs);

	visualize(area,true,(float)-1,"D8 Upslope Area w/ Flats Resolved");

//	output_ascii_data("zout",area);

	return 0;
}
