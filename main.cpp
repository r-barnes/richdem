#include "utility.h"
#include "data_structures.h"
#include "data_io.h"
#include "d8_methods.h"
//#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"
#include "flat_resolution.h"
//#include "visualize.h"
#include <sys/time.h>

int main(int argc, char **argv){
	if(argc!=2){ //TODO
		printf("RichDEM was built by Richard Barnes (rbarnes@umn.edu, http://finog.org)\nIt was designed to:\n\t*Use the fastest available algorithms.\n\t*Run in parallel whenever possible.\n\t*Use straight-forward, easy-to-debug code.\n\t*Advance the state-of-the-art of DEM processing.\n\nIt is suggested you edit main.cpp to suit your needs.\nSyntax: ./richdem <INPUT FILE>\n");
		return -1;
	}

	float_2d elevations;
	load_ascii_data(argv[1],elevations);

	print2d("%2d ",elevations);
//	visualize(diff, false, 0);
//	visualize(elevations, false, 0);

//	pit_fill_wang(elevations);
//	print2d("%2d ",elevations);

//	visualize(elevations, false, 0);

	char_2d flowdirs(elevations);
	d8_flow_directions(elevations,flowdirs);

	timeval startTime;
	gettimeofday(&startTime, NULL);

	resolve_flats(elevations,flowdirs);
	printf("\033[96mResolve time: %lf\033[39m\n",timediff(startTime));

	print2d("%2d ",elevations);

//	print2d("%4.1f ",elevations);

//	visualize(elevations,false,0);
//	d8_flow_directions(elevations,flowdirs);
//	visualize(flowdirs,true,NO_FLOW);

//	print2d("%4.2f ",elevations);
//	print2d("%2d ",flowdirs);

//	d8_flow_directions(elevations,flowdirs,false);

//	float_2d flowdirs(elevations);
//	pit_fill_wang(elevations);

//	float_2d area(elevations);

//	dinf_flow_directions(elevations,flowdirs);
///	dinf_upslope_area(flowdirs, area);

//	output_ascii_data("zout",area);


	return 0;
}
