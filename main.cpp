#include "utility.h"
#include "data_structures.h"
#include "load_data.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"

int main(int argc, char **argv){
	float_2d elevations,elevations2;
	uint_2d area;
	try{
		load_ascii_data(argv[1],elevations);
		elevations2=elevations;
		pit_fill_wang(elevations2);
		pit_fill_barnes(elevations);

//		char_2d flowdirs(elevations);

	//	pit_fill_yonghe2009(elevations);
	//	pit_fill_wang(elevations);
//	dinf_flow_directions(elevations,flowdirs);
	//	dinf_upslope_area(flowdirs);

//		d8_flow_directions(elevations,flowdirs);

//		uint_2d area(elevations);
//		print_flow(flowdirs);
//		d8_upslope_area(flowdirs,area);
		return 0;
	} catch (int e) {
		diagnostic("Unfortunately, I was unable to continue.\nClosing...\n");
		return -1;
	}
}
