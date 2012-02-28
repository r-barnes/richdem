#include "utility.h"
#include "data_structures.h"
#include "load_data.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"

int main(int argc, char **argv){
	float_2d elevations;
	try{
		load_ascii_data(argv[1],elevations);
		float_2d flowdirs(elevations);

	//	pit_fill_yonghe2009(elevations);
	//	pit_fill_wang(elevations);
	//	dinf_flow_directions(elevations,flowdirs);
	//	dinf_upslope_area(flowdirs);

		dinf_flow_directions(elevations,flowdirs);

		float_2d area(elevations);
//		print_flow(flowdirs);

		dinf_upslope_area(flowdirs,area);
		return 0;
	} catch (int e) {
		diagnostic("Unfortunately, I was unable to continue.\nClosing...\n");
		return -1;
	}
}
