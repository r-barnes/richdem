#include "utility.h"
#include "data_structures.h"
#include "data_io.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"
#include "flat_resolution.h"

int main(int argc, char **argv){
	float_2d elevations;
	try{
		load_ascii_data(argv[1],elevations);
//		pit_fill_barnes3(elevations);
//		print2d("%2.0f ", elevations);

		char_2d flowdirs(elevations);
		d8_flow_directions(elevations,flowdirs);
//		print2d("%2d ",flowdirs);
		resolve_flats(elevations,flowdirs);

		d8_flow_directions(elevations,flowdirs);

		print2d("%4.2f ",elevations);
		print2d("%2d ",flowdirs);

//		d8_flow_directions(elevations,flowdirs,false);

//		float_2d flowdirs(elevations);
//		pit_fill_wang(elevations);

//		float_2d area(elevations);

//		dinf_flow_directions(elevations,flowdirs);
///		dinf_upslope_area(flowdirs, area);

//		output_ascii_data("zout",area);


		return 0;
	} catch (int e) {
		diagnostic("Unfortunately, I was unable to continue.\nClosing...\n");
		return -1;
	}
}
