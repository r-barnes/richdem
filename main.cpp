#include "utility.h"
#include "data_structures.h"
#include "load_data.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"

int main(int argc, char **argv){
	float_2d elevations;
	float_2d flowdirs;
	float no_data;
	int data_cells;
	try{
		data_cells=load_ascii_data(argv[1],elevations,no_data);

	//	pit_fill_yonghe2009(elevations);
	//	pit_fill_wang(elevations);
	//	dinf_flow_directions(elevations,flowdirs);
	//	dinf_upslope_area(flowdirs);

		dinf_flow_directions(elevations,flowdirs,no_data);

//		print_flow(flowdirs);

		dinf_upslope_area(flowdirs,data_cells);
		return 0;
	} catch (int e) {
		diagnostic("Unfortunately, I was unable to continue.\nClosing...\n");
		return -1;
	}
}
