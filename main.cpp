#include "utility.h"
#include "data_structures.h"
#include "data_io.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"
#include "flat_resolution.h"
#include "watershed.h"
//#include "visualize.h"
#include "debug.h"
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <map>

int main(int argc, char **argv){
	float_2d elevations;
	load_ascii_data(argv[1],elevations);

	barnes_flood(elevations);

	float_2d flowdirs(elevations);
	dinf_flow_directions(elevations,flowdirs);

	int_2d flat_resolution_mask(elevations), groups;
	resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);

	dinf_flow_flats(flat_resolution_mask,groups,flowdirs);

	float_2d area(elevations);
	dinf_upslope_area(flowdirs, area);

	output_ascii_data(argv[2],area,8);

	return 0;
}
