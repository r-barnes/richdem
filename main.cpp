#include "utility.h"
#include "data_structures.h"
#include "load_data.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"

int main(int argc, char **argv){
	float_2d elevations;
	float_2d flowdirs;
	load_ascii_data(argv[1],elevations);

//	pit_fill_yonghe2009(elevations);
//	pit_fill_wang(elevations);
	dinf_flow_directions(elevations,flowdirs);
}
