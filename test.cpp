#include "utility.h"
#include "data_structures.h"
#include "data_io.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"
#include "flat_resolution.h"
#include "debug.h"
#include <string>
#include <sys/time.h>

int main(int argc, char **argv){
	float_2d elevations,elevations2;
	load_ascii_data(argv[1],elevations);

	elevations2=elevations;

	int_2d labels,labels2;
	find_watersheds(elevations,labels);
	find_watersheds_test(elevations2,labels2);

	diagnostic_arg("Diff was: %f\n",avg_diff(labels,labels2));
}
