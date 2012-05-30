#include "utility.h"
#include "data_structures.h"
#include "data_io.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"
#include "flat_resolution.h"
#include "watershed.h"
#include "debug.h"
#include <string>
#include <sys/time.h>
#include "unit_test.h"

int main(){
	float_2d elevations;
	load_ascii_data("unit_test/bf03.dem", elevations);

	{
		float_2d slope_percent, unit_slope_percent;
		d8_slope(elevations, slope_percent, TATTRIB_SLOPE_PERCENT);
		load_ascii_data("unit_test/bf03_slope_percent.txt", unit_slope_percent);
		printf("Average SLOPE PERCENT difference: %lf\n",unit_avg_diff(slope_percent, unit_slope_percent));
	}
	
}
