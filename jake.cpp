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

void jacob_wetland_metric(const float_2d &smooted_cti, const int_2d &hydric_soils, const float_2d &smoothed_percent_slope, const float_2d &smoothed_profile_curvature, float_2d &result){
	Timer timer;

	diagnostic("\n###Jacob's Wetland Metric\n");

	diagnostic_arg("Jacob's Wetland Metric matrix will require approximately %ldMB of RAM.\n", flow_accumulation.width()*flow_accumulation.height()*((long)sizeof(float))/1024/1024);
	diagnostic("Setting up Jacob's Wetland Metric matrix...");
	result.copyprops(smoothed_cit);
	result.no_data=-1;	//TODO
	diagnostic("succeeded.\n");

	diagnostic("Calculating Jacob's Wetland Metric...\n");
	timer.start();
	#pragma omp parallel for collapse(2)
	for(int x=0;x<result.width();x++)
	for(int y=0;y<result.height();y++){
		if(smoothed_cti(x,y)==smoothed_cti.no_data){
			result(x,y)=result.no_data;
			continue;
		}
		result(x,y)=(-4)+(0.5*smoothed_cti(x,y))+(hydric_soils)-(0.2*smoothed_percent_slope(x,y))-(0.5*smoothed_profile_curvature(x,y));
	}
	diagnostic_arg("succeeded in %lfs.\n",timer.lap());
}

int main(int argc, char **argv){
	Timer running_calc_time,running_io_time,total_time;

	total_time.start();

	float_2d elevations;
	running_io_time.start();
	load_ascii_data(argv[1],elevations);
	running_io_time.stop();

	running_calc_time.start();
	barnes_flood(elevations);

	float_2d flowdirs(elevations);
	dinf_flow_directions(elevations,flowdirs);

	int_2d flat_resolution_mask(elevations), groups;
	resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);
	dinf_flow_flats(flat_resolution_mask,groups,flowdirs);
	flat_resolution_mask.clear();
	groups.clear();

	float_2d area;
	dinf_upslope_area(flowdirs, area);
	flowdirs.clear();

	float_2d percent_slope;
	d8_slope(elevations, percent_slope, TATTRIB_SLOPE_PERCENT);

	float_2d CTI;
	d8_CTI(area, percent_slope, CTI);
	area.clear();

	float_2d profile_curvature;
	d8_profile_curvature(elevations, profile_curvature);

	elevations.clear();

	running_calc_time.stop();
	total_time.stop();

	diagnostic_arg("Total time was: %lfs\n", total_time.accumulated());
	diagnostic_arg("Time spent in computation: %lfs\n",running_calc_time.accumulated());
	diagnostic_arg("Time spent in I/O: %lfs\n",running_io_time.accumulated());

	return 0;
}

