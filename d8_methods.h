#ifndef _d8_methods_included
#define _d8_methods_included

#include "data_structures.h"
#include "interface.h"
#ifdef _OPENMP
	#include <omp.h>
#endif

#define SLOPE_RISERUN	1
#define SLOPE_PERCENT	2
#define SLOPE_RADIAN	3
#define SLOPE_DEGREE	4

#define TATTRIB_SLOPE				1
#define TATTRIB_CURVATURE			2
#define TATTRIB_PLANFORM_CURVATURE	3
#define TATTRIB_PROFILE_CURVATURE	4
#define TATTRIB_ASPECT				5
#define TATTRIB_SLOPE_RISERUN		6
#define TATTRIB_SLOPE_PERCENT		7
#define TATTRIB_SLOPE_RADIAN		8
#define TATTRIB_SLOPE_DEGREE		9

void d8_flow_flats(const int_2d &flat_resolution_mask, const int_2d &groups, char_2d &flowdirs);
void d8_upslope_area(const char_2d &flowdirs, int_2d &area);
void d8_slope(const float_2d &elevations, float_2d &slopes, int slope_type=TATTRIB_SLOPE_RISERUN);
void d8_aspect(const float_2d &elevations, float_2d &aspects);
void d8_curvature(const float_2d &elevations, float_2d &curvatures);
void d8_profile_curvature(const float_2d &elevations, float_2d &profile_curvatures);
void d8_planform_curvature(const float_2d &elevations, float_2d &planform_curvatures);

void find_watersheds(float_2d &elevations, int_2d &labels, bool alter_elevations=false);
void watershed_area(const int_2d &labels);

void d8_SPI(const float_2d &flow_accumulation, const float_2d &percent_slope, float_2d &result);
void d8_CTI(const float_2d &flow_accumulation, const float_2d &percent_slope, float_2d &result);

//234
//105
//876
template<class T>
int d8_FlowDir(const array2d<T> &elevations, const int x, const int y){
	T minimum_elevation=elevations(x,y);
	int flowdir=NO_FLOW;

	if (elevations.edge_grid(x,y)){
		if(x==0 && y==0)
			return 2;
		else if(x==0 && y==elevations.height()-1)
			return 8;
		else if(x==elevations.width()-1 && y==0)
			return 4;
		else if(x==elevations.width()-1 && y==elevations.height()-1)
			return 6;
		else if(x==0)
			return 1;
		else if(x==elevations.width()-1)
			return 5;
		else if(y==0)
			return 3;
		else if(y==elevations.height()-1)
			return 7;
	}

	//NOTE: Since the very edges of the DEM are defined to always flow outwards, if they have defined elevations, it is not necessary to check of a neighbour is IN_GRID in the following
	//NOTE: It is assumed that the no_data datum is an extremely negative number, such that all water which makes it to the edge of the DEM's region of defined elevations is sucked directly off the grid, rather than piling up on the edges.
	for(int n=1;n<=8;n++)
		if(	elevations(x+dx[n],y+dy[n])<minimum_elevation || (elevations(x+dx[n],y+dy[n])==minimum_elevation && flowdir>0 && flowdir%2==0 && n%2==1) ){
			minimum_elevation=elevations(x+dx[n],y+dy[n]);
			flowdir=n;
		}

	return flowdir;
}

template<class T>
void d8_flow_directions(const array2d<T> &elevations, char_2d &flowdirs){
	diagnostic_arg("The D8 flow directions will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*((long)sizeof(char))/1024/1024);
	diagnostic("Resizing flow directions matrix...");
	flowdirs.resize(elevations.width(),elevations.height(),false);
	diagnostic("succeeded.\n");

	diagnostic("Setting no_data value on flowdirs matrix...");
	flowdirs.no_data=d8_NO_DATA;
	diagnostic("succeeded.\n");

	diagnostic("Initializing D8 flow directions...");
	flowdirs.init(NO_FLOW);
	diagnostic("succeeded.\n");

	diagnostic("%%Calculating D8 flow directions...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++){
		progress_bar(x*elevations.height()*100/(elevations.width()*elevations.height()));
		for(int y=0;y<elevations.height();y++)
			if(elevations(x,y)==elevations.no_data)
				flowdirs(x,y)=flowdirs.no_data;
			else
				flowdirs(x,y)=d8_FlowDir(elevations,x,y);
	}
	diagnostic_arg(SUCCEEDED_IN,progress_bar(-1));
}

#endif
