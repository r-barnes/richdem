#ifndef _d8_methods_included
#define _d8_methods_included

#include "data_structures.h"
#include "interface.h"
#include <omp.h>

#define SLOPE_RISERUN	1
#define SLOPE_PERCENT	2
#define SLOPE_RADIAN	3
#define SLOPE_DEGREE	4

void d8_upslope_area(const char_2d &flowdirs, uint_2d &area);
void d8_slope(const float_2d &elevations, float_2d &slopes, int slope_type=SLOPE_RISERUN);

//234
//105
//876
template<class T>
int d8_FlowDir(const array2d<T> &elevations, const int x, const int y){
	float minimum_elevation=elevations(x,y);
	int flowdir=NO_FLOW;

	if (EDGE_GRID(x,y,elevations.width(),elevations.height())){
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
void d8_flow_directions(const array2d<T> &elevations, char_2d &flowdirs, bool init=true){
	diagnostic_arg("The D8 flow directions will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Resizing flow directions matrix...");
	flowdirs.resize(elevations.width(),elevations.height(),!init);
	diagnostic("succeeded.\n");

	if(init){
		diagnostic("Setting no_data value on flowdirs matrix...");
		flowdirs.no_data=d8_NO_DATA;
		diagnostic("succeeded.\n");

		diagnostic("Initializing D8 flow directions...\n");
		flowdirs.init(NO_FLOW);
		diagnostic("\tsucceeded.\n");
	}

	diagnostic("Calculating D8 flow directions...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++){
		progress_bar(x*omp_get_num_threads()*elevations.height()*100/(elevations.width()*elevations.height()));
		for(int y=0;y<elevations.height();y++){
			if(elevations(x,y)==elevations.no_data){
				if(init)
					flowdirs(x,y)=d8_NO_DATA;
				else
					continue;	//Masking
			} else if(flowdirs(x,y)==NO_FLOW)
				flowdirs(x,y)=d8_FlowDir(elevations,x,y);
		}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}

#endif
