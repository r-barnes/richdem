#include "data_structures.h"
#include "utility.h"
#include "d8_methods.h"
#include "interface.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <omp.h>

//D8 Directions
static int const dx[9]={0,-1,-1,0,1,1,1,0,-1};
static int const dy[9]={0,0,-1,-1,-1,0,1,1,1};
static int const inverse_flow[9]={0,5,6,7,8,1,2,3,4};
//std::string fd[9]={"·","←","↖","↑","↗","→","↘","↓","↙"};

typedef struct grid_cell_type {
	int x;
	int y;
	grid_cell_type(int x0, int y0){
		x=x0;
		y=y0;
	}
} grid_cell;

int d8_FlowDir(float_2d &elevations, int x, int y){
	float minimum_elevation=elevations(x,y);
	int flowdir=NO_FLOW;
//	flat tflat={x,y};
//	bool has_different=false,has_same=false;

	if (elevations(x,y)==elev_NO_DATA) return d8_NO_DATA; //Missing data

	for(int n=1;n<=8;n++){
		if(!IN_GRID(x+dx[n],y+dy[n],elevations.size1(),elevations.size2())) continue;
		if(elevations(x+dx[n],y+dy[n])==elev_NO_DATA) continue;
//		if(elevations(x+dx[n],y+dy[n])!=elevations(x,y))
//			has_different=true;
//		else
//			has_same=true;
		if(	elevations(x+dx[n],y+dy[n])<minimum_elevation ||
			(elevations(x+dx[n],y+dy[n])==minimum_elevation && flowdir>0 && 				flowdir%2==0 && n%2==1)
		)
		{
			minimum_elevation=elevations(x+dx[n],y+dy[n]);
			flowdir=n;
		}
	}

//	if(has_same && has_different){
//		tflat.x=x;
//		tflat.y=y;
//		flats.push(tflat);
//	}

	return flowdir;
}

int d8_flow_directions(float_2d elevations, char_2d flowdirs){
	diagnostic_arg("The D8 flow directions will require approximately %ldMB of RAM.\n",elevations.size1()*elevations.size2()*sizeof(char)/1024/1024);
	diagnostic("Resizing flow directions matrix...");
	try{
		flowdirs.resize(elevations.size1(),elevations.size2());
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Calculating D8 flow directions...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.size1();x++){
		progress_bar(x*omp_get_num_threads()*elevations.size2()/(elevations.size1()*elevations.size2()));
		for(int y=0;y<elevations.size2();y++){
			if(EDGE_GRID(x,y,elevations.size1(),elevations.size2())){ //Edges of grid undefined
				flowdirs(x,y)=NO_FLOW;
				continue;
			}
			if(flowdirs(x,y)!=0) continue; //Flow direction already set
				flowdirs(x,y)=d8_FlowDir(elevations,x,y);
		}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}

/*int d8_upslope_area(float_2d &elevations, char_2d &flowdirs){
	char_2d dependents;
	uint_2d area;
	vector<grid_cell*> sources;

	diagnostic_arg("Dependents count matrix and area matrix will require approximately %ldMB of RAM.\n",elevations.size1()*elevations.size2()*(sizeof(char)+sizeof(unsigned int))/1024/1024);
	diagnostic("Resizing count & area matricies...");
	try{
		dependents.resize(elevations.size1(),elevations.size2());
		area.resize(elevations.size1(),elevations.size2());
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");	

	diagnostic("Searching for area sources...");
	for(int x=0;x<elevations.size1();x++)
	for(int y=0;y<elevations.size2();y++){
		int d=0;
		for(int n=1;n<=8;n++)
			if(inverse_flow[flowdirs(x+dx[n],y+dy[n])]==n)
				d++;
		dependents(x,y)=d;
		if(d==0)
			sources.push(
	}
}*/
