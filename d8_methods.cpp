#include "data_structures.h"
#include "utility.h"
#include "d8_methods.h"
#include "interface.h"
#include <omp.h>

int dx[9]={0,-1,-1,0,1,1,1,0,-1};
int dy[9]={0,0,-1,-1,-1,0,1,1,1};
//std::string fd[9]={"·","←","↖","↑","↗","→","↘","↓","↙"};

int d8_FlowDir(float_2d &elevations, int x, int y){
	float minimum_elevation=elevations(x,y);
	int flowdir=NO_FLOW;
//	flat tflat={x,y};
//	bool has_different=false,has_same=false;

	if (elevations(x,y)==elev_NO_DATA) return d8_NO_DATA; //Missing data

	for(int n=1;n<=8;n++){
		if(elevations(x+dx[n],y+dy[n])==elev_NO_DATA) continue;
		if(!IN_GRID(x+dx[n],y+dy[n],elevations.size1(),elevations.size2())) continue;
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
	for(long int x=0;x<elevations.size1();x++){
		#pragma omp master
		{
		progress_bar(100*omp_get_num_threads()*x*elevations.size2()/(elevations.size1()*elevations.size2())); //Todo: Should I check to see if ftell fails here?
		}
		for(long int y=0;y<elevations.size2();y++){
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
