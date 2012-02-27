#include "data_structures.h"
#include "utility.h"
#include "d8_methods.h"
#include "interface.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <omp.h>
#include <queue>

//D8 Directions
static int const dx[9]={0,-1,-1,0,1,1,1,0,-1};
static int const dy[9]={0,0,-1,-1,-1,0,1,1,1};
static int const inverse_flow[9]={0,5,6,7,8,1,2,3,4};
//std::string fd[9]={"·","←","↖","↑","↗","→","↘","↓","↙"};

int d8_FlowDir(const float_2d &elevations, const int x, const int y){
	float minimum_elevation=elevations(x,y);
	int flowdir=NO_FLOW;

	if (EDGE_GRID(x,y,elevations.size1(),elevations.size2())) return 0;
	if (elevations(x,y)==elevations.no_data) return d8_NO_DATA; //No data for this cell

	for(int n=1;n<=8;n++){
		if(!IN_GRID(x+dx[n],y+dy[n],elevations.size1(),elevations.size2())) continue;
		if(elevations(x+dx[n],y+dy[n])==elev_NO_DATA) continue;

		if(	elevations(x+dx[n],y+dy[n])<minimum_elevation || (elevations(x+dx[n],y+dy[n])==minimum_elevation && flowdir>0 && flowdir%2==0 && n%2==1) ){
			minimum_elevation=elevations(x+dx[n],y+dy[n]);
			flowdir=n;
		}
	}

	return flowdir;
}

int d8_flow_directions(const float_2d &elevations, char_2d &flowdirs){
	diagnostic_arg("The D8 flow directions will require approximately %ldMB of RAM.\n",elevations.size1()*elevations.size2()*sizeof(char)/1024/1024);
	diagnostic("Resizing flow directions matrix...");
	try{
		flowdirs.resize(elevations.size1(),elevations.size2(),false);
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Calculating D8 flow directions...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.size1();x++){
		progress_bar(x*omp_get_num_threads()*elevations.size2()*100/(elevations.size1()*elevations.size2()));
		for(int y=0;y<elevations.size2();y++)
			flowdirs(x,y)=d8_FlowDir(elevations,x,y);
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}








//234
//105
//876
bool does_cell_flow_into_me(const int x, const int y, int n, const char_2d &flowdirs){
	return (flowdirs(x,y)!=d8_NO_DATA && n==inverse_flow[flowdirs(x,y)]);
}

int d8_upslope_area(const char_2d &flowdirs){
	char_2d dependency;
	uint_2d area;
	std::queue<grid_cell*> sources;

	diagnostic_arg("The sources queue will require at most approximately %ldMB of RAM.\n",flowdirs.size1()*flowdirs.size2()*sizeof(grid_cell)/1024/1024);

	diagnostic_arg("The dependency matrix will require approximately %ldMB of RAM.\n",flowdirs.size1()*flowdirs.size2()*sizeof(char)/1024/1024);
	diagnostic("Resizing dependency matrix...");
	try{
		dependency.resize(flowdirs.size1(),flowdirs.size2(),false);
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic_arg("The area matrix will require approximately %ldMB of RAM.\n",flowdirs.size1()*flowdirs.size2()*sizeof(unsigned int)/1024/1024);
	diagnostic("Resizing the area matrix...");
	try{
		area.resize(flowdirs.size1(),flowdirs.size2(),false);
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Calculating dependency matrix...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<flowdirs.size1();x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.size2()*100/(flowdirs.size1()*flowdirs.size2()));
		for(int y=0;y<flowdirs.size2();y++){
			dependency(x,y)=0;
			for(int n=1;n<=8;n++){
				if(!IN_GRID(x+dx[n],y+dy[n],flowdirs.size1(),flowdirs.size2())) continue;
				if(does_cell_flow_into_me(x+dx[n],y+dy[n],n,flowdirs))
					dependency(x,y)++;
			}
		}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	diagnostic("Locating source cells...\n");
	progress_bar(-1);
	for(int x=0;x<flowdirs.size1();x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.size2()*100/(flowdirs.size1()*flowdirs.size2()));
		for(int y=0;y<flowdirs.size2();y++)
			if(flowdirs(x,y)!=d8_NO_DATA && dependency(x,y)==0)
				sources.push(new grid_cell(x,y));
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	diagnostic("Calculating up-slope areas...\n");
	progress_bar(-1);
	long int ccount=0;
	while(sources.size()>0){
		grid_cell *c=sources.front();
		sources.pop();

		ccount++;
		progress_bar(ccount*100/flowdirs.data_cells);

		area(c->x,c->y)=1;
		for(int n=0;n<8;n++){
			if(!IN_GRID(c->x+dx[n],c->y+dy[n],flowdirs.size1(),flowdirs.size2())) continue;
			if(does_cell_flow_into_me(c->x+dx[n],c->y+dy[n],n,flowdirs))
				area(c->x,c->y)+=area(c->x+dx[n],c->y+dy[n]);
		}
		if(flowdirs(c->x,c->y)!=d8_NO_DATA && flowdirs(c->x,c->y)!=0)
			if( (--dependency(c->x+dx[flowdirs(c->x,c->y)],c->y+dy[flowdirs(c->x,c->y)]))==0)
				sources.push(new grid_cell(c->x+dx[flowdirs(c->x,c->y)],c->y+dy[flowdirs(c->x,c->y)]));

		delete c;
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}
