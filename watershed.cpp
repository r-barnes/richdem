#include "data_structures.h"
#include "interface.h"
#include "utility.h"
#include <vector>
#include <map>
#include <stack>
#include <queue>

//FindWatersheds works in the same way as BarnesFlood, save that it labels watersheds, working inwards from the edges of the DEM. This is helpful for checking flow accumulation
//Procedure:	FindWatersheds
//Description:
//		Same as BarnesFlood. Labels starts out as no_data. If it is found that a
//		no_data labels cell coincides with a non-no_data elevations cell, then this
//		is the beginning of a new watershed. Cells which are flooded from a labeled
//		cell take on that cell's label
//Inputs:
//		elevations		A 2D array of cell elevations
//		labels			A 2D array to hold the watershed labels
//Requirements:
//		None
//Effects:
//		"labels" will be altered to contain the labels of the watersheds
//Returns:
//		None
void find_watersheds(float_2d &elevations, int_2d &labels){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cellz_compare> open;
	std::stack<grid_cellz, std::vector<grid_cellz> > meander;
	bool_2d closed;
	unsigned long processed_cells=0;
	unsigned long pitc=0,openc=0;
	int clabel=1;	//TODO: Thought this was more clear than zero in the results.

	diagnostic("\n###Barnes Flood\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*((long)sizeof(bool))/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(false);
	diagnostic("succeeded.\n");

	diagnostic_arg("The labels matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*((long)sizeof(bool))/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	labels.copyprops(elevations);
	labels.no_data=-1;
	labels.init(labels.no_data);
	diagnostic("succeeded.\n");

	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*((long)sizeof(grid_cellz))/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.width();x++){
		open.push(grid_cellz(x,0,elevations(x,0) ));
		open.push(grid_cellz(x,elevations.height()-1,elevations(x,elevations.height()-1) ));
		closed(x,0)=true;
		closed(x,elevations.height()-1)=true;
	}
	for(int y=1;y<elevations.height()-1;y++){
		open.push(grid_cellz(0,y,elevations(0,y)	));
		open.push(grid_cellz(elevations.width()-1,y,elevations(elevations.width()-1,y) ));
		closed(0,y)=true;
		closed(elevations.width()-1,y)=true;
	}
	diagnostic("succeeded.\n");

	diagnostic("Performing the Barnes Flood...\n");
	progress_bar(-1);
	while(open.size()>0 || meander.size()>0){
		grid_cellz c;
		if(meander.size()>0){
			c=meander.top();
			meander.pop();
			pitc++;
		} else {
			c=open.top();
			open.pop();
			openc++;
		}
		processed_cells++;

		if(labels(c.x,c.y)==labels.no_data && elevations(c.x,c.y)!=elevations.no_data)	//Implies a cell without a label which borders the edge of the DEM or a region of no_data
			labels(c.x,c.y)=clabel++;

		for(int n=1;n<=8;n++){
			int nx=c.x+dx[n];
			int ny=c.y+dy[n];
			if(!elevations.in_grid(nx,ny)) continue;
			if(closed(nx,ny)) 
				continue;

			if(labels(c.x,c.y)!=labels.no_data)
				labels(nx,ny)=labels(c.x,c.y);

			closed(nx,ny)=true;
			if(elevations(nx,ny)<=c.z){
				elevations(nx,ny)=c.z;
				meander.push(grid_cellz(nx,ny,c.z));
			} else
				open.push(grid_cellz(nx,ny,elevations(nx,ny)));
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed. %ld in pits, %ld not in pits.\n",processed_cells,pitc,openc);
}


void watershed_area(const int_2d &labels){
	std::map<int, int> wsheds;
	for(int x=0;x<labels.width();x++)
	for(int y=0;y<labels.height();y++)
		if(labels(x,y)!=labels.no_data)
			wsheds[labels(x,y)]++;

	for(std::map<int, int>::iterator i=wsheds.begin();i!=wsheds.end();i++)
		printf("Watershed %d has area %d\n",(*i).first,(*i).second);
}
