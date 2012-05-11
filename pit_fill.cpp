#include "utility.h"
#include "data_structures.h"
#include "interface.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <queue>
#include <stack>

//Procedure:	BarnesFlood
//Description:
//		The BarnesFlood starts on the edges of the DEM and then works its way
//      inwards using a priority queue to determine the lowest cell which has
//      a path to the edge. The neighbours of this cell are added to the priority
//      queue if they are higher. If they are lower, they are added to a "pit"
//      queue which is used to flood pits. Cells which are higher than a pit being
//      filled are added to the priority queue. In this way, pits are filled without
//      incurring the expense of the priority queue.
//Inputs:
//		elevations		A 2D array of cell elevations
//Requirements:
//		None
//Effects:
//		"elevations" will be altered to contain a pit-filled version of the original
//Returns:
//		None
void barnes_flood(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cellz_compare> open;
//	std::queue<grid_cellz> meander;
	std::stack<grid_cellz, std::vector<grid_cellz> > meander;
	bool_2d closed;
	unsigned long processed_cells=0;
	unsigned long pitc=0,openc=0;

	diagnostic("\n###Barnes Flood\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*((long)sizeof(bool))/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(false);
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

	diagnostic("%%Performing the Barnes Flood...\n");
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

		for(int n=1;n<=8;n++){
			int nx=c.x+dx[n];
			int ny=c.y+dy[n];
			if(!elevations.in_grid(nx,ny)) continue;
			if(closed(nx,ny)) 
				continue;

			closed(nx,ny)=true;
			if(elevations(nx,ny)<=c.z){
				elevations(nx,ny)=c.z;
				meander.push(grid_cellz(nx,ny,c.z));
			} else
				open.push(grid_cellz(nx,ny,elevations(nx,ny)));
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg(SUCCEEDED_IN,progress_bar(-1));
	diagnostic_arg("%ld cells processed. %ld in pits, %ld not in pits.\n",processed_cells,pitc,openc);
}
