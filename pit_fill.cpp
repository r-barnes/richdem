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
	ProgressBar progress;

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
	progress.start( elevations.width()*elevations.height() );
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
		progress.update(processed_cells);
	}
	diagnostic_arg(SUCCEEDED_IN,progress.stop());
	diagnostic_arg("%ld cells processed. %ld in pits, %ld not in pits.\n",processed_cells,pitc,openc);
}







//Procedure:	d8_edge_flow
//Description:
//		Determines what flow direction a cell on the edge of a DEM should have
//		that flow direction should always be outward. This does that mapping
//Inputs:
//		elevations		A 2D array of cell elevations
//		flowdirs		A 2D array of D8 flow directions
//Requirements:
//		None
//Effects:
//		None
//Returns:
//		flowdirs.no_data if the elevation of the cell is elevations.no_data
//		Otherwise, a D8 direction pointing off the grid
char d8_edge_flow(int x, int y, const float_2d &elevations, const char_2d &flowdirs){
	if(!elevations.edge_grid(x,y))	//NOTE: Shouldn't happen
		throw "Barnes Flood+Flow Directions tried to initialize with a non-edge cell!";
	else if(elevations(x,y)==elevations.no_data)
		return flowdirs.no_data;
	else if(x==0 && y==0)
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
	else	//NOTE: Avoids a compiler warning. Control can now never reach end of a non-void function.
		throw "Barnes Flood+Flow Directions tried to initialize with a non-edge cell!";
}

//Procedure:	BarnesFloodFlowDirs
//Description:
//		Floods DEM inwards from its edges building up flow directions as it
//		goes. The result is a D8 flowdirs array in which pits have been carved
//		towards a drainage point. Based on Metz 2011.
//Inputs:
//		elevations		A 2D array of cell elevations
//		flowdirs		A 2D array for storing D8 flow directions
//Requirements:
//		None
//Effects:
//		"flowdirs" will be altered to contain D8 flow directions for each cell
//Returns:
//		None
void barnes_flood_flowdirs(const float_2d &elevations, char_2d &flowdirs){
	std::priority_queue<grid_cellzk, std::vector<grid_cellzk>, grid_cellzk_compare> open;
	bool_2d closed;
	unsigned long processed_cells=0;
	int cell_num=0;
	ProgressBar progress;

	diagnostic("\n###Barnes Flood+Flow Directions\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n", elevations.width()*elevations.height()*((long)sizeof(bool))/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(false);
	diagnostic("succeeded.\n");

	diagnostic_arg("The flowdirs matrix will require approximately %ldMB of RAM.\n", elevations.width()*elevations.height()*((long)sizeof(bool))/1024/1024);
	diagnostic("Setting up the flowdirs matrix...");
	flowdirs.copyprops(elevations);
	flowdirs.no_data=NO_FLOW;
	diagnostic("succeeded.\n");

	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*((long)sizeof(grid_cellz))/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.width();x++){
		open.push(grid_cellzk(x,0,elevations(x,0),++cell_num ));
		open.push(grid_cellzk(x,elevations.height()-1,elevations(x,elevations.height()-1), ++cell_num ));
		flowdirs(x,0)=d8_edge_flow(x,0,elevations,flowdirs);
		flowdirs(x,elevations.height()-1)=d8_edge_flow(x, elevations.height()-1, elevations, flowdirs);
		closed(x,0)=true;
		closed(x,elevations.height()-1)=true;
	}
	for(int y=1;y<elevations.height()-1;y++){
		open.push(grid_cellzk(0,y,elevations(0,y),++cell_num	));
		open.push(grid_cellzk(elevations.width()-1,y,elevations(elevations.width()-1,y), ++cell_num ));
		flowdirs(0,y)=d8_edge_flow(0,y,elevations,flowdirs);
		flowdirs(elevations.width()-1,y)=d8_edge_flow(elevations.width()-1, y, elevations, flowdirs);
		closed(0,y)=true;
		closed(elevations.width()-1,y)=true;
	}
	diagnostic("succeeded.\n");

	const int d8_order[9]={0,1,3,5,7,2,4,6,8};
	diagnostic("%%Performing the Barnes Flood+Flow Directions...\n");
	progress.start( elevations.width()*elevations.height() );
	while(open.size()>0){
		grid_cellzk c=open.top();
		open.pop();
		processed_cells++;

		for(int no=1;no<=8;no++){
			int n=d8_order[no];
			int nx=c.x+dx[n];
			int ny=c.y+dy[n];
			if(!elevations.in_grid(nx,ny)) continue;
			if(closed(nx,ny)) 
				continue;

			closed(nx,ny)=true;

			if(elevations(nx,ny)==elevations.no_data)
				flowdirs(nx,ny)=flowdirs.no_data;
			else
				flowdirs(nx,ny)=inverse_flow[n];

			open.push(grid_cellzk(nx,ny,elevations(nx,ny),++cell_num));
		}
		progress.update(processed_cells);
	}
	diagnostic_arg(SUCCEEDED_IN,progress.stop());
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}
