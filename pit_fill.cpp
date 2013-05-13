#include "utility.hpp"
#include "data_structures.hpp"
#include "interface.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <queue>
#include <stack>

//barnes_flood
/**
  @brief  Floods DEM inwards from the edges, filling all pits and removing all digital dams.
  @author Richard Barnes

      The BarnesFlood starts on the edges of the DEM and then works its way
      inwards using a priority queue to determine the lowest cell which has
      a path to the edge. The neighbours of this cell are added to the priority
      queue if they are higher. If they are lower, they are added to a "pit"
      queue which is used to flood pits. Cells which are higher than a pit being
      filled are added to the priority queue. In this way, pits are filled without
      incurring the expense of the priority queue.

  @param[in,out]  &elevations
    A grid of cell elevations
*/
void barnes_flood(float_2d &elevations){
  grid_cellz_pq open;
//  std::queue<grid_cellz> meander;
  std::stack<grid_cellz, std::vector<grid_cellz> > meander;
  bool_2d closed;
  unsigned long processed_cells=0;
  unsigned long pitc=0;
  ProgressBar progress;

  diagnostic("\n###Barnes Flood\n");
  diagnostic("Setting up boolean flood array matrix...");
  closed.copyprops(elevations);
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
    open.push(grid_cellz(0,y,elevations(0,y)  ));
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
    } else {
      c=open.top();
      open.pop();
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
        if(elevations(nx,ny)<c.z){
          ++pitc;
          elevations(nx,ny)=c.z;
        }
        meander.push(grid_cellz(nx,ny,c.z));
      } else
        open.push(grid_cellz(nx,ny,elevations(nx,ny)));
    }
    progress.update(processed_cells);
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
  diagnostic_arg("%ld cells processed. %ld in pits.\n",processed_cells,pitc);
}







//d8_edge_flow
/**
  @brief  Helper function which returns a flow direction pointing to an edge if (x,y) is a non-no_data edge cell
  @author Richard Barnes

  @param[in]  x             x-coordinate of the cell
  @param[in]  y             y-coordinate of the cell
  @param[in]  &elevations   A grid of cell elevations
  @param[in]  &flowdirs     A grid of D8 flow directions

  @pre  (x,y) must be an edge cell or an error will be thrown

  @return If (x,y) is a non-no_data edge cell, then a D8 direction is returned which points to the edge
*/
char d8_edge_flow(int x, int y, const float_2d &elevations, const char_2d &flowdirs){
  if(!elevations.edge_grid(x,y))  //NOTE: Shouldn't happen
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
  else  //NOTE: Avoids a compiler warning. Control can now never reach end of a non-void function.
    throw "Barnes Flood+Flow Directions tried to initialize with a non-edge cell!";
}




//barnes_flood_flowdirs
/**
  @brief  Determines D8 flow directions by flooding inwards, pits are implicitly carved to drainage points. Based on Metz 2011.
  @author Richard Barnes

  @param[in]   &elevations  A grid of cell elevations
  @param[out]  &flowdirs    A grid of D8 flow directions

  @post \pname{flowdirs} takes the properties and dimensions of \pname{elevations}
*/
void barnes_flood_flowdirs(const float_2d &elevations, char_2d &flowdirs){
  grid_cellzk_pq open;
  bool_2d closed;
  unsigned long processed_cells=0;
  int cell_num=0;
  ProgressBar progress;

  diagnostic("\n###Barnes Flood+Flow Directions\n");
  diagnostic("Setting up boolean flood array matrix...");
  closed.copyprops(elevations);
  closed.init(false);
  diagnostic("succeeded.\n");

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
    open.push(grid_cellzk(0,y,elevations(0,y),++cell_num  ));
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
