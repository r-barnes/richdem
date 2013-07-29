#ifndef pit_fill_include
#define pit_fill_include
#include "data_structures.hpp"

void barnes_flood(float_2d &elevations);




//priority_flood_flowdirs
/**
  @brief  Determines D8 flow directions and implicitly fills pits.
  @author Richard Barnes (rbarnes@umn.edu)

    This version of Priority-Flood starts on the edges of the DEM and then
    works its way inwards using a priority queue to determine the lowest cell
    which has a path to the edge. The neighbours of this cell are given D8 flow
    directions which point to it. All depressions are implicitly filled and
    digital dams removed.

    Based on Metz 2011.

  @param[in]   &elevations  A grid of cell elevations
  @param[out]  &flowdirs    A grid of D8 flow directions

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM.

  @post
    1. **flowdirs** contains a D8 flow direction of each cell or a value
       _NO_FLOW_ for those cells which are not part of the DEM.
    2. **flowdirs** has no cells which are not part of a continuous flow
       path leading to the edge of the DEM.
*/
template <class T>
void priority_flood_flowdirs(const array2d<T> &elevations, char_2d &flowdirs){
  grid_cellzk_pq open;
  bool_2d closed;
  unsigned long processed_cells=0;
  ProgressBar progress;

  diagnostic("\n###Priority-Flood+Flow Directions\n");
  diagnostic("Setting up boolean flood array matrix...");
  closed.copyprops(elevations);
  closed.init(false);
  diagnostic("succeeded.\n");

  diagnostic("Setting up the flowdirs matrix...");
  flowdirs.copyprops(elevations);
  flowdirs.no_data=NO_FLOW;
  diagnostic("succeeded.\n");

  diagnostic_arg("The priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*((long)sizeof(grid_cellz))/1024/1024);
  diagnostic("Adding cells to the priority queue...");
  for(int x=0;x<elevations.width();x++){
    open.push_cell(x,0,elevations(x,0));
    open.push_cell(x,elevations.height()-1,elevations(x,elevations.height()-1));
    flowdirs(x,0)=3;
    flowdirs(x,elevations.height()-1)=7;
    closed(x,0)=true;
    closed(x,elevations.height()-1)=true;
  }
  for(int y=1;y<elevations.height()-1;y++){
    open.push_cell(0,y,elevations(0,y) );
    open.push_cell(elevations.width()-1,y,elevations(elevations.width()-1,y) );
    flowdirs(0,y)=1;
    flowdirs(elevations.width()-1,y)=5;
    closed(0,y)=true;
    closed(elevations.width()-1,y)=true;
  }
  diagnostic("succeeded.\n");

  flowdirs(0,0)=2;
  flowdirs(flowdirs.width()-1,0)=4;
  flowdirs(0,flowdirs.height()-1)=8;
  flowdirs(flowdirs.width()-1,flowdirs.height()-1)=6;

  const int d8_order[9]={0,1,3,5,7,2,4,6,8};
  diagnostic("%%Performing Priority-Flood+Flow Directions...\n");
  progress.start( elevations.width()*elevations.height() );
  while(open.size()>0){
    grid_cellz c=open.top();
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

      open.push_cell(nx,ny,elevations(nx,ny));
    }
    progress.update(processed_cells);
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
  diagnostic_arg("%ld cells processed.\n",processed_cells);
}



#endif
