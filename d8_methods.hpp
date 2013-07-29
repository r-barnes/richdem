#ifndef _d8_methods_included
#define _d8_methods_included

#include "data_structures.hpp"
#include <queue>
#include "interface.hpp"
#include <algorithm>
#ifdef _OPENMP
  #include <omp.h>
#endif

/// Used with #d8_terrain_attribute to get an ill-defined curvature thing
/// (TODO)
#define TATTRIB_CURVATURE           2

/// Used with #d8_terrain_attribute to get planform curvature as per
/// Zevenbergen and Thorne 1987 */
#define TATTRIB_PLANFORM_CURVATURE  3

/// Used with #d8_terrain_attribute to get profile curvature as per
/// Zevenbergen and Thorne 1987 */
#define TATTRIB_PROFILE_CURVATURE   4

/// Used with #d8_terrain_attribute to get aspect as per Horn 1981
#define TATTRIB_ASPECT              5

/// Used with #d8_terrain_attribute to get slope as per Horn 1981
#define TATTRIB_SLOPE_RISERUN       6

/// Used with #d8_terrain_attribute to get slope as per Horn 1981, the value
/// is multiplied by 100
#define TATTRIB_SLOPE_PERCENT       7

/// Used with #d8_terrain_attribute to get slope as per Horn 1981, an arc
/// tangent of the value is taken
#define TATTRIB_SLOPE_RADIAN        8

/// Used with #d8_terrain_attribute to get slope as per Horn 1981, an arc 
/// tangent of the value is taken and converted to degrees
#define TATTRIB_SLOPE_DEGREE        9

void d8_slope(
  const float_2d &elevations,
  float_2d &slopes,
  int slope_type=TATTRIB_SLOPE_RISERUN
);

void d8_profile_curvature(
  const float_2d &elevations,
  float_2d &profile_curvatures
);

void d8_planform_curvature(
  const float_2d &elevations,
  float_2d &planform_curvatures
);

void find_watersheds(
  float_2d &elevations,
  int_2d &labels,
  bool alter_elevations=false
);

void d8_aspect(const float_2d &elevations, float_2d &aspects);
void d8_curvature(const float_2d &elevations, float_2d &curvatures);
void watershed_area(const int_2d &labels);

void d8_SPI(
  const float_2d &flow_accumulation,
  const float_2d &percent_slope,
  float_2d &result
);

void d8_CTI(
  const float_2d &flow_accumulation,
  const float_2d &percent_slope,
  float_2d &result
);

//234
//105
//876
//d8_FlowDir
/**
  @brief  Calculates the D8 flow direction of a cell
  @author Richard Barnes (rbarnes@umn.edu)

  This calculates the D8 flow direction of a cell using the D8
  neighbour system, as defined in utility.h. Cells on the edge
  of the grid flow off the nearest edge.

  Helper function for d8_flow_directions().

  @param[in]  &elevations  A DEM
  @param[in]  x            x coordinate of cell
  @param[in]  y            y coordinate of cell

  @returns The D8 flow direction of the cell
*/
template<class T>
static int d8_FlowDir(const array2d<T> &elevations, const int x, const int y){
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

  /*NOTE: Since the very edges of the DEM are defined to always flow outwards,
  if they have defined elevations, it is not necessary to check if a neighbour
  is IN_GRID in the following
  NOTE: It is assumed that the no_data datum is an extremely negative
  number, such that all water which makes it to the edge of the DEM's region
  of defined elevations is sucked directly off the grid, rather than piling up
  on the edges.*/
  for(int n=1;n<=8;n++)
    if(
      elevations(x+dx[n],y+dy[n])<minimum_elevation
      || (elevations(x+dx[n],y+dy[n])==minimum_elevation
            && flowdir>0 && flowdir%2==0 && n%2==1)
    ){
      minimum_elevation=elevations(x+dx[n],y+dy[n]);
      flowdir=n;
    }

  return flowdir;
}



//d8_flow_directions
/**
  @brief  Calculates the D8 flow directions of a DEM
  @author Richard Barnes (rbarnes@umn.edu)

  This calculates the D8 flow directions of a DEM. Its argument
  'flowdirs' will return a grid with flow directions using the D8
  neighbour system, as defined in utility.h. The choice of data type
  for array2d must be able to hold exact values for all neighbour
  identifiers (usually [-1,7]).

  Uses d8_FlowDir() as a helper function.

  @todo                    Combine dinf and d8 neighbour systems

  @param[in]  &elevations  A DEM
  @param[out] &flowdirs    Returns the flow direction of each cell
*/
template<class T, class U>
void d8_flow_directions(
  const array2d<T> &elevations,
  array2d<U> &flowdirs
){
  ProgressBar progress;

  diagnostic("Setting up the flow directions matrix...");
  flowdirs.copyprops(elevations);
  flowdirs.init(NO_FLOW);
  flowdirs.no_data=d8_NO_DATA;
  diagnostic("succeeded.\n");

  diagnostic("%%Calculating D8 flow directions...\n");
  progress.start( elevations.width()*elevations.height() );
  #pragma omp parallel for
  for(int x=0;x<elevations.width();x++){
    progress.update( x*elevations.height() );
    for(int y=0;y<elevations.height();y++)
      if(elevations(x,y)==elevations.no_data)
        flowdirs(x,y)=flowdirs.no_data;
      else
        flowdirs(x,y)=d8_FlowDir(elevations,x,y);
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
}




//234
//105
//876
//d8_masked_FlowDir
/**
  @brief  Helper function to d8_flow_flats()
  @author Richard Barnes (rbarnes@umn.edu)

  This determines a cell's flow direction, taking into account flat membership.
  It is a helper function to d8_flow_flats()

  @param[in]  &flat_mask      A mask from resolve_flats_barnes()
  @param[in]  &groups         A grouping from resolve_flats_barnes()
  @param[in]  x               x coordinate of cell
  @param[in]  y               y coordinate of cell

  @returns    The flow direction of the cell
*/
static int d8_masked_FlowDir(
  const int_2d &flat_mask,
  const int_2d &groups,
  const int x,
  const int y
){
  int minimum_elevation=flat_mask(x,y);
  int flowdir=NO_FLOW;

  //It is safe to do this without checking to see that (nx,ny) is within
  //the grid because we only call this function on interior cells
  for(int n=1;n<=8;n++){
    int nx=x+dx[n];
    int ny=y+dy[n];
    if( groups(nx,ny)!=groups(x,y))
      continue;
    if(  flat_mask(nx,ny)<minimum_elevation || (flat_mask(nx,ny)==minimum_elevation && flowdir>0 && flowdir%2==0 && n%2==1) ){
      minimum_elevation=flat_mask(nx,ny);
      flowdir=n;
    }
  }

  return flowdir;
}

//d8_flow_flats
/**
  @brief  Calculates flow directions in flats
  @author Richard Barnes (rbarnes@umn.edu)

  This determines flow directions within flats which have been resolved
  using resolve_flats_barnes().

  Uses the helper function d8_masked_FlowDir()

  @param[in]  &flat_mask      A mask from resolve_flats_barnes()
  @param[in]  &groups         A grouping from resolve_flats_barnes()
  @param[out] &flowdirs       Returns flat-resolved flow directions

  @pre
    1. **flat_mask** contains the number of increments to be applied to each
       cell to form a gradient which will drain the flat it is a part of.
    2. Any cell without a local gradient has a value of #NO_FLOW in
       **flowdirs**; all other cells have defined flow directions.
    3. If a cell is part of a flat, it has a value greater than zero in
       **groups** indicating which flat it is a member of; otherwise, it has a
       value of 0.

  @post
    1. Every cell whose flow direction could be resolved by this algorithm
       (all drainable flats) will have a defined flow direction in
       **flowdirs**. Any cells which could not be resolved (non-drainable
       flats) will still be marked #NO_FLOW.
*/
template<class U>
void d8_flow_flats(
  const int_2d &flat_mask,
  const int_2d &groups,
  array2d<U> &flowdirs
){
  ProgressBar progress;

  diagnostic("%%Calculating D8 flow directions using flat mask...\n");
  progress.start( flat_mask.width()*flat_mask.height() );
  #pragma omp parallel for
  for(int x=1;x<flat_mask.width()-1;x++){
    progress.update( x*flat_mask.height() );
    for(int y=1;y<flat_mask.height()-1;y++)
      if(flat_mask(x,y)==flat_mask.no_data)
        continue;
      else if (flowdirs(x,y)==NO_FLOW)
        flowdirs(x,y)=d8_masked_FlowDir(flat_mask,groups,x,y);
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
}






//d8_upslope_area
/**
  @brief  Calculates the D8 up-slope area, given the D8 flow directions
  @author Richard Barnes (rbarnes@umn.edu)

  This calculates the D8 up-slope area of a grid of D8 flow directions using
  by calculating each cell's dependency on its neighbours and then using
  a priority-queue to process cells in a top-of-the-watershed-down fashion

  @param[in]  &flowdirs  A D8 flowdir grid from d8_flow_directions()
  @param[out] &area      Returns the up-slope area of each cell
*/
template<class T, class U>
void d8_upslope_area(const array2d<T> &flowdirs, array2d<U> &area){
  char_2d dependency;
  std::queue<grid_cell> sources;
  ProgressBar progress;

  diagnostic("\n###D8 Upslope Area\n");

  diagnostic_arg(
    "The sources queue will require at most approximately %ldMB of RAM.\n",
    flowdirs.width()*flowdirs.height()*((long)sizeof(grid_cell))/1024/1024
  );

  diagnostic("Resizing dependency matrix...");
  dependency.copyprops(flowdirs);
  diagnostic("succeeded.\n");

  diagnostic("Setting up the area matrix...");
  area.copyprops(flowdirs);
  area.init(0);
  area.no_data=d8_NO_DATA;
  diagnostic("succeeded.\n");

  diagnostic("%%Calculating dependency matrix & setting no_data cells...\n");
  progress.start( flowdirs.width()*flowdirs.height() );
  #pragma omp parallel for
  for(int x=0;x<flowdirs.width();x++){
    progress.update( x*flowdirs.height() );
    for(int y=0;y<flowdirs.height();y++){
      dependency(x,y)=0;
      if(flowdirs(x,y)==flowdirs.no_data){
        area(x,y)=area.no_data;
        continue;
      }
      for(int n=1;n<=8;n++)
        if(!flowdirs.in_grid(x+dx[n],y+dy[n]))
          continue;
        else if(flowdirs(x+dx[n],y+dy[n])==NO_FLOW)
          continue;
        else if(flowdirs(x+dx[n],y+dy[n])==flowdirs.no_data)
          continue;
        else if(n==inverse_flow[(int)flowdirs(x+dx[n],y+dy[n])])
          ++dependency(x,y);
    }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());

  diagnostic("%%Locating source cells...\n");
  progress.start( flowdirs.width()*flowdirs.height() );
  for(int x=0;x<flowdirs.width();x++){
    progress.update( x*flowdirs.height() );
    for(int y=0;y<flowdirs.height();y++)
      if(flowdirs(x,y)==flowdirs.no_data)
        continue;
      else if(dependency(x,y)==0)
        sources.push(grid_cell(x,y));
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());

  diagnostic("%%Calculating up-slope areas...\n");
  progress.start(flowdirs.data_cells);
  long int ccount=0;
  while(sources.size()>0){
    grid_cell c=sources.front();
    sources.pop();

    ccount++;
    progress.update(ccount);

    area(c.x,c.y)+=1;

    if(flowdirs(c.x,c.y)==NO_FLOW)
      continue;

    int nx=c.x+dx[(int)flowdirs(c.x,c.y)];
    int ny=c.y+dy[(int)flowdirs(c.x,c.y)];
    if(flowdirs.in_grid(nx,ny) && area(nx,ny)!=area.no_data){
      area(nx,ny)+=area(c.x,c.y);
      if((--dependency(nx,ny))==0)
        sources.push(grid_cell(nx,ny));
    }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
}




//d8_upslope_cells
/**
  @brief  Calculates which cells ultimately D8-flow through a given cell
  @author Richard Barnes (rbarnes@umn.edu)

  Given the coordinates x, y of a cell, this returns a grid indicating
  which cells ultimately flow into the indicated cell.
  1=Upslope cell
  2=Member of initializing line
  All other cells have a no_data value

  @param[in]  &flowdirs  A D8 flowdir grid from d8_flow_directions()
  @param[out] &area      Returns the up-slope area of each cell
*/
template<class T, class U>
void d8_upslope_cells(
  int x0, int y0, int x1, int y1,
  const array2d<T> &flowdirs,array2d<U> &upslope_cells
){
  diagnostic("Setting up the upslope_cells matrix...");
  upslope_cells.copyprops(flowdirs);
  upslope_cells.init(d8_NO_DATA);
  upslope_cells.no_data=d8_NO_DATA;
  diagnostic("succeeded.\n");
  ProgressBar progress;

  std::queue<grid_cell> expansion;

  if(x0>x1){
    std::swap(x0,x1);
    std::swap(y0,y1);
  }

  //Modified Bresenham Line-Drawing Algorithm
  int deltax=x1-x0;
  int deltay=y1-y0;
  float error=0;
  float deltaerr=(float)deltay/(float)deltax;
  if (deltaerr<0)
    deltaerr=-deltaerr;
  diagnostic_arg("Line slope is %f\n",deltaerr);
  int y=y0;
  for(int x=x0;x<=x1;x++){
    expansion.push(grid_cell(x,y));
    upslope_cells(x,y)=2;
    error+=deltaerr;
    if (error>=0.5) {
      expansion.push(grid_cell(x+1,y));
      upslope_cells(x+1,y)=2;
      y+=sgn(deltay);
      error-=1;
    }
  }

  progress.start(flowdirs.data_cells);
  long int ccount=0;
  while(expansion.size()>0){
    grid_cell c=expansion.front();
    expansion.pop();

    progress.update(ccount++);

    for(int n=1;n<=8;n++)
      if(!flowdirs.in_grid(c.x+dx[n],c.y+dy[n]))
        continue;
      else if(flowdirs(c.x+dx[n],c.y+dy[n])==NO_FLOW)
        continue;
      else if(flowdirs(c.x+dx[n],c.y+dy[n])==flowdirs.no_data)
        continue;
      else if(upslope_cells(c.x+dx[n],c.y+dy[n])==upslope_cells.no_data && n==inverse_flow[flowdirs(c.x+dx[n],c.y+dy[n])]){
        expansion.push(grid_cell(c.x+dx[n],c.y+dy[n]));
        upslope_cells(c.x+dx[n],c.y+dy[n])=1;
      }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
  diagnostic_arg("Found %ld up-slope cells.\n",ccount);
}




//d8_flats_alter_dem
/**
  @brief  Alters the elevations of the DEM so that all flats drain
  @author Richard Barnes (rbarnes@umn.edu)

  This alters elevations within the DEM so that flats which have been
  resolved using resolve_flats_barnes() will drain.

  @param[in]     &flat_mask   A mask from resolve_flats_barnes()
  @param[in]     &groups      A grouping from resolve_flats_barnes()
  @param[in,out] &elevations  2D array of elevations

  @pre
    1. **flat_mask** contains the number of increments to be applied to each
       cell to form a gradient which will drain the flat it is a part of.
    2. If a cell is part of a flat, it has a value greater than zero in
       **labels** indicating which flat it is a member of; otherwise, it has a
       value of 0.

  @post
    1. Every cell whose part of a flat which could be drained will have its
       elevation altered in such a way as to guarantee that its flat does
       drain.
*/
template<class U>
void d8_flats_alter_dem(
  const int_2d &flat_mask,
  const int_2d &labels,
  array2d<U> &elevations
){
  ProgressBar progress;

  diagnostic("%%Calculating D8 flow directions using flat mask...\n");
  progress.start( flat_mask.width()*flat_mask.height() );
  #pragma omp parallel for
  for(int x=1;x<flat_mask.width()-1;x++){
    progress.update( x*flat_mask.height() );
    for(int y=1;y<flat_mask.height()-1;y++){
      if(labels(x,y)==0)
        continue;

      bool higher[9];
      for(int n=1;n<=8;++n)
        higher[n]=elevations(x,y)>elevations(x+dx[n],y+dy[n]);
      //TODO: nextafterf is the floating point version; should use an
      //overloaded version instead to be able to handle both double and float
      for(int i=0;i<flat_mask(x,y);++i)
        elevations(x,y)=nextafterf(elevations(x,y),std::numeric_limits<U>::infinity());
      for(int n=1;n<=8;++n){
        int nx=x+dx[n];
        int ny=y+dy[n];      
        if(labels(nx,ny)==labels(x,y))
          continue;
        if(elevations(x,y)<elevations(nx,ny))
          continue;
        if(!higher[n])
          diagnostic_arg("Attempting to raise (%d,%d) resulted in an invalid alteration of the DEM!\n",x,y);
      }
    }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
}


#endif
