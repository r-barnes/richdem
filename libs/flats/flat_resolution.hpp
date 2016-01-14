/**
  @file
  Flat Resolution

  Richard Barnes (rbarnes@umn.edu), 2012

  Contains code to generate an elevation mask which is guaranteed to drain
  a flat using a convergent flow pattern (unless it's a mesa)
*/
#ifndef _flat_resolution_included
#define _flat_resolution_included

#include "utility.hpp"
#include "data_structures.hpp"
#include <deque>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#ifdef _OPENMP
  #include <omp.h>
#endif

//Procedure: BuildAwayGradient
/**
  @brief Build a gradient away from the high edges of the flats
  @author Richard Barnes (rbarnes@umn.edu)

  The queue of high-edge cells developed in find_flat_edges() is copied
  into the procedure. A breadth-first expansion labels cells by their
  distance away from terrain of higher elevation. The maximal distance
  encountered is noted.

  @param[in]  &elevations   A 2D array of cell elevations
  @param[in]  &flowdirs     A 2D array indicating each cell's flow direction
  @param[out] &flat_mask    A 2D array for storing flat_mask
  @param[in]  &edges        The high-edge FIFO queue from find_flat_edges()
  @param[out] &flat_height  Vector with length equal to max number of labels
  @param[in]  &labels       2D array storing labels developed in label_this()

  @pre
    1. Every cell in **labels** is marked either 0, indicating that the cell is
       not part of a flat, or a number greater than zero which identifies the
       flat to which the cell belongs.
    2. Any cell without a local gradient is marked #NO_FLOW in **flowdirs**.
    3. Every cell in **flat_mask** is initialized to 0.
    4. **edges** contains, in no particular order, all the high edge cells of
       the DEM (those flat cells adjacent to higher terrain) which are part of
       drainable flats.

  @post
    1. **flat_height** will have an entry for each label value of **labels**
       indicating the maximal number of increments to be applied to the flat
       identified by that label.
    2. **flat_mask** will contain the number of increments to be applied
       to each cell to form a gradient away from higher terrain; cells not in a
       flat will have a value of 0.
*/
template <class T>
static void BuildAwayGradient(
  const array2d<T> &flowdirs, int_2d &flat_mask, std::deque<grid_cell> edges,
  std::vector<int> &flat_height, const int_2d &labels
){
  Timer timer;
  timer.start();

  int x,y,nx,ny;
  int loops=1;
  grid_cell iteration_marker(-1,-1);


  diagnostic("Performing Barnes flat resolution's away gradient...");

  //Incrementation
  edges.push_back(iteration_marker);
  while(edges.size()!=1){  //Only iteration marker is left in the end
    x=edges.front().x;
    y=edges.front().y;
    edges.pop_front();

    if(x==-1){  //I'm an iteration marker
      loops++;
      edges.push_back(iteration_marker);
      continue;
    }

    if(flat_mask(x,y)>0) continue;  //I've already been incremented!

    //If I incremented, maybe my neighbours should too
    flat_mask(x,y)=loops;
    flat_height[labels(x,y)]=loops;
    for(int n=1;n<=8;n++){
      nx=x+dx[n];
      ny=y+dy[n];
      if(labels.in_grid(nx,ny)
          && labels(nx,ny)==labels(x,y)
          && flowdirs(nx,ny)==NO_FLOW)
        edges.push_back(grid_cell(nx,ny));
    }
  }

  timer.stop();
  diagnostic_arg("succeeded in %fs!\n",timer.accumulated());
}



//Procedure: BuildTowardsCombinedGradient
/**
  @brief Builds gradient away from the low edges of flats, combines gradients
  @author Richard Barnes (rbarnes@umn.edu)

  The queue of low-edge cells developed in find_flat_edges() is copied
  into the procedure. A breadth-first expansion labels cells by their
  distance away from terrain of lower elevation. This is combined with
  the gradient from BuildAwayGradient() to give the final increments of
  each cell in forming the flat mask.

  @param[in]  &elevations   A 2D array of cell elevations
  @param[in]  &flowdirs     A 2D array indicating each cell's flow direction
  @param[in,out] &flat_mask A 2D array for storing flat_mask
  @param[in]  &edges        The low-edge FIFO queue from find_flat_edges()
  @param[in]  &flat_height  Vector with length equal to max number of labels
  @param[in]  &labels       2D array storing labels developed in label_this()

  @pre
    1. Every cell in **labels** is marked either 0, indicating that the cell
       is not part of a flat, or a number greater than zero which identifies
       the flat to which the cell belongs.
    2. Any cell without a local gradient is marked #NO_FLOW in **flowdirs**.
    3. Every cell in **flat_mask** has either a value of 0, indicating
       that the cell is not part of a flat, or a value greater than zero
       indicating the number of increments which must be added to it to form a
       gradient away from higher terrain.
    4. **flat_height** has an entry for each label value of **labels**
       indicating the maximal number of increments to be applied to the flat
       identified by that label in order to form the gradient away from higher
       terrain.
    5. **edges** contains, in no particular order, all the low edge cells of
       the DEM (those flat cells adjacent to lower terrain).

  @post
    1. **flat_mask** will contain the number of increments to be applied
       to each cell to form a superposition of the gradient away from higher
       terrain with the gradient towards lower terrain; cells not in a flat
       have a value of 0.
*/
template <class T>
static void BuildTowardsCombinedGradient(
  const array2d<T> &flowdirs, int_2d &flat_mask, std::deque<grid_cell> edges,
  std::vector<int> &flat_height, const int_2d &labels
){
  Timer timer;
  timer.start();

  int x,y,nx,ny;
  int loops=1;
  grid_cell iteration_marker(-1,-1);


  diagnostic("Barnes flat resolution: toward and combined gradients...");

  //Make previous flat_mask negative so that we can keep track of where we are
  #pragma omp parallel for collapse(2)
  for(int x=0;x<flat_mask.width();x++)
  for(int y=0;y<flat_mask.height();y++)
    flat_mask(x,y)*=-1;


  //Incrementation
  edges.push_back(iteration_marker);
  while(edges.size()!=1){  //Only iteration marker is left in the end
    x=edges.front().x;
    y=edges.front().y;
    edges.pop_front();

    if(x==-1){  //I'm an iteration marker
      loops++;
      edges.push_back(iteration_marker);
      continue;
    }

    if(flat_mask(x,y)>0) continue;  //I've already been incremented!

    //If I incremented, maybe my neighbours should too
    if(flat_mask(x,y)!=0)  //If !=0, it _will_ be less than 0.
      flat_mask(x,y)=(flat_height[labels(x,y)]+flat_mask(x,y))+2*loops;
    else
      flat_mask(x,y)=2*loops;

    for(int n=1;n<=8;n++){
      nx=x+dx[n];
      ny=y+dy[n];
      if(labels.in_grid(nx,ny)
          && labels(nx,ny)==labels(x,y)
          && flowdirs(nx,ny)==NO_FLOW)
        edges.push_back(grid_cell(nx,ny));
    }
  }

  timer.stop();
  diagnostic_arg("succeeded in %fs!\n",timer.accumulated());
}


//Procedure: label_this
/**
  @brief Labels all the cells of a flat with a common label.
  @author Richard Barnes (rbarnes@umn.edu)

  Performs a flood fill operation which labels all the cells of a flat
  with a common label. Each flat will have a unique label

  @param[in]  x           x-coordinate of flood fill seed
  @param[in]  y           y-coordinate of flood-fill seed
  @param[in]  label       Label to apply to the cells
  @param[out] &labels     2D array which will contain the labels
  @param[in]  &elevations 2D array of cell elevations

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM.
    2. **labels** has the same dimensions as **elevations**.
    3. **(x0,y0)** belongs to the flat which is to be labeled.
    4. **label** is a unique label which has not been previously applied to a
       flat.
    5. **labels** is initialized to zero prior to the first call to this
       function.

  @post
    1. **(x0,y0)** and every cell reachable from it by passing over only cells
       of the same elevation as it (all the cells in the flat to which c
       belongs) will be marked as **label** in **labels**.
*/
template<class T>
static void label_this(
  int x0, int y0, const int label, int_2d &labels,
  const array2d<T> &elevations
){
  std::queue<grid_cell> to_fill;
  to_fill.push(grid_cell(x0,y0));
  const T target_elevation=elevations(x0,y0);

  while(to_fill.size()>0){
    grid_cell c=to_fill.front();
    to_fill.pop();
    if(elevations(c.x,c.y)!=target_elevation)
      continue;
    if(labels(c.x,c.y)>0)
      continue;
    labels(c.x,c.y)=label;
    for(int n=1;n<=8;n++)
      if(labels.in_grid(c.x+dx[n],c.y+dy[n]))
        to_fill.push(grid_cell(c.x+dx[n],c.y+dy[n]));
  }
}

//Procedure: find_flat_edges
/**
  @brief Identifies cells adjacent to higher and lower terrain
  @author Richard Barnes (rbarnes@umn.edu)

  Cells adjacent to lower and higher terrain are identified and
  added to the appropriate queue

  @param[out] &low_edges  Queue for storing cells adjacent to lower terrain
  @param[out] &high_edges Queue for storing cells adjacent to higher terrain
  @param[in]  &flowdirs   2D array indicating flow direction for each cell
  @param[in]  &elevations 2D array of cell elevations

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM.
    2. Any cell without a local gradient is marked #NO_FLOW in **flowdirs**.

  @post
    1. **high_edges** will contain, in no particular order, all the high edge
       cells of the DEM: those flat cells adjacent to higher terrain.
    2. **low_edges** will contain, in no particular order, all the low edge
       cells of the DEM: those flat cells adjacent to lower terrain.
*/
template <class T, class U>
static void find_flat_edges(
  std::deque<grid_cell> &low_edges,
  std::deque<grid_cell> &high_edges,
  const array2d<U> &flowdirs,
  const array2d<T> &elevations
){
  int nx,ny;
  int cells_without_flow=0;
  ProgressBar progress;
  diagnostic("%%Searching for flats...\n");
  progress.start( flowdirs.width()*flowdirs.height() );
  for(int x=0;x<flowdirs.width();x++){
    progress.update( x*flowdirs.height() );
    for(int y=0;y<flowdirs.height();y++){
      if(flowdirs(x,y)==flowdirs.no_data)
        continue;
      if(flowdirs(x,y)==NO_FLOW)
        cells_without_flow++;
      for(int n=1;n<=8;n++){
        nx=x+dx[n];
        ny=y+dy[n];

        if(!flowdirs.in_grid(nx,ny)) continue;
        if(flowdirs(nx,ny)==flowdirs.no_data) continue;

        if(flowdirs(x,y)!=NO_FLOW && flowdirs(nx,ny)==NO_FLOW && elevations(nx,ny)==elevations(x,y)){
          low_edges.push_back(grid_cell(x,y));
          break;
        } else if(flowdirs(x,y)==NO_FLOW && elevations(x,y)<elevations(nx,ny)){
          high_edges.push_back(grid_cell(x,y));
          break;
        }
      }
    }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
  diagnostic_arg("%d cells had no flow direction.\n",cells_without_flow);
}


//Procedure: resolve_flats_barnes
/**
  @brief  Performs the flat resolution by Barnes, Lehman, and Mulla.
  @author Richard Barnes (rbarnes@umn.edu)

  TODO

  @param[in]  &elevations 2D array of cell elevations
  @param[in]  &flowdirs   2D array indicating flow direction of each cell
  @param[in]  &flat_mask  2D array which will hold incremental elevation mask
  @param[in]  &labels     2D array indicating flat membership

  @pre
    1. **elevations** contains the elevations of every cell or the _NoData_
        value for cells not part of the DEM.
    2. Any cell without a local gradient is marked #NO_FLOW in **flowdirs**.

  @post
    1. **flat_mask** will have a value greater than or equal to zero for every
       cell, indicating its number of increments. These can be used be used
       in conjunction with **labels** to determine flow directions without
       altering the DEM, or to alter the DEM in subtle ways to direct flow.
    2. **labels** will have values greater than or equal to 1 for every cell
       which is in a flat. Each flat's cells will bear a label unique to that
       flat.
*/
template <class T, class U>
void resolve_flats_barnes(
  const array2d<T> &elevations,
  const array2d<U> &flowdirs,
  int_2d &flat_mask,
  int_2d &labels
){
  Timer timer;
  timer.start();

  std::deque<grid_cell> low_edges,high_edges;  //TODO: Need estimate of size

  diagnostic("\n###Barnes Flat Resolution\n");

  diagnostic("Setting up labels matrix...");
  labels.copyprops(flowdirs);
  labels.init(0);
  diagnostic("succeeded.\n");

  diagnostic("Setting up flat resolution mask...");
  flat_mask.copyprops(elevations);
  flat_mask.init(0);
  flat_mask.no_data=-1;
  diagnostic("succeeded!\n");

  find_flat_edges(low_edges, high_edges, flowdirs, elevations);

  if(low_edges.size()==0){
    if(high_edges.size()>0)
      diagnostic("There were flats, but none of them had outlets!\n");
    else
      diagnostic("There were no flats!\n");
    return;
  }

  diagnostic("Labeling flats...");
  int group_number=1;
  for(std::deque<grid_cell>::iterator i=low_edges.begin();i!=low_edges.end();++i)
    if(labels(i->x,i->y)==0)
      label_this(i->x, i->y, group_number++, labels, elevations);
  diagnostic("succeeded!\n");

  diagnostic_arg("Found %d unique flats.\n",group_number);

  diagnostic("Removing flats without outlets from the queue...");
  std::deque<grid_cell> temp;
  for(std::deque<grid_cell>::iterator i=high_edges.begin();i!=high_edges.end();++i)
    if(labels(i->x,i->y)!=0)
      temp.push_back(*i);
  diagnostic("succeeded.\n");

  if(temp.size()<high_edges.size())  //TODO: Prompt for intervention?
    diagnostic("\033[91mNot all flats have outlets; the DEM contains sinks/pits/depressions!\033[39m\n");
  high_edges=temp;
  temp.clear();

  diagnostic_arg("The flat height vector will require approximately %ldMB of RAM.\n",
        group_number*((long)sizeof(int))/1024/1024);
  diagnostic("Creating flat height vector...");
  std::vector<int> flat_height(group_number);
  diagnostic("succeeded!\n");

  BuildAwayGradient(
    flowdirs, flat_mask, high_edges, flat_height, labels
  );
  BuildTowardsCombinedGradient(
    flowdirs, flat_mask, low_edges, flat_height, labels
  );

  timer.stop();
  diagnostic_arg("Calculation time for Barnes algorithm: %fs\n", timer.accumulated());
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
  @param[in]  &labels         The labels output from resolve_flats_barnes()
  @param[in]  x               x coordinate of cell
  @param[in]  y               y coordinate of cell

  @returns    The flow direction of the cell
*/
static int d8_masked_FlowDir(
  const int_2d &flat_mask,
  const int_2d &labels,
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
    if( labels(nx,ny)!=labels(x,y))
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
  @param[in]  &labels         The labels output from resolve_flats_barnes()
  @param[out] &flowdirs       Returns flat-resolved flow directions

  @pre
    1. **flat_mask** contains the number of increments to be applied to each
       cell to form a gradient which will drain the flat it is a part of.
    2. Any cell without a local gradient has a value of #NO_FLOW in
       **flowdirs**; all other cells have defined flow directions.
    3. If a cell is part of a flat, it has a value greater than zero in
       **labels** indicating which flat it is a member of; otherwise, it has a
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
  const int_2d &labels,
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
        flowdirs(x,y)=d8_masked_FlowDir(flat_mask,labels,x,y);
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
}


//d8_flats_alter_dem
/**
  @brief  Alters the elevations of the DEM so that all flats drain
  @author Richard Barnes (rbarnes@umn.edu)

  This alters elevations within the DEM so that flats which have been
  resolved using resolve_flats_barnes() will drain.

  @param[in]     &flat_mask   A mask from resolve_flats_barnes()
  @param[in]     &labels      A grouping from resolve_flats_barnes()
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
