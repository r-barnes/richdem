/**
  @file
  @brief Resolve flats according to Barnes (2014)
  @author Richard Barnes (rbarnes@umn.edu), 2012

  Contains code to generate an elevation mask which is guaranteed to drain
  a flat using a convergent flow pattern (unless it's a mesa)
*/
#ifndef _richdem_flat_resolution_hpp_
#define _richdem_flat_resolution_hpp_

#include "richdem/common/logger.hpp"
#include "richdem/common/ProgressBar.hpp"
#include "richdem/common/grid_cell.hpp"
#include "richdem/flowmet/d8_flowdirs.hpp"
#include <deque>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>

namespace richdem {

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
  const Array2D<int32_t> &flat_mask,
  const Array2D<int32_t> &labels,
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
    2. Any cell without a local gradient has a value of NO_FLOW in
       **flowdirs**; all other cells have defined flow directions.
    3. If a cell is part of a flat, it has a value greater than zero in
       **labels** indicating which flat it is a member of; otherwise, it has a
       value of 0.

  @post
    1. Every cell whose flow direction could be resolved by this algorithm
       (all drainable flats) will have a defined flow direction in
       **flowdirs**. Any cells which could not be resolved (non-drainable
       flats) will still be marked NO_FLOW.
*/
template<class U>
void d8_flow_flats(
  const Array2D<int32_t> &flat_mask,
  const Array2D<int32_t> &labels,
  Array2D<U> &flowdirs
){
  ProgressBar progress;

  RDLOG_ALG_NAME<<"Calculating D8 flow directions using flat mask...";
  progress.start( flat_mask.width()*flat_mask.height() );
  #pragma omp parallel for
  for(int x=1;x<flat_mask.width()-1;x++){
    progress.update( x*flat_mask.height() );
    for(int y=1;y<flat_mask.height()-1;y++)
      if(flat_mask(x,y)==flat_mask.noData())
        continue;
      else if (flowdirs(x,y)==NO_FLOW)
        flowdirs(x,y)=d8_masked_FlowDir(flat_mask,labels,x,y);
  }
  RDLOG_TIME_USE<<"Succeeded in = "<<progress.stop()<<" s";
}

//Procedure: BuildAwayGradient
/**
  @brief Build a gradient away from the high edges of the flats
  @author Richard Barnes (rbarnes@umn.edu)

  The queue of high-edge cells developed in find_flat_edges() is copied
  into the procedure. A breadth-first expansion labels cells by their
  distance away from terrain of higher elevation. The maximal distance
  encountered is noted.

  @param[in]  &flowdirs     A 2D array indicating each cell's flow direction
  @param[out] &flat_mask    A 2D array for storing flat_mask
  @param[in]  &edges        The high-edge FIFO queue from find_flat_edges()
  @param[out] &flat_height  Vector with length equal to max number of labels
  @param[in]  &labels       2D array storing labels developed in label_this()

  @pre
    1. Every cell in **labels** is marked either 0, indicating that the cell is
       not part of a flat, or a number greater than zero which identifies the
       flat to which the cell belongs.
    2. Any cell without a local gradient is marked NO_FLOW in **flowdirs**.
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
template<class U>
static void BuildAwayGradient(
  const Array2D<U>       &flowdirs,
  Array2D<int32_t>       &flat_mask,
  std::deque<GridCell>  edges,
  std::vector<int>       &flat_height,
  const Array2D<int32_t> &labels
){
  Timer timer;
  timer.start();

  int loops = 1;
  GridCell iteration_marker(-1,-1);

  RDLOG_PROGRESS<<"Performing Barnes flat resolution's away gradient...";

  //Incrementation
  edges.push_back(iteration_marker);
  while(edges.size()!=1){  //Only iteration marker is left in the end
    int x = edges.front().x;
    int y = edges.front().y;
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
      int nx = x+dx[n];
      int ny = y+dy[n];
      if(labels.inGrid(nx,ny)
          && labels(nx,ny)==labels(x,y)
          && flowdirs(nx,ny)==NO_FLOW)
        edges.push_back(GridCell(nx,ny));
    }
  }

  timer.stop();
  RDLOG_TIME_USE<<"Succeeded in = "<<timer.accumulated()<<" s";
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

  @param[in]  &flowdirs     A 2D array indicating each cell's flow direction
  @param[in,out] &flat_mask A 2D array for storing flat_mask
  @param[in]  &edges        The low-edge FIFO queue from find_flat_edges()
  @param[in]  &flat_height  Vector with length equal to max number of labels
  @param[in]  &labels       2D array storing labels developed in label_this()

  @pre
    1. Every cell in **labels** is marked either 0, indicating that the cell
       is not part of a flat, or a number greater than zero which identifies
       the flat to which the cell belongs.
    2. Any cell without a local gradient is marked NO_FLOW in **flowdirs**.
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
template<class U>
static void BuildTowardsCombinedGradient(
  const Array2D<U>       &flowdirs,
  Array2D<int32_t>       &flat_mask,
  std::deque<GridCell>  edges,
  std::vector<int>       &flat_height,
  const Array2D<int32_t> &labels
){
  Timer timer;
  timer.start();

  int loops = 1;
  GridCell iteration_marker(-1,-1);


  RDLOG_PROGRESS<<"Barnes flat resolution: toward and combined gradients...";

  //Make previous flat_mask negative so that we can keep track of where we are
  #pragma omp parallel for collapse(2)
  for(int x=0;x<flat_mask.width();x++)
  for(int y=0;y<flat_mask.height();y++)
    flat_mask(x,y)*=-1;


  //Incrementation
  edges.push_back(iteration_marker);
  while(edges.size()!=1){  //Only iteration marker is left in the end
    int x = edges.front().x;
    int y = edges.front().y;
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
      int nx = x+dx[n];
      int ny = y+dy[n];
      if(labels.inGrid(nx,ny)
          && labels(nx,ny)==labels(x,y)
          && flowdirs(nx,ny)==NO_FLOW)
        edges.push_back(GridCell(nx,ny));
    }
  }

  timer.stop();
  RDLOG_TIME_USE<<"Succeeded in = "<<timer.accumulated()<<" s";
}


//Procedure: label_this
/**
  @brief Labels all the cells of a flat with a common label.
  @author Richard Barnes (rbarnes@umn.edu)

  Performs a flood fill operation which labels all the cells of a flat
  with a common label. Each flat will have a unique label

  @param[in]  x0
            x-coordinate of flood fill seed
  @param[in]  y0          y-coordinate of flood-fill seed
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
  int x0,
  int y0,
  const int label,
  Array2D<int32_t> &labels,
  const Array2D<T> &elevations
){
  std::queue<GridCell> to_fill;
  to_fill.push(GridCell(x0,y0));
  const T target_elevation = elevations(x0,y0);

  while(to_fill.size()>0){
    GridCell c = to_fill.front();
    to_fill.pop();
    if(elevations(c.x,c.y)!=target_elevation)
      continue;
    if(labels(c.x,c.y)>0)
      continue;
    labels(c.x,c.y)=label;
    for(int n=1;n<=8;n++)
      if(labels.inGrid(c.x+dx[n],c.y+dy[n]))
        to_fill.push(GridCell(c.x+dx[n],c.y+dy[n]));
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
    2. Any cell without a local gradient is marked NO_FLOW in **flowdirs**.

  @post
    1. **high_edges** will contain, in no particular order, all the high edge
       cells of the DEM: those flat cells adjacent to higher terrain.
    2. **low_edges** will contain, in no particular order, all the low edge
       cells of the DEM: those flat cells adjacent to lower terrain.
*/
template <class T, class U>
static void find_flat_edges(
  std::deque<GridCell> &low_edges,
  std::deque<GridCell> &high_edges,
  const Array2D<U>      &flowdirs,
  const Array2D<T>      &elevations
){
  int cells_without_flow=0;
  ProgressBar progress;
  RDLOG_PROGRESS<<"Searching for flats...";
  progress.start( flowdirs.width()*flowdirs.height() );
  for(int x=0;x<flowdirs.width();x++){
    progress.update( x*flowdirs.height() );
    for(int y=0;y<flowdirs.height();y++){
      if(flowdirs(x,y)==flowdirs.noData())
        continue;
      if(flowdirs(x,y)==NO_FLOW)
        cells_without_flow++;
      for(int n=1;n<=8;n++){
        int nx = x+dx[n];
        int ny = y+dy[n];

        if(!flowdirs.inGrid(nx,ny)) continue;
        if(flowdirs(nx,ny)==flowdirs.noData()) continue;

        if(flowdirs(x,y)!=NO_FLOW && flowdirs(nx,ny)==NO_FLOW && elevations(nx,ny)==elevations(x,y)){
          low_edges.push_back(GridCell(x,y));
          break;
        } else if(flowdirs(x,y)==NO_FLOW && elevations(x,y)<elevations(nx,ny)){
          high_edges.push_back(GridCell(x,y));
          break;
        }
      }
    }
  }
  RDLOG_TIME_USE<<"Succeeded in = "<<progress.stop()<<" s";
  RDLOG_MISC<<"Cells with no flow direction = "<<cells_without_flow;
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
    2. Any cell without a local gradient is marked NO_FLOW in **flowdirs**.

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
  const Array2D<T> &elevations,
  const Array2D<U> &flowdirs,
  Array2D<int32_t> &flat_mask,
  Array2D<int32_t> &labels
){
  Timer timer;
  timer.start();

  std::deque<GridCell> low_edges,high_edges;  //TODO: Need estimate of size

  RDLOG_ALG_NAME<<"Flat Resolution (Barnes 2014)";
  RDLOG_CITATION<<"Barnes, R., Lehman, C., Mulla, D., 2014a. An efficient assignment of drainage direction over flat surfaces in raster digital elevation models. Computers & Geosciences 62, 128–135. doi:10.1016/j.cageo.2013.01.009";

  RDLOG_PROGRESS<<"Setting up labels matrix...";
  labels.templateCopy(elevations);
  labels.resize(flowdirs);
  labels.setAll(0);

  RDLOG_PROGRESS<<"Setting up flat resolution mask...";
  flat_mask.templateCopy(elevations);
  flat_mask.resize(elevations);
  flat_mask.setAll(0);
  flat_mask.setNoData(-1);

  find_flat_edges(low_edges, high_edges, flowdirs, elevations);

  if(low_edges.size()==0){
    if(high_edges.size()>0)
      RDLOG_WARN<<"There were flats, but none of them had outlets!";
    else
      RDLOG_WARN<<"There were no flats!";
    return;
  }

  RDLOG_PROGRESS<<"Labeling flats...";
  int group_number=1;
  for(auto i=low_edges.begin();i!=low_edges.end();++i)
    if(labels(i->x,i->y)==0)
      label_this(i->x, i->y, group_number++, labels, elevations);

  RDLOG_MISC<<"Unique flats = "<<group_number;

  RDLOG_PROGRESS<<"Removing flats without outlets from the queue...";
  std::deque<GridCell> temp;
  for(auto i=high_edges.begin();i!=high_edges.end();++i)
    if(labels(i->x,i->y)!=0)
      temp.push_back(*i);

  if(temp.size()<high_edges.size())  //TODO: Prompt for intervention?
    RDLOG_WARN<<"Not all flats have outlets; the DEM contains sinks/pits/depressions!";
  high_edges=temp;
  temp.clear();

  RDLOG_MEM_USE<<"The flat height vector will require approximately "
               <<(group_number*((long)sizeof(int))/1024/1024)
               <<"MB of RAM.";
        
  RDLOG_PROGRESS<<"Creating flat height vector...";
  std::vector<int> flat_height(group_number);

  BuildAwayGradient(
    flowdirs, flat_mask, high_edges, flat_height, labels
  );
  BuildTowardsCombinedGradient(
    flowdirs, flat_mask, low_edges, flat_height, labels
  );

  RDLOG_TIME_USE<<"Wall-time = "<<timer.stop()<<" s";
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
  const Array2D<int32_t> &flat_mask,
  const Array2D<int32_t> &labels,
  Array2D<U> &elevations
){
  ProgressBar progress;

  RDLOG_ALG_NAME<<"Calculating D8 flow directions using flat mask...";
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
          RDLOG_WARN<<"Attempting to raise ("<<x<<","<<y<<") resulted in an invalid alteration of the DEM!";
      }
    }
  }
  RDLOG_TIME_USE<<"Succeeded in = "<<progress.stop()<<" s";
}



//TODO: Documentation
template<class T, class U>
void barnes_flat_resolution_d8(Array2D<T> &elevations, Array2D<U> &flowdirs, bool alter){
  d8_flow_directions(elevations,flowdirs);

  Array2D<int32_t> flat_mask, labels;

  resolve_flats_barnes(elevations,flowdirs,flat_mask,labels);

  if(alter){  
    //NOTE: If this value appears anywhere an error's occurred
    flowdirs.setAll(155); //TODO
    d8_flats_alter_dem(flat_mask, labels, elevations);
    d8_flow_directions(elevations,flowdirs);
  } else {
    d8_flow_flats(flat_mask,labels,flowdirs);
  }

  flowdirs.templateCopy(elevations);
}

}

#endif
