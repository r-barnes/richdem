/**
  @file
  @brief Resolve flats according to Barnes (2014)
  @author Richard Barnes (rbarnes@umn.edu), 2012

  Contains code to generate an elevation mask which is guaranteed to drain
  a flat using a convergent flow pattern (unless it's a mesa)
*/
#ifndef _richdem_Barnes2014_flat_resolution_hpp_
#define _richdem_Barnes2014_flat_resolution_hpp_

#include <richdem/common/logger.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <richdem/common/grid_cell.hpp>
#include <richdem/common/Array2D.hpp>

#include <richdem/flats/find_flats.hpp>

#include <deque>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>

namespace richdem {



/**
  @brief Build a gradient away from the high edges of the flats
  @author Richard Barnes (rbarnes@umn.edu)

  The queue of high-edge cells developed in FindFlatEdges() is copied
  into the procedure. A breadth-first expansion labels cells by their
  distance away from terrain of higher elevation. The maximal distance
  encountered is noted.

  @param[in]  &flats        2D array indicating flat membership from FindFlats()
  @param[out] &flat_mask    A 2D array for storing flat_mask
  @param[in]  &edges        The high-edge FIFO queue from FindFlatEdges()
  @param[out] &flat_height  Vector with length equal to max number of labels
  @param[in]  &labels       2D array storing labels developed in LabelFlat()

  @pre
    1. Every cell in **labels** is marked either 0, indicating that the cell is
       not part of a flat, or a number greater than zero which identifies the
       flat to which the cell belongs.
    2. Any cell without a local gradient is marked IS_A_FLAT in **flats**.
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
static void BuildAwayGradient(
  const Array2D<int8_t>    &flats,
  Array2D<int32_t>         &flat_mask,
  std::deque<GridCell>     &high_edges,
  std::vector<int>         &flat_height,
  const Array2D<int32_t>   &labels
){
  Timer timer;
  timer.start();

  int loops = 1;
  GridCell iteration_marker(-1,-1);

  RDLOG_PROGRESS<<"Performing Barnes flat resolution's away gradient...";

  //Incrementation
  high_edges.push_back(iteration_marker);
  while(high_edges.size()!=1){  //Only iteration marker is left in the end
    const auto c = high_edges.front();
    high_edges.pop_front();

    if(c.x==-1){  //I'm an iteration marker
      loops++;
      high_edges.push_back(iteration_marker);
      continue;
    }

    if(flat_mask(c.x,c.y)>0)
      continue;  //I've already been incremented!

    //If I incremented, maybe my neighbours should too
    flat_mask(c.x,c.y)           = loops;
    flat_height[labels(c.x,c.y)] = loops;
    for(int n=1;n<=8;n++){
      const int nx = c.x+dx[n];
      const int ny = c.y+dy[n];
      if(
           labels.inGrid(nx,ny)
        && labels(nx,ny)==labels(c.x,c.y)
        && flats(nx,ny)==IS_A_FLAT
      ){
        high_edges.emplace_back(nx,ny);
      }
    }
  }

  timer.stop();
  RDLOG_TIME_USE<<"Succeeded in = "<<timer.accumulated()<<" s";
}



/**
  @brief Builds gradient away from the low edges of flats, combines gradients
  @author Richard Barnes (rbarnes@umn.edu)

  The queue of low-edge cells developed in FindFlatEdges() is copied
  into the procedure. A breadth-first expansion labels cells by their
  distance away from terrain of lower elevation. This is combined with
  the gradient from BuildAwayGradient() to give the final increments of
  each cell in forming the flat mask.

  @param[in]  &flats        2D array indicating flat membership from FindFlats()
  @param[in,out] &flat_mask A 2D array for storing flat_mask
  @param[in]  &edges        The low-edge FIFO queue from FindFlatEdges()
  @param[in]  &flat_height  Vector with length equal to max number of labels
  @param[in]  &labels       2D array storing labels developed in LabelFlat()

  @pre
    1. Every cell in **labels** is marked either 0, indicating that the cell
       is not part of a flat, or a number greater than zero which identifies
       the flat to which the cell belongs.
    2. Any cell without a local gradient is marked IS_A_FLAT in **flats**.
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
static void BuildTowardsCombinedGradient(
  Array2D<int8_t>        &flats,
  Array2D<int32_t>       &flat_mask,
  std::deque<GridCell>   &low_edges,
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
    flat_mask(x,y) *= -1;


  //Incrementation
  low_edges.push_back(iteration_marker);
  while(low_edges.size()!=1){  //Only iteration marker is left in the end
    const auto c = low_edges.front();
    low_edges.pop_front();

    if(c.x==-1){  //I'm an iteration marker
      loops++;
      low_edges.push_back(iteration_marker);
      continue;
    }

    if(flat_mask(c.x,c.y)>0)
      continue;  //I've already been incremented!

    //If I incremented, maybe my neighbours should too
    if(flat_mask(c.x,c.y)!=0)  //If !=0, it _will_ be less than 0.
      flat_mask(c.x,c.y)=(flat_height[labels(c.x,c.y)]+flat_mask(c.x,c.y))+2*loops;
    else
      flat_mask(c.x,c.y)=2*loops;

    for(int n=1;n<=8;n++){
      const int nx = c.x+dx[n];
      const int ny = c.y+dy[n];
      if(
           labels.inGrid(nx,ny)
        && labels(nx,ny)==labels(c.x,c.y)
        && flats (nx,ny)==IS_A_FLAT
      ){
        low_edges.emplace_back(nx,ny);
      }
    }
  }

  timer.stop();
  RDLOG_TIME_USE<<"Succeeded in = "<<timer.accumulated()<<" s";
}



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
static void LabelFlat(
  const int x0,                 //Cell at the edge of the flat
  const int y0,                 //Cell at the edge of the flat
  const int label,              //Label to give the flat
  Array2D<int32_t> &labels,     //Labels array to which label is applied
  const Array2D<T> &elevations  //Elevations array, to determine if we are still in the flat
){
  //We'll perform a breadth-first traversal on the flat, so we need a queue
  std::queue<GridCell> to_fill;
  to_fill.emplace(x0,y0);
  const T target_elevation = elevations(x0,y0);

  while(to_fill.size()>0){
    const auto c = to_fill.front();
    to_fill.pop();

    //Higher and lower cells are not part of the flat
    if(elevations(c.x,c.y)!=target_elevation) 
      continue;

    //Already labeled cells are in the same flat, but have been dealt with
    if(labels(c.x,c.y)>0)                     
      continue;

    //Okay, this cell is part of the flat. Label it.
    labels(c.x,c.y) = label;

    //Consider this cell's neighbours
    for(int n=1;n<=8;n++){
      const int nx = c.x+dx[n];
      const int ny = c.y+dy[n];
      if(labels.inGrid(nx,ny)) //TODO: Should probably avoid adding cells that cannot be part of the flat. This'll likely speed things up.
        to_fill.emplace(nx,ny);
    }
  }
}



/**
  @brief Identifies cells adjacent to higher and lower terrain
  @author Richard Barnes (rbarnes@umn.edu)

  Cells adjacent to lower and higher terrain are identified and
  added to the appropriate queue

  @param[out] &low_edges  Queue for storing cells adjacent to lower terrain
  @param[out] &high_edges Queue for storing cells adjacent to higher terrain
  @param[in]  &flats      2D array indicating flat membership from FindFlats()
  @param[in]  &elevations 2D array of cell elevations

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM.
    2. Any cell without a local gradient is marked IS_A_FLAT in **flats**.

  @post
    1. **high_edges** will contain, in no particular order, all the high edge
       cells of the DEM: those flat cells adjacent to higher terrain.
    2. **low_edges** will contain, in no particular order, all the low edge
       cells of the DEM: those flat cells adjacent to lower terrain.
*/
template <class T>
static void FindFlatEdges(
  std::deque<GridCell>  &low_edges,
  std::deque<GridCell>  &high_edges,
  const Array2D<int8_t> &flats,
  const Array2D<T>      &elevations
){
  int cells_without_flow=0;
  ProgressBar progress;
  RDLOG_PROGRESS<<"Searching for flats...";

  progress.start( flats.size() );

  #pragma omp parallel for collapse(2)
  for(int y=0;y<flats.height();y++)
  for(int x=0;x<flats.width();x++){
    ++progress;

    if(flats(x,y)==IS_A_FLAT)
      cells_without_flow++;

    if(flats.isNoData(x,y))
      continue;

    for(int n=1;n<=8;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];

      if(!flats.inGrid(nx,ny))
        continue;

      //If the focal cell has flow, but has a neighbour without flow at the same
      //elevation, then the focal cell is at the edge of a flat containing that
      //neighbour
      if(
           flats( x, y)==NOT_A_FLAT
        && flats(nx,ny)==IS_A_FLAT
        && elevations(nx,ny)==elevations(x,y)
      ){
        #pragma omp critical
        low_edges.emplace_back(x,y);
        break;

      //If the focal cell has no flow and has a neighbour which is at a higher
      //elevation, then the focal cell is at the edge of a flat
      } else if(
           flats(x,y)==IS_A_FLAT     
        && elevations(x,y)<elevations(nx,ny)
      ){
        #pragma omp critical
        high_edges.emplace_back(x,y);
        break;
      }
    }
  }

  RDLOG_TIME_USE<<"Succeeded in = "<<progress.stop()<<" s";
  RDLOG_MISC<<"Cells with no flow direction = "<<cells_without_flow;
  RDLOG_MISC<<"Low edge cells               = "<<low_edges.size();
  RDLOG_MISC<<"High edge cells              = "<<high_edges.size();
}



/**
  @brief  Generates a flat resolution mask in the style of Barnes (2014)
  @author Richard Barnes (rbarnes@umn.edu)

  This algorithm is a modification of that presented by Barnes (2014). It has
  been rejiggered so that a knowledge of flow directions is not necessary to run
  it.

  @param[in]  &elevations 2D array of cell elevations
  @param[in]  &flat_mask  2D array which will hold incremental elevation mask
  @param[in]  &labels     2D array indicating flat membership

  @pre
    1. **elevations** contains the elevations of every cell or the _NoData_
        value for cells not part of the DEM.

  @post
    1. **flat_mask** will have a value greater than or equal to zero for every
       cell, indicating its number of increments. These can be used be used
       in conjunction with **labels** to determine flow directions without
       altering the DEM, or to alter the DEM in subtle ways to direct flow.
    2. **labels** will have values greater than or equal to 1 for every cell
       which is in a flat. Each flat's cells will bear a label unique to that
       flat. A value of 0 means the cell is not part of a flat.
*/
template <class T>
void GetFlatMask(
  const Array2D<T>         &elevations,
  Array2D<int32_t>         &flat_mask,
  Array2D<int32_t>         &labels
){
  Timer timer;
  timer.start();

  std::deque<GridCell> low_edges;
  std::deque<GridCell> high_edges;  //TODO: Need estimate of size

  RDLOG_ALG_NAME<<"Barnes (2014) Flat Resolution Flat Mask Generation";
  RDLOG_CITATION<<"Barnes, R., Lehman, C., Mulla, D., 2014a. An efficient assignment of drainage direction over flat surfaces in raster digital elevation models. Computers & Geosciences 62, 128–135. doi:10.1016/j.cageo.2013.01.009";
  
  Array2D<int8_t> flats;
  FindFlats(elevations, flats);

  RDLOG_PROGRESS<<"Setting up labels matrix...";
  labels.templateCopy(elevations);
  labels.resize(elevations);
  labels.setAll(0);

  RDLOG_PROGRESS<<"Setting up flat resolution mask...";
  flat_mask.templateCopy(elevations);
  flat_mask.resize(elevations);
  flat_mask.setAll(0);
  flat_mask.setNoData(-1);

  FindFlatEdges(low_edges, high_edges, flats, elevations);

  if(low_edges.size()==0){
    if(high_edges.size()>0)
      RDLOG_WARN<<"There were flats, but none of them had outlets! Quitting flat resolution.";
    else
      RDLOG_WARN<<"There were no flats! Quitting flat resolution.";
    return;
  }

  RDLOG_PROGRESS<<"Labeling flats...";
  int group_number=1;
  for(const auto &i: low_edges)
    if(labels(i.x,i.y)==0) //If the cell has not already been labeled
      LabelFlat(i.x, i.y, group_number++, labels, elevations);

  RDLOG_MISC<<"Unique flats = "<<group_number;

  RDLOG_PROGRESS<<"Removing flats without outlets from the queue...";
  std::deque<GridCell> temp;
  for(const auto &i: high_edges)
    if(labels(i.x,i.y)!=0)
      temp.push_back(i);

  if(temp.size()<high_edges.size())  //TODO: Prompt for intervention?
    RDLOG_WARN<<"Not all flats have outlets; the DEM contains sinks/pits/depressions!";
  high_edges = temp;
  temp.clear();

  RDLOG_MEM_USE<<"The flat height vector will require approximately "
               <<(group_number*((long)sizeof(int))/1024/1024)
               <<"MB of RAM.";
        
  RDLOG_PROGRESS<<"Creating flat height vector...";
  std::vector<int> flat_height(group_number);

  BuildAwayGradient           (flats, flat_mask, high_edges, flat_height, labels);
  BuildTowardsCombinedGradient(flats, flat_mask, low_edges,  flat_height, labels);

  RDLOG_TIME_USE<<"Wall-time = "<<timer.stop()<<" s";
}



/**
  @brief  Alters the elevations of the DEM so that all flats drain
  @author Richard Barnes (rbarnes@umn.edu)

  This alters elevations within the DEM so that flats which have been
  resolved using GetFlatMask() will drain.

  @param[in]     &flat_mask   A mask from GetFlatMask()
  @param[in]     &labels      A grouping from GetFlatMask()
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
void ResolveFlatsEpsilon_Barnes2014(
  const Array2D<int32_t> &flat_mask,
  const Array2D<int32_t> &labels,
  Array2D<U>             &elevations
){
  ProgressBar progress;

  RDLOG_ALG_NAME<<"Barnes (2014) Flat Resolution (DEM modification)...";
  RDLOG_CITATION<<"Barnes, R., Lehman, C., Mulla, D., 2014a. An efficient assignment of drainage direction over flat surfaces in raster digital elevation models. Computers & Geosciences 62, 128–135. doi:10.1016/j.cageo.2013.01.009";
  progress.start( flat_mask.size() );

  int raise_warn = 0;

  #pragma omp parallel for collapse(2)
  for(int y=1;y<flat_mask.height()-1;y++)
  for(int x=1;x<flat_mask.width()-1;x++){
    ++progress;

    //This cell is not part of a flat, so it does not need to be considered
    if(labels(x,y)==0)
      continue;

    //This array holds information about which cells the focal cell was lower
    //than prior to the focal cell being raised. This is used to check to see if
    //the focal cell has been raised too high.
    bool lower[9];
    for(int n=1;n<=8;++n)
      lower[n] = elevations(x,y)<elevations(x+dx[n],y+dy[n]);

    //Raise the focal cell by the appropriate number of increments
    for(int i=0;i<flat_mask(x,y);++i)
      elevations(x,y) = std::nextafter(elevations(x,y),std::numeric_limits<U>::infinity());

    //Check the surrounding cells to see if we are inappropriately higher than
    //any of them
    for(int n=1;n<=8;++n){
      const int nx = x+dx[n];
      const int ny = y+dy[n];
      //This neighbour is part of the flat, so it does not matter if we are
      //higher than it
      if(labels(nx,ny)==labels(x,y))
        continue;
      //This neighbour is higher than me, so it does not matter
      if(elevations(x,y)<elevations(nx,ny))
        continue;
      //I am higher than or of equal elevation to my neighbour, but I was
      //originally lower than my neighbour. This is a problem.
      if(lower[n])
        raise_warn++;
    }
  }
  RDLOG_WARN<<"Cells inappropriately raised above surrounding terrain = "<<raise_warn;
  RDLOG_TIME_USE<<"Succeeded in = "<<progress.stop()<<" s";
}













/**
  @brief  Calculates flow directions in flats
  @author Richard Barnes (rbarnes@umn.edu)

  This determines flow directions within flats which have been resolved
  using GetFlatMask().

  Uses the helper function D8MaskedFlowdir()

  @param[in]  &flat_mask      A mask from GetFlatMask()
  @param[in]  &labels         The labels output from GetFlatMask()
  @param[out] &flowdirs       Returns flat-resolved flow directions

  @pre
    1. **flat_mask** contains the number of increments to be applied to each
       cell to form a gradient which will drain the flat it is a part of.
    2. Any cell without a local gradient has a value of NO_FLOW_GEN in
       **flowdirs**; all other cells have defined flow directions.
    3. If a cell is part of a flat, it has a value greater than zero in
       **labels** indicating which flat it is a member of; otherwise, it has a
       value of 0.

  @post
    1. Every cell whose flow direction could be resolved by this algorithm
       (all drainable flats) will have a defined flow direction in
       **flowdirs**. Any cells which could not be resolved (non-drainable
       flats) will still be marked NO_FLOW_GEN.
*/
template<class U>
void ResolveFlatsFlowdirs_Barnes2014(
  const Array2D<int32_t>   &flat_mask,
  const Array2D<int32_t>   &labels,
  Array2D<U>               &flowdirs
){
  ProgressBar progress;

  RDLOG_ALG_NAME<<"Barnes (2014) Flat Resolution Flow Direction Modification";
  RDLOG_CITATION<<"Barnes, R., Lehman, C., Mulla, D., 2014a. An efficient assignment of drainage direction over flat surfaces in raster digital elevation models. Computers & Geosciences 62, 128–135. doi:10.1016/j.cageo.2013.01.009";

  progress.start( flat_mask.size() );

  #pragma omp parallel for collapse(2)
  for(int y=1;y<flat_mask.width()-1;y++)
  for(int x=1;x<flat_mask.height()-1;x++){
    const int ci = y*flat_mask.width()+x;

    if(flat_mask.isNoData(x,y))
      continue;
    if (flowdirs.at(9*ci)!=NO_FLOW_GEN)
      continue;

    int minimum_elevation = flat_mask(ci);
    int flowdir           = NO_FLOW_GEN;

    //It is safe to do this without checking to see that (nx,ny) is within
    //the grid because we only call this function on interior cells
    for(int n=1;n<=8;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];
      if( labels(nx,ny)!=labels(ci))
        continue;
      if( flat_mask(nx,ny)<minimum_elevation || (flat_mask(nx,ny)==minimum_elevation && flowdir>0 && flowdir%2==0 && n%2==1) ){
        minimum_elevation=flat_mask(nx,ny);
        flowdir=n;
      }
    }

    flowdirs.at(9*ci+0)       = HAS_FLOW_GEN;
    flowdirs.at(9*ci+flowdir) = 1;
  }

  RDLOG_TIME_USE<<"Succeeded in = "<<progress.stop()<<" s";
}

}

#endif
