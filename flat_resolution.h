//Flat Resolution
//Richard Barnes (rbarnes@umn.edu), 2012
//Develops an elevation mask which is guaranteed to drain
//a flat using a convergent flow pattern (unless it's a mesa)

#ifndef _flat_resolution_included
#define _flat_resolution_included

#include "utility.h"
#include "interface.h"
#include "data_structures.h"
#include <deque>
#include <vector>
#include <queue>
#include "debug.h"
#ifdef _OPENMP
  #include <omp.h>
#endif



//Procedure: BuildAwayGradient
/**
  @brief Build a gradient away from the high edges of the flats
  @author Richard Barnes

  The queue of high-edge cells developed in find_flat_edges() is copied
  into the procedure. A breadth-first expansion labels cells by their
  distance away from terrain of higher elevation. The maximal distance
  encountered is noted.

  @param[in]  &elevations   A 2D array of cell elevations
  @param[in]  &flowdirs     A 2D array indicating each cell's flow direction
  @param[out] &incrementations  A 2D array for storing incrementations
  @param[in]  &edges        The high-edge FIFO queue from find_flat_edges()
  @param[out] &flat_height  Vector with length equal to max number of labels
  @param[in]  &labels       2D array storing labels developed in label_this()

  @pre
    "incrementations" is initiliazed to "0"

  @post
    "incrementations" will contain the D8 distance of every flat cell from
    terrain of differing elevation
    "flat_height" will contain, for each flat, the maximal distance any
    of its cells are from terrain of differing elevation
*/
template <class T, class U>
static void BuildAwayGradient(
  const array2d<T> &elevations, const array2d<U> &flowdirs,
  int_2d &flat_mask, std::deque<grid_cell> edges, 
  std::vector<int> &flat_height, const int_2d &labels
){
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
      if(elevations.in_grid(nx,ny) 
          && elevations(nx,ny)==elevations(x,y) 
          && flowdirs(nx,ny)==NO_FLOW)
        edges.push_back(grid_cell(nx,ny));
    }
  }

  diagnostic("succeeded!\n");
}



//Procedure: BuildTowardsCombinedGradient
/**
  @brief Builds gradient away from the low edges of flats, combines gradients
  @author Richard Barnes

  The queue of low-edge cells developed in find_flat_edges() is copied
  into the procedure. A breadth-first expansion labels cells by their
  distance away from terrain of lower elevation. This is combined with
  the gradient from BuildAwayGradient() to give the final increments of
  each cell in forming the flat mask.

  @param[in]  &elevations   A 2D array of cell elevations
  @param[in]  &flowdirs     A 2D array indicating each cell's flow direction
  @param[in,out]  &incrementations  A 2D array for storing incrementations
  @param[in]  &edges        The low-edge FIFO queue from find_flat_edges()
  @param[in]  &flat_height  Vector with length equal to max number of labels
  @param[in]  &labels       2D array storing labels developed in label_this()

  @pre
    "incrementations" has been adjusted by BuildAwayGradient()

  @post
    "incrementations" will contain a convergent flow pattern which drains
    every cell of the flat.
*/
template <class T, class U>
static void BuildTowardsCombinedGradient(
  const array2d<T> &elevations, const array2d<U> &flowdirs,
  int_2d &flat_mask, std::deque<grid_cell> edges, 
  std::vector<int> &flat_height, const int_2d &labels
){
  int x,y,nx,ny;
  int loops=1;
  grid_cell iteration_marker(-1,-1);

  diagnostic("Barnes flat resolution: toward and combined gradients...");

  //Make previous incrementations negative so that we can keep track of where we are
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
      if(elevations.in_grid(nx,ny) 
          && elevations(nx,ny)==elevations(x,y) 
          && flowdirs(nx,ny)==NO_FLOW)
        edges.push_back(grid_cell(nx,ny));
    }
  }

  diagnostic("succeeded!\n");
}


//Procedure: label_this
/**
  @brief Labels all the cells of a flat with a common label.
  @author Richard Barnes

  Performs a flood fill operation which labels all the cells of a flat
  with a common label. Each flat will have a unique label

  @param[in]  x           x-coordinate of flood fill seed
  @param[in]  y           y-coordinate of flood-fill seed
  @param[in]  label       Label to apply to the cells
  @param[out] &labels     2D array which will contain the labels
  @param[in]  &elevations 2D array of cell elevations

  @pre
    "labels" is initialized to "-1"
  @post
    The 2D array "labels" is modified such that each cell
    which can be reached from (x,y) while traversing only
    cells which share the same elevation as (x,y)
    is labeled with "label"
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
    if(elevations(c.x,c.y)!=target_elevation || labels(c.x,c.y)>0) continue;
    labels(c.x,c.y)=label;
    for(int n=1;n<=8;n++)
      if(labels.in_grid(c.x+dx[n],c.y+dy[n]))
        to_fill.push(grid_cell(c.x+dx[n],c.y+dy[n]));
  }
}

//Procedure: find_flat_edges
/**
  @brief Identifies cells adjacent to higher and lower terrain
  @author Richard Barnes

  Cells adjacent to lower and higher terrain are identified and
  added to the appropriate queue

  @param[out] &low_edges  Queue for storing cells adjacent to lower terrain
  @param[out] &high_edges Queue for storing cells adjacent to higher terrain
  @param[in]  &flowdirs   2D array indicating flow direction for each cell
  @param[in]  &elevations 2D array of cell elevations

  @pre
    "flowdirs" should be set such that cells without a defined flow direction
    have a value #NO_FLOW.
    "low_edges" and "high_edges" should be traversible
*/
template <class T, class U>
static void find_flat_edges(
  std::deque<grid_cell> &low_edges,
  std::deque<grid_cell> &high_edges,
  const array2d<U> &flowdirs,
  const array2d<T> &elevations
){
  int nx,ny;
  ProgressBar progress;
  diagnostic("%%Searching for flats...\n");
  progress.start( flowdirs.width()*flowdirs.height() );
  for(int x=0;x<flowdirs.width();x++){
    progress.update( x*flowdirs.height() );
    for(int y=0;y<flowdirs.height();y++){
      if(flowdirs(x,y)==flowdirs.no_data)
        continue;
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
}


//Procedure: resolve_flats_barnes
/**
  @brief  Performs the flat resolution by Barnes, Lehman, and Mulla.
  @author Richard Barnes

  TODO

  @param[in]  &elevations 2D array of cell elevations
  @param[in]  &flowdirs   2D array indicating flow direction of each cell
  @param[in]  &flat_mask  2D array which will hold incremental elevation mask
  @param[in]  &labels     2D array indicating flat membership
*/
template <class T, class U>
void resolve_flats_barnes(
  const array2d<T> &elevations,
  const array2d<U> &flowdirs,
  int_2d &flat_mask,
  int_2d &labels
){
  std::deque<grid_cell> low_edges,high_edges;  //TODO: Need estimate of size

  diagnostic("\n###Barnes Flat Resolution\n");

  diagnostic("Setting up labels matrix...");
  labels.copyprops(flowdirs);
  labels.init(-1);
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
  int group_number=0;
  for(std::deque<grid_cell>::iterator i=low_edges.begin();i!=low_edges.end();++i)
    if(labels(i->x,i->y)==-1)
      label_this(i->x, i->y, group_number++, labels, elevations);
  diagnostic("succeeded!\n");

  diagnostic_arg("Found %d unique flats.\n",group_number);

  diagnostic("Removing flats without outlets from the queue...");
  std::deque<grid_cell> temp;
  for(std::deque<grid_cell>::iterator i=high_edges.begin();i!=high_edges.end();++i)
    if(labels(i->x,i->y)!=-1)
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
    elevations, flowdirs, flat_mask,
    high_edges, flat_height, labels
  );
  BuildTowardsCombinedGradient(
    elevations, flowdirs, flat_mask,
    low_edges, flat_height, labels
  );
}





//Procedure: flat_mask
/**
  @brief  Creates binary mask indicating cells with undefined flow direction
  @author Richard Barnes

  Primarily for use in debugging (TODO).

  @param[in]  &flowdirs 2D array indicating flow direction for each cell
  @param[out] &fmask    2D binary array indicating whether cell is in a flat

  @post
  The 2D array fmask is modified such that the value 3 represents
  cells with a NO_DATA value in flowdirs, 1 represents cells without
  a defined flow direction, and 0 represents cells with a defined
  flow direction.
*/
template <class T>
void flat_mask(const array2d<T> &flowdirs, uint_2d &fmask){
  diagnostic("Locating pits based on undefined flow directions...");
  fmask.copyprops(flowdirs);
  fmask.init(0);
  fmask.no_data=3;
  #pragma omp parallel for
  for(int x=0;x<flowdirs.width();x++)
  for(int y=0;y<flowdirs.height();y++)
    if(flowdirs(x,y)==flowdirs.no_data)
      fmask(x,y)=fmask.no_data;
    else
      fmask(x,y)=(flowdirs(x,y)==NO_FLOW);
  diagnostic("succeeded!\n");
}
#endif
