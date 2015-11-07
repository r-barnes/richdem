#ifndef __distpdf_barnes2014pf_hp__
#define __distpdf_barnes2014pf_hp__

#include "Array2D.hpp"
#include "common.hpp"
#include <queue>
#include <vector>

template<class elev_t, class label_t>
label_t GetNewLabel(
  int x,
  int y,
  label_t &current_label,
  uint8_t edge,
  const Array2D<elev_t> &dem,
  const Array2D<label_t> &labels
){
  if(labels(x,y)!=0)
    return labels(x,y);

  if((edge & GRID_TOP) && y==0)
    return 1;

  if((edge & GRID_BOTTOM) && y==dem.viewHeight()-1)
    return 1;

  if((edge & GRID_LEFT) && x==0)
    return 1;

  if((edge & GRID_RIGHT) && x==dem.viewWidth()-1)
    return 1;

  for(int n=1;n<=8;n++){
    int nx = x+dx[n];
    int ny = y+dy[n];
    if(!dem.in_grid(nx,ny))
      continue;
    if(labels(nx,ny)!=0 && dem(nx,ny)<=dem(x,y))
      return labels(nx,ny);
  }

  return current_label++;
}

template<class elev_t, class label_t>
void PriorityFlood(
  Array2D<elev_t>                         &dem,
  Array2D<label_t>                        &labels,
  std::vector<std::map<label_t, elev_t> > &my_graph,
  uint8_t edge
){
  GridCellZ_pq<elev_t> pq;
  std::queue< GridCellZ<elev_t> > pit;

  label_t current_label = 2;

  labels.init(0);

  for(int x=0;x<dem.viewWidth();x++){
    const int the_y = dem.viewHeight()-1;
    pq.emplace(x,0,    dem(x,0    ));
    pq.emplace(x,the_y,dem(x,the_y));
  }

  for(int y=1;y<dem.viewHeight()-1;y++){
    const int the_x = dem.viewWidth()-1;
    pq.emplace(0,    y,dem(0,    y));
    pq.emplace(the_x,y,dem(the_x,y));
  }

  while(!pq.empty() || !pit.empty()){
    GridCellZ<elev_t> c;
    if(pit.size()>0){
      c=pit.front();
      pit.pop();
    } else {
      c=pq.top();
      pq.pop();
    }

    //At this point the cell's label is guaranteed to be positive and in the
    //range [1,MAX_INTEGER] (unless we overflow).
    auto my_label = labels(c.x,c.y) = GetNewLabel(c.x,c.y,current_label,edge,dem,labels);

    for(int n=1;n<=8;n++){
      auto nx = c.x+dx[n]; //Neighbour's x-coordinate
      auto ny = c.y+dy[n]; //Neighbour's y-coordinate

      //Check to see if the neighbour coordinates are valid
      if(!dem.in_grid(nx,ny))
        continue;

      auto n_label = labels(nx,ny); //Neighbour's label
      //Does the neighbour have a label? If so, it is part of the edge, has
      //already been assigned a label by a parent cell which must be of lower or
      //equal elevation to the current cell, or has already been processed, in
      //which case its elevation is lower or equal to this cell.
      if(n_label!=0){
        //If the neighbour's label were the same as the current cell's, then the
        //current cell's flow and the neighbour's flow eventually comingle. If
        //the neighbour's label is different it has been added by a cell whose
        //flow drains the opposite side of a watershed from this cell. Here, we
        //make a note of the height of that watershed.
        if(n_label!=my_label){
          auto elev_over = std::max(dem(nx,ny),dem(c.x,c.y)); //TODO: I think this should always be the neighbour.
          //If count()==0, then we haven't seen this watershed before.
          //Otherwise, only make a note of the spill-over elevation if it is
          //lower than what we've seen before.
          if(my_graph[my_label].count(n_label)==0 || elev_over<my_graph[my_label][n_label]){
            my_graph[my_label][n_label] = elev_over;
            my_graph[n_label][my_label] = elev_over;
          }
        }
      } else {
        //The neighbour is not one we've seen before, so mark it as being part of
        //our watershed and add it as an unprocessed item to the queue.
        labels(nx,ny) = labels(c.x,c.y);

        //If the neighbour is lower than this cell, elevate it to the level of
        //this cell so that a depression is not formed. The flow directions will
        //be fixed later, after all the depressions have been filled.
        if(dem(nx,ny)<=c.z){
          dem(nx,ny) = c.z;
          pit.emplace(nx,ny,c.z);
        //Otherwise, if the neighbour is higher, do not adjust its elevation.
        } else {
          pq.emplace(nx,ny,dem(nx,ny));
        }
      }
    }
  }

  my_graph.resize(current_label);
}

#endif