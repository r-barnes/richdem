#include "Array2D.hpp"
#include "common.hpp"
#include <iostream>
#include <queue>
#include <time.h>
#include <cstdint>

const uint8_t GRID_LEFT   = 1;
const uint8_t GRID_TOP    = 2;
const uint8_t GRID_RIGHT  = 4;
const uint8_t GRID_BOTTOM = 8;

class Flag{
 public:
  std::vector< std::vector<bool> > flagArray;

  Flag(int width, int height){
    flagArray.resize(height, std::vector<bool>());
    for(int y=0;y<height;y++)
      flagArray[y].resize(width,false);
  }

  void set(int x, int y){
    flagArray[y][x] = true;
  }

  int isSet(int x, int y){
    return flagArray[y][x];
  }
};

template<class elev_t, class label_t>
void WatershedsMeet(
  label_t my_label, label_t n_label, elev_t my_elev, elev_t n_elev,
  std::map<label_t, std::map<label_t, elev_t> > &my_graph
){
  if(n_label==0)
    return;
  if(my_label==n_label)
    return;

  auto elev_over = std::max(my_elev,n_elev); //TODO: I think this should always be the neighbour.
  //If count()==0, then we haven't seen this watershed before.
  //Otherwise, only make a note of the spill-over elevation if it is
  //lower than what we've seen before.
  if(my_graph[my_label].count(n_label)==0 || elev_over<my_graph[my_label][n_label]){
    my_graph[my_label][n_label] = elev_over;
    my_graph[n_label][my_label] = elev_over;
  }
}

template<class elev_t, class label_t>
void ProcessTraceQue_onepass(Array2D<elev_t> &dem, Array2D<label_t> &labels, Flag &flag, std::queue<GridCellZ<elev_t> > &traceQueue, GridCellZ_pq<elev_t> &priorityQueue, std::map<label_t, std::map<label_t, elev_t> > &my_graph){
  while (!traceQueue.empty()){
    GridCellZ<elev_t> c = traceQueue.front();
    traceQueue.pop();

    bool bInPQ = false;
    for(int n=1;n<=8;n++){
      int nx = c.x+dx[n];
      int ny = c.y+dy[n];

      WatershedsMeet(labels(c.x,c.y),labels(nx,ny),dem(c.x,c.y),dem(nx,ny),my_graph);

      if (flag.isSet(nx,ny))
        continue;    
      
      
      //The neighbour is unprocessed and higher than the central cell
      if(c.z<dem(nx,ny)){
        traceQueue.emplace(nx,ny,dem(nx,ny));
        labels(nx,ny) = labels(c.x,c.y);
        flag.set(nx,ny);
        continue;
      }

      //Decide  whether (nx, ny) is a true border cell
      if (!bInPQ) {
        bool isBoundary = true;
        for(int nn=1;nn<=8;nn++){
          int nnx = nx+dx[n];
          int nny = ny+dy[n];
          if (flag.isSet(nnx,nny) && dem(nnx,nny)<dem(nx,ny)){
            isBoundary = false;
            break;
          }
        }
        if(isBoundary){
          priorityQueue.push(c);
          bInPQ = true;
        }
      }
    }
  }
}

template<class elev_t, class label_t>
void ProcessPit_onepass(Array2D<elev_t> &dem, Array2D<label_t> &labels, Flag &flag, std::queue<GridCellZ<elev_t> > &depressionQue, std::queue<GridCellZ<elev_t> > &traceQueue, GridCellZ_pq<elev_t> &priorityQueue, std::map<label_t, std::map<label_t, elev_t> > &my_graph){
  while (!depressionQue.empty()){
    GridCellZ<elev_t> c = depressionQue.front();
    depressionQue.pop();

    for(int n=1;n<=8;n++){
      int nx = c.x+dx[n];
      int ny = c.y+dy[n];

      WatershedsMeet(labels(c.x,c.y),labels(nx,ny),dem(c.x,c.y),dem(nx,ny),my_graph);

      if (flag.isSet(nx,ny))
        continue;    

      labels(nx,ny) = labels(c.x,c.y);
      flag.set(nx,ny);

      if (dem(nx,ny) > c.z) { //Slope cell
        traceQueue.emplace(nx,ny,dem(nx,ny));
      } else {                //Depression cell
        dem(nx,ny) = c.z;
        depressionQue.emplace(nx,ny,c.z);
      }
    }
  }
}

template<class elev_t, class label_t>
void Zhou2015Labels(
  Array2D<elev_t>                               &dem,
  Array2D<label_t>                              &labels,
  label_t                                        current_label, //NOTE: Should start at at least 2 (TODO: Explain why)
  std::map<label_t, std::map<label_t, elev_t> > &my_graph,
  uint8_t edge
){
  std::queue<GridCellZ<elev_t> > traceQueue;
  std::queue<GridCellZ<elev_t> > depressionQue;

  Flag flag(dem.viewWidth(), dem.viewHeight());

  labels.init(0);

  GridCellZ_pq<elev_t> priorityQueue;

  for(int x=0;x<dem.viewWidth();x++){
    const int height = dem.viewHeight()-1;
    priorityQueue.emplace(x,0,     dem(x,0     ));
    priorityQueue.emplace(x,height,dem(x,height));
    labels(x,     0) = current_label++;
    labels(x,height) = current_label++;
    flag.set(x,     0);
    flag.set(x,height);
  }

  for(int y=1;y<dem.viewHeight()-1;y++){
    const int width = dem.viewWidth()-1;
    priorityQueue.emplace(0,    y,dem(0,    y));
    priorityQueue.emplace(width,y,dem(width,y));
    labels(0,    y) = current_label++;
    labels(width,y) = current_label++;
    flag.set(0,    y);
    flag.set(width,y);
  }

  if(edge & GRID_TOP)
    labels.setRow(0,1);

  if(edge & GRID_BOTTOM)
    labels.setRow(dem.viewHeight()-1,1);

  if(edge & GRID_LEFT)
    labels.setCol(0,1);

  if(edge & GRID_RIGHT)
    labels.setCol(dem.viewWidth()-1,1);

  while (!priorityQueue.empty()){
    GridCellZ<elev_t> c = priorityQueue.top();
    priorityQueue.pop();

    for(int n=1;n<=8;n++){
      int nx = c.x+dx[n];
      int ny = c.y+dy[n];

      if (!dem.in_grid(nx,ny))
        continue;

      WatershedsMeet(labels(c.x,c.y),labels(nx,ny),dem(c.x,c.y),dem(nx,ny),my_graph);

      if(flag.isSet(nx,ny))
        continue;

      labels(nx,ny) = labels(c.x,c.y);
      flag.set(nx,ny);

      if(dem(nx,ny)<=c.z){ //Depression cell
        dem(nx,ny) = c.z;
        depressionQue.emplace(nx,ny,c.z);
        ProcessPit_onepass(dem,labels,flag,depressionQue,traceQueue,priorityQueue,my_graph);
      } else {          //Slope cell
        traceQueue.emplace(nx,ny,dem(nx,ny));
      }     
      ProcessTraceQue_onepass(dem,labels,flag,traceQueue,priorityQueue,my_graph);
    }
  }
}