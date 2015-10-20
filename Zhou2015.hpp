#include "Array2D.hpp"
#include "data_structures.h"
#include <iostream>
#include <queue>
#include <time.h>
#include <cstdint>

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

  int isProcessed(int x, int y){
    //if the cell is outside the DEM, it is regared as processed
    if (y<0 || y>=flagArray.size() || x<0 || x>=flagArray[0].size())
      return true;
    return flagArray[y][x];
  }

  int isProcessedDirect(int x, int y){
    return flagArray[y][x];
  }
};

template<class elev_t>
void ProcessTraceQue_onepass(Array2D<elev_t> &dem, Flag &flag, std::queue<grid_cellz<elev_t> > &traceQueue, grid_cellz_pq<elev_t> &priorityQueue, int &count){
  int  total  = 0;
  int  nPSC   = 0;
  bool bInPQ  = false;
  while (!traceQueue.empty()){
    grid_cellz<elev_t> c = traceQueue.front();
    traceQueue.pop();
    total++;

    bInPQ=false;
    for(int n=1;n<=8;n++){
      int nx = c.x+dx[n];
      int ny = c.y+dy[n];

      if (flag.isProcessedDirect(nx,ny))
        continue;    
      
      elev_t n_spill = dem(nx,ny);
      
      if (n_spill <= c.z)   {
        if (!bInPQ) {
          //decide  whether (nx, ny) is a true border cell
          bool isBoundary=true;
          for(int nn=1;nn<=8;nn++){
            int nnx = nx+dx[n];
            int nny = ny+dy[n];
            if (flag.isProcessedDirect(nnx,nny) && dem(nnx,nny)<n_spill){
              isBoundary = false;
              break;
            }
          }
          if(isBoundary){
            priorityQueue.push(c);
            bInPQ = true;
            nPSC++;
          }
        }
        continue; 
      }
      //otherwise
      //N is unprocessed and N is higher than C
      traceQueue.emplace(nx,ny,n_spill);
      flag.set(nx,ny);    
    }
  }
  count+=total-nPSC;
}

template<class elev_t>
void ProcessPit_onepass(Array2D<elev_t> &dem, Flag& flag, std::queue<grid_cellz<elev_t> > &depressionQue, std::queue<grid_cellz<elev_t> > &traceQueue, grid_cellz_pq<elev_t> &priorityQueue,int &count){
  while (!depressionQue.empty()){
    grid_cellz<elev_t> c = depressionQue.front();
    depressionQue.pop();
    count++;

    for(int n=1;n<=8;n++){
      int nx = c.x+dx[n];
      int ny = c.y+dy[n];

      if (flag.isProcessedDirect(nx,ny))
        continue;    

      elev_t n_spill = dem(nx,ny);
      if (n_spill > c.z) { //slope cell
        traceQueue.emplace(nx,ny,n_spill);
        flag.set(nx,ny);
        continue;
      }

      //depression cell
      flag.set(nx,ny);
      dem(nx,ny) = c.z;
      depressionQue.emplace(nx,ny,c.z);
    }
  }
}

template<class elev_t>
void Zhou2015_PriorityFlood(Array2D<elev_t> &dem){
  std::queue<grid_cellz<elev_t> > traceQueue;
  std::queue<grid_cellz<elev_t> > depressionQue;

  Flag flag(dem.viewWidth(), dem.viewHeight());

  grid_cellz_pq<elev_t> priorityQueue;
  int count               = 0;

  //Load border cells and cells adjacent to NoData cells into the priority queue
  int validElementsCount = 0;
  for (int y=0;y<dem.viewHeight();y++)  //row
  for (int x=0;x<dem.viewWidth ();x++){ //col
    if(dem.isNoData(x,y)){
      flag.set(x,y);
      continue;
    }

    validElementsCount++;
    for(int n=1;n<=8;n++){
      int nx = x+dx[n];
      int ny = y+dy[n];
      if(!dem.in_grid(nx,ny) || dem.isNoData(nx,ny)){
        priorityQueue.emplace(x,y,dem(x,y));
        flag.set(x,y);
        break;
      }
    }
  }

  while (!priorityQueue.empty()){
    grid_cellz<elev_t> c = priorityQueue.top();
    priorityQueue.pop();

    count++;

    for(int n=1;n<=8;n++){
      int nx = c.x+dx[n];
      int ny = c.y+dy[n];

      if (flag.isProcessed(nx,ny))
        continue;

      elev_t n_spill = dem(nx,ny);
      if(n_spill<=c.z){
        //depression cell
        dem(nx,ny) = c.z;
        flag.set(nx,ny);
        depressionQue.emplace(nx,ny,c.z);
        ProcessPit_onepass(dem,flag,depressionQue,traceQueue,priorityQueue,count);
      } else {
        //slope cell
        flag.set(nx,ny);
        traceQueue.emplace(nx,ny,n_spill);
      }     
      ProcessTraceQue_onepass(dem,flag,traceQueue,priorityQueue,count);
    }
  }
}