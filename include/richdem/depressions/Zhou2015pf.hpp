#ifndef __distpdf_zhou2015pf_hp__
#define __distpdf_zhou2015pf_hp__

#include "richdem/common/Array2D.hpp"
#include <queue>
#include <vector>
#include <map>
#include <iostream> //TODO

typedef char label_t;

template<class elev_t>
void ProcessTraceQue_onepass(Array2D<elev_t> &dem, Array2D<label_t> &labels, std::queue<int> &traceQueue, std::priority_queue<std::pair<elev_t, int>, std::vector< std::pair<elev_t, int> >, std::greater< std::pair<elev_t, int> > > &priorityQueue){
  while (!traceQueue.empty()){
    int c = traceQueue.front();
    traceQueue.pop();

    bool bInPQ = false;
    for(int n=1;n<=8;n++){
      int ni = dem.nToI(c, dx[n], dy[n]);
      if(ni==-1)
        continue;

      if(labels(ni)!=0)
        continue;    
      
      //The neighbour is unprocessed and higher than the central cell
      if(dem(c)<dem(ni)){
        traceQueue.emplace(ni);
        labels(ni) = labels(c);
        continue;
      }

      //Decide  whether (nx, ny) is a true border cell
      if (!bInPQ) {
        bool isBoundary = true;
        for(int nn=1;nn<=8;nn++){
          int nni = dem.nToI(ni, dx[n], dy[n]);
          if(nni==-1)
            continue;

          if (labels(nni)!=0 && dem(nni)<dem(ni)){
            isBoundary = false;
            break;
          }
        }
        if(isBoundary){
          priorityQueue.emplace(dem(c),c);
          bInPQ = true;
        }
      }
    }
  }
}

template<class elev_t>
void ProcessPit_onepass(elev_t c_elev, Array2D<elev_t> &dem, Array2D<label_t> &labels, std::queue<int> &depressionQue, std::queue<int> &traceQueue, std::priority_queue<std::pair<elev_t, int>, std::vector< std::pair<elev_t, int> >, std::greater< std::pair<elev_t, int> > > &priorityQueue){
  while (!depressionQue.empty()){
    int c = depressionQue.front();
    depressionQue.pop();

    for(int n=1;n<=8;n++){
      int ni = dem.nToI(c, dx[n], dy[n]);
      if(ni==-1)
        continue;

      if(labels(ni)!=0)
        continue;    

      labels(ni) = labels(c);

      if (dem(ni) > c_elev) { //Slope cell
        traceQueue.emplace(ni);
      } else {                //Depression cell
        dem(ni) = c_elev;
        depressionQue.emplace(ni);
      }
    }
  }
}

template<class elev_t>
void Zhou2015Labels(
  Array2D<elev_t>                         &dem
){
  std::queue<int> traceQueue;
  std::queue<int> depressionQue;

  Timer timer;
  timer.start();

  Array2D<label_t> labels;
  labels.resize(dem);

  labels.setAll(0);

  std::priority_queue<std::pair<elev_t, int>, std::vector< std::pair<elev_t, int> >, std::greater< std::pair<elev_t, int> > > priorityQueue;

  auto PlaceCell = [&](int x, int y){
    int i = dem.xyToI(x,y);
    priorityQueue.emplace(dem(i),i);
  };

  for(int x=0;x<dem.width();x++)     //Top Row
    PlaceCell(x,0);

  for(int x=0;x<dem.width();x++)     //Bottom Row
    PlaceCell(x,dem.height()-1);

  for(int y=1;y<dem.height()-1;y++)  //Left Column
    PlaceCell(0,y);

  for(int y=1;y<dem.height()-1;y++)  //Right Column
    PlaceCell(dem.width()-1,y);

  while (!priorityQueue.empty()){
    std::pair<elev_t, int> c = priorityQueue.top();
    priorityQueue.pop();

    labels(c.second) = 10;

    for(int n=1;n<=8;n++){
      int ni = dem.nToI(c.second, dx[n], dy[n]);
      if(ni==-1)
        continue;

      if(labels(ni)!=0)
        continue;

      labels(ni) = labels(c.second);

      if(dem(ni)<=c.first){ //Depression cell
        dem(ni) = c.first;
        depressionQue.emplace(ni);
        ProcessPit_onepass(c.first,dem,labels,depressionQue,traceQueue,priorityQueue);
      } else {          //Slope cell
        traceQueue.emplace(ni);
      }     
      ProcessTraceQue_onepass(dem,labels,traceQueue,priorityQueue);
    }
  }

  timer.stop();
  std::cerr<<"Zhou2015 completed in "<<timer.accumulated()<<"s."<<std::endl;
}

#endif