/**
  @file
  @brief Defines the Priority-Flood algorithm described by Zhou, G., Sun, Z., Fu, S., 2016. An efficient variant of the Priority-Flood algorithm for filling depressions in raster digital elevation models. Computers & Geosciences 90, Part A, 87 – 96. doi:http://dx.doi.org/10.1016/j.cageo.2016.02.021.

    The code herein has been extensive modified by Richard Barnes (rbarnes@umn.edu) for inclusion with RichDEM.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_zhou2016pf_hpp_
#define _richdem_zhou2016pf_hpp_

#include <richdem/common/logger.hpp>
#include <richdem/common/Array2D.hpp>
#include <richdem/common/timer.hpp>
#include <queue>
#include <vector>
#include <map>
#include <iostream>

namespace richdem {

typedef char label_t;

template<class elev_t>
void ProcessTraceQue_onepass(
  Array2D<elev_t> &dem,
  Array2D<label_t> &labels,
  std::queue<int> &traceQueue,
  std::priority_queue<std::pair<elev_t, int>, std::vector< std::pair<elev_t, int> >, std::greater< std::pair<elev_t, int> > > &priorityQueue
){
  while (!traceQueue.empty()){
    auto c = traceQueue.front();
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
void ProcessPit_onepass(
  elev_t c_elev,
  Array2D<elev_t> &dem,
  Array2D<label_t> &labels,
  std::queue<int> &depressionQue,
  std::queue<int> &traceQueue,
  std::priority_queue<std::pair<elev_t, int>, std::vector< std::pair<elev_t, int> >, std::greater< std::pair<elev_t, int> > > &priorityQueue
){
  while (!depressionQue.empty()){
    auto c = depressionQue.front();
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

/**
  @brief  Fills all pits and removes all digital dams from a DEM, quickly
  @author G. Zhou, Z. Sun, S. Fu, Richard Barnes (this implementation)

    Works similarly to the Priority-Flood described by Barnes et al. (2014), but
    reduces the number of items which must pass through the priority queue, thus
    achieving greater efficiencies.

  @param[in,out]  &dem   A grid of cell elevations

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM. Note that the _NoData_ value is assumed to
       be a negative number less than any actual data value.

  @post
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM.
    2. **elevations** contains no landscape depressions or digital dams.
*/
template<class elev_t>
void PriorityFlood_Zhou2016(
  Array2D<elev_t> &dem
){
  std::queue<int> traceQueue;
  std::queue<int> depressionQue;

  RDLOG_ALG_NAME<<"Priority-Flood (Zhou2016 version)";
  RDLOG_CITATION<<"Zhou, G., Sun, Z., Fu, S., 2016. An efficient variant of the Priority-Flood algorithm for filling depressions in raster digital elevation models. Computers & Geosciences 90, Part A, 87 – 96. doi:http://dx.doi.org/10.1016/j.cageo.2016.02.021";

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
  RDLOG_TIME_USE<<"Zhou2016 wall-time = "<<timer.accumulated()<<" s";
}

}

#endif
