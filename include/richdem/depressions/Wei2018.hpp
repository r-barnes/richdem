#ifndef _richdem_wei2008_hpp_
#define _richdem_wei2008_hpp_

#include <richdem/common/Array2D.hpp>
#include <richdem/common/logger.hpp>
#include <richdem/common/grid_cell.hpp>
#include <richdem/common/timer.hpp>
#include <iostream>
#include <queue>


namespace richdem {

template<class T>
static void InitPriorityQue(
  Array2D<T>& dem,
  Array2D<bool>& flag,
  GridCellZ_pq<T>& priorityQueue
){
  std::queue<GridCellZ<T> > depressionQue;

  // push border cells into the PQ
  for(int y = 0; y < dem.height(); y++)
  for(int x = 0; x < dem.width(); x++){
    if (flag(x,y)) continue;

    if (dem.isNoData(x,y)) {
      flag(x,y)=true;
      for (int n=1;n<=8; n++){
        const auto nx = x+dx[n];
        const auto ny = y+dy[n];

        if(!dem.inGrid(nx,ny))
          continue;

        if (flag(nx,ny)) 
          continue;

        if (!dem.isNoData(nx, ny)){
          priorityQueue.emplace(nx,ny,dem(nx, ny));
          flag(nx,ny)=true;
        }
      }
    } else if(dem.isEdgeCell(x,y)){
      //on the DEM border
      priorityQueue.emplace(x,y,dem(x,y));
      flag(x,y)=true;          
    }
  }
}



template<class T>
static void ProcessTraceQue(
  Array2D<T>& dem,
  Array2D<bool>& flag,
  std::queue<GridCellZ<T> >& traceQueue,
  GridCellZ_pq<T>& priorityQueue
){
  std::queue<GridCellZ<T>  > potentialQueue;
  int indexThreshold=2;  //index threshold, default to 2
  while (!traceQueue.empty()){
    const auto node = traceQueue.front();
    traceQueue.pop();
    bool Mask[5][5]={{false},{false},{false},{false},{false}};
    for (int n=1;n<=8; n++){
      const auto nx = node.x+dx[n];
      const auto ny = node.y+dy[n];
      if(flag(nx,ny))
        continue;

      if (dem(nx,ny)>node.z){
        traceQueue.emplace(nx,ny, dem(nx,ny));
        flag(nx,ny)=true;
      } else {
        //initialize all masks as false   
        bool have_spill_path_or_lower_spill_outlet=false; //whether cell n has a spill path or a lower spill outlet than node if n is a depression cell
        for(int k=1;k<=8; k++){
          const auto nny = ny+dy[k];
          const auto nnx = nx+dx[k];
          if((Mask[nny-node.y+2][nnx-node.x+2]) ||
            (flag(nnx,nny) && dem(nnx,nny)<node.z)
            )
          {
            Mask[ny-node.y+2][nx-node.x+2]=true;
            have_spill_path_or_lower_spill_outlet=true;
            break;
          }
        }
        
        if(!have_spill_path_or_lower_spill_outlet){
          if (n<indexThreshold) potentialQueue.push(node);
          else
            priorityQueue.push(node);
          break; // make sure node is not pushed twice into PQ
        }
      }
    }//end of for loop
  }

  while (!potentialQueue.empty()){
    const auto node = potentialQueue.front();
    potentialQueue.pop();

    //first case
    for (int n=1;n<=8; n++){
      const auto nx = node.x+dx[n];
      const auto ny = node.y+dy[n];
      if(flag(nx,ny))
        continue;

      priorityQueue.push(node);
      break;
    }   
  }
}



template<class T>
static void ProcessPit(
  Array2D<T>& dem, 
  Array2D<bool>& flag, 
  std::queue<GridCellZ<T> >& depressionQue,
  std::queue<GridCellZ<T> >& traceQueue,
  GridCellZ_pq<T>& priorityQueue
){
  while (!depressionQue.empty()){
    auto node = depressionQue.front();
    depressionQue.pop();
    for (int n=1;n<=8; n++){
      const auto nx = node.x+dx[n];
      const auto ny = node.y+dy[n];
      if (flag(nx,ny))
        continue;    

      const auto iSpill = dem(nx,ny);
      if (iSpill > node.z){ //slope cell
        flag(nx,ny)=true;
        traceQueue.emplace(nx,ny,iSpill);
        continue;
      }

      //depression cell
      flag(nx,ny) = true;
      dem(nx, ny) = node.z;
      depressionQue.emplace(nx,ny,node.z);
    }
  }
}



template<class T>
void PriorityFlood_Wei2018(Array2D<T> &dem){
  std::queue<GridCellZ<T> > traceQueue;
  std::queue<GridCellZ<T> > depressionQue;
  
  RDLOG_ALG_NAME<<"Priority-Flood (Wei2018 version)";
  RDLOG_CITATION<<"Wei, H., Zhou, G., Fu, S., 2018. Efficient Priority-Flood depression filling in raster digital elevation models. International Journal of Digital Earth 0, 1–13. https://doi.org/10.1080/17538947.2018.1429503";

  Timer timer;
  timer.start();

  Array2D<bool> flag(dem.width(),dem.height(),false);

  GridCellZ_pq<T> priorityQueue;

  int numberofall   = 0;
  int numberofright = 0;

  InitPriorityQue(dem,flag,priorityQueue); 
  while (!priorityQueue.empty()){
    const auto tmpNode = priorityQueue.top();
    priorityQueue.pop();

    for (int n=1;n<=8; n++){
      auto ny = tmpNode.y+dy[n];
      auto nx = tmpNode.x+dx[n];

      if(!dem.inGrid(nx,ny))
        continue;

      if(flag(nx,ny))
        continue;

      auto iSpill = dem(nx,ny);
      if (iSpill <= tmpNode.z){
        //depression cell
        dem(nx,ny) = tmpNode.z;
        flag(nx,ny) = true;
        depressionQue.emplace(nx,ny,tmpNode.z);
        ProcessPit(dem,flag,depressionQue,traceQueue,priorityQueue);
      } else {
        //slope cell
        flag(nx,ny) = true;
        traceQueue.emplace(nx,ny,iSpill);
      }     
      ProcessTraceQue(dem,flag,traceQueue,priorityQueue); 
    }
  }

  timer.stop();
  RDLOG_TIME_USE<<"Wei2018 wall-time = "<<timer.accumulated()<<" s";  
}

}

#endif
