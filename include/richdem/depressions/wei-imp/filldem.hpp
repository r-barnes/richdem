#ifndef _richdem_wei2008_hpp_
#define _richdem_wei2008_hpp_

#include <richdem/common/Array2D.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <queue>
#include <algorithm>
#include "utils.h"
#include <time.h>
#include <list>
#include <stack>
#include <vector>


namespace richdem {


#include <functional>
class Node {
 public:
  int row;
  int col;
  float spill;
  int N;
  Node(){
    row   = 0;
    col   = 0;
    spill = -9999.0;
    N     = -1;
  }

  Node(int row0, int col0, float spill0) : row(row0), col(col0), spill(spill0) {}

  struct Greater : public std::binary_function< Node, Node, bool > {
    bool operator()(const Node n1, const Node n2) const {
      return n1.spill > n2.spill;
    }
  };

  bool operator>(const Node& a)  {
    return this->spill > a.spill;
  }
};




typedef std::priority_queue<Node, std::vector<Node>, Node::Greater> PriorityQueue;

template<class T>
void InitPriorityQue(
  Array2D<T>& dem,
  Flag& flag,
  PriorityQueue& priorityQueue
){
  std::queue<Node> depressionQue;

  // push border cells into the PQ
  for(int y = 0; y < dem.height(); y++)
  for(int x = 0; x < dem.width(); x++){
    if (flag.IsProcessedDirect(y,x)) continue;

    if (dem.isNoData(x,y)) {
      flag.SetFlag(y,x);
      for (int i = 0; i < 8; i++){
        auto iRow = Get_rowTo(i, y);
        auto iCol = Get_colTo(i, x);
        if (flag.IsProcessed(iRow,iCol)) 
          continue;

        if (!dem.isNoData(iCol, iRow)){
          priorityQueue.emplace(iRow,iCol,dem(iCol, iRow));
          flag.SetFlag(iRow,iCol);
        }
      }
    } else if(dem.isEdgeCell(x,y)){
      //on the DEM border
      priorityQueue.emplace(y,x,dem(x,y));
      flag.SetFlag(y,x);          
    }
  }
}



template<class T>
void ProcessTraceQue(
  Array2D<T>& dem,
  Flag& flag,
  std::queue<Node>& traceQueue,
  PriorityQueue& priorityQueue
){
  std::queue<Node > potentialQueue;
  int indexThreshold=2;  //index threshold, default to 2
  while (!traceQueue.empty()){
    const auto node = traceQueue.front();
    traceQueue.pop();
    bool Mask[5][5]={{false},{false},{false},{false},{false}};
    for (int i = 0; i < 8; i++){
      auto iRow = Get_rowTo(i,node.row);
      auto iCol = Get_colTo(i,node.col);
      if(flag.IsProcessedDirect(iRow,iCol))
        continue;

      if (dem(iCol,iRow)>node.spill){
        traceQueue.emplace(iRow,iCol, dem(iCol,iRow));
        flag.SetFlag(iRow,iCol);
      } else {
        //initialize all masks as false   
        bool have_spill_path_or_lower_spill_outlet=false; //whether cell i has a spill path or a lower spill outlet than node if i is a depression cell
        for(int k = 0; k < 8; k++){
          auto kRow = Get_rowTo(k,iRow);
          auto kCol = Get_colTo(k,iCol);
          if((Mask[kRow-node.row+2][kCol-node.col+2]) ||
            (flag.IsProcessedDirect(kRow,kCol)&&dem(kCol,kRow)<node.spill)
            )
          {
            Mask[iRow-node.row+2][iCol-node.col+2]=true;
            have_spill_path_or_lower_spill_outlet=true;
            break;
          }
        }
        
        if(!have_spill_path_or_lower_spill_outlet){
          if (i<indexThreshold) potentialQueue.push(node);
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
    for (int i = 0; i < 8; i++){
      auto iRow = Get_rowTo(i,node.row);
      auto iCol = Get_colTo(i,node.col);
      if(flag.IsProcessedDirect(iRow,iCol))
        continue;

      priorityQueue.push(node);
      break;
    }   
  }
}



template<class T>
void ProcessPit(
  Array2D<T>& dem, 
  Flag& flag, 
  std::queue<Node>& depressionQue,
  std::queue<Node>& traceQueue,
  PriorityQueue& priorityQueue
){
  while (!depressionQue.empty()){
    auto node = depressionQue.front();
    depressionQue.pop();
    for (int i = 0; i < 8; i++){
      auto iRow = Get_rowTo(i, node.row);
      auto iCol = Get_colTo(i,  node.col);
      if (flag.IsProcessedDirect(iRow,iCol))
        continue;    

      auto iSpill = dem(iCol,iRow);
      if (iSpill > node.spill){ //slope cell
        flag.SetFlag(iRow,iCol);
        traceQueue.emplace(iRow,iCol,iSpill);
        continue;
      }

      //depression cell
      flag.SetFlag(iRow,iCol);
      dem(iCol, iRow) = node.spill;
      depressionQue.emplace(iRow,iCol,node.spill);
    }
  }
}



template<class T>
void fillDEM(Array2D<T> &dem){
  std::queue<Node> traceQueue;
  std::queue<Node> depressionQue;
  
  time_t timeStart, timeEnd;
  std::cout<<"Using our proposed variant to fill DEM"<<std::endl;
  timeStart = time(NULL);
  Flag flag;
  if (!flag.Init(dem.width(),dem.height())) {
    printf("Failed to allocate memory!\n");
    return;
  }

  PriorityQueue priorityQueue;

  int numberofall   = 0;
  int numberofright = 0;

  InitPriorityQue(dem,flag,priorityQueue); 
  while (!priorityQueue.empty()){
    auto tmpNode = priorityQueue.top();
    priorityQueue.pop();
    auto row   = tmpNode.row;
    auto col   = tmpNode.col;
    auto spill = tmpNode.spill;

    for (int i = 0; i < 8; i++){
      auto iRow = Get_rowTo(i, row);
      auto iCol = Get_colTo(i, col);

      if (flag.IsProcessed(iRow,iCol))
        continue;

      auto iSpill = dem(iCol,iRow);
      if (iSpill <= spill){
        //depression cell
        dem(iCol,iRow) = spill;
        flag.SetFlag(iRow,iCol);
        depressionQue.emplace(iRow,iCol,spill);
        ProcessPit(dem,flag,depressionQue,traceQueue,priorityQueue);
      } else {
        //slope cell
        flag.SetFlag(iRow,iCol);
        traceQueue.emplace(iRow,iCol,iSpill);
      }     
      ProcessTraceQue(dem,flag,traceQueue,priorityQueue); 
    }
  }
  timeEnd = time(NULL);
  double consumeTime = difftime(timeEnd, timeStart);
  std::cout<<"Time used:"<<consumeTime<<" seconds"<<std::endl;
}

}

#endif
