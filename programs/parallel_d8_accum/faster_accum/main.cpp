#include "richdem/common/Array2D.hpp"
#include "richdem/common/grid_cell.hpp"
#include "richdem/common/timer.hpp"
#include <cstdint>
#include <iostream>
#include <queue>

typedef uint8_t flowdir_t;
typedef int64_t accum_t;
typedef uint8_t c_dependency_t;

void FlowAccumulation(
  const Array2D<flowdir_t> &flowdirs,
  Array2D<accum_t>         &accum
){
  //Each cell that flows points to a neighbouring cell. But not every cell
  //is pointed at. Cells which are not pointed at are the peaks from which
  //flow originates. In order to calculate the flow accumulation we begin at
  //peaks and pass flow downwards. Once flow has been passed downwards, the
  //cell receiving the flow no longer needs to wait for the cell which
  //passed the flow. When the receiving cell has no cells on which it is
  //waiting, it then becomes a peak itself. The number of cells pointing at
  //a cell is its "dependency count". In this section of the code we find
  //each cell's dependency count.

  accum.resize(flowdirs,0);
  accum.setNoData(ACCUM_NO_DATA);

  Array2D<c_dependency_t> dependencies(flowdirs,0);
  for(int32_t y=0;y<flowdirs.height();y++)
  for(int32_t x=0;x<flowdirs.width();x++){
    if(flowdirs.isNoData(x,y)){  //This cell is a no_data cell
      accum(x,y) = ACCUM_NO_DATA;
      continue;                
    }         

    int n = flowdirs(x,y);       //The neighbour this cell flows into
    if(n==NO_FLOW)               //This cell does not flow into a neighbour
      continue;

    int nx = x+dx[n];            //x-coordinate of the neighbour
    int ny = y+dy[n];            //y-coordinate of the neighbour
      
    //Neighbour is not on the grid
    if(!flowdirs.inGrid(nx,ny))
      continue;

    //Neighbour is valid and is part of the grid. The neighbour depends on this
    //cell, so increment its dependency count.
    dependencies(nx,ny)++;
  }

  //Now that we know how many dependencies each cell has, we can determine which
  //cells are the peaks: the sources of flow. We make a note of where the peaks
  //are for later use.
  std::queue<GridCell> sources;
  for(int32_t y=0;y<dependencies.height();y++)
  for(int32_t x=0;x<dependencies.width();x++)
    //Valid cell with no dependencies: a peak!
    if(dependencies(x,y)==0 && !flowdirs.isNoData(x,y))
      sources.emplace(x,y);

  //Now that we know where the sources of flow are, we can start at this. Each
  //cell will have at least an accumulation of 1: itself. It then passes this
  //and any other accumulation it has gathered along its flow path to its
  //neighbour and decrements the neighbours dependency count. When a neighbour
  //has no more dependencies, it becomes a source.
  while(!sources.empty()){         //There are sources remaining
    GridCell c = sources.front();  //Grab a source. Order is not important here.
    sources.pop();                 //We've visited this source. Discard it.

    if(flowdirs.isNoData(c.x,c.y)) //Oh snap! This isn't a real cell!
      continue;

    accum(c.x,c.y)++;              //This is a real cell, and it accumulates
                                   //one cell's worth of flow automatically.

    int n = flowdirs(c.x,c.y);     //Who is this source's neighbour?

    if(n==NO_FLOW)                 //This cell doesn't flow anywhere.
      continue;                    //Move on to the next source.

    int nx = c.x+dx[n];            //Okay, this cell is going somewhere.
    int ny = c.y+dy[n];            //Make a note of where

    //This cell flows of the edge of the grid. Move on to next source.
    if(!flowdirs.inGrid(nx,ny))
      continue;
    //This cell flows into a no_data cell. Move on to next source.
    if(flowdirs.isNoData(nx,ny))
      continue;

    //This cell has a neighbour it flows into. Add to its accumulation.
    accum(nx,ny) += accum(c.x,c.y);
    //Decrement the neighbour's dependencies.
    dependencies(nx,ny)--;

    //The neighbour has no more dependencies, so it has become a source
    if(dependencies(nx,ny)==0)
      sources.emplace(nx,ny);
  }
}


int main(int argc, char **argv){
  if(argc!=2){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input File>"<<std::endl;
    return -1;
  }

  Array2D<flowdir_t> fds(argv[1], false);
  Array2D<accum_t> accum;

  Timer timer_calc;
  timer_calc.start();
  FlowAccumulation(fds,accum);
  timer_calc.stop();

  std::cout<<"Calc time: "<<timer_calc.accumulated()<<std::endl;
}