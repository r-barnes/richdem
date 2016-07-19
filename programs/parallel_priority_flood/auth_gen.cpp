#include <iostream>
#include <string>
#include <cstdlib>
#include "Array2D.hpp"
#include "common.hpp"

//improved_priority_flood
/**
  @brief  Fills all pits and removes all digital dams from a DEM
  @author Richard Barnes (rbarnes@umn.edu)

    Priority-Flood starts on the edges of the DEM and then works its way
    inwards using a priority queue to determine the lowest cell which has
    a path to the edge. The neighbours of this cell are added to the priority
    queue if they are higher. If they are lower, they are raised to the
    elevation of the cell adding them, thereby filling in pits. The neighbors
    are then added to a "pit" queue which is used to flood pits. Cells which
    are higher than a pit being filled are added to the priority queue. In this
    way, pits are filled without incurring the expense of the priority queue.

  @param[in,out]  &elevations   A grid of cell elevations

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
void PriorityFlood(Array2D<elev_t> &dem){
  GridCellZ_pq<elev_t> pq;
  std::queue< GridCellZ<elev_t> > pit;

  Array2D<char> closed(dem.viewWidth(), dem.viewHeight(), false);

  auto PlaceCell = [&](int x, int y){
    pq.emplace(x,y,dem(x,y));
    closed(x,y) = true;
  };

  for(int x=0;x<dem.viewWidth();x++){
    PlaceCell(x,0);
    PlaceCell(x,dem.viewHeight()-1);
  }

  for(int y=1;y<dem.viewHeight()-1;y++){
    PlaceCell(0,y);
    PlaceCell(dem.viewWidth()-1,y);
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

    for(int n=1;n<=8;n++){
      auto nx = c.x+dx[n]; //Neighbour's x-coordinate
      auto ny = c.y+dy[n]; //Neighbour's y-coordinate

      //Check to see if the neighbour coordinates are valid
      if(!dem.in_grid(nx,ny))
        continue;
      if(closed(nx,ny))
        continue;

      closed(nx,ny) = true;

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

template<class T>
int PerformAlgorithm(std::string filename, std::string output){
  Array2D<T> dem(filename,false);

  PriorityFlood(dem);

  dem.saveGDAL(output,filename,0,0);

  return 0;
}

int main(int argc, char **argv){
  if(argc!=3){
    std::cerr<<argv[0]<<" <INPUT> <OUTPUT_PREFIX>"<<std::endl;
    return -1;
  }

  switch(peekGDALType(argv[1])){
    case GDT_Unknown:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[1]))<<std::endl;
      return -1;
    case GDT_Byte:
      return PerformAlgorithm<uint8_t >(argv[1],argv[2]);
    case GDT_UInt16:
      return PerformAlgorithm<uint16_t>(argv[1],argv[2]);
    case GDT_Int16:
      return PerformAlgorithm<int16_t >(argv[1],argv[2]);
    case GDT_UInt32:
      return PerformAlgorithm<uint32_t>(argv[1],argv[2]);
    case GDT_Int32:
      return PerformAlgorithm<int32_t >(argv[1],argv[2]);
    case GDT_Float32:
      return PerformAlgorithm<float   >(argv[1],argv[2]);
    case GDT_Float64:
      return PerformAlgorithm<double  >(argv[1],argv[2]);
    case GDT_CInt16:
    case GDT_CInt32:
    case GDT_CFloat32:
    case GDT_CFloat64:
      std::cerr<<"Complex types are unsupported. Sorry!"<<std::endl;
      return -1;
    default:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[1]))<<std::endl;
      return -1;
  }

  return 0;
}