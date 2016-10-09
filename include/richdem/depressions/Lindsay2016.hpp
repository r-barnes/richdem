#ifndef _richdem_lindsay2016_hpp_
#define _richdem_lindsay2016_hpp_

#include "richdem/common/Array2D.hpp"
#include "richdem/common/grid_cell.hpp"
#include <limits>

enum LindsayMode {
  COMPLETE_BREACHING,
  SELECTIVE_BREACHING,
  CONSTRAINED_BREACHING
};

enum LindsayCellType {
  UNVISITED,
  VISITED,
  EDGE
};

template<class T>
void Lindsay2016(
  Array2D<T>  &dem,
  LindsayMode mode,
  uint32_t    maxpathlen,
  T           maxdepth
){
  Array2D<uint32_t>     backlinks(dem);
  Array2D<uint8_t>      visited(dem);
  Array2D<uint8_t>      pits(dem);
  std::vector<uint32_t> flood_array;
  GridCellZk_pq<T>      pq;

  dem.setNoData(-9999);

  uint32_t total_pits = 0;

  visited.setAll(LindsayCellType::UNVISITED);

  //Seed the priority queue
  for(int y=0;y<dem.height();y++)
  for(int x=0;x<dem.width();x++){
    if(dem.isNoData(x,y))             //Don't evaluate NoData cells
      continue;

    if(dem.isEdgeCell(x,y)){          //Valid edge cells go on priority-queue
      pq.emplace(x,y,dem(x,y));
      visited(x,y) = LindsayCellType::EDGE;
      continue;
    }

    //For identifying the lowest neighbour of a pit cell
    T lowest_neighbour = std::numeric_limits<T>::max();    
    for(int n=1;n<=8;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];
      if(!dem.inGrid(nx,ny))
        continue;
      if(dem.isNoData(nx,ny)){        //Cells which can drain into NoData go on priority-queue
        pq.emplace(x,y, dem(x,y));
        visited(x,y) = LindsayCellType::EDGE;
        goto nextcell; //VELOCIRAPTOR
      }
      lowest_neighbour = std::min(dem(nx,ny),lowest_neighbour);
    }

    //Was this a pit cell? If so: elevate it to a level just lower than its
    //lowest neighbour. This makes the breaching/tunneling procedures work
    //better.
    if(dem(x,y)<lowest_neighbour)
      dem(x,y) = std::nextafter(lowest_neighbour, std::numeric_limits<T>::lowest());

    if(dem(x,y)<=lowest_neighbour){
      pits(x,y) = true;
      total_pits++; //TODO: May not need this
    }

    nextcell:;
  }

  while(!pq.empty()){
    const auto c = pq.top();
    pq.pop();

    for(int n=1;n<=8;n++){
      const int nx = c.x+dx[n];
      const int ny = c.y+dy[n];

      if(!dem.inGrid(nx,ny))
        continue;
      if(dem.isNoData(nx,ny))
        continue;
      if(visited(nx,ny))
        continue;

      const auto my_e = dem(nx,ny);

      pq.emplace(nx,ny,my_e);
      flood_array.emplace_back(dem.xyToI(nx,ny));
      visited(nx,ny)   = LindsayCellType::VISITED;
      backlinks(nx,ny) = dem.xyToI(c.x,c.y);

      if(!pits(nx,ny))
        continue;

      //Locate a cell that is lower than the pit cell, or an edge cell
      auto     cc                      = dem.xyToI(nx,ny);
      uint32_t pathlen                 = 1;
      T        maxheight               = dem(cc);
      typename Array2D<T>::i_t pathend = cc;

      while(dem(cc)>=my_e && visited(cc)!=LindsayCellType::EDGE){
        cc = backlinks(cc);
        pathlen++;
        maxheight = std::max(maxheight, dem(cc));
        pathend   = cc;
      }
      auto end_e = dem(pathend);

      cc = dem.xyToI(nx,ny);

      if(mode==COMPLETE_BREACHING || (pathlen<=maxpathlen && maxheight-end_e<=maxdepth)){
        while(cc!=pathend){ //TODO: Use nextafter to enforce drainage along the path. Probs not necessary bc of flood_array
          dem(cc) = end_e;
          cc      = backlinks(cc);
        }
      } else if(mode==CONSTRAINED_BREACHING){ //TODO: Refine this with regards to the paper
        while(cc!=pathend){
          if(maxheight-end_e<=maxdepth)
            dem(cc) = end_e;
          else
            dem(cc) -= maxdepth;
          cc = backlinks(cc);
        }
      }
    }
  }

  for(const auto f: flood_array){
    auto parent = backlinks(f);
    if(dem(f)<=dem(parent))
      dem(f)=std::nextafter(dem(parent),std::numeric_limits<T>::max());
  }
}

#endif