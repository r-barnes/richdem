/**
  @file
  @brief Defines all the Priority-Flood algorithms described by Barnes (2014) "Priority-Flood: An Optimal Depression-Filling and Watershed-Labeling Algorithm for Digital Elevation Models".

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_priority_flood_hpp_
#define _richdem_priority_flood_hpp_

#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/grid_cell.hpp"
#include "richdem/flowmet/d8_flowdirs.hpp"
#include <queue>
#include <limits>
#include <iostream>
#include <cstdlib> //Used for exit

namespace richdem {

/**
  @brief  Determine if a DEM has depressions
  @author Richard Barnes (rbarnes@umn.edu)

    Priority-Flood starts on the edges of the DEM and then works its way
    inwards using a priority queue to determine the lowest cell which has
    a path to the edge. The neighbours of this cell are added to the priority
    queue. If the neighbours are lower than the cell which is adding them, then
    they are part of a depression and the question is answered.

  @param[in]  &elevations   A grid of cell elevations

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM. Note that the _NoData_ value is assumed to
       be a negative number less than any actual data value.

  @return True if the DEM contains depressions; otherwise, false.

  @correctness
    The correctness of this command is determined by inspection. (TODO)
*/
template <Topology topo, class elev_t>
bool HasDepressions(const Array2D<elev_t> &elevations){
  GridCellZ_pq<elev_t> open;
  ProgressBar progress;

  RDLOG_ALG_NAME<<"HasDepressions (Based on Priority-Flood)";
  RDLOG_CITATION<<"Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024";
  RDLOG_CONFIG  <<"topology = "<<TopologyName(topo);

  const int *const dx   = topo == Topology::D8?d8x:topo==Topology::D4?d4x:NULL;
  const int *const dy   = topo == Topology::D8?d8y:topo==Topology::D4?d4y:NULL;
  const int        nmax = topo == Topology::D8?  8:topo==Topology::D4?  4:   0;

  RDLOG_PROGRESS<<"Setting up boolean flood array matrix...";
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);

  RDLOG_MEM_USE<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
               <<"MB of RAM.";

  RDLOG_PROGRESS<<"Adding perimeter cells to the priority queue...";
  for(int x=0;x<elevations.width();x++){
    open.emplace(x,0,elevations(x,0) );
    open.emplace(x,elevations.height()-1,elevations(x,elevations.height()-1) );
    closed(x,0)=true;
    closed(x,elevations.height()-1)=true;
  }
  for(int y=1;y<elevations.height()-1;y++){
    open.emplace(0,y,elevations(0,y)  );
    open.emplace(elevations.width()-1,y,elevations(elevations.width()-1,y) );
    closed(0,y)=true;
    closed(elevations.width()-1,y)=true;
  }

  RDLOG_PROGRESS<<"Searching for depressions...";
  progress.start( elevations.size() );
  while(open.size()>0){
    GridCellZ<elev_t> c=open.top();
    open.pop();
    ++progress;

    for(int n=1;n<=nmax;n++){
      int nx = c.x+dx[n];
      int ny = c.y+dy[n];
      if(!elevations.inGrid(nx,ny)) continue;
      if(closed(nx,ny))
        continue;

      closed(nx,ny) = true;
      if(elevations(nx,ny)<elevations(c.x,c.y)){
        RDLOG_TIME_USE<<"Succeeded in    = "<<progress.stop() <<" s";
        RDLOG_MISC<<"Depression found.";
        return true;
      }
      open.emplace(nx,ny,elevations(nx,ny));
    }
  }
  RDLOG_TIME_USE<<"t Succeeded in    = "<<progress.stop() <<" s";
  RDLOG_MISC<<"No depressions found.";
  return false;
}






/**
  @brief  Fills all pits and removes all digital dams from a DEM
  @author Richard Barnes (rbarnes@umn.edu)

    Priority-Flood starts on the edges of the DEM and then works its way
    inwards using a priority queue to determine the lowest cell which has
    a path to the edge. The neighbours of this cell are added to the priority
    queue. If the neighbours are lower than the cell which is adding them, then
    they are raised to match its elevation; this fills depressions.

  @param[in,out]  &elevations   A grid of cell elevations

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM. Note that the _NoData_ value is assumed to
       be a negative number less than any actual data value.

  @post
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM.
    2. **elevations** contains no landscape depressions or digital dams.

  @correctness
    The correctness of this command is determined by inspection. (TODO)
*/
template <Topology topo, class elev_t>
void PriorityFlood_Original(Array2D<elev_t> &elevations){
  GridCellZ_pq<elev_t> open;
  uint64_t processed_cells = 0;
  uint64_t pitc            = 0;
  ProgressBar progress;

  RDLOG_ALG_NAME<<"A Priority-Flood (Original)";
  RDLOG_CITATION<<"Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024";
  RDLOG_PROGRESS<<"Setting up boolean flood array matrix...";
  RDLOG_CONFIG  <<"topology = "<<TopologyName(topo);

  const int *const dx   = topo == Topology::D8?d8x:topo==Topology::D4?d4x:NULL;
  const int *const dy   = topo == Topology::D8?d8y:topo==Topology::D4?d4y:NULL;
  const int        nmax = topo == Topology::D8?  8:topo==Topology::D4?  4:   0;

  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);

  RDLOG_MEM_USE<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
               <<"MB of RAM.";

  RDLOG_PROGRESS<<"Adding perimeter cells to the priority queue...";
  for(int x=0;x<elevations.width();x++){
    open.emplace(x,0,elevations(x,0) );
    open.emplace(x,elevations.height()-1,elevations(x,elevations.height()-1) );
    closed(x,0)=true;
    closed(x,elevations.height()-1)=true;
  }
  for(int y=1;y<elevations.height()-1;y++){
    open.emplace(0,y,elevations(0,y)  );
    open.emplace(elevations.width()-1,y,elevations(elevations.width()-1,y) );
    closed(0,y)=true;
    closed(elevations.width()-1,y)=true;
  }

  RDLOG_PROGRESS<<"Performing the original Priority Flood...";
  progress.start( elevations.size() );
  while(open.size()>0){
    GridCellZ<elev_t> c=open.top();
    open.pop();
    processed_cells++;

    for(int n=1;n<=nmax;n++){
      int nx = c.x+dx[n];
      int ny = c.y+dy[n];
      if(!elevations.inGrid(nx,ny)) continue;
      if(closed(nx,ny))
        continue;

      closed(nx,ny) = true;
      if(elevations(nx,ny)<elevations(c.x,c.y))
        ++pitc;
      elevations(nx,ny) = std::max(elevations(nx,ny),elevations(c.x,c.y));
      open.emplace(nx,ny,elevations(nx,ny));
    }
    progress.update(processed_cells);
  }
  RDLOG_TIME_USE<<"Succeeded in    = "<<progress.stop() <<" s";
  RDLOG_MISC    <<"Cells processed = "<<processed_cells;
  RDLOG_MISC    <<"Cells in pits   = "<<pitc;
}



/**
  @brief  Fills all pits and removes all digital dams from a DEM, but faster
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

  @correctness
    The correctness of this command is determined by inspection. (TODO)
*/
template <Topology topo, class elev_t>
void PriorityFlood_Barnes2014(Array2D<elev_t> &elevations){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  uint64_t processed_cells = 0;
  uint64_t pitc            = 0;
  ProgressBar progress;

  RDLOG_ALG_NAME << "Priority-Flood (Improved)";
  RDLOG_CITATION << "Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024";
  RDLOG_CONFIG   <<"topology = "<<TopologyName(topo);

  const int *const dx   = topo == Topology::D8?d8x:topo==Topology::D4?d4x:NULL;
  const int *const dy   = topo == Topology::D8?d8y:topo==Topology::D4?d4y:NULL;
  const int        nmax = topo == Topology::D8?  8:topo==Topology::D4?  4:   0;

  RDLOG_PROGRESS << "Setting up boolean flood array matrix...";
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);

  RDLOG_MEM_USE<<"Priority queue requires approx = "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
               <<"MB of RAM.";

  RDLOG_PROGRESS<<"Adding cells to the priority queue...";

  for(int x=0;x<elevations.width();x++){
    open.emplace(x,0,elevations(x,0) );
    open.emplace(x,elevations.height()-1,elevations(x,elevations.height()-1) );
    closed(x,0)=true;
    closed(x,elevations.height()-1)=true;
  }
  for(int y=1;y<elevations.height()-1;y++){
    open.emplace(0,y,elevations(0,y)  );
    open.emplace(elevations.width()-1,y,elevations(elevations.width()-1,y) );
    closed(0,y)=true;
    closed(elevations.width()-1,y)=true;
  }

  RDLOG_PROGRESS<<"Performing the improved Priority-Flood...";
  progress.start( elevations.size() );
  while(open.size()>0 || pit.size()>0){
    GridCellZ<elev_t> c;
    if(pit.size()>0){
      c=pit.front();
      pit.pop();
    } else {
      c=open.top();
      open.pop();
    }
    processed_cells++;

    for(int n=1;n<=nmax;n++){
      int nx=c.x+dx[n];
      int ny=c.y+dy[n];
      if(!elevations.inGrid(nx,ny)) continue;
      if(closed(nx,ny))
        continue;

      closed(nx,ny)=true;
      if(elevations(nx,ny)<=c.z){
        if(elevations(nx,ny)<c.z){
          ++pitc;
          elevations(nx,ny)=c.z;
        }
        pit.push(GridCellZ<elev_t>(nx,ny,c.z));
      } else
        open.emplace(nx,ny,elevations(nx,ny));
    }
    progress.update(processed_cells);
  }
  RDLOG_TIME_USE<<"Succeeded in "<<std::fixed<<std::setprecision(1)<<progress.stop()<<" s";
  RDLOG_MISC    <<"Cells processed = "<<processed_cells;
  RDLOG_MISC    <<"Cells in pits = "  <<pitc;
}


/**
  @brief  Modifies floating-point cell elevations to guarantee drainage.
  @author Richard Barnes (rbarnes@umn.edu)

    This version of Priority-Flood starts on the edges of the DEM and then
    works its way inwards using a priority queue to determine the lowest cell
    which has a path to the edge. The neighbours of this cell are added to the
    priority queue if they are higher. If they are lower, then their elevation
    is increased by a small amount to ensure that they have a drainage path and
    they are added to a "pit" queue which is used to flood pits. Cells which
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
    2. **elevations** has no landscape depressions, digital dams, or flats.

  @correctness
    The correctness of this command is determined by inspection. (TODO)
*/
template <Topology topo, class elev_t>
void PriorityFloodEpsilon_Barnes2014(Array2D<elev_t> &elevations){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  ProgressBar progress;
  uint64_t processed_cells = 0;
  uint64_t pitc            = 0;
  auto     PitTop          = elevations.noData();
  int      false_pit_cells = 0;

  RDLOG_ALG_NAME<<"Priority-Flood+Epsilon";
  RDLOG_CITATION<<"Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024";
  RDLOG_CONFIG  <<"topology = "<<TopologyName(topo);

  const int *const dx   = topo == Topology::D8?d8x:topo==Topology::D4?d4x:NULL;
  const int *const dy   = topo == Topology::D8?d8y:topo==Topology::D4?d4y:NULL;
  const int        nmax = topo == Topology::D8?  8:topo==Topology::D4?  4:   0;

  RDLOG_PROGRESS<<"Setting up boolean flood array matrix...";
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);

  RDLOG_PROGRESS<<"Adding cells to the priority queue...";
  for(int x=0;x<elevations.width();x++){
    open.emplace(x,0,elevations(x,0) );
    open.emplace(x,elevations.height()-1,elevations(x,elevations.height()-1) );
    closed(x,0)=true;
    closed(x,elevations.height()-1)=true;
  }
  for(int y=1;y<elevations.height()-1;y++){
    open.emplace(0,y,elevations(0,y)  );
    open.emplace(elevations.width()-1,y,elevations(elevations.width()-1,y) );
    closed(0,y)=true;
    closed(elevations.width()-1,y)=true;
  }

  RDLOG_PROGRESS<<"Performing Priority-Flood+Epsilon...";
  progress.start( elevations.size() );
  while(open.size()>0 || pit.size()>0){
    GridCellZ<elev_t> c;
    if(pit.size()>0 && open.size()>0 && open.top().z==pit.front().z){
      c=open.top();
      open.pop();
      PitTop=elevations.noData();
    } else if(pit.size()>0){
      c=pit.front();
      pit.pop();
      if(PitTop==elevations.noData())
        PitTop=elevations(c.x,c.y);
    } else {
      c=open.top();
      open.pop();
      PitTop=elevations.noData();
    }
    processed_cells++;

    for(int n=1;n<=nmax;n++){
      int nx=c.x+dx[n];
      int ny=c.y+dy[n];

      if(!elevations.inGrid(nx,ny)) continue;

      if(closed(nx,ny))
        continue;
      closed(nx,ny)=true;

      if(elevations(nx,ny)==elevations.noData())
        pit.push(GridCellZ<elev_t>(nx,ny,elevations.noData()));

      else if(elevations(nx,ny)<=std::nextafter(c.z,std::numeric_limits<elev_t>::infinity())){
        if(PitTop!=elevations.noData() && PitTop<elevations(nx,ny) && std::nextafter(c.z,std::numeric_limits<elev_t>::infinity())>=elevations(nx,ny))
          ++false_pit_cells;
        ++pitc;
        elevations(nx,ny)=std::nextafter(c.z,std::numeric_limits<elev_t>::infinity());
        pit.emplace(nx,ny,elevations(nx,ny));
      } else
        open.emplace(nx,ny,elevations(nx,ny));
    }
    progress.update(processed_cells);
  }
  RDLOG_TIME_USE<<"succeeded in "<<progress.stop()<<" s";
  RDLOG_MISC<<"Cells processed = "<<processed_cells;
  RDLOG_MISC<<"Cells in pits = "  <<pitc           ;
  if(false_pit_cells)
    RDLOG_WARN<<"\033[91mW In assigning negligible gradients to depressions, some depressions rose above the surrounding cells. This implies that a larger storage type should be used. The problem occured for "<<false_pit_cells<<" of "<<elevations.numDataCells()<<".\033[39m";
}


///Priority-Flood+Epsilon is not available for integer data types
template<Topology topo>
void PriorityFloodEpsilon_Barnes2014(Array2D<uint8_t> &elevations){
  throw std::runtime_error("Priority-Flood+Epsilon is only available for floating-point data types!");
}

///Priority-Flood+Epsilon is not available for integer data types
template<Topology topo>
void PriorityFloodEpsilon_Barnes2014(Array2D<uint16_t> &elevations){
  throw std::runtime_error("Priority-Flood+Epsilon is only available for floating-point data types!");
}

///Priority-Flood+Epsilon is not available for integer data types
template<Topology topo>
void PriorityFloodEpsilon_Barnes2014(Array2D<int16_t> &elevations){
  throw std::runtime_error("Priority-Flood+Epsilon is only available for floating-point data types!");
}

///Priority-Flood+Epsilon is not available for integer data types
template<Topology topo>
void PriorityFloodEpsilon_Barnes2014(Array2D<uint32_t> &elevations){
  throw std::runtime_error("Priority-Flood+Epsilon is only available for floating-point data types!");
}

///Priority-Flood+Epsilon is not available for integer data types
template<Topology topo>
void PriorityFloodEpsilon_Barnes2014(Array2D<int32_t> &elevations){
  throw std::runtime_error("Priority-Flood+Epsilon is only available for floating-point data types!");
}


/**
  @brief  Determines D8 flow directions and implicitly fills pits.
  @author Richard Barnes (rbarnes@umn.edu)

    This version of Priority-Flood starts on the edges of the DEM and then
    works its way inwards using a priority queue to determine the lowest cell
    which has a path to the edge. The neighbours of this cell are given D8 flow
    directions which point to it. All depressions are implicitly filled and
    digital dams removed.

    Based on Metz 2011.

  @param[in]   &elevations  A grid of cell elevations
  @param[out]  &flowdirs    A grid of D8 flow directions

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM. Note that the _NoData_ value is assumed to
       be a negative number less than any actual data value.

  @post
    1. **flowdirs** contains a D8 flow direction of each cell or a value
       _NO_FLOW_ for those cells which are not part of the DEM.
    2. **flowdirs** has no cells which are not part of a continuous flow
       path leading to the edge of the DEM.

  @correctness
    The correctness of this command is determined by inspection. (TODO)
*/
template <class elev_t>
void PriorityFloodFlowdirs_Barnes2014(const Array2D<elev_t> &elevations, Array2D<d8_flowdir_t> &flowdirs){
  GridCellZk_pq<elev_t> open;
  uint64_t processed_cells = 0;
  ProgressBar progress;

  RDLOG_ALG_NAME<<"Priority-Flood+Flow Directions";
  RDLOG_CITATION<<"Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024";
  RDLOG_PROGRESS<<"Setting up boolean flood array matrix...";
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);

  RDLOG_PROGRESS<<"Setting up the flowdirs matrix...";
  flowdirs.resize(elevations.width(),elevations.height());
  flowdirs.setNoData(NO_FLOW);

  RDLOG_PROGRESS<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
                <<"MB of RAM.";

  RDLOG_PROGRESS<<"Adding cells to the priority queue...";
  for(int x=0;x<elevations.width();x++){
    open.emplace(x,0,elevations(x,0));
    open.emplace(x,elevations.height()-1,elevations(x,elevations.height()-1));
    flowdirs(x,0)=3;
    flowdirs(x,elevations.height()-1)=7;
    closed(x,0)=true;
    closed(x,elevations.height()-1)=true;
  }
  for(int y=1;y<elevations.height()-1;y++){
    open.emplace(0,y,elevations(0,y) );
    open.emplace(elevations.width()-1,y,elevations(elevations.width()-1,y) );
    flowdirs(0,y)=1;
    flowdirs(elevations.width()-1,y)=5;
    closed(0,y)=true;
    closed(elevations.width()-1,y)=true;
  }

  flowdirs(0,0)=2;
  flowdirs(flowdirs.width()-1,0)=4;
  flowdirs(0,flowdirs.height()-1)=8;
  flowdirs(flowdirs.width()-1,flowdirs.height()-1)=6;

  const int d8_order[9]={0,1,3,5,7,2,4,6,8};
  RDLOG_PROGRESS<<"Performing Priority-Flood+Flow Directions...";
  progress.start( elevations.size() );
  while(open.size()>0){
    auto c=open.top();
    open.pop();
    processed_cells++;

    for(int no=1;no<=8;no++){
      int n=d8_order[no];
      int nx=c.x+dx[n];
      int ny=c.y+dy[n];
      if(!elevations.inGrid(nx,ny))
        continue;
      if(closed(nx,ny))
        continue;

      closed(nx,ny)=true;

      if(elevations(nx,ny)==elevations.noData())
        flowdirs(nx,ny)=flowdirs.noData();
      else
        flowdirs(nx,ny)=d8_inverse[n];

      open.emplace(nx,ny,elevations(nx,ny));
    }
    progress.update(processed_cells);
  }
  RDLOG_TIME_USE<<"succeeded in "<<progress.stop()<<" s";
  RDLOG_MISC<<"Cells processed = "<<processed_cells;
}











/**
  @brief  Indicates which cells are in pits
  @author Richard Barnes (rbarnes@umn.edu)

    This version of Priority-Flood starts on the edges of the DEM and then
    works its way inwards using a priority queue to determine the lowest cell
    which has a path to the edge. If a cell is lower than this cell then it is
    part of a pit and is given a value 1 to indicate this. The result is a grid
    where every cell which is in a pit is labeled.

  @param[in]   &elevations   A grid of cell elevations
  @param[out]  &pit_mask     A grid of indicating which cells are in pits

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM. Note that the _NoData_ value is assumed to
       be a negative number less than any actual data value.

  @post
    1. **pit_mask** contains a 1 for each cell which is in a pit and a 0 for
       each cell which is not. The value 3 indicates NoData

  @correctness
    The correctness of this command is determined by inspection. (TODO)
*/
//TODO: Can I use a smaller data type?
template <Topology topo, class elev_t>
void pit_mask(const Array2D<elev_t> &elevations, Array2D<uint8_t> &pit_mask){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  uint64_t processed_cells = 0;
  uint64_t pitc            = 0;
  ProgressBar progress;

  RDLOG_ALG_NAME<<"Pit Mask";
  RDLOG_CITATION<<"Barnes, R. 2016. RichDEM: Terrain Analysis Software. http://github.com/r-barnes/richdem"; //TODO
  RDLOG_CONFIG  <<"topology = "<<TopologyName(topo);  

  const int *const dx   = topo == Topology::D8?d8x:topo==Topology::D4?d4x:NULL;
  const int *const dy   = topo == Topology::D8?d8y:topo==Topology::D4?d4y:NULL;
  const int        nmax = topo == Topology::D8?  8:topo==Topology::D4?  4:   0;

  RDLOG_PROGRESS<<"Setting up boolean flood array matrix...";
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);

  RDLOG_PROGRESS<<"Setting up the pit mask matrix...";
  pit_mask.resize(elevations.width(),elevations.height());
  pit_mask.setNoData(3);

  RDLOG_MEM_USE<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
               <<"MB of RAM.";

  RDLOG_PROGRESS<<"Adding cells to the priority queue...";
  for(int x=0;x<elevations.width();x++){
    open.emplace(x,0,elevations(x,0) );
    open.emplace(x,elevations.height()-1,elevations(x,elevations.height()-1) );
    closed(x,0)=true;
    closed(x,elevations.height()-1)=true;
  }
  for(int y=1;y<elevations.height()-1;y++){
    open.emplace(0,y,elevations(0,y)  );
    open.emplace(elevations.width()-1,y,elevations(elevations.width()-1,y) );
    closed(0,y)=true;
    closed(elevations.width()-1,y)=true;
  }

  RDLOG_PROGRESS<<"Performing the pit mask...";
  progress.start( elevations.size() );
  while(open.size()>0 || pit.size()>0){
    GridCellZ<elev_t> c;
    if(pit.size()>0){
      c=pit.front();
      pit.pop();
    } else {
      c=open.top();
      open.pop();
    }
    processed_cells++;

    for(int n=1;n<=nmax;n++){
      int nx = c.x+dx[n];
      int ny = c.y+dy[n];
      if(!elevations.inGrid(nx,ny)) continue;
      if(closed(nx,ny))
        continue;

      closed(nx,ny)=true;
      if(elevations(nx,ny)<=c.z){
        if(elevations(nx,ny)<c.z){
          pitc++;
          pit_mask(nx,ny)=1;
        }
        pit.emplace(nx,ny,c.z);
      } else{
        pit_mask(nx,ny)=0;
        open.emplace(nx,ny,elevations(nx,ny));
      }
    }

    if(elevations.isNoData(c.x,c.y))
      pit_mask(c.x,c.y) = pit_mask.noData();

    progress.update(processed_cells);
  }
  RDLOG_TIME_USE<<"Succeeded in "<<progress.stop()<<" s";
  RDLOG_MISC<<"Cells processed = "<<processed_cells;
  RDLOG_MISC<<"Cells in depressions = "  <<pitc;
}



/**
  @brief  Gives a common label to all cells which drain to a common point
  @author Richard Barnes (rbarnes@umn.edu)

    All the edge cells of the DEM are given unique labels. This version of
    Priority-Flood starts on the edges of the DEM and then works its way
    inwards using a priority queue to determine the lowest cell which has a
    path to the edge. The neighbours of this cell are then given its label. All
    depressions are implicitly filled and digital dams removed. The result is
    a grid of cells where all cells with a common label drain to a common
    point.

  @param[in,out] elevations        A grid of cell elevations
  @param[out]    labels            A grid to hold the watershed labels
  @param[in]     alter_elevations
    If true, then **elevations** is altered as though PriorityFlood_Barnes2014()
    had been applied. The result is that all cells drain to the edges of the
    DEM. Otherwise, **elevations** is not altered.

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM. Note that the _NoData_ value is assumed to
       be a negative number less than any actual data value.

  @post
    1. **elevations** contains no depressions or digital dams, if
       **alter_elevations** was set.
    2. **labels** contains a label for each cell indicating its membership in a
       given watershed. Cells bearing common labels drain to common points.

  @correctness
    The correctness of this command is determined by inspection. (TODO)
*/
template <Topology topo, class elev_t>
void PriorityFloodWatersheds_Barnes2014(
  Array2D<elev_t> &elevations, Array2D<int32_t> &labels, bool alter_elevations
){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  unsigned long processed_cells=0;
  unsigned long pitc=0,openc=0;
  int clabel=1;  //TODO: Thought this was more clear than zero in the results.
  ProgressBar progress;

  RDLOG_ALG_NAME<<"Priority-Flood+Watershed Labels";
  RDLOG_CITATION<<"Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024";
  RDLOG_CONFIG  <<"topology = "<<TopologyName(topo);  

  const int *const dx   = topo == Topology::D8?d8x:topo==Topology::D4?d4x:NULL;
  const int *const dy   = topo == Topology::D8?d8y:topo==Topology::D4?d4y:NULL;
  const int        nmax = topo == Topology::D8?  8:topo==Topology::D4?  4:   0;

  RDLOG_PROGRESS<<"Setting up boolean flood array matrix...";
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);

  RDLOG_PROGRESS<<"Setting up watershed label matrix...";
  labels.resize(elevations.width(),elevations.height(),-1);
  labels.setNoData(-1);

  RDLOG_MEM_USE<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
               <<"MB of RAM.";

  RDLOG_PROGRESS<<"Adding cells to the priority queue...";
  for(int x=0;x<elevations.width();x++){
    open.emplace(x,0,elevations(x,0) );
    open.emplace(x,elevations.height()-1,elevations(x,elevations.height()-1) );
    closed(x,0)=true;
    closed(x,elevations.height()-1)=true;
  }
  for(int y=1;y<elevations.height()-1;y++){
    open.emplace(0,y,elevations(0,y)  );
    open.emplace(elevations.width()-1,y,elevations(elevations.width()-1,y) );
    closed(0,y)=true;
    closed(elevations.width()-1,y)=true;
  }

  RDLOG_PROGRESS<<"Performing Priority-Flood+Watershed Labels...";
  progress.start( elevations.size() );
  while(open.size()>0 || pit.size()>0){
    GridCellZ<elev_t> c;
    if(pit.size()>0){
      c=pit.front();
      pit.pop();
      pitc++;
    } else {
      c=open.top();
      open.pop();
      openc++;
    }
    processed_cells++;

    //Since all interior cells will be flowing into a cell which has already
    //been processed, the following line identifies only the edge cells of the
    //DEM. Each edge cell seeds its own watershed/basin. The result of this will
    //be many small watersheds/basins around the edge of the DEM.
    if(labels(c.x,c.y)==labels.noData() && elevations(c.x,c.y)!=elevations.noData())  //Implies a cell without a label which borders the edge of the DEM or a region of no_data
      labels(c.x,c.y)=clabel++;

    for(int n=1;n<=nmax;n++){
      int nx=c.x+dx[n];
      int ny=c.y+dy[n];
      if(!elevations.inGrid(nx,ny)) continue;
      if(closed(nx,ny))
        continue;

      //Since the neighbouring cell is not closed, its flow is directed to this
      //cell. Therefore, it is part of the same watershed/basin as this cell.
      labels(nx,ny)=labels(c.x,c.y);

      closed(nx,ny)=true;
      if(elevations(nx,ny)<=c.z){
        if(alter_elevations)
          elevations(nx,ny)=c.z;
        pit.push(GridCellZ<elev_t>(nx,ny,c.z));
      } else
        open.push(GridCellZ<elev_t>(nx,ny,elevations(nx,ny)));
    }
    progress.update(processed_cells);
  }

  RDLOG_TIME_USE<<"succeeded in "<<progress.stop()<<" s";

  RDLOG_MISC<<"Cells processed   = "  <<processed_cells;
  RDLOG_MISC<<"Cells in pits     = "  <<pitc;
  RDLOG_MISC<<"Cells not in pits = "  <<openc;
}



/**
  @brief  Fill depressions, but only if they're small
  @author Richard Barnes (rbarnes@umn.edu)

    Priority-Flood starts on the edges of the DEM and then works its way
    inwards using a priority queue to determine the lowest cell which has
    a path to the edge. The neighbours of this cell are added to the priority
    queue if they are higher. If they are lower, they are raised to the
    elevation of the cell adding them, thereby filling in pits. The neighbors
    are then added to a "pit" queue which is used to flood pits. Cells which
    are higher than a pit being filled are added to the priority queue. In this
    way, pits are filled without incurring the expense of the priority queue.

    When a depression is encountered this command measures its size before 
    filling it. Only small depressions are filled.

  @param[in,out]  &elevations   A grid of cell elevations
  @param[in]      max_dep_size  Depression must have <=max_dep_size cells to be
                                filled

  @pre
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM. Note that the _NoData_ value is assumed to
       be a negative number less than any actual data value.

  @post
    1. **elevations** contains the elevations of every cell or a value _NoData_
       for cells not part of the DEM.
    2. **elevations** all landscape depressions <=max_dep_size are filled.

  @correctness
    The correctness of this command is determined by inspection. (TODO)
*/
template <Topology topo, class elev_t>
void PriorityFlood_Barnes2014_max_dep(
  Array2D<elev_t> &elevations,
  uint64_t max_dep_size
){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  uint64_t processed_cells = 0;
  uint64_t pitc            = 0;
  ProgressBar progress;

  RDLOG_ALG_NAME<<"\nPriority-Flood (Improved) Limited to Maximum Depression Area";
  RDLOG_CITATION<<"Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024";
  RDLOG_CONFIG  <<"topology = "<<TopologyName(topo);

  const int *const dx   = topo == Topology::D8?d8x:topo==Topology::D4?d4x:NULL;
  const int *const dy   = topo == Topology::D8?d8y:topo==Topology::D4?d4y:NULL;
  const int        nmax = topo == Topology::D8?  8:topo==Topology::D4?  4:   0;

  RDLOG_PROGRESS<<"Setting up boolean flood array matrix...";
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);

  RDLOG_MEM_USE<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
               <<"MB of RAM.";

  RDLOG_PROGRESS<<"Adding cells to the priority queue...";

  for(int x=0;x<elevations.width();x++){
    open.emplace(x,0,elevations(x,0) );
    open.emplace(x,elevations.height()-1,elevations(x,elevations.height()-1) );
    closed(x,0)=true;
    closed(x,elevations.height()-1)=true;
  }
  for(int y=1;y<elevations.height()-1;y++){
    open.emplace(0,y,elevations(0,y)  );
    open.emplace(elevations.width()-1,y,elevations(elevations.width()-1,y) );
    closed(0,y)=true;
    closed(elevations.width()-1,y)=true;
  }

  elev_t dep_elev = 0;             //Elevation of the rim/spill point of the depression we're in
  std::vector<GridCell> dep_cells; //Cells comprising the depression we're in

  RDLOG_PROGRESS<<"Performing the improved Priority-Flood...";
  progress.start( elevations.size() );
  while(open.size()>0 || pit.size()>0){
    GridCellZ<elev_t> c;
    if(pit.size()>0){                       //There are cells inside a depression which should be processed
      c=pit.front();                        //Get next cell of the depression
      pit.pop();                      
      dep_cells.push_back(c);               //Add cell to list of cells in depression
    } else {                                //We're not in a depression
      c=open.top();                         //Get next highest cell
      open.pop();
      if(dep_cells.size()<=max_dep_size){   //Have we just crawled out of a small depression?
        for(const auto &pc: dep_cells)      //Loop through cells in depression
          elevations(pc.x,pc.y) = dep_elev; //Raise each cell to the level of depression's rim/spill point
        dep_cells.clear();                  //We're done with this depression now
      } else if(dep_cells.size()>0) {       //Have we just crawled out of a depression?
        dep_cells.clear();                  //Not interested in small depressions
      }
    }
    processed_cells++;

    for(int n=1;n<=nmax;n++){
      int nx=c.x+dx[n];
      int ny=c.y+dy[n];
      if(!elevations.inGrid(nx,ny)) continue;
      if(closed(nx,ny))
        continue;

      closed(nx,ny)=true;
      if(elevations(nx,ny)<c.z){                //Cell <= current elev can be processed quickly with a queue
          ++pitc;                               //Count it
        pit.push(GridCellZ<elev_t>(nx,ny,c.z)); //Add this cell to depression-prcessing queue
        dep_elev = c.z;                         //Note rim/spill-point elevation
      } else{
        open.emplace(nx,ny,elevations(nx,ny));  //Not a depression, carry on
      }
    }
    progress.update(processed_cells);
  }
  RDLOG_TIME_USE<<"Succeeded in "<<std::fixed<<std::setprecision(1)<<progress.stop()<<" s";
  RDLOG_MISC<<"Cells processed = "<<processed_cells;
  RDLOG_MISC<<"Cells in pits = "  <<pitc;
}

}

#endif
