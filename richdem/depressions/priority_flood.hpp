#ifndef pit_fill_include
#define pit_fill_include
#include "../common/Array2D.hpp"
#include "../common/grid_cell.hpp"
#include "../flowdirs/d8_flowdirs.hpp"
#include <queue>
#include <limits>
#include <iostream>
#include <cstdlib> //Used for exit



//original_priority_flood
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
*/
template <class elev_t>
void original_priority_flood(Array2D<elev_t> &elevations){
  GridCellZ_pq<elev_t> open;
  unsigned long processed_cells=0;
  unsigned long pitc=0;
  ProgressBar progress;

  std::cerr<<"\n###Barnes Flood"<<std::endl;
  std::cerr<<"Setting up boolean flood array matrix..."<<std::flush;
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
           <<"MB of RAM."
           <<std::endl;
  std::cerr<<"Adding cells to the priority queue..."<<std::endl;
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
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"%%Performing the original Priority Flood..."<<std::endl;
  progress.start( elevations.width()*elevations.height() );
  while(open.size()>0){
    GridCellZ<elev_t> c=open.top();
    open.pop();
    processed_cells++;

    for(int n=1;n<=8;n++){
      int nx=c.x+dx[n];
      int ny=c.y+dy[n];
      if(!elevations.inGrid(nx,ny)) continue;
      if(closed(nx,ny))
        continue;

      closed(nx,ny)=true;
      if(elevations(nx,ny)<elevations(c.x,c.y)) ++pitc;
      elevations(nx,ny)=std::max(elevations(nx,ny),elevations(c.x,c.y));
      open.emplace(nx,ny,elevations(nx,ny));
    }
    progress.update(processed_cells);
  }
  std::cerr<<"\t\033[96msucceeded in "<<progress.stop()<<"s.\033[39m"<<std::endl;
  std::cerr<<processed_cells<<" cells processed. "<<pitc<<" in pits."<<std::endl;
}



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
template <class elev_t>
void improved_priority_flood(Array2D<elev_t> &elevations){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  unsigned long processed_cells=0;
  unsigned long pitc=0;
  ProgressBar progress;

  std::cerr<<"\n###Improved Priority-Flood"<<std::endl;
  std::cerr<<"Citation: Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024"<<std::endl;
  std::cerr<<"Setting up boolean flood array matrix..."<<std::flush;
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
           <<"MB of RAM."
           <<std::endl;
  std::cerr<<"Adding cells to the priority queue..."<<std::flush;
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
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"%%Performing the improved Priority-Flood..."<<std::endl;
  progress.start( elevations.width()*elevations.height() );
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

    for(int n=1;n<=8;n++){
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
  std::cerr<<"\t\033[96msucceeded in "<<std::fixed<<std::setprecision(1)<<progress.stop()<<"s.\033[39m"<<std::endl;
  std::cerr<<processed_cells<<" cells processed. "<<pitc<<" in pits."<<std::endl;
}


//priority_flood_epsilon
/**
  @brief  Assigns every cell an elevation which guarantees drainage
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
*/
template <class elev_t>
void priority_flood_epsilon(Array2D<elev_t> &elevations){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  ProgressBar progress;
  unsigned long processed_cells = 0;
  unsigned long pitc            = 0;
  auto PitTop                   = elevations.noData();
  int false_pit_cells           = 0;

  std::cerr<<"\n###Priority-Flood+Epsilon"<<std::endl;
  std::cerr<<"Citation: Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024"<<std::endl;
  std::cerr<<"Setting up boolean flood array matrix..."<<std::flush;
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
           <<"MB of RAM."
           <<std::endl;
  std::cerr<<"Adding cells to the priority queue..."<<std::flush;
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
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"%%Performing Priority-Flood+Epsilon..."<<std::endl;
  progress.start( elevations.width()*elevations.height() );
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

    for(int n=1;n<=8;n++){
      int nx=c.x+dx[n];
      int ny=c.y+dy[n];

      if(!elevations.inGrid(nx,ny)) continue;

      if(closed(nx,ny))
        continue;
      closed(nx,ny)=true;

      if(elevations(nx,ny)==elevations.noData())
        pit.push(GridCellZ<elev_t>(nx,ny,elevations.noData()));

      else if(elevations(nx,ny)<=nextafterf(c.z,std::numeric_limits<float>::infinity())){
        if(PitTop!=elevations.noData() && PitTop<elevations(nx,ny) && nextafterf(c.z,std::numeric_limits<float>::infinity())>=elevations(nx,ny))
          ++false_pit_cells;
        ++pitc;
        elevations(nx,ny)=nextafterf(c.z,std::numeric_limits<float>::infinity());
        pit.push(GridCellZ<elev_t>(nx,ny,elevations(nx,ny)));
      } else
        open.emplace(nx,ny,elevations(nx,ny));
    }
    progress.update(processed_cells);
  }
  std::cerr<<"\t\033[96msucceeded in "<<progress.stop()<<"s.\033[39m"<<std::endl;
  std::cerr<<processed_cells<<" cells processed. "<<pitc<<" in pits."<<std::endl;
  if(false_pit_cells)
    std::cerr<<"\033[91mIn assigning negligible gradients to depressions, some depressions rose above the surrounding cells. This implies that a larger storage type should be used. The problem occured for "<<false_pit_cells<<" of "<<elevations.numDataCells()<<std::endl;
}


template<>
void priority_flood_epsilon(Array2D<uint8_t> &elevations){
  std::cerr<<"Priority-Flood+Epsilon is only available for floating-point data types!"<<std::endl;
  exit(-1);
}

template<>
void priority_flood_epsilon(Array2D<uint16_t> &elevations){
  std::cerr<<"Priority-Flood+Epsilon is only available for floating-point data types!"<<std::endl;
  exit(-1);
}

template<>
void priority_flood_epsilon(Array2D<int16_t> &elevations){
  std::cerr<<"Priority-Flood+Epsilon is only available for floating-point data types!"<<std::endl;
  exit(-1);
}

template<>
void priority_flood_epsilon(Array2D<uint32_t> &elevations){
  std::cerr<<"Priority-Flood+Epsilon is only available for floating-point data types!"<<std::endl;
  exit(-1);
}

template<>
void priority_flood_epsilon(Array2D<int32_t> &elevations){
  std::cerr<<"Priority-Flood+Epsilon is only available for floating-point data types!"<<std::endl;
  exit(-1);
}


//priority_flood_flowdirs
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
*/
template <class elev_t>
void priority_flood_flowdirs(const Array2D<elev_t> &elevations, Array2D<int8_t> &flowdirs){
  GridCellZk_pq<elev_t> open;
  unsigned long processed_cells=0;
  ProgressBar progress;

  std::cerr<<"\n###Priority-Flood+Flow Directions"<<std::endl;
  std::cerr<<"Citation: Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024"<<std::endl;
  std::cerr<<"Setting up boolean flood array matrix..."<<std::flush;
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"Setting up the flowdirs matrix..."<<std::flush;  
  flowdirs.resize(elevations.width(),elevations.height());
  flowdirs.setNoData(NO_FLOW);
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
           <<"MB of RAM."
           <<std::endl;

  std::cerr<<"Adding cells to the priority queue..."<<std::endl;
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
  std::cerr<<"succeeded."<<std::endl;

  flowdirs(0,0)=2;
  flowdirs(flowdirs.width()-1,0)=4;
  flowdirs(0,flowdirs.height()-1)=8;
  flowdirs(flowdirs.width()-1,flowdirs.height()-1)=6;

  const int d8_order[9]={0,1,3,5,7,2,4,6,8};
  std::cerr<<"%%Performing Priority-Flood+Flow Directions..."<<std::endl;
  progress.start( elevations.width()*elevations.height() );
  while(open.size()>0){
    GridCellZ<elev_t> c=open.top();
    open.pop();
    processed_cells++;

    for(int no=1;no<=8;no++){
      int n=d8_order[no];
      int nx=c.x+dx[n];
      int ny=c.y+dy[n];
      if(!elevations.inGrid(nx,ny)) continue;
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
  std::cerr<<"\t\033[96msucceeded in "<<progress.stop()<<"s.\033[39m"<<std::endl;
  std::cerr<<processed_cells<<" cells processed."<<std::endl;
}











//pit_mask
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
       each cell which is not.
*/
template <class elev_t>
void pit_mask(const Array2D<elev_t> &elevations, Array2D<int32_t> &pit_mask){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  unsigned long processed_cells=0;
  unsigned long pitc=0;
  ProgressBar progress;

  std::cerr<<"\n###Pit Mask"<<std::endl;
  std::cerr<<"Setting up boolean flood array matrix..."<<std::flush;
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"Setting up the pit mask matrix..."<<std::endl;
  pit_mask.resize(elevations.width(),elevations.height());
  pit_mask.setNoData(3);
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
           <<"MB of RAM."
           <<std::endl;
  std::cerr<<"Adding cells to the priority queue..."<<std::flush;
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
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"%%Performing the pit mask..."<<std::endl;
  progress.start( elevations.width()*elevations.height() );
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

    for(int n=1;n<=8;n++){
      int nx=c.x+dx[n];
      int ny=c.y+dy[n];
      if(!elevations.inGrid(nx,ny)) continue;
      if(closed(nx,ny))
        continue;

      closed(nx,ny)=true;
      if(elevations(nx,ny)<=c.z){
        if(elevations(nx,ny)<c.z)
          pit_mask(nx,ny)=1;
        pit.push(GridCellZ<elev_t>(nx,ny,c.z));
      } else{
        pit_mask(nx,ny)=0;
        open.emplace(nx,ny,elevations(nx,ny));
      }
    }

    if(elevations(c.x,c.y)==elevations.noData())
      pit_mask(c.x,c.y)=pit_mask.noData();

    progress.update(processed_cells);
  }
  std::cerr<<"\t\033[96msucceeded in "<<progress.stop()<<"s.\033[39m"<<std::endl;
  std::cerr<<processed_cells<<" cells processed. "<<pitc<<" in pits."<<std::endl;
}



//priority_flood_watersheds
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
    If true, then **elevations** is altered as though improved_priority_flood()
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
*/
template<class elev_t>
void priority_flood_watersheds(
  Array2D<elev_t> &elevations, Array2D<int32_t> &labels, bool alter_elevations
){
  GridCellZ_pq<elev_t> open;
  std::queue<GridCellZ<elev_t> > pit;
  unsigned long processed_cells=0;
  unsigned long pitc=0,openc=0;
  int clabel=1;  //TODO: Thought this was more clear than zero in the results.
  ProgressBar progress;

  std::cerr<<"\n###Priority-Flood+Watershed Labels"<<std::endl;
  std::cerr<<"Citation: Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024"<<std::endl;
  std::cerr<<"Setting up boolean flood array matrix..."<<std::flush;
  Array2D<int8_t> closed(elevations.width(),elevations.height(),false);
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"Setting up watershed label matrix..."<<std::flush;
  labels.resize(elevations.width(),elevations.height(),-1);
  labels.setNoData(-1);
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"The priority queue will require approximately "
           <<(elevations.width()*2+elevations.height()*2)*((long)sizeof(GridCellZ<elev_t>))/1024/1024
           <<"MB of RAM."
           <<std::endl;
  std::cerr<<"Adding cells to the priority queue..."<<std::endl;
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
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"%%Performing Priority-Flood+Watershed Labels..."<<std::endl;
  progress.start( elevations.width()*elevations.height() );
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

    for(int n=1;n<=8;n++){
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

  std::cerr<<"\t\033[96msucceeded in "<<progress.stop()<<"s.\033[39m"<<std::endl;

  std::cerr<<processed_cells<<" cells processed. "
           <<pitc           <<" in pits, "
           <<openc          <<" not in pits."
           <<std::endl;
}



#endif
