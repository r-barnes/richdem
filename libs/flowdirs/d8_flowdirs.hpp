#ifndef _richdem_d8_flowdirs_hpp_
#define _richdem_d8_flowdirs_hpp_

#include "../common/Array2D.hpp"
#include "../common/interface.hpp"

///Value used to indicate that a flow direction cell has no data
#define d8_NO_DATA    113

///Value used to indicate that a cell does not have a defined flow direction
//(i.e. that it has no local gradient)
#define NO_FLOW       0

//Neighbour directions
//345
//102
//876

//Inverse neighbour directions
//678
//201
//543

///sqrt(2), used to generate distances from a central cell to its neighbours
#define SQRT2         1.414213562373095048801688724209698078569671875376948

//D8 Directions
#ifndef d8flowdirs_dxdy
#define d8flowdirs_dxdy
///x offsets of D8 neighbours, from a central cell
const int dx[9]={0,-1,1,-1,0,1,1,0,-1};  //TODO: These should be merged with my new dinf_d8 to reflect a uniform and intelligent directional system
///y offsets of D8 neighbours, from a central cell
const int dy[9]={0,0,0,-1,-1,-1,1,1,1};
#endif
///Arrows indicating flow directions
const wchar_t fd[9]={L'·',L'←',L'↖',L'↑',L'↗',L'→',L'↘',L'↓',L'↙'};
///Distances from a central cell to each of its 8 neighbours
const double dr[9]={0,1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2}; //TODO: Check these for new D8 directions
///dx[] and dy[] offsets are labeled 0-8. This maps the inverse path from each of those cells.
const int inverse_flow[9]={0,2,6,7,8,1,3,4,5}; //Inverse of a given n from chart below
//derived from the following neighbour directions
//234
//105
//876






//234
//105
//876
//d8_FlowDir
/**
  @brief  Calculates the D8 flow direction of a cell
  @author Richard Barnes (rbarnes@umn.edu)

  This calculates the D8 flow direction of a cell using the D8
  neighbour system, as defined in utility.h. Cells on the edge
  of the grid flow off the nearest edge.

  Helper function for d8_flow_directions().

  @param[in]  &elevations  A DEM
  @param[in]  x            x coordinate of cell
  @param[in]  y            y coordinate of cell

  @returns The D8 flow direction of the cell
*/
template<class T>
static int d8_FlowDir(const Array2D<T> &elevations, const int x, const int y){
  T minimum_elevation = elevations(x,y);
  int flowdir         = NO_FLOW;

  if (elevations.isEdgeCell(x,y)){
    if(x==0 && y==0)
      return 2;
    else if(x==0 && y==elevations.viewHeight()-1)
      return 8;
    else if(x==elevations.viewWidth()-1 && y==0)
      return 4;
    else if(x==elevations.viewWidth()-1 && y==elevations.viewHeight()-1)
      return 6;
    else if(x==0)
      return 1;
    else if(x==elevations.viewWidth()-1)
      return 5;
    else if(y==0)
      return 3;
    else if(y==elevations.viewHeight()-1)
      return 7;
  }

  /*NOTE: Since the very edges of the DEM are defined to always flow outwards,
  if they have defined elevations, it is not necessary to check if a neighbour
  is IN_GRID in the following
  NOTE: It is assumed that the no_data datum is an extremely negative
  number, such that all water which makes it to the edge of the DEM's region
  of defined elevations is sucked directly off the grid, rather than piling up
  on the edges.*/
  for(int n=1;n<=8;n++)
    if(
      elevations(x+dx[n],y+dy[n])<minimum_elevation
      || (elevations(x+dx[n],y+dy[n])==minimum_elevation
            && flowdir>0 && flowdir%2==0 && n%2==1) //TODO: What is this modulus stuff for?
    ){
      minimum_elevation=elevations(x+dx[n],y+dy[n]);
      flowdir=n;
    }

  return flowdir;
}



//d8_flow_directions
/**
  @brief  Calculates the D8 flow directions of a DEM
  @author Richard Barnes (rbarnes@umn.edu)

  This calculates the D8 flow directions of a DEM. Its argument
  'flowdirs' will return a grid with flow directions using the D8
  neighbour system, as defined in utility.h. The choice of data type
  for array2d must be able to hold exact values for all neighbour
  identifiers (usually [-1,7]).

  Uses d8_FlowDir() as a helper function.

  @todo                    Combine dinf and d8 neighbour systems

  @param[in]  &elevations  A DEM
  @param[out] &flowdirs    Returns the flow direction of each cell
*/
template<class T, class U>
void d8_flow_directions(
  const Array2D<T> &elevations,
        Array2D<U> &flowdirs
){
  ProgressBar progress;

  std::cerr<<"Setting up the flow directions matrix..."<<std::flush;
  flowdirs.resize(elevations);
  flowdirs.init(NO_FLOW);
  flowdirs.setNoData(d8_NO_DATA);
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"%%Calculating D8 flow directions..."<<std::endl;
  progress.start( elevations.viewWidth()*elevations.viewHeight() );
  #pragma omp parallel for
  for(int x=0;x<elevations.viewWidth();x++){
    progress.update( x*elevations.viewHeight() );
    for(int y=0;y<elevations.viewHeight();y++)
      if(elevations(x,y)==elevations.noData())
        flowdirs(x,y)=flowdirs.noData();
      else
        flowdirs(x,y)=d8_FlowDir(elevations,x,y);
  }
  std::cerr<<"Succeeded in "<<progress.stop()<<"s."<<std::endl;
}

#endif