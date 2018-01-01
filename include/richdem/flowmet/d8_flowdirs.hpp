/**
  @file
  @brief Functions for calculating D8 flow directions

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_d8_flowdirs_hpp_
#define _richdem_d8_flowdirs_hpp_

#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/ProgressBar.hpp"

namespace richdem {

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
    if(elevations.isTopLeft(x,y))
      return 2;
    else if(elevations.isBottomLeft(x,y))
      return 8;
    else if(elevations.isTopRight(x,y))
      return 4;
    else if(elevations.isBottomRight(x,y))
      return 6;
    else if(elevations.isLeftCol(x,y))
      return 1;
    else if(elevations.isRightCol(x,y))
      return 5;
    else if(elevations.isTopRow(x,y))
      return 3;
    else if(elevations.isBottomRow(x,y))
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

  RDLOG_ALG_NAME<<"D8 Flow Directions";
  RDLOG_CITATION<<"TODO";

  RDLOG_PROGRESS<<"Setting up the flow directions matrix...";
  flowdirs.resize(elevations);
  flowdirs.setAll(NO_FLOW);
  flowdirs.setNoData(FLOWDIR_NO_DATA);

  RDLOG_PROGRESS<<"Calculating D8 flow directions...";
  progress.start( elevations.width()*elevations.height() );
  #pragma omp parallel for
  for(int y=0;y<elevations.height();y++){
    progress.update( y*elevations.width() );
    for(int x=0;x<elevations.width();x++)
      if(elevations(x,y)==elevations.noData())
        flowdirs(x,y) = flowdirs.noData();
      else
        flowdirs(x,y) = d8_FlowDir(elevations,x,y);
  }
  RDLOG_TIME_USE<<"Succeeded in = "<<progress.stop()<<" s";
}

}

#endif
