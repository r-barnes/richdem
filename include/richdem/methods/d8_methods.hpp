/**
  @file
  @brief Defines a number of functions for calculating terrain attributes

  Richard Barnes (rbarnes@umn.edu), 2015
*/
//TODO: Curvature may be "Least Squares Fitted Plane" per http://gis.stackexchange.com/questions/37066/how-to-calculate-terrain-curvature
#ifndef _richdem_d8_methods_hpp_
#define _richdem_d8_methods_hpp_

#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/constants.hpp"
#include "richdem/common/grid_cell.hpp"
#include "richdem/common/ProgressBar.hpp"
#include <queue>
#include <stdexcept>

namespace richdem {

/**
  @brief  Returns the sign (+1, -1, 0) of a number. Branchless.
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  val  Input value

  @return
    -1 for a negative input, +1 for a positive input, and 0 for a zero input
*/
template <class T>
static inline T sgn(T val){
  return (T(0) < val) - (val < T(0));
}


/**
  @brief  Calculates the D8 flow accumulation, given the D8 flow directions
  @author Richard Barnes (rbarnes@umn.edu)

  This calculates the D8 flow accumulation of a grid of D8 flow directions by
  calculating each cell's dependency on its neighbours and then using a
  priority-queue to process cells in a top-of-the-watershed-down fashion

  @param[in]  &flowdirs  A D8 flowdir grid from d8_flow_directions()
  @param[out] &area      Returns the up-slope area of each cell
*/
template<class T, class U>
void d8_flow_accum(const Array2D<T> &flowdirs, Array2D<U> &area){
  std::queue<GridCell> sources;
  ProgressBar progress;

  RDLOG_ALG_NAME<<"D8 Flow Accumulation";
  RDLOG_CITATION<<"TODO";

  RDLOG_MEM_USE<<"The sources queue will require at most approximately "
               <<(flowdirs.size()*((long)sizeof(GridCell))/1024/1024)
               <<"MB of RAM.";

  RDLOG_PROGRESS<<"Resizing dependency matrix...";
  Array2D<int8_t> dependency(flowdirs,0);

  RDLOG_PROGRESS<<"Setting up the area matrix...";
  area.resize(flowdirs,0);
  area.setNoData(-1);

  RDLOG_PROGRESS<<"Calculating dependency matrix & setting noData() cells...";
  progress.start( flowdirs.size() );
  #pragma omp parallel for
  for(int y=0;y<flowdirs.height();y++){
    progress.update( y*flowdirs.width() );
    for(int x=0;x<flowdirs.width();x++){
      if(flowdirs.isNoData(x,y)){
        area(x,y) = area.noData();
        continue;
      }

      int n = flowdirs(x,y); //The neighbour this cell flows into
      if(n==NO_FLOW)         //This cell does not flow into a neighbour
        continue;

      int nx = x+dx[n];      //x-coordinate of the neighbour
      int ny = y+dy[n];      //y-coordinate of the neighbour

      //Neighbour is not on the grid
      if(!flowdirs.inGrid(nx,ny))
        continue;

      //Neighbour is valid and is part of the grid. The neighbour depends on this
      //cell, so increment its dependency count.
      ++dependency(nx,ny);
    }
  }
  RDLOG_TIME_USE<<"Dependency calculation time = "<<progress.stop()<<" s";

  RDLOG_PROGRESS<<"Locating source cells...";
  for(int y=0;y<flowdirs.height();y++)
  for(int x=0;x<flowdirs.width();x++)
    if(dependency(x,y)==0 && !flowdirs.isNoData(x,y))
      sources.emplace(x,y);

  RDLOG_PROGRESS<<"Calculating flow accumulation areas...";
  progress.start(flowdirs.numDataCells());
  long int ccount=0;
  while(sources.size()>0){
    GridCell c=sources.front();
    sources.pop();

    ccount++;
    progress.update(ccount);

    area(c.x,c.y)++;

    int n = flowdirs(c.x,c.y);

    if(n==NO_FLOW)
      continue;

    int nx=c.x+dx[n];
    int ny=c.y+dy[n];

    if(!flowdirs.inGrid(nx,ny))
      continue;
    if(flowdirs.isNoData(nx,ny))
      continue;

    area(nx,ny)+=area(c.x,c.y);
    --dependency(nx,ny);

    if(dependency(nx,ny)==0)
      sources.emplace(nx,ny);
  }
  RDLOG_TIME_USE<<"Flow accumulation calculation time = "<<progress.stop()<<" s";

  //TODO: Explain this better
  int loops=0;
  for(int i=-1;i>=-8;i--)
    loops+=dependency.countval(-1);
  RDLOG_MISC<<"Input contained at least = "<<loops<<" loops";
}




//d8_upslope_cells
/**
  @brief  Calculates which cells ultimately D8-flow through a given cell
  @author Richard Barnes (rbarnes@umn.edu)

  Given the coordinates x0,y0 of a cell and x1,y1 of another, possibly distinct,
  cell this draws a line between the two using the Bresenham Line-Drawing
  Algorithm and returns a grid showing all the cells whose flow ultimately
  passes through the indicated cells.

  The grid has the values:

  * 1=Upslope cell
  * 2=Member of initializing line
  * All other cells have a noData() value

  @param[in]  x0              x-coordinate of start of line
  @param[in]  y0              y-coordinate of start of line
  @param[in]  x1              x-coordinate of end of line
  @param[in]  y1              y-coordinate of end of line
  @param[in]  &flowdirs       A D8 flowdir grid from d8_flow_directions()
  @param[out] &upslope_cells  A grid of 1/2/NoData, as in the description
*/
template<class T, class U>
void d8_upslope_cells(
  int x0,
  int y0,
  int x1,
  int y1,
  const Array2D<T> &flowdirs,
  Array2D<U>       &upslope_cells
){
  //TODO: ALG NAME?
  RDLOG_PROGRESS<<"Setting up the upslope_cells matrix..."<<std::flush;
  upslope_cells.resize(flowdirs);
  upslope_cells.setAll(FLOWDIR_NO_DATA);
  upslope_cells.setNoData(FLOWDIR_NO_DATA);

  ProgressBar progress;

  std::queue<GridCell> expansion;

  if(x0>x1){
    std::swap(x0,x1);
    std::swap(y0,y1);
  }

  //Modified Bresenham Line-Drawing Algorithm
  int deltax     = x1-x0;
  int deltay     = y1-y0;
  float error    = 0;
  float deltaerr = (float)deltay/(float)deltax;

  if (deltaerr<0)
    deltaerr = -deltaerr;

  RDLOG_MISC<<"Line slope is "<<deltaerr;
  int y=y0;
  for(int x=x0;x<=x1;x++){
    expansion.push(GridCell(x,y));
    upslope_cells(x,y)=2;
    error+=deltaerr;
    if (error>=0.5) {
      expansion.push(GridCell(x+1,y));
      upslope_cells(x+1,y) = 2;
      y                   += sgn(deltay);
      error               -= 1;
    }
  }

  progress.start(flowdirs.data_cells);
  long int ccount=0;
  while(expansion.size()>0){
    GridCell c=expansion.front();
    expansion.pop();

    progress.update(ccount++);

    for(int n=1;n<=8;n++)
      if(!flowdirs.inGrid(c.x+dx[n],c.y+dy[n]))
        continue;
      else if(flowdirs(c.x+dx[n],c.y+dy[n])==NO_FLOW)
        continue;
      else if(flowdirs(c.x+dx[n],c.y+dy[n])==flowdirs.noData())
        continue;
      else if(upslope_cells(c.x+dx[n],c.y+dy[n])==upslope_cells.noData() && n==d8_inverse[flowdirs(c.x+dx[n],c.y+dy[n])]){
        expansion.push(GridCell(c.x+dx[n],c.y+dy[n]));
        upslope_cells(c.x+dx[n],c.y+dy[n])=1;
      }
  }
  RDLOG_TIME_USE<<"Succeeded in "<<progress.stop();
  RDLOG_MISC<<"Found "<<ccount<<" up-slope cells."; //TODO
}






}

#endif
