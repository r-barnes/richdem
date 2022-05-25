#ifndef _richdem_upslope_cells_functions_hpp_
#define _richdem_upslope_cells_functions_hpp_

#include <richdem/common/Array2D.hpp>
#include <richdem/common/Array3D.hpp>
#include <richdem/common/logger.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <queue>

namespace richdem {

/**
 * @brief Add all non-nodata cells to the return queue and sets the
 * ``upslope_cells`` value to 2 there.
 * 
 * @tparam T upslope cells dtype
 * @tparam U mask dtype
 * @param[in] mask The mask
 * @param[out] upslope_cells grid for indicating upslope cells
 * @return std::queue<GridCell> 
 */
template<class T, class U>
std::queue<GridCell> queue_from_mask(
  const Array2D<U>  &mask,
  Array2D<T>              &upslope_cells
){
  std::queue<GridCell> expansion;
  // this cannot be parallelized, because it appends everything to the same queue
  for (auto x=0;x<mask.width();x++){
    for (auto y=0;y<mask.height();y++){
      if(mask(x,y)!=mask.noData()){
        upslope_cells(x,y)=2;
        expansion.emplace(x,y);
      }
    }
  }

  return expansion;
}

/**
 * @brief Draws a line between the points (x0,y0) and (x1,y1). Adds visited
 * cells to the return queue and sets them to 2 in ``upslope_cells``.
 * 
 * @tparam T upslope cells dtype
 * @param[in] x0 x coordinate of first point
 * @param[in] y0 y coordinate of first point
 * @param[in] x1 x coordinate of second point
 * @param[in] y1 y coordinate of second point
 * @param[out] upslope_cells grid for indicating upslope cells
 * @return std::queue<GridCell> A queue with gridcells that were visited
 */
template<class T>
std::queue<GridCell> queue_from_linepoints(
  int x0,
  int y0,
  int x1,
  int y1,
  Array2D<T>       &upslope_cells
){
  std::queue<GridCell> expansion;

  // if the slope>1, swap x and y coordinates.
  // It works exatly the same, but points
  // are drawn with swapped coordinates
  bool swapxy = false;
  if(std::abs(float(y0-y1)/float(x0-x1))>1){
    std::cout << "swapped x and y" << std::endl;
    std::swap(x0,y0);
    std::swap(x1,y1);
    swapxy= true;
  }

  if(x0>x1){
    std::swap(x0,x1);
    std::swap(y0,y1);
  }

  // modified Bresenham Line-Drawing Algorithm from
  // https://www.geeksforgeeks.org/bresenhams-line-generation-algorithm/
  int deltaerr = 2*(y1-y0);
  int error = deltaerr - (x1-x0);

  RDLOG_MISC<<"Line slope is "<<deltaerr;
  int y=y0;
  for(int x=x0;x<=x1;x++){
    if (!swapxy){
      expansion.emplace(x,y);
      upslope_cells(x,y)=2;
    } else {
      expansion.emplace(y,x);
      upslope_cells(y,x)=2;
    }
    error += deltaerr;
    if (error>=0 && x!=x1) { //safeguard for segfault if we draw to the edge of the grid
      // Draw the point above if we move to te next pixel
      // the difference between a "solid" and non-solid line:
      // o o X x x <--the X
      // x x x 0 0 
      if(!swapxy){
        expansion.emplace(x+1,y);
        upslope_cells(x+1,y) = 2;
      } else {
        expansion.emplace(y,x+1);
        upslope_cells(y,x+1) = 2;
      }
      y++;
      error -= 2*(x1-x0);
    }
  }
  return expansion;
}

//upslope_from_props
/**
 * @brief Calculates uplsope cells for a given set of input cells, given a
 * proportions matrix of flow. All cells that have flow into input cells will be added
 * 
 * @param[in]  &expansion      queue with starting cells
 * @param[in]  &elevation      DEM
 * @param[out] &upslope_cells  A grid of 1/2/NoData
 */
template<Topology topo, class T, class U>
void upslope_cells_props(
  std::queue<GridCell> &expansion,
  const Array3D<T>     &props,
  Array2D<U>           &upslope_cells
){
  //TODO: ALG NAME?
  RDLOG_PROGRESS<<"Setting up the upslope_cells matrix..."<<std::flush;
  // upslope_cells.resize(props.width(),props.height());
  // upslope_cells.setAll(0);
  // upslope_cells.setNoData(0);

  static_assert(topo==Topology::D8 || topo==Topology::D4);
  constexpr auto dx = get_dx_for_topology<topo>();
  constexpr auto dy = get_dy_for_topology<topo>();
  constexpr auto nmax = get_nmax_for_topology<topo>();

  ProgressBar progress;

  progress.start(props.size());
  long int ccount=0;
  while(expansion.size()>0){
    GridCell c=expansion.front();
    expansion.pop();

    progress.update(ccount++);

    for(int n=1;n<=8;n++)
      if(!props.inGrid(c.x+dx[n],c.y+dy[n]))
        continue;
      else if(upslope_cells.isNoData(c.x+dx[n],c.y+dy[n]) && props(c.x+dx[n],c.y+dy[n],d8_inverse[n])>0) {
        expansion.emplace(c.x+dx[n],c.y+dy[n]);
        upslope_cells(c.x+dx[n],c.y+dy[n])=1;
      }
  }
  RDLOG_TIME_USE<<"Succeeded in "<<progress.stop();
  RDLOG_MISC<<"Found "<<ccount<<" up-slope cells."; //TODO
}

//upslope_cells multiple flow implementation
/**
 * @brief Calculates uplsope cells for a given set of input cells, assuming
 * multiple flow. That is, all cells that are higher than the cell being
 * processed will be added, regardless of whether the current cell is the lowest
 * neighbour.
 * 
 * @param[in]  &expansion      queue with starting cells
 * @param[in]  &elevation      DEM
 * @param[out] &upslope_cells  A grid of 1/2/NoData
 */
template<Topology topo, class T, class U>
void upslope_cells_mflow(
  std::queue<GridCell> &expansion,
  const Array2D<T>     &elevation,
  Array2D<U>           &upslope_cells
){
  //TODO: ALG NAME?
  RDLOG_PROGRESS<<"Setting up the upslope_cells matrix..."<<std::flush;
  // upslope_cells.resize(elevation);
  // upslope_cells.setAll(0);
  // upslope_cells.setNoData(0);

  static_assert(topo==Topology::D8 || topo==Topology::D4);
  constexpr auto dx = get_dx_for_topology<topo>();
  constexpr auto dy = get_dy_for_topology<topo>();
  constexpr auto nmax = get_nmax_for_topology<topo>();

  ProgressBar progress;

  progress.start(elevation.size());
  long int ccount=0;
  while(expansion.size()>0){
    GridCell c=expansion.front();
    expansion.pop();

    progress.update(ccount++);

    for(int n=1;n<=8;n++)
      if(!elevation.inGrid(c.x+dx[n],c.y+dy[n]))
        continue;
      else if(elevation.isNoData(c.x+dx[n],c.y+dy[n]))
        continue;
      else if(upslope_cells.isNoData(c.x+dx[n],c.y+dy[n]) && elevation(c.x,c.y)<elevation(c.x+dx[n],c.y+dy[n])) {
        expansion.emplace(c.x+dx[n],c.y+dy[n]);
        upslope_cells(c.x+dx[n],c.y+dy[n])=1;
      }
  }
  RDLOG_TIME_USE<<"Succeeded in "<<progress.stop();
  RDLOG_MISC<<"Found "<<ccount<<" up-slope cells."; //TODO
}

}
#endif