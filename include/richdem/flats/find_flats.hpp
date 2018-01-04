#ifndef _richdem_find_flats_hpp_
#define _richdem_find_flats_hpp_

#include <richdem/common/logger.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <richdem/common/Array2D.hpp>

namespace richdem {

#define FLAT_NO_DATA    -1
#define NOT_A_FLAT       0
#define IS_A_FLAT        1

/**
  @brief  Finds flats: cells with no local gradient
  @author Richard Barnes (rbarnes@umn.edu)

  TODO

  @param[in]  &elevations     An elevation field
  @param[in]  &flats          An array of flat indicators (post-conditions)

  @post
    1. flats contains the value 1 for each cell in a flat (all cells without a
       local gradient), 0 for each cell not in a flat (all cells with a local
       gradient), and -1 for each no data cell.
*/
template<class T>
void FindFlats(
  const Array2D<T>   &elevations,
  Array2D<int8_t>    &flats
){
  flats.resize(elevations);
  flats.setNoData(FLAT_NO_DATA);

  ProgressBar progress;

  progress.start( elevations.size() );

  #pragma omp parallel for
  for(int y=0;y<elevations.height();y++)
  for(int x=0;x<elevations.width();x++){
    if(elevations.isNoData(x,y)){
      flats(x,y) = FLAT_NO_DATA;
      continue;
    }

    if(elevations.isEdgeCell(x,y)){
      flats(x,y) = NOT_A_FLAT;
      continue;
    }

    //We'll now assume that the cell is a flat unless proven otherwise
    flats(x,y) = IS_A_FLAT;

    for(int n=1;n<=8;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];
      if(elevations(nx,ny)<elevations(x,y) || elevations.isNoData(nx,ny)){
        flats(x,y) = NOT_A_FLAT;
        break;
      }
    }

    //We handled the base case just above the for loop
  }

  RDLOG_TIME_USE<<"Succeeded in = "<<progress.stop()<<" s";
}

}

#endif
