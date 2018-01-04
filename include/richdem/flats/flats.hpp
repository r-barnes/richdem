#ifndef _richdem_flats_hpp_
#define _richdem_flats_hpp_

#include <richdem/flats/Barnes2014.hpp>

namespace richdem {

/**
  @brief  Alters the elevations of the DEM so that all flats drain
  @author Richard Barnes (rbarnes@umn.edu)

  This alters elevations within the DEM so that all cells will have a drainage
  path.

  @param[in]     &elevations  An elevations field

  @post
    1. Every cell which is part of a flat that can be drained will have its
       elevation altered in such a way as to guarantee that it does drain.
*/
template<class T>
void ResolveFlatsEpsilon(
  Array2D<T> &elevations
){
  Array2D<int32_t> flat_mask, labels;
  GetFlatMask(elevations, flat_mask, labels);
  ResolveFlatsEpsilon_Barnes2014(flat_mask, labels, elevations);
}

}

#endif
