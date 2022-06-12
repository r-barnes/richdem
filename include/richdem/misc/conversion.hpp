#pragma once

#include <richdem/common/constants.hpp>
#include <richdem/common/Array2D.hpp>

#include <stdexcept>
#include <string>

namespace richdem {

void convert_arc_flowdirs_to_richdem_d8(
  const Array2D<d8_flowdir_t> &arc_flowdirs,
  Array2D<d8_flowdir_t> &rd_flowdirs
){
  assert(arc_flowdirs.width() == rd_flowdirs.width());
  assert(arc_flowdirs.height() == rd_flowdirs.height());

  for(auto i=arc_flowdirs.i0();i<arc_flowdirs.size();i++){
    switch(arc_flowdirs(i)){
      case 0:   rd_flowdirs(i) = 0; break;
      case 1:   rd_flowdirs(i) = 5; break;
      case 2:   rd_flowdirs(i) = 6; break;
      case 4:   rd_flowdirs(i) = 7; break;
      case 8:   rd_flowdirs(i) = 8; break;
      case 16:  rd_flowdirs(i) = 1; break;
      case 32:  rd_flowdirs(i) = 2; break;
      case 64:  rd_flowdirs(i) = 3; break;
      case 128: rd_flowdirs(i) = 4; break;
      default:
        throw std::runtime_error("Unknown flow direction " + std::to_string(arc_flowdirs(i)));
    }
  }
}

}