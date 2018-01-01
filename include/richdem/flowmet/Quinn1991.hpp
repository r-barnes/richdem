#ifndef _richdem_Quinn1991_hpp_
#define _richdem_Quinn1991_hpp_

#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/ProgressBar.hpp"

namespace richdem {

template<class E>
std::vector<float> FM_Quinn(const Array2D<E> &elevations){
  RDLOG_ALG_NAME<<"Quinn (1991) Flow Accumulation (aka MFD, MD8)";
  RDLOG_CITATION<<"Quinn, P., Beven, K., Chevallier, P., Planchon, O., 1991. The Prediction Of Hillslope Flow Paths For Distributed Hydrological Modelling Using Digital Terrain Models. Hydrological Processes 5, 59–79."; 
  return FM_Holmgren(elevations, 1.0);
}

}

#endif
