#ifndef _richdem_depressions_
#define _richdem_depressions_

#include <richdem/common/constants.hpp>
#include <richdem/depressions/Barnes2014.hpp>
#include <richdem/depressions/Lindsay2016.hpp>
#include <richdem/depressions/Wei2018.hpp>
#include <richdem/depressions/Zhou2016.hpp>
#include <stdexcept>

namespace richdem {

template<Topology topo, class T>
void FillDepressions(Array2D<T> &dem){ 
  if(topo==Topology::D8)
    PriorityFlood_Zhou2016(dem);
  else if(topo==Topology::D4)
    PriorityFlood_Barnes2014<Topology::D4>(dem);
  else
    throw std::runtime_error("Unknown topology!");
}

template<Topology topo, class T> void FillDepressionsEpsilon(Array2D<T> &dem){ PriorityFloodEpsilon_Barnes2014<topo>(dem); }
template<Topology topo, class T> void BreachDepressions     (Array2D<T> &dem){ CompleteBreaching_Lindsay2016<topo>  (dem); }

}

#endif
