#ifndef _richdem_depressions_
#define _richdem_depressions_

#include <richdem/common/constants.hpp>
#include <richdem/depressions/Barnes2014.hpp>
#include <richdem/depressions/Lindsay2016.hpp>
#include <richdem/depressions/Wei2018.hpp>
#include <richdem/depressions/Zhou2016.hpp>

namespace richdem {

template<class T> void FillDepressionsD8       (Array2D<T> &dem){ PriorityFlood_Zhou2016<Topology::D8>         (dem); }
template<class T> void FillDepressionsEpsilonD8(Array2D<T> &dem){ PriorityFloodEpsilon_Barnes2014<Topology::D8>(dem); }
template<class T> void BreachDepressionsD8     (Array2D<T> &dem){ CompleteBreaching_Lindsay2016<Topology::D8>  (dem); }

template<class T> void FillDepressionsD4       (Array2D<T> &dem){ PriorityFlood_Zhou2016<Topology::D4>         (dem); }
template<class T> void FillDepressionsEpsilonD4(Array2D<T> &dem){ PriorityFloodEpsilon_Barnes2014<Topology::D4>(dem); }
template<class T> void BreachDepressionsD4     (Array2D<T> &dem){ CompleteBreaching_Lindsay2016<Topology::D4>  (dem); }


}

#endif
