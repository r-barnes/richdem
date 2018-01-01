#ifndef _richdem_depressions_
#define _richdem_depressions_

#include <richdem/depressions/Barnes2014.hpp>
#include <richdem/depressions/Zhou2016.hpp>

namespace richdem {

template<class T> void FillDepressions       (Array2D<T> &dem){ Zhou2016              (dem); }
template<class T> void FillDepressionsEpsilon(Array2D<T> &dem){ PriorityFloodEpsilon_Barnes2014(dem); }

}

#endif
