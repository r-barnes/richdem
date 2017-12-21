#ifndef _richdem_depressions_
#define _richdem_depressions_

#include <richdem/depressions/priority_flood.hpp>
#include <richdem/depressions/Zhou2016pf.hpp>

namespace richdem {

template<class T> void FillDepressions       (Array2D<T> &dem){ Zhou2016              (dem); }
template<class T> void FillDepressionsEpsilon(Array2D<T> &dem){ priority_flood_epsilon(dem); }

}

#endif
