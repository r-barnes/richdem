#ifndef _richdem_hpp_
#define _richdem_hpp_

#include "common/Array2D.hpp"
#include "common/constants.hpp"
#include "common/grid_cell.hpp"
#include "common/memory.hpp"
#include "common/ProgressBar.hpp"
#include "common/timer.hpp"
#include "common/version.hpp"

#include "depressions/Lindsay2016.hpp"
#include "depressions/priority_flood.hpp"
#include "depressions/Zhou2016pf.hpp"
#include "depressions/depressions.hpp"

#include "flats/flat_resolution_dinf.hpp"
#include "flats/flat_resolution.hpp"

#include "flowdirs/d8_flowdirs.hpp"
#include "flowdirs/dinf_flowdirs.hpp"

#include "methods/d8_methods.hpp"
#include "methods/dall_methods.hpp"
#include "methods/dinf_methods.hpp"

#ifdef USEGDAL
#include "common/gdal.hpp"
#endif

#endif
