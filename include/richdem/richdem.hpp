#ifndef _richdem_hpp_
#define _richdem_hpp_

#include "common/Array2D.hpp"
#include "common/constants.hpp"
#include "common/grid_cell.hpp"
#include "common/ManagedVector.hpp"
#include "common/memory.hpp"
#include "common/ProgressBar.hpp"
#include "common/random.hpp"
#include "common/timer.hpp"
#include "common/version.hpp"

#include "depressions/Barnes2014.hpp"
#include "depressions/depressions.hpp"
#include "depressions/Lindsay2016.hpp"
#include "depressions/Zhou2016.hpp"

#include "flats/flat_resolution.hpp"
#include "flats/flat_resolution_dinf.hpp"

#include "flowmet/d8_flowdirs.hpp"
#include "flowmet/dinf_flowdirs.hpp"
#include "flowmet/Fairfield1991.hpp"
#include "flowmet/Freeman1991.hpp"
#include "flowmet/Holmgren1994.hpp"
#include "flowmet/OCallaghan1984.hpp"
#include "flowmet/Orlandini2003.hpp"
#include "flowmet/Quinn1991.hpp"
#include "flowmet/Seibert2007.hpp"
#include "flowmet/Tarboton1997.hpp"

#include "methods/d8_methods.hpp"
#include "methods/dinf_methods.hpp"
#include "methods/flow_accumulation.hpp"
#include "methods/flow_accumulation_generic.hpp"
#include "methods/strahler.hpp"
#include "methods/terrain_attributes.hpp"

#ifdef USEGDAL
#include "common/gdal.hpp"
#endif

#endif
