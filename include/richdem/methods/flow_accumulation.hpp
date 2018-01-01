#ifndef _richdem_flow_accumulation_hpp_
#define _richdem_flow_accumulation_hpp_

#include <richdem/flowmet/Fairfield1991.hpp>
#include <richdem/flowmet/Freeman1991.hpp>
#include <richdem/flowmet/Holmgren1994.hpp>
#include <richdem/flowmet/OCallaghan1984.hpp>
#include <richdem/flowmet/Orlandini2003.hpp>
#include <richdem/flowmet/Quinn1991.hpp>
#include <richdem/flowmet/Seibert2007.hpp>
#include <richdem/flowmet/Tarboton1997.hpp>
#include <richdem/methods/flow_accumulation_generic.hpp>

namespace richdem {

template<class elev_t, class accum_t> void FA_Tarboton          (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FM_Tarboton<elev_t>         , elevations, accum); }
template<class elev_t, class accum_t> void FA_Holmgren          (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum, double xparam) {FlowAccumulation(FM_Holmgren<elev_t>         , elevations, accum, xparam); }
template<class elev_t, class accum_t> void FA_Quinn             (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FM_Quinn<elev_t>            , elevations, accum); }
template<class elev_t, class accum_t> void FA_Freeman           (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum, double xparam) {FlowAccumulation(FM_Freeman<elev_t>          , elevations, accum, xparam); }
template<class elev_t, class accum_t> void FA_FairfieldLeymarie (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FM_FairfieldLeymarie<elev_t>, elevations, accum); }
template<class elev_t, class accum_t> void FA_Rho8              (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FM_Rho8<elev_t>             , elevations, accum); }
template<class elev_t, class accum_t> void FA_D8                (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FM_D8<elev_t>               , elevations, accum); }
template<class elev_t, class accum_t> void FA_OCallaghan        (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FM_D8<elev_t>               , elevations, accum); }

}

#endif
