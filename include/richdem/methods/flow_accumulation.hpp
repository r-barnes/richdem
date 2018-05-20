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

template<class elev_t, class accum_t> void FA_Tarboton           (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_Tarboton                       (elevations, props         );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_Dinfinity          (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_Dinfinity                      (elevations, props         );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_Holmgren           (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum, double xparam) { Array3D<float> props(elevations); FM_Holmgren                       (elevations, props, xparam );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_Quinn              (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_Quinn                          (elevations, props         );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_Freeman            (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum, double xparam) { Array3D<float> props(elevations); FM_Freeman                        (elevations, props, xparam );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_FairfieldLeymarieD8(const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_FairfieldLeymarie<Topology::D8>(elevations, props         );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_FairfieldLeymarieD4(const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_FairfieldLeymarie<Topology::D4>(elevations, props         );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_Rho8               (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_Rho8                           (elevations, props         );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_Rho4               (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_Rho4                           (elevations, props         );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_OCallaghanD8       (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_OCallaghan<Topology::D8>       (elevations, props         );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_OCallaghanD4       (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_OCallaghan<Topology::D4>       (elevations, props         );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_D8                 (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_D8                             (elevations, props         );  FlowAccumulation(props, accum); }
template<class elev_t, class accum_t> void FA_D4                 (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                { Array3D<float> props(elevations); FM_D4                             (elevations, props         );  FlowAccumulation(props, accum); }

}

#endif
