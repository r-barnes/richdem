#ifndef _richdem_Holmgren1994_
#define _richdem_Holmgren1994_

#include "richdem/common/constants.hpp"
#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/Array3D.hpp"
#include "richdem/common/ProgressBar.hpp"

namespace richdem {

template<class E>
void FM_Holmgren(
  const Array2D<E> &elevations,
  Array3D<float> &props,
  const double xparam
){
  RDLOG_ALG_NAME<<"Holmgren (1994) Flow Accumulation (aka MFD, MD8)";
  RDLOG_CITATION<<"Holmgren, P., 1994. Multiple flow direction algorithms for runoff modelling in grid based elevation models: an empirical evaluation. Hydrological processes 8, 327–334.";
  RDLOG_CONFIG<<"x = "<<xparam;

  props.setAll(NO_FLOW_GEN);
  props.setNoData(NO_DATA_GEN);

  constexpr double L1   = 0.5;
  constexpr double L2   = 0.354; //TODO: More decimal places
  constexpr double L[9] = {0,L1,L2,L1,L2,L1,L2,L1,L2};

  ProgressBar progress;
  progress.start(elevations.size());

  #pragma omp parallel for collapse(2)
  for(int y=0;y<elevations.height();y++)
  for(int x=0;x<elevations.width();x++){
    ++progress;

    if(elevations.isNoData(x,y)){
      props(x,y,0) = NO_DATA_GEN;
      continue;
    }

    if(elevations.isEdgeCell(x,y))
      continue;
    
    const E e = elevations(x,y);

    double C = 0;

    for(int n=1;n<=8;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];

      if(!elevations.inGrid(nx,ny))
        continue;
      if(elevations.isNoData(nx,ny)) //TODO: Don't I want water to drain this way?
        continue;

      const E ne = elevations(nx,ny);

      if(ne<e){
        const double rise = e-ne;
        const double run  = dr[n];
        const double grad = rise/run;
        props(x,y,n)      = std::pow(grad * L[n],xparam);
        C                += props(x,y,n);
      }
    }

    if(C>0){
      props(x,y,0) = HAS_FLOW_GEN;

      C = 1/C;

      for(int n=1;n<=8;n++){
        if(props(x,y,n)>0)
          props(x,y,n) *= C;
        else
          props(x,y,n) = 0;
      }
    }
  }
  progress.stop();
}

}

#endif
