#ifndef _richdem_Freeman1991_hpp_
#define _richdem_Freeman1991_hpp_

#include "richdem/common/constants.hpp"
#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/Array3D.hpp"
#include "richdem/common/ProgressBar.hpp"

namespace richdem {

template<class E>
void FM_Freeman(
  const Array2D<E> &elevations,
  Array3D<float> &props,
  const double xparam
){
  RDLOG_ALG_NAME<<"Freeman (1991) Flow Accumulation (aka MFD, MD8)";
  RDLOG_CITATION<<"Freeman, T.G., 1991. Calculating catchment area with divergent flow based on a regular grid. Computers & Geosciences 17, 413–422.";
  RDLOG_CONFIG<<"p = "<<xparam;

  props.setAll(NO_FLOW_GEN);
  props.setNoData(NO_DATA_GEN);

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

    const E e    = elevations(x,y);

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
        const auto cval   = std::pow(grad,xparam);
        props(x,y,n)      = cval;
        C                += cval;
      }
    }

    if(C>0){
      props(x,y,0) = HAS_FLOW_GEN;

      C = 1/C; //TODO

      for(int n=1;n<=8;n++){
        auto &this_por = props(x,y,n);
        if(this_por>0)
          this_por *= C;
        else
          this_por = 0;
      }
    }
  }
  progress.stop();
}

}

#endif
