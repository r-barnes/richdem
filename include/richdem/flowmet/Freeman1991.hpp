#ifndef _richdem_Freeman1991_hpp_
#define _richdem_Freeman1991_hpp_

#include "richdem/common/constants.hpp"
#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/ProgressBar.hpp"

namespace richdem {

template<class E>
std::vector<float> FM_Freeman(
  const Array2D<E> &elevations,
  const double xparam
){
  RDLOG_ALG_NAME<<"Freeman (1991) Flow Accumulation (aka MFD, MD8)";
  RDLOG_CITATION<<"Freeman, T.G., 1991. Calculating catchment area with divergent flow based on a regular grid. Computers & Geosciences 17, 413–422.";
  RDLOG_CONFIG<<"p = "<<xparam;

  std::vector<float> props(9*elevations.size(),NO_FLOW_GEN);

  ProgressBar progress;
  progress.start(elevations.size());

  #pragma omp parallel for collapse(2)
  for(int y=1;y<elevations.height()-1;y++)
  for(int x=1;x<elevations.width()-1;x++){
    ++progress;

    const E e    = elevations(x,y);
    const int ci = elevations.xyToI(x,y);

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
        props[9*ci+n]     = cval;
        C                += cval;
      }
    }

    if(C>0){
      props.at(9*ci+0) = HAS_FLOW_GEN;

      C = 1/C; //TODO

      for(int n=1;n<=8;n++){
        auto &this_por = props.at(9*ci+n);
        if(this_por>0)
          this_por *= C;
        else
          this_por = 0;
      }
    }
  }
  progress.stop();

  return props;
}

}

#endif
