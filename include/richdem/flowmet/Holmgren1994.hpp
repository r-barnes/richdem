#ifndef _richdem_Holmgren1994_
#define _richdem_Holmgren1994_

#include "richdem/common/constants.hpp"
#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/ProgressBar.hpp"

namespace richdem {

template<class E>
std::vector<float> FM_Holmgren(const Array2D<E> &elevations, const double xparam){
  RDLOG_ALG_NAME<<"Holmgren (1994) Flow Accumulation (aka MFD, MD8)";
  RDLOG_CITATION<<"Holmgren, P., 1994. Multiple flow direction algorithms for runoff modelling in grid based elevation models: an empirical evaluation. Hydrological processes 8, 327–334.";
  RDLOG_CONFIG<<"x = "<<xparam;

  std::vector<float> props(9*elevations.size(),NO_FLOW_GEN);

  constexpr double L1   = 0.5;
  constexpr double L2   = 0.354; //TODO: More decimal places
  constexpr double L[9] = {0,L1,L2,L1,L2,L1,L2,L1,L2};

  ProgressBar progress;
  progress.start(elevations.size());

  #pragma omp parallel for collapse(2)
  for(int y=1;y<elevations.height()-1;y++)
  for(int x=1;x<elevations.width()-1;x++){
    ++progress;
    const E e = elevations(x,y);

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
        props.at(9*ci+n)  = std::pow(grad * L[n],xparam);
        C                += props.at(9*ci+n);
      }
    }

    if(C>0){
      props.at(9*ci+0) = HAS_FLOW_GEN;

      C = 1/C;

      for(int n=1;n<=8;n++){
        if(props[9*ci+n]>0)
          props.at(9*ci+n) *= C;
        else
          props.at(9*ci+n) = 0;
      }
    }
  }
  progress.stop();

  return props;
}

}

#endif
