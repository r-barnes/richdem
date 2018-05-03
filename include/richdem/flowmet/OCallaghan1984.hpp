#ifndef _richdem_OCallaghan1984_hpp_
#define _richdem_OCallaghan1984_hpp_

#include "richdem/common/constants.hpp"
#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/Array3D.hpp"
#include "richdem/common/ProgressBar.hpp"

namespace richdem {

//TODO: Add Marks et al (1984)
template<class E>
void FM_OCallaghanD4(
  const Array2D<E> &elevations,
  Array3D<float> &props
){
  RDLOG_ALG_NAME<<"O'Callaghan (1984)/Marks (1984) D4 Flow Accumulation";
  RDLOG_CITATION<<"O'Callaghan, J.F., Mark, D.M., 1984. The Extraction of Drainage Networks from Digital Elevation Data. Computer vision, graphics, and image processing 28, 323--344.";

  props.setAll(NO_FLOW_GEN);

  ProgressBar progress;
  progress.start(elevations.size());

  #pragma omp parallel for collapse(2)
  for(int y=1;y<elevations.height()-1;y++)
  for(int x=1;x<elevations.width()-1;x++){
    ++progress;

    const int ci = elevations.xyToI(x,y);
    const E   e  = elevations(x,y);

    int lowest_n      = 0;
    E   lowest_n_elev = std::numeric_limits<E>::max();
    for(int n=1;n<=8;n++){
      if(n_diag[n])         //Skip diagonals
        continue;

      const int ni = ci + elevations.nshift(n);

      if(elevations.isNoData(ni)) //TODO: Don't I want water to drain this way?
        continue;

      const E ne = elevations(ni);

      if(ne>=e)
        continue;

      if(ne<lowest_n_elev){
        lowest_n_elev = ne;
        lowest_n      = n;
      }
    }

    if(lowest_n==0)
      continue;

    props(x,y,0) = HAS_FLOW_GEN;

    assert(elevations(ci)>=elevations(ci+elevations.nshift(lowest_n))); //Ensure flow goes downhill

    props(x,y,lowest_n) = 1;
  }
  progress.stop();
}

//TODO: Add Marks et al (1984)
template<class E>
void FM_OCallaghanD8(
  const Array2D<E> &elevations,
  Array3D<float> &props
){
  RDLOG_ALG_NAME<<"O'Callaghan (1984)/Marks (1984) D8 Flow Accumulation";
  RDLOG_CITATION<<"O'Callaghan, J.F., Mark, D.M., 1984. The Extraction of Drainage Networks from Digital Elevation Data. Computer vision, graphics, and image processing 28, 323--344.";

  props.setAll(NO_FLOW_GEN);

  ProgressBar progress;
  progress.start(elevations.size());

  #pragma omp parallel for collapse(2)
  for(int y=1;y<elevations.height()-1;y++)
  for(int x=1;x<elevations.width()-1;x++){
    ++progress;

    const int ci = elevations.xyToI(x,y);
    const E   e  = elevations(x,y);

    int lowest_n      = 0;
    E   lowest_n_elev = std::numeric_limits<E>::max();
    for(int n=1;n<=8;n++){
      const int ni = ci + elevations.nshift(n);

      if(elevations.isNoData(ni)) //TODO: Don't I want water to drain this way?
        continue;

      const E ne = elevations(ni);

      if(ne>=e)
        continue;

      if(ne<lowest_n_elev){
        lowest_n_elev = ne;
        lowest_n      = n;
      }
    }

    if(lowest_n==0)
      continue;

    props(x,y,0) = HAS_FLOW_GEN;

    assert(elevations(ci)>=elevations(ci+elevations.nshift(lowest_n))); //Ensure flow goes downhill

    props(x,y,lowest_n) = 1;
  }
  progress.stop();
}


template<class E>
void FM_D8(const Array2D<E> &elevations, Array3D<float> &props){
  FM_OCallaghanD8(elevations, props);
}

template<class E>
void FM_D4(const Array2D<E> &elevations, Array3D<float> &props){
  FM_OCallaghanD4(elevations, props);
}

}

#endif
