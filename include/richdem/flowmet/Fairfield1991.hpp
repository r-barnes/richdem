#ifndef _richdem_Fairfield1991_hpp_
#define _richdem_Fairfield1991_hpp_

#include "richdem/common/constants.hpp"
#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/Array3D.hpp"
#include "richdem/common/ProgressBar.hpp"
#include "richdem/common/random.hpp"

namespace richdem {

template <Topology topo, class elev_t>
void FM_FairfieldLeymarie(const Array2D<elev_t> &elevations, Array3D<float> &props){
  RDLOG_ALG_NAME<<"Fairfield (1991) Rho8/Rho4 Flow Accumulation";
  RDLOG_CITATION<<"Fairfield, J., Leymarie, P., 1991. Drainage networks from grid digital elevation models. Water resources research 27, 709–717.";

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

    const elev_t e    = elevations(x,y);

    int    greatest_n     = 0; //TODO: Use a constant
    double greatest_slope = 0;
    for(int n=1;n<=8;n++){
      if(topo==Topology::D4 && n_diag[n]) //Skip diagonals
        continue;

      const int nx = x+dx[n];
      const int ny = y+dy[n];

      if(!elevations.inGrid(nx,ny))
        continue;
      if(elevations.isNoData(nx,ny)) //TODO: Don't I want water to drain this way?
        continue;

      const elev_t ne = elevations(nx,ny);

      if(ne>=e)
        continue;

      double rho_slope = (e-ne);
      if(topo==Topology::D8 && n_diag[n])
        rho_slope *= 1/(2-uniform_rand_real(0,1));
      else if(topo==Topology::D4 && (n==D8_NORTH || n==D8_SOUTH))
        rho_slope *= 1/(1/uniform_rand_real(0,1)-1);

      if(rho_slope>greatest_slope){
        greatest_n     = n;
        greatest_slope = rho_slope;
      }
    }

    if(greatest_n==0)
      continue;

    props(x,y,0)          = HAS_FLOW_GEN;
    props(x,y,greatest_n) = 1;

    assert(elevations(x,y)>=elevations(x+dx[greatest_n],y+dy[greatest_n])); //Ensure flow goes downhill
  }
  progress.stop();
}



template<class elev_t>
void FM_Rho8(const Array2D<elev_t> &elevations, Array3D<float> &props){
  //Algorithm headers are taken care of in FM_FairfieldLeymarie()
  FM_FairfieldLeymarie<Topology::D8>(elevations, props);
}



template<class elev_t>
void FM_Rho4(const Array2D<elev_t> &elevations, Array3D<float> &props){
  //Algorithm headers are taken care of in FM_FairfieldLeymarie()
  FM_FairfieldLeymarie<Topology::D4>(elevations, props);
}

}

#endif
