#ifndef _flowdirs_Tarboton1997_hpp_
#define _flowdirs_Tarboton1997_hpp_

#include "richdem/common/constants.hpp"
#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/Array3D.hpp"
#include "richdem/common/ProgressBar.hpp"

namespace richdem {

template<class elev_t>
void FM_Tarboton(
  const Array2D<elev_t> &elevations,
  Array3D<float> &props
){
  RDLOG_ALG_NAME<<"Tarboton (1997) Flow Accumulation (aka D-Infinity, D∞)";
  RDLOG_CITATION<<"Tarboton, D.G., 1997. A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water resources research 33, 309–319.";

  props.setAll(NO_FLOW_GEN);
  props.setNoData(NO_DATA_GEN);

  //TODO: Assumes that the width and height of grid cells are equal and scaled
  //to 1.
  constexpr double d1   = 1;
  constexpr double d2   = 1;
  const     float  dang = std::atan2(d2,d1);

  const auto nwrap = [](int8_t n){ return (n==9)?1:n; };

  //Table 1 of Tarboton (1997)
  //          Column #  =   0    1    2    3    4    5   6    7
  // const int    dy_e1[8] = { 0 , -1 , -1 ,  0 ,  0 ,  1 , 1 ,  0 };
  // const int    dx_e1[8] = { 1 ,  0 ,  0 , -1 , -1 ,  0 , 0 ,  1 };
  // const int    dy_e2[8] = {-1 , -1 , -1 , -1 ,  1 ,  1 , 1 ,  1 };
  // const int    dx_e2[8] = { 1 ,  1 , -1 , -1 , -1 , -1 , 1 ,  1 };
  // const double ac[8]    = { 0.,  1.,  1.,  2.,  2.,  3., 3.,  4.};
  // const double af[8]    = { 1., -1.,  1., -1.,  1., -1., 1., -1.};

  //I remapped the foregoing table for ease of use with RichDEM. The facets
  //are renumbered as follows:
  //    3->1    2->2    1->3    0->4    7->5    6->6    5->7    4->8
  //This gives the following table
  //  Remapped Facet #  =  -   1    2     3    4    5   6    7    8  
  //  Tarboton Facet #  =  -   3    2     1    0    7   6    5    4  
  const int    dy_e1[9] = {0,  0 , -1 ,  -1 ,  0 ,  0 , 1 ,  1 ,  0  };
  const int    dx_e1[9] = {0, -1 ,  0 ,   0 ,  1 ,  1 , 0 ,  0 , -1  };
  const int    dy_e2[9] = {0, -1 , -1 ,  -1 , -1 ,  1 , 1 ,  1 ,  1  };
  const int    dx_e2[9] = {0, -1 , -1 ,   1 ,  1 ,  1 , 1 , -1 , -1  };
  //const double ac[9]    = {0,  2.,  1.,   1.,  0.,  4., 3.,  3.,  2. };
  const double af[9]    = {0, -1.,  1.,  -1.,  1., -1., 1., -1.,  1. };

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

    int8_t nmax = -1;
    double smax = 0;
    float  rmax = 0;

    for(int n=1;n<=8;n++){
      if(!elevations.inGrid (x+dx_e1[n],y+dy_e1[n]))
        continue;
      if(elevations.isNoData(x+dx_e1[n],y+dy_e1[n]))
        continue;
      if(!elevations.inGrid (x+dx_e2[n],y+dy_e2[n]))
        continue;
      if(elevations.isNoData(x+dx_e2[n],y+dy_e2[n]))
        continue;

      //Is is assumed that cells with a value of NoData have very negative
      //elevations with the result that they draw flow off of the grid.

      //Choose elevations based on Table 1 of Tarboton (1997), Barnes TODO
      const double e0 = elevations(x,y);
      const double e1 = elevations(x+dx_e1[n],y+dy_e1[n]);
      const double e2 = elevations(x+dx_e2[n],y+dy_e2[n]);

      const double s1 = (e0-e1)/d1;
      const double s2 = (e1-e2)/d2;

      double r = std::atan2(s2,s1);
      double s;

      if(r<1e-7){
        r = 0;
        s = s1;
      } else if(r>dang-1e-7){
        r = dang;
        s = (e0-e2)/sqrt(d1*d1+d2*d2);
      } else {
        s = sqrt(s1*s1+s2*s2);
      }

      if(s>smax){
        smax = s;
        nmax = n;
        rmax = r;
      }
    }

    if(nmax==-1)
      continue;

    props(x,y,0) = HAS_FLOW_GEN;

    if(af[nmax]==1 && rmax==0)
      rmax = dang;
    else if(af[nmax]==1 && rmax==dang)
      rmax = 0;
    else if(af[nmax]==1)
      rmax = M_PI/4-rmax;

    //Code used by Tarboton to calculate the angle Rg. This should give the same
    //result despite the rearranged table
    // double rg = NO_FLOW;
    // if(nmax!=-1)
    //   rg = (af[nmax]*rmax+ac[nmax]*M_PI/2);

    if(rmax==0){
      props(x,y,nmax) = 1;
    } else if(rmax==dang){
      props(x,y,nwrap(nmax+1)) = 1;
    } else {
      props(x,y,nmax)          = rmax/(M_PI/4.);
      props(x,y,nwrap(nmax+1)) = 1-rmax/(M_PI/4.);      
    }
  }
  progress.stop();
}

template<class E>
void FM_Dinfinity(const Array2D<E> &elevations, Array3D<float> &props){
  FM_Tarboton(elevations, props);
}

}

#endif
