/**
  @file
  @brief Functions for calculating flow directions according to a variety of authors

  Richard Barnes (rbarnes@umn.edu), 2016
*/
#ifndef _richdem_dall_flowdirs_hpp_
#define _richdem_dall_flowdirs_hpp_

#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/ProgressBar.hpp"
#include "richdem/common/grid_cell.hpp"
#include "richdem/common/random.hpp"

#include <iomanip> //TODO: Cut

namespace richdem {

// template<class A>
// static inline void StrahlerNumber(
//   const int x,
//   const int y,
//   const int n,
//   dep_t                &dep,
//   Array2D<A>           &strahler,
//   std::queue<GridCell> &q,
//   A                     flow
// ){
//   constexpr A msb  = 1 << (std::numeric_limits<A>::digits-1); //First time this Strahler # was seen
//   constexpr A smsb = 1 << (std::numeric_limits<A>::digits-2); //This Strahler number arose from combining two other neighbours and should be overriden if appropriate

//   const auto s = strahler(x,y) & (~msb) & (~smsb);

//   const int   nx = x+dx[n];
//   const int   ny = y+dy[n];
//         auto& ns = strahler(nx,ny);

//   const bool first_time = ns & msb;
//   const bool emergent   = ns & smsb;

//   //Clear info bits
//   ns &= ~msb;
//   ns &= ~smsb;

//   if(ns==0){
//     ns = s | msb;
//   } else if(ns==s && first_time){
//     ns = (s+1) | smsb;
//   } else if(ns==s && emergent){
//     ns = s | msb;
//   }

//   // if(--dep(nx,ny)<=0)
//   //   q.emplace(nx,ny);
// }




// template<class AccumF, class E, class A>
// static void KernelSeibertMcGlynn(
//   const FDMode mode,
//   AccumF accumf,
//   const Array2D<E> &elevations,
//   Array2D<A> &accum,
//   dep_t &dep,
//   std::queue<GridCell> &q,
//   const int x,
//   const int y,
//   const double xparam
// ){
//   //TODO: Assumes that the width and height of grid cells are equal and scaled
//   //to 1.
//   constexpr double d1   = 1;
//   constexpr double d2   = 1;
//   constexpr float  dang = std::atan2(d2,d1);

//   const auto nwrap = [](const int n){ return (n==9)?1:n; };

//   constexpr double L1   = 0.5;
//   constexpr double L2   = 0.354; //TODO: More decimal places
//   constexpr double L[9] = {0,L1,L2,L1,L2,L1,L2,L1,L2};

//   std::array<double,9> svals    = {{0,0,0,0,0,0,0,0,0}};
//   std::array<double,9> rvals    = {{0,0,0,0,0,0,0,0,0}};

//   //Table 1 of Tarboton (1997)
//   //          Column #  =   0    1    2    3    4    5   6    7
//   // const int    dy_e1[8] = { 0 , -1 , -1 ,  0 ,  0 ,  1 , 1 ,  0 };
//   // const int    dx_e1[8] = { 1 ,  0 ,  0 , -1 , -1 ,  0 , 0 ,  1 };
//   // const int    dy_e2[8] = {-1 , -1 , -1 , -1 ,  1 ,  1 , 1 ,  1 };
//   // const int    dx_e2[8] = { 1 ,  1 , -1 , -1 , -1 , -1 , 1 ,  1 };
//   // const double ac[8]    = { 0.,  1.,  1.,  2.,  2.,  3., 3.,  4.};
//   // const double af[8]    = { 1., -1.,  1., -1.,  1., -1., 1., -1.};

//   //I remapped the foregoing table for ease of use with RichDEM. The facets
//   //are renumbered as follows:
//   //    3->1    2->2    1->3    0->4    7->5    6->6    5->7    4->8
//   //This gives the following table
//   //  Remapped Facet #  =  -   1    2     3    4    5   6    7    8  
//   //  Tarboton Facet #  =  -   3    2     1    0    7   6    5    4  
//   const int    dy_e1[9] = {0,  0 , -1 ,  -1 ,  0 ,  0 , 1 ,  1 ,  0  };
//   const int    dx_e1[9] = {0, -1 ,  0 ,   0 ,  1 ,  1 , 0 ,  0 , -1  };
//   const int    dy_e2[9] = {0, -1 , -1 ,  -1 , -1 ,  1 , 1 ,  1 ,  1  };
//   const int    dx_e2[9] = {0, -1 , -1 ,   1 ,  1 ,  1 , 1 , -1 , -1  };
//   //const double ac[9]    = {0,  2.,  1.,   1.,  0.,  4., 3.,  3.,  2. };
//   const double af[9]    = {0, -1.,  1.,  -1.,  1., -1., 1., -1.,  1. };

//   for(int n=1;n<=8;n++){
//     if(!elevations.inGrid (x+dx_e1[n],y+dy_e1[n]))
//       continue;
//     if(elevations.isNoData(x+dx_e1[n],y+dy_e1[n]))
//       continue;
//     if(!elevations.inGrid (x+dx_e2[n],y+dy_e2[n]))
//       continue;
//     if(elevations.isNoData(x+dx_e2[n],y+dy_e2[n]))
//       continue;

//     //Is is assumed that cells with a value of NoData have very negative
//     //elevations with the result that they draw flow off of the grid.

//     //Choose elevations based on Table 1 of Tarboton (1997), Barnes TODO
//     const double e0 = elevations(x,y);
//     const double e1 = elevations(x+dx_e1[n],y+dy_e1[n]);
//     const double e2 = elevations(x+dx_e2[n],y+dy_e2[n]);

//     const double s1 = (e0-e1)/d1;
//     const double s2 = (e1-e2)/d2;

//     rvals[n] = std::atan2(s2,s1);

//     if(rvals[n]<1e-7){
//       rvals[n] = 0;
//       svals[n] = s1;
//     } else if(rvals[n]>dang-1e-7){
//       rvals[n] = dang;
//       svals[n] = (e0-e2)/sqrt(d1*d1+d2*d2);
//     } else {
//       svals[n] = sqrt(s1*s1+s2*s2);
//     }

//     if(svals[n]<0){
//       svals[n] = 0;
//       continue;
//     }

//     if(af[n]==1 && rvals[n]==0)
//       rvals[n] = dang;
//     else if(af[n]==1 && rvals[n]==dang)
//       rvals[n] = 0;
//     else if(af[n]==1)
//       rvals[n] = M_PI/4-rvals[n];

//     assert(0<=rvals[n] && rvals[n]<=dang);
//   }

//   double C = 0;
//   for(int n=1;n<=8;n++){
//     svals[n] = std::pow(svals[n]*L[n], xparam);
//     C       += svals[n];
//   }

//   if(C==0)
//     return;

//   assert(C>0);
//   C = accum(x,y)/C;
//   assert(C>0);

//   for(int n=1;n<=8;n++){
//     if(svals[n]<=0)
//       continue;

//     if(rvals[n]==0){
//       accumf(x,y,n,          dep,accum,q, C*svals[n]);
//     } else if(rvals[n]==dang){
//       accumf(x,y,nwrap(n+1), dep,accum,q, C*svals[n]);
//     } else {
//       assert(0<=rvals[n] && rvals[n]<=dang);
//       assert(C>0);
//       assert(svals[n]>0);
//       accumf(x,y,n,          dep,accum,q, C*svals[n]*(  rvals[n]/(M_PI/4.)));
//       accumf(x,y,nwrap(n+1), dep,accum,q, C*svals[n]*(1-rvals[n]/(M_PI/4.)));
//     }
//   }
// }



template<class elev_t>
std::vector<float> FP_Tarboton(const Array2D<elev_t> &elevations){
  RDLOG_ALG_NAME<<"Tarboton (1997) Flow Accumulation (aka D-Infinity, D∞)";
  RDLOG_CITATION<<"Tarboton, D.G., 1997. A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water resources research 33, 309–319.";

  std::vector<float> props(9*elevations.size(),0); //TODO: NO_FLOW

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
  for(int y=1;y<elevations.height()-1;y++)
  for(int x=1;x<elevations.width()-1;x++){
    ++progress;

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
      props.at(9*elevations.xyToI(x,y)+nmax) = 1;
    } else if(rmax==dang){
      props.at(9*elevations.xyToI(x,y)+nwrap(nmax+1)) = 1;
    } else {
      props.at(9*elevations.xyToI(x,y)+nmax)          = rmax/(M_PI/4.);
      props.at(9*elevations.xyToI(x,y)+nwrap(nmax+1)) = 1-rmax/(M_PI/4.);      
    }
  }
  progress.stop();

  return props;
}



template<class E>
std::vector<float> FP_Holmgren(const Array2D<E> &elevations, const double xparam){
  RDLOG_ALG_NAME<<"Holmgren (1994) Flow Accumulation (aka MFD, MD8)";
  RDLOG_CITATION<<"Holmgren, P., 1994. Multiple flow direction algorithms for runoff modelling in grid based elevation models: an empirical evaluation. Hydrological processes 8, 327–334.";
  RDLOG_CONFIG<<"x = "<<xparam;

  std::vector<float> props(9*elevations.size(),0);

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

    C = 1/C;

    for(int n=1;n<=8;n++){
      if(props[9*ci+n]>0)
        props.at(9*ci+n) *= C;
      else
        props.at(9*ci+n) = 0;
    }
  }
  progress.stop();

  return props;
}

template<class E>
std::vector<float> FP_Quinn(const Array2D<E> &elevations){
  RDLOG_ALG_NAME<<"Quinn (1991) Flow Accumulation (aka MFD, MD8)";
  RDLOG_CITATION<<"Quinn, P., Beven, K., Chevallier, P., Planchon, O., 1991. The Prediction Of Hillslope Flow Paths For Distributed Hydrological Modelling Using Digital Terrain Models. Hydrological Processes 5, 59–79."; 
  return FP_Holmgren(elevations, 1.0);
}



template<class E>
std::vector<float> FP_Freeman(
  const Array2D<E> &elevations,
  const double xparam
){
  RDLOG_ALG_NAME<<"Freeman (1991) Flow Accumulation (aka MFD, MD8)";
  RDLOG_CITATION<<"Freeman, T.G., 1991. Calculating catchment area with divergent flow based on a regular grid. Computers & Geosciences 17, 413–422.";
  RDLOG_CONFIG<<"p = "<<xparam;

  std::vector<float> props(9*elevations.size(),0);

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
        props[9*ci+n]     = std::pow(grad,xparam);
        C                += props.at(9*ci+n);
      }
    }

    C = 1/C; //TODO

    for(int n=1;n<=8;n++){
      auto &this_por = props.at(9*ci+n);
      if(this_por>0)
        this_por *= C;
      else
        this_por = 0;
    }
  }
  progress.stop();

  return props;
}

template<class E>
std::vector<float> FP_FairfieldLeymarie(const Array2D<E> &elevations){
  RDLOG_ALG_NAME<<"Fairfield (1991) \"Rho8\" Flow Accumulation";
  RDLOG_CITATION<<"Fairfield, J., Leymarie, P., 1991. Drainage networks from grid digital elevation models. Water resources research 27, 709–717.";

  std::vector<float> props(9*elevations.size(),0);

  ProgressBar progress;
  progress.start(elevations.size());

  #pragma omp parallel for collapse(2)
  for(int y=1;y<elevations.height()-1;y++)
  for(int x=1;x<elevations.width()-1;x++){
    ++progress;

    const int ci = elevations.xyToI(x,y);
    const E e    = elevations(x,y);

    int    greatest_n     = 0; //TODO: Use a constant
    double greatest_slope = 0;
    for(int n=1;n<=8;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];

      if(!elevations.inGrid(nx,ny))
        continue;
      if(elevations.isNoData(nx,ny)) //TODO: Don't I want water to drain this way?
        continue;

      const E ne = elevations(nx,ny);

      if(ne>=e)
        continue;

      double rho_slope = (e-ne);
      if(n_diag[n])
        rho_slope *= 1/(2-uniform_rand_real(0,1));

      if(rho_slope>greatest_slope){
        greatest_n     = n;
        greatest_slope = rho_slope;
      }
    }

    if(greatest_n==0)
      continue;

    props.at(9*ci+greatest_n) = 1;

    assert(elevations(x,y)>=elevations(x+dx[greatest_n],y+dy[greatest_n])); //Ensure flow goes downhill
  }
  progress.stop();

  return props;
}

template<class E>
std::vector<float> FP_Rho8(const Array2D<E> &elevations){
  //Algorithm headers are taken care of in FP_FairfieldLeymarie()
  return FP_FairfieldLeymarie(elevations);
}



//TODO: Add Marks et al (1984)
template<class E>
static std::vector<float> FP_OCallaghan(const Array2D<E> &elevations){
  RDLOG_ALG_NAME<<"O'Callaghan (1984)/Marks (1984) Flow Accumulation (aka D8)";
  RDLOG_CITATION<<"O'Callaghan, J.F., Mark, D.M., 1984. The Extraction of Drainage Networks from Digital Elevation Data. Computer vision, graphics, and image processing 28, 323--344.";

  std::vector<float> props(9*elevations.size(),0);

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

    assert(elevations(ci)>=elevations(ci+elevations.nshift(lowest_n))); //Ensure flow goes downhill

    props.at(9*ci+lowest_n) = 1;
  }
  progress.stop();

  return props;
}


template<class E>
std::vector<float> FP_D8(const Array2D<E> &elevations){
  return FP_OCallaghan(elevations);
}





template<class F, class E, class A, typename... Args>
static void FlowAccumulation(F func, const Array2D<E> &elevations, Array2D<A> &accum, Args... args ){
  Timer overall;
  overall.start();

  RDLOG_ALG_NAME<<"Generic Flow Accumulation Algorithm";

  accum.resize(elevations,1);
  accum.setNoData(ACCUM_NO_DATA);

  const auto props = func(elevations, args...);

  //Create dependencies array
  RDLOG_PROGRESS<<"Creating dependencies array..."<<std::endl;
  Array2D<int8_t> deps(elevations, 0);
  for(int y=1;y<elevations.height()-1;y++)
  for(int x=1;x<elevations.width()-1;x++){
    const int ci = elevations.xyToI(x,y);
    for(int n=1;n<=8;n++)
      if(props[9*ci+n]>0){
        const int ni = ci + elevations.nshift(n);
        deps(ni)++;
      }
  }

  RDLOG_DEBUG<<"Source cells found = "<<deps.size(); //TODO: Switch log target

  //Find sources
  std::queue<int> q;
  for(auto i=deps.i0();i<deps.size();i++)
    if(deps(i)==0 && !elevations.isNoData(i))
      q.emplace(i);

  RDLOG_PROGRESS<<"Calculating flow accumulation...";
  ProgressBar progress;
  progress.start(elevations.size());
  while(!q.empty()){
    ++progress;

    const auto ci = q.front();
    q.pop();

    assert(!elevations.isNoData(ci));

    const auto c_accum = accum(ci);

    for(int n=1;n<=8;n++){
      if(props[9*ci+n]<=0)
        continue;
      const int ni = ci+elevations.nshift(n);
      if(elevations.isNoData(ni))
        continue;
      accum(ni) += props[9*ci+n]*c_accum;
      if(--deps(ni)==0)
        q.emplace(ni);
      assert(deps(ni)>=0);
    }
  }
  progress.stop();

  for(auto i=elevations.i0();i<elevations.size();i++)
    if(elevations.isNoData(i))
      accum(i)=accum.noData();

  RDLOG_TIME_USE<<"Wall-time       = "<<overall.stop()<<" s"     ;
}



template<class elev_t, class accum_t> void FA_Tarboton          (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FP_Tarboton<elev_t>         , elevations, accum); }
template<class elev_t, class accum_t> void FA_Holmgren          (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum, double xparam) {FlowAccumulation(FP_Holmgren<elev_t>         , elevations, accum, xparam); }
template<class elev_t, class accum_t> void FA_Quinn             (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FP_Quinn<elev_t>            , elevations, accum); }
template<class elev_t, class accum_t> void FA_Freeman           (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum, double xparam) {FlowAccumulation(FP_Freeman<elev_t>          , elevations, accum, xparam); }
template<class elev_t, class accum_t> void FA_FairfieldLeymarie (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FP_FairfieldLeymarie<elev_t>, elevations, accum); }
template<class elev_t, class accum_t> void FA_Rho8              (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FP_Rho8<elev_t>             , elevations, accum); }
template<class elev_t, class accum_t> void FA_D8                (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FP_D8<elev_t>               , elevations, accum); }
template<class elev_t, class accum_t> void FA_OCallaghan        (const Array2D<elev_t> &elevations, Array2D<accum_t> &accum)                {FlowAccumulation(FP_D8<elev_t>               , elevations, accum); }










// enum OrlandiniMode {
//   LAD,
//   LTD
// };


// template<class AccumF, class E, class A>
// static void KernelOrlandini(
//   const FDMode mode,
//   AccumF accumf,
//   const Array2D<E> &elevations,
//   Array2D<A> &accum,
//   dep_t &dep,
//   std::queue<GridCell> &q,
//   const int x,
//   const int y,
//   Array2D<double> &delta,
//   const OrlandiniMode omode,
//   const double lambda
// ){
//   //TODO: Assumes that the width and height of grid cells are equal and scaled
//   //to 1.
//   constexpr double d1   = 1;
//   constexpr double d2   = 1;
//   constexpr float  dang = std::atan2(d2,d1);

//   //Table 1 of Tarboton (1997)
//   //          Column #  =   0    1    2    3    4    5   6    7
//   // const int    dy_e1[8] = { 0 , -1 , -1 ,  0 ,  0 ,  1 , 1 ,  0 };
//   // const int    dx_e1[8] = { 1 ,  0 ,  0 , -1 , -1 ,  0 , 0 ,  1 };
//   // const int    dy_e2[8] = {-1 , -1 , -1 , -1 ,  1 ,  1 , 1 ,  1 };
//   // const int    dx_e2[8] = { 1 ,  1 , -1 , -1 , -1 , -1 , 1 ,  1 };

//   //I remapped the foregoing table for ease of use with RichDEM. The facets
//   //are renumbered as follows:
//   //    3->1    2->2    1->3    0->4    7->5    6->6    5->7    4->8
//   //This gives the following table
//   //  Remapped Facet #  =  -   1    2     3    4    5   6    7    8  
//   //  Tarboton Facet #  =  -   3    2     1    0    7   6    5    4  
//   const int    dy_e1[9] = {0,  0 , -1 ,  -1 ,  0 ,  0 , 1 ,  1 ,  0  };
//   const int    dx_e1[9] = {0, -1 ,  0 ,   0 ,  1 ,  1 , 0 ,  0 , -1  };
//   const int    dy_e2[9] = {0, -1 , -1 ,  -1 , -1 ,  1 , 1 ,  1 ,  1  };
//   const int    dx_e2[9] = {0, -1 , -1 ,   1 ,  1 ,  1 , 1 , -1 , -1  };

//   //Table 1 of Orlandini et al (2003)
//   //                        021 023 063 069 089 087 047 041
//   //                p1  =     2   2   6   6   8   8   4   4
//   //                p2  =     1   3   3   9   9   7   7   1
//   //             sigma  =     1  -1   1  -1   1  -1   1  -1
//   //Orlandini et al number neighbour cells like so:
//   //    369
//   //    208
//   //    147
//   //Therefore, facets correspond as follows:
//   //RichDEM   Orlandini
//   //    012 = 023
//   //    023 = 036
//   //    034 = 069
//   //    045 = 089
//   //    056 = 078
//   //    067 = 047
//   //    078 = 014
//   //    081 = 012
//   //We also convert Orlandini's coordinate system to the RichDEM system like so:
//   //    2->1, 3->2, 6->3, 9->4, 8->5, 7->6, 4->7, 1->8
//   //The remapped Table 1 is then
//   //                         023  063  069   089  087   047  041  021 
//   //                p1  =      1    3    3     5    5     7    7    1 
//   //                p2  =      2    2    4     4    6     6    8    8 
//   //             sigma  =     -1    1   -1     1   -1     1   -1    1 
//   const int       p1[9] = {0,  1 ,  3 ,  3 ,   5,   5,    7,   7,   1  };
//   const int       p2[9] = {0,  2 ,  2 ,  4 ,   4,   6,    6,   8,   8  };
//   const int    sigma[9] = {0,  -1,  1 , -1 ,   1,  -1,    1,  -1,   1  };

//   int8_t nmax = -1;
//   double smax = 0;
//   float  rmax = 0;

//   for(int n=1;n<=8;n++){
//     if(!elevations.inGrid (x+dx_e1[n],y+dy_e1[n]))
//       continue;
//     if(elevations.isNoData(x+dx_e1[n],y+dy_e1[n]))
//       continue;
//     if(!elevations.inGrid (x+dx_e2[n],y+dy_e2[n]))
//       continue;
//     if(elevations.isNoData(x+dx_e2[n],y+dy_e2[n]))
//       continue;

//     //Is is assumed that cells with a value of NoData have very negative
//     //elevations with the result that they draw flow off of the grid.

//     //Choose elevations based on Table 1 of Tarboton (1997), Barnes TODO
//     const double e0 = elevations(x,y);
//     const double e1 = elevations(x+dx_e1[n],y+dy_e1[n]);
//     const double e2 = elevations(x+dx_e2[n],y+dy_e2[n]);

//     const double s1 = (e0-e1)/d1;
//     const double s2 = (e1-e2)/d2;

//     double r = std::atan2(s2,s1);
//     double s;

//     if(r<1e-7){
//       r = 0;
//       s = s1;
//     } else if(r>dang-1e-7){
//       r = dang;
//       s = (e0-e2)/sqrt(d1*d1+d2*d2);
//     } else {
//       s = sqrt(s1*s1+s2*s2);
//     }

//     if(s>smax){
//       smax = s;
//       nmax = n;
//       rmax = r;
//     }
//   }

//   if(nmax==-1)
//     return;

//   double d1me;
//   double d2me;

//   if(omode==OrlandiniMode::LTD){
//     d1me = lambda*delta(x,y) + sigma[nmax]*d1*std::sin(rmax);
//     d2me = lambda*delta(x,y) - sigma[nmax]*std::sqrt(d1*d1+d2*d2)*std::sin(dang-rmax);
//   } else if(omode==OrlandiniMode::LAD){
//     d1me = lambda*delta(x,y) + sigma[nmax]*rmax;
//     d2me = lambda*delta(x,y) - sigma[nmax]*(dang-rmax);
//   }

//   int p;

//   if(std::abs(d1me)<=std::abs(d2me)){
//     delta(x,y) = d1me;
//     p = p1[nmax];
//   } else {
//     delta(x,y) = d2me;
//     p = p2[nmax];
//   }

//   //TODO: This fixes issues, but makes me think I have one or more constants wrong
//   if(elevations(x,y)<elevations(x+dx[p],y+dy[p])){
//     if(std::abs(d1me)<=std::abs(d2me)){
//       delta(x,y) = d2me;
//       p = p2[nmax];
//     } else {
//       delta(x,y) = d1me;
//       p = p1[nmax];
//     }    
//   }

//   const int nx = x+dx[p];
//   const int ny = y+dy[p];

//   if(elevations(x,y)<elevations(nx,ny)){
//     RDLOG_DEBUG<<"Uh oh";
//     RDLOG_DEBUG<<"p="<<(int)p;
//     RDLOG_DEBUG<<"nmax="<<(int)nmax;
//     for(int iy=std::max(0,y-1);iy<=std::min(elevations.height()-1,y+1);iy++){
//       for(int ix=std::max(0,x-1);ix<=std::min(elevations.width()-1,x+1);ix++)
//         RDLOG_DEBUG<<std::setprecision(15)<<elevations(ix,iy)<<" ";
//       RDLOG_DEBUG;
//     }
//   }
//   RDLOG_DEBUG<<std::flush;
//   assert(elevations(x,y)>=elevations(nx,ny)); //Ensure flow goes downhill

//   delta(nx,ny) = delta(x,y);

//   accumf(x,y,p,          dep,accum,q, accum(x,y));
// }








/*


template<class F, class E, class A, typename... Args>
static void DistanceDispersionEstimate(
  F f,
  const Array2D<E> &elevations,
  Array2D<A> &accum,
  Args&&... args
){
  ProgressBar progress;
  std::queue<GridCell> q;
  dep_t dep,depbk;

  Timer overall;
  overall.start();

  accum.setAll(0);
  accum.setNoData(ACCUM_NO_DATA);
  for(typename Array2D<E>::i_t i=0;i<elevations.size();i++)
    if(elevations.isNoData(i))
      accum(i) = ACCUM_NO_DATA;

  dep.resize(elevations);
  dep.setAll(0);

  RDLOG_PROGRESS<<"Calculating dependencies...";
  progress.start(elevations.size());
  //#pragma omp parallel for collapse(2)
  for(int y=0;y<elevations.height();y++)
  for(int x=0;x<elevations.width();x++){
    ++progress;
    if(!elevations.isNoData(x,y))
      f(FDMode::CALC_DEPENDENCIES,elevations,accum,dep,q,x,y,std::forward<Args>(args)...);
  }
  progress.stop();

  std::vector<int>  dispersion(elevations.size(),0);
  std::vector<int>  count     (elevations.size(),0);

  progress.start(2000);
  for(int i=0;i<2000;i++){
    ++progress;

    q.empty();

    std::vector<bool> visited   (elevations.size(),false);
    
    int rx,ry;
    do {
      rx = uniform_rand_int(0,elevations.width ()-1);
      ry = uniform_rand_int(0,elevations.height()-1);
    } while(elevations.isNoData(rx,ry));
    q.emplace(rx,ry);

    dep.setAll(0);

    int cells_in_level = 1;
    int level          = 0;

    while(!q.empty()){
      ++progress;

      const auto c = q.front();
      q.pop();

      if(visited[elevations.xyToI(c.x,c.y)])
        continue;
      visited[elevations.xyToI(c.x,c.y)] = true;
      dispersion[level]++;


      accum(c.x,c.y) += 1;
      f(FDMode::CALC_ACCUM,elevations,accum,dep,q,c.x,c.y,std::forward<Args>(args)...);

      cells_in_level--;
      if(cells_in_level==0){
        level++;
        cells_in_level = q.size();
      }
    }

    for(int l=0;l<level;l++)
      count[l]++;
  }
  progress.stop();

  std::cout<<"Lvl,Dispersion,Count\n";
  for(int i=0;i<500;i++)
    std::cout<<i<<","<<dispersion[i]<<","<<count[i]<<"\n";
}


*/



// template<class E, class A>
// void FP_SeibertMcGlynn(const Array2D<E> &elevations, Array2D<A> &accum, double x){
//   RDLOG_ALG_NAME<<"Seibert and McGlynn (2007) Flow Accumulation (aka MD-Infinity, MD∞)";
//   RDLOG_WARN<<"TODO: This flow accumulation method is not yet functional.";
//   RDLOG_CONFIG<<"x = "<<x;
//   KernelFlowdir(KernelSeibertMcGlynn<decltype(PassAccumulation<A>),E,A>,PassAccumulation<A>,elevations,accum,x);
// }

// template<class E, class A>
// void FP_Orlandini(const Array2D<E> &elevations, Array2D<A> &accum, OrlandiniMode mode, double lambda){
//   RDLOG_ALG_NAME<<"Orlandini et al. (2003) Flow Accumulation (aka D8-LTD, D8-LAD)";
//   RDLOG_CITATION<<"Orlandini, S., Moretti, G., Franchini, M., Aldighieri, B., Testa, B., 2003. Path-based methods for the determination of nondispersive drainage directions in grid-based digital elevation models: TECHNICAL NOTE. Water Resources Research 39(6). doi:10.1029/2002WR001639.";
//   RDLOG_CONFIG<<"lambda = "<<lambda;
//   Array2D<double> delta(elevations,0);
//   KernelFlowdir(KernelOrlandini<decltype(PassAccumulation<A>),E,A>,PassAccumulation<A>,elevations,accum,delta,mode,lambda);
// }






/*


template<class A>
static inline void CleanseStrahler(Array2D<A> &accum){
  constexpr A msb  = 1 << (std::numeric_limits<A>::digits-1); //First time this Strahler # was seen
  constexpr A smsb = 1 << (std::numeric_limits<A>::digits-2); //This Strahler number arose from combining two other neighbours and should be overriden if appropriate
  for(typename Array2D<A>::i_t i=0;i<accum.size();i++)
    accum(i) = accum(i) & (~msb) & (~smsb);
}


template<class E, class A>
void Strahler_FairfieldLeymarie(const Array2D<E> &elevations, Array2D<A> &accum){
  RDLOG_ALG_NAME<<"Fairfield (1991) \"Rho8\" Strahler";
  RDLOG_CITATION<<"Fairfield, J., Leymarie, P., 1991. Drainage networks from grid digital elevation models. Water resources research 27, 709–717.";
  Array2D<d8_flowdir_t> fd(elevations);
  KernelFlowdir(KernelFairfieldLeymarie<decltype(StrahlerNumber<A>),E,A>,StrahlerNumber<A>,elevations,accum,fd);
  CleanseStrahler(accum);
}

template<class E, class A>
void Strahler_Rho8(const Array2D<E> &elevations, Array2D<A> &accum){
  //Algorithm headers are taken care of in Strahler_FairfieldLeymarie()
  Strahler_FairfieldLeymarie(elevations,accum);
}

template<class E, class A>
void Strahler_Quinn(const Array2D<E> &elevations, Array2D<A> &accum){
  RDLOG_ALG_NAME<<"Quinn (1991) Strahler";
  RDLOG_CITATION<<"Quinn, P., Beven, K., Chevallier, P., Planchon, O., 1991. The Prediction Of Hillslope Flow Paths For Distributed Hydrological Modelling Using Digital Terrain Models. Hydrological Processes 5, 59–79."; 
  KernelFlowdir(KernelHolmgren<decltype(StrahlerNumber<A>),E,A>,StrahlerNumber<A>,elevations,accum,(double)1.0);
  CleanseStrahler(accum);
}

template<class E, class A>
void Strahler_Holmgren(const Array2D<E> &elevations, Array2D<A> &accum, double x){
  RDLOG_ALG_NAME<<"Holmgren (1994) Strahler";
  RDLOG_CITATION<<"Holmgren, P., 1994. Multiple flow direction algorithms for runoff modelling in grid based elevation models: an empirical evaluation. Hydrological processes 8, 327–334.";
  KernelFlowdir(KernelHolmgren<decltype(StrahlerNumber<A>),E,A>,StrahlerNumber<A>,elevations,accum,x);
  CleanseStrahler(accum);
}

template<class E, class A>
void Strahler_Tarboton(const Array2D<E> &elevations, Array2D<A> &accum){
  RDLOG_ALG_NAME<<"Tarboton (1997) \"D-Infinity\" Strahler";
  RDLOG_CITATION<<"Tarboton, D.G., 1997. A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water resources research 33, 309–319.";
  Array2D< std::pair<float,int8_t> > fd(elevations);
  KernelFlowdir(KernelTarboton<decltype(StrahlerNumber<A>),E,A>,StrahlerNumber<A>,elevations,accum,fd);
  CleanseStrahler(accum);
}

// template<class E, class A>
// void Strahler_SeibertMcGlynn(const Array2D<E> &elevations, Array2D<A> &accum, double xparam){
//   RDLOG_ALG_NAME<<"Seibert and McGlynn Strahler (TODO)";
//   RDLOG_WARN<<"TODO: This flow accumulation method is not yet functional.";
//   KernelFlowdir(KernelSeibertMcGlynn<decltype(StrahlerNumber<A>),E,A>,StrahlerNumber<A>,elevations,accum,xparam);
//   CleanseStrahler(accum);
// }

*/



}

#endif
