#ifndef _richdem_Seibert2007_hpp_
#define _richdem_Seibert2007_hpp_

namespace richdem {

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

// template<class E, class A>
// void FM_SeibertMcGlynn(const Array2D<E> &elevations, Array2D<A> &accum, double x){
//   RDLOG_ALG_NAME<<"Seibert and McGlynn (2007) Flow Accumulation (aka MD-Infinity, MD∞)";
//   RDLOG_WARN<<"TODO: This flow accumulation method is not yet functional.";
//   RDLOG_CONFIG<<"x = "<<x;
//   KernelFlowdir(KernelSeibertMcGlynn<decltype(PassAccumulation<A>),E,A>,PassAccumulation<A>,elevations,accum,x);
// }

}

#endif
