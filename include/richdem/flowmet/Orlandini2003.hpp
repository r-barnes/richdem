#ifndef _richdem_Orlandini2003_hpp_
#define _richdem_Orlandini2003_hpp_

namespace richdem {

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

// template<class E, class A>
// void FM_Orlandini(const Array2D<E> &elevations, Array2D<A> &accum, OrlandiniMode mode, double lambda){
//   RDLOG_ALG_NAME<<"Orlandini et al. (2003) Flow Accumulation (aka D8-LTD, D8-LAD)";
//   RDLOG_CITATION<<"Orlandini, S., Moretti, G., Franchini, M., Aldighieri, B., Testa, B., 2003. Path-based methods for the determination of nondispersive drainage directions in grid-based digital elevation models: TECHNICAL NOTE. Water Resources Research 39(6). doi:10.1029/2002WR001639.";
//   RDLOG_CONFIG<<"lambda = "<<lambda;
//   Array2D<double> delta(elevations,0);
//   KernelFlowdir(KernelOrlandini<decltype(PassAccumulation<A>),E,A>,PassAccumulation<A>,elevations,accum,delta,mode,lambda);
// }

}

#endif
