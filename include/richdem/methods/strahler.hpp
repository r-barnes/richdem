#ifndef _richdem_strahler_hpp_
#define _richdem_strahler_hpp_

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

#endif
