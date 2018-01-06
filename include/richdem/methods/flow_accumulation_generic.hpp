#ifndef _flow_accumulatin_generic_
#define _flow_accumulatin_generic_

#include <richdem/common/Array2D.hpp>
#include <richdem/common/logger.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <queue>

namespace richdem {


/**
  @brief  Calculate flow accumulation from a flow metric array
  @author Richard Barnes (rbarnes@umn.edu)

  Given a flow metric function \p func, this calculations the flow accumulation.
   
  @param[in]     func         The flow metric to use
  @param[in]     &elevations  An elevation field
  @param[in,out] &accum       Accumulation matrix: must be already initialized
  @param[in]     args         Arguments passed to the flow metric (e.g. exponent)

  @pre
    1. The accumulation matrix must already be initialized to the amount of flow
       each cell will generate. A good default value is 1, in which case the
       accumulation matrix will be modified to show how many cells' flow
       ultimately passes through each cell.

  @post
    1. \p accum is modified so that each cell indicates how much upstrema flow
       passes through it (in addition to flow generated within the cell itself).
*/
template<class F, class E, class A, typename... Args>
void FlowAccumulation(F func, const Array2D<E> &elevations, Array2D<A> &accum, Args... args ){
  Timer overall;
  overall.start();

  RDLOG_ALG_NAME<<"Generic Flow Accumulation Algorithm";

  accum.setNoData(ACCUM_NO_DATA);

  if(accum.width()!=elevations.width() || accum.height()!=elevations.height())
    throw std::runtime_error("Accumulation array must have same dimensions as elevations!");

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

}

#endif
