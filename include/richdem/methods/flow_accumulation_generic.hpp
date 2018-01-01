#ifndef _flow_accumulatin_generic_
#define _flow_accumulatin_generic_

#include <richdem/common/Array2D.hpp>
#include <richdem/common/logger.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <queue>

namespace richdem {

template<class F, class E, class A, typename... Args>
void FlowAccumulation(F func, const Array2D<E> &elevations, Array2D<A> &accum, Args... args ){
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

}

#endif
