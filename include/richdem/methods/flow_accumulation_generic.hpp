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
template<class A>
void FlowAccumulation(const Array3D<float> &props, Array2D<A> &accum){
  Timer overall;
  overall.start();

  RDLOG_ALG_NAME<<"Generic Flow Accumulation Algorithm";

  accum.setNoData(ACCUM_NO_DATA);

  if(accum.width()!=props.width() || accum.height()!=props.height())
    throw std::runtime_error("Accumulation array must have same dimensions as proportions array!");

  //Create dependencies array
  RDLOG_PROGRESS<<"Creating dependencies array..."<<std::endl;
  Array2D<int8_t> deps(props, 0);
  for(int y=1;y<props.height()-1;y++)
  for(int x=1;x<props.width()-1;x++){
    const int ci = accum.xyToI(x,y);
    if(props.isNoData(ci))
      continue;
    for(int n=1;n<=8;n++)
      if(props(x,y,n)>0){
        const int ni = ci + accum.nshift(n);
        deps(ni)++;
      }
  }

  //Find sources
  std::queue<int> q;
  for(auto i=deps.i0();i<deps.size();i++)
    if(deps(i)==0 && !props.isNoData(i))
      q.emplace(i);

  RDLOG_DEBUG<<"Source cells found = "<<q.size(); //TODO: Switch log target

  RDLOG_PROGRESS<<"Calculating flow accumulation...";
  ProgressBar progress;
  progress.start(props.size());
  while(!q.empty()){
    ++progress;

    const auto ci = q.front();
    q.pop();

    assert(!props.isNoData(ci));

    const auto c_accum = accum(ci);

    for(int n=1;n<=8;n++){
      if(props.getIN(ci,n)<=0) //No Flow in this direction or other flags
        continue;
      const int ni = ci+accum.nshift(n);
      if(props.isNoData(ni))
        continue;
      accum(ni) += props.getIN(ci,n)*c_accum;
      if(--deps(ni)==0)
        q.emplace(ni);
      assert(deps(ni)>=0);
    }
  }
  progress.stop();

  for(auto i=props.i0();i<props.size();i++)
    if(props.isNoData(i))
      accum(i) = accum.noData();

  RDLOG_TIME_USE<<"Wall-time       = "<<overall.stop()<<" s"     ;
}

}

#endif
