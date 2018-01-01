/**
  @file
  @brief Couples the Barnes (2014) flat resolution algorithm with the Tarboton (1997) D-infinity flow metric
  @author Richard Barnes
*/
#ifndef _richdem_flat_resolution_dinf_hpp_
#define _richdem_flat_resolution_dinf_hpp_

#include "richdem/flats/flat_resolution.hpp"
#include "richdem/flowmet/dinf_flowdirs.hpp"
#include "richdem/common/logger.hpp"

namespace richdem {

static const float d8_to_dinf[9]={-1, 4*M_PI/4, 3*M_PI/4, 2*M_PI/4, 1*M_PI/4, 0, 7*M_PI/4, 6*M_PI/4, 5*M_PI/4};

static float dinf_masked_FlowDir(
  const Array2D<int32_t> &flat_resolution_mask,
  const Array2D<int32_t> &groups,
  const int x,
  const int y
){
  double smax = 0;
  int    nmax = -1;
  double rmax = 0;

  double e0,e1,e2,d1,d2,s1,s2,r,s;

  //Yes, this should be 0-8, this is the Tarboton neighbour system
  for(int n=0;n<8;n++){
    //TODO: Can these ever give !IN_GRID errors?
    if(groups(x+dx_e1[n],y+dy_e1[n])!=groups(x,y)) continue;
    if(groups(x+dx_e2[n],y+dy_e2[n])!=groups(x,y)) continue;

    e0 = flat_resolution_mask(x,y);
    e1 = flat_resolution_mask(x+dx_e1[n],y+dy_e1[n]);
    e2 = flat_resolution_mask(x+dx_e2[n],y+dy_e2[n]);
    d1 = 1;
    d2 = 1;
    s1 = (e0-e1)/d1;
    s2 = (e1-e2)/d2;
    r  = atan2(s2,s1);
    s  = sqrt(s1*s1+s2*s2);
    if(r<0){
      r = 0;
      s = s1;
    } else if(r>atan2(d2,d1)){
      r = atan2(d2,d1);
      s = (e0-e2)/sqrt(d1*d1+d2*d2);
    }
    if(s>smax){
      smax = s;
      nmax = n;
      rmax = r;
    }
  }

  double rg=NO_FLOW;
  if(nmax!=-1){
    rg=(af[nmax]*rmax+ac[nmax]*M_PI/2);
  } else {
    for(int n=1;n<=8;n++){  //TODO: I have a feeling this is potentially unsafe as it may create dependency loops. Does it? TODO: Switch this to dinf_dx
      if(groups(x+dx[n],y+dy[n])==groups(x,y) && flat_resolution_mask(x+dx[n],y+dy[n])<flat_resolution_mask(x,y)){
        rg=d8_to_dinf[n];
        break;
      }
    }
  }

  return rg;
}

void dinf_flow_flats(
  const Array2D<int32_t> &flat_resolution_mask, 
  const Array2D<int32_t> &groups, 
  Array2D<float> &flowdirs
){
  ProgressBar progress;

  RDLOG_ALG_NAME<<"Dinf Flow Flats";
  RDLOG_CITATION<<"TODO";

  RDLOG_PROGRESS<<"Calculating Dinf flow directions using flat mask...";
  progress.start( flat_resolution_mask.width()*flat_resolution_mask.height() );
  #pragma omp parallel for
  for(int x=1;x<flat_resolution_mask.width()-1;x++){
    progress.update( x*flat_resolution_mask.height() );
    for(int y=1;y<flat_resolution_mask.height()-1;y++)
      if(flat_resolution_mask(x,y)==flat_resolution_mask.noData())
        continue;
      else if(flowdirs(x,y)==NO_FLOW)
        flowdirs(x,y)=dinf_masked_FlowDir(flat_resolution_mask,groups,x,y);
  }
  RDLOG_TIME_USE<<"Succeeded in "<<progress.stop()<<" s";
}

template<class T>
void resolve_flats_barnes_dinf(
  const Array2D<T> &elevations,
  Array2D<float> &flowdirs
){
  dinf_flow_directions(elevations,flowdirs);

  Array2D<int32_t> flat_mask,labels;
  resolve_flats_barnes(elevations,flowdirs,flat_mask,labels);

  dinf_flow_flats(flat_mask,labels,flowdirs);
}

}

#endif
