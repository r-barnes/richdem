#ifndef _richdem_dinf_flowdirs_hpp_
#define _richdem_dinf_flowdirs_hpp_

#include "richdem/common/Array2D.hpp"

#define dinf_NO_DATA 113

//Table 1 of Tarboton
static const int    dy_e1[8] = {0,-1,-1,0,0,1,1,0};
static const int    dx_e1[8] = {1,0,0,-1,-1,0,0,1};
static const int    dy_e2[8] = {-1,-1,-1,-1,1,1,1,1};
static const int    dx_e2[8] = {1,1,-1,-1,-1,-1,1,1};
static const double ac   [8] = {0.,1.,1.,2.,2.,3.,3.,4.};
static const double af   [8] = {1.,-1.,1.,-1.,1.,-1.,1.,-1.};

template <class T>
float dinf_FlowDir(const Array2D<T> &elevations, const int x, const int y){
  double smax = 0;
  int    nmax = -1;
  double rmax = 0;

  double e0,e1,e2,d1,d2,s1,s2,r,s;

  if (elevations.isEdgeCell(x,y)){
    if(x==0 && y==0)
      return 3*M_PI/4;  //D8: 2
    else if(x==0 && y==elevations.height()-1)
      return 5*M_PI/4;  //D8: 8
    else if(x==elevations.width()-1 && y==0)
      return 1*M_PI/4;  //D8: 4
    else if(x==elevations.width()-1 && y==elevations.height()-1)
      return 7*M_PI/4;  //D8: 6
    else if(x==0)
      return 4*M_PI/4;  //D8: 1
    else if(x==elevations.width()-1)
      return 0*M_PI/4;  //D8: 5
    else if(y==0)
      return 2*M_PI/4;  //D8: 3
    else if(y==elevations.height()-1)
      return 6*M_PI/4;  //D8: 7
  }
  
  //Since I am not on the edge of the grid if I've made it this far, may neighbours cannot be off the grid
  //Yes, this should be 0-8, this is the Tarboton neighbour system
  for(int n=0;n<8;n++){
    //Very negative no_data's should be acceptable, and suck water of the grid.
    //if(elevations(x+dx_e1[n],y+dy_e1[n])==elevations.no_data) continue;
    //if(elevations(x+dx_e2[n],y+dy_e2[n])==elevations.no_data) continue;
    //Therefore, these lines are not really necessary.
    //I leave them here to make it very clear that they are not necessary.

    e0 = elevations(x,y);
    e1 = elevations(x+dx_e1[n],y+dy_e1[n]);
    e2 = elevations(x+dx_e2[n],y+dy_e2[n]);
    d1 = 1;
    d2 = 1;
    s1 = (e0-e1)/d1;
    s2 = (e1-e2)/d2;
    r  = atan2(s2,s1);

    if(r<0){
      r = 0;
      s = s1;
    } else if(r>atan2(d2,d1)){
      r = atan2(d2,d1);
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

  double rg=NO_FLOW;
  if(nmax!=-1)
    rg=(af[nmax]*rmax+ac[nmax]*M_PI/2);

  return rg;
}

template <class T>
void dinf_flow_directions(const Array2D<T> &elevations, Array2D<float> &flowdirs){
  ProgressBar progress;

  std::cerr<<"\n###Dinf Flow Directions"<<std::endl;

  std::cerr<<"Setting up the Dinf flow directions matrix..."<<std::flush;
  flowdirs.resize(elevations);
  flowdirs.setNoData(dinf_NO_DATA);
  flowdirs.setAll(NO_FLOW);
  std::cerr<<"succeeded.\n"<<std::endl;

  std::cerr<<"%%Calculating Dinf flow directions..."<<std::flush;
  progress.start( elevations.width()*elevations.height() );
  #pragma omp parallel for
  for(int x=0;x<elevations.width();x++){
    progress.update( x*elevations.height() );
    for(int y=0;y<elevations.height();y++)
      if(elevations(x,y)==elevations.noData())
        flowdirs(x,y) = flowdirs.noData();
      else
        flowdirs(x,y) = dinf_FlowDir(elevations,x,y);
  }
  std::cerr<<"succeeded in "<<progress.stop()<<"s."<<std::endl;
}

#endif
