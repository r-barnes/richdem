/**
  @file
  @brief Defines the D-infinite flow routing method described by Tarboton (1997)

  This file implements the D-infinite flow routing method originally described by
Tarboton (1997). It incorporates minor alterations and additional safe-guards
described in Barnes (TODO).

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_dinf_flowdirs_hpp_
#define _richdem_dinf_flowdirs_hpp_

#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/ProgressBar.hpp"

namespace richdem {

///Value used to indicate that a flow direction cell has no data
#define dinf_NO_DATA -1

//Table 1 of Tarboton (1997), Barnes TODO
//              Column #  =   0    1    2    3    4    5   6    7
static const int dy_e1[8] = { 0 , -1 , -1 ,  0 ,  0 ,  1 , 1 ,  0 };
static const int dx_e1[8] = { 1 ,  0 ,  0 , -1 , -1 ,  0 , 0 ,  1 };
static const int dy_e2[8] = {-1 , -1 , -1 , -1 ,  1 ,  1 , 1 ,  1 };
static const int dx_e2[8] = { 1 ,  1 , -1 , -1 , -1 , -1 , 1 ,  1 };
static const double ac[8] = { 0.,  1.,  1.,  2.,  2.,  3., 3.,  4.};
static const double af[8] = { 1., -1.,  1., -1.,  1., -1., 1., -1.};

/**
  @brief  Determine the D-infinite flow direction of a cell
  @author Implementation by Richard Barnes (rbarnes@umn.edu)

    This function determines the D-infinite flow direction of a cell, as
    described by Tarboton (1997) and Barnes (2013, TODO). TODO

  @param[in] elevations   A 2D grid of elevation data
  @param[in] x            x-coordinate of cell to determine flow direction for
  @param[in] y            y-coordinate of cell to determine flow direction for

  @return A floating-point value between [0,2*Pi) indicating flow direction
*/
template <class T>
static float dinf_FlowDir(const Array2D<T> &elevations, const int x, const int y){
  //Ensure that flow is pulled off the edge of the grid
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

  int    nmax = -1;
  double smax = 0;
  double rmax = 0;

  //I am not on the edge of the grid. All my neighbours can be examined.

  for(int n=0;n<8;n++){
    //Is is assumed that cells with a value of NoData have very negative
    //elevations with the result that they draw flow off of the grid.

    //Choose elevations based on Table 1 of Tarboton (1997), Barnes TODO
    const double e0 = elevations(x,y);
    const double e1 = elevations(x+dx_e1[n],y+dy_e1[n]);
    const double e2 = elevations(x+dx_e2[n],y+dy_e2[n]);

    //TODO: Assumes that the width and height of grid cells are equal and scaled
    //to 1.
    const double d1 = 1;
    const double d2 = 1;

    const double s1 = (e0-e1)/d1;
    const double s2 = (e1-e2)/d2;
    double r        = atan2(s2,s1);

    double s;

    if(r<0){
      r = 0;
      s = s1;
    } else if(r>atan2(d2,d1)){
      r = atan2(d2,d1); //TODO: This is a constant
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

  double rg = NO_FLOW;
  if(nmax!=-1)
    rg = (af[nmax]*rmax+ac[nmax]*M_PI/2);

  return rg;
}


/**
  @brief  Determine the D-infinite flow direction of every cell in a grid
  @author Richard Barnes (rbarnes@umn.edu)

    This function runs dinf_FlowDir() on every cell in a grid which has a data
    value.

  @param[in]  &elevations  A 2D grid of elevation data
  @param[out] &flowdirs    A 2D grid which will contain the flow directions
*/
template <class T>
void dinf_flow_directions(const Array2D<T> &elevations, Array2D<float> &flowdirs){
  ProgressBar progress;

  RDLOG_ALG_NAME<<"Dinf Flow Directions";
  RDLOG_CITATION<<"Tarboton, D.G. 1997. A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water Resources Research. Vol. 33. pp 309-319.";

  RDLOG_PROGRESS<<"Setting up the Dinf flow directions matrix...";
  flowdirs.resize(elevations);
  flowdirs.setNoData(dinf_NO_DATA);
  flowdirs.setAll(NO_FLOW);

  RDLOG_PROGRESS<<"Calculating Dinf flow directions...";
  progress.start( elevations.size() );
  #pragma omp parallel for
  for(int y=0;y<elevations.height();y++){
    progress.update( y*elevations.width() );
    for(int x=0;x<elevations.width();x++)
      if(elevations(x,y)==elevations.noData())
        flowdirs(x,y) = flowdirs.noData();
      else
        flowdirs(x,y) = dinf_FlowDir(elevations,x,y);
  }
  RDLOG_TIME_USE<<"Succeeded in = "<<progress.stop()<<" s";
}

}

#endif
