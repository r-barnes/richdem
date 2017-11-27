/**
  @file
  @brief Terrain attributes that can only be calculated with Tarboton's D-infinity flow metric
  @author Richard Barnes (rbarnes@umn.edu), 2015

  This file implements the D-infinite flow routing method originally described by
  Tarboton (1997). It incorporates minor alterations and additional safe-guards
  described in Barnes (2013, TODO).
*/
#ifndef _richdem_dinf_methods_hpp_
#define _richdem_dinf_methods_hpp_

#include <cmath>
#include <queue>
#include "richdem/common/logger.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/constants.hpp"
#include "richdem/common/ProgressBar.hpp"
#include "richdem/common/grid_cell.hpp"

namespace richdem {

//TODO: Can these be merged with the regular D8 directions?
//X- and Y-offests of D-inf neighbours (TODO: More explanation, and why there are 9)
static const int dinf_dx[9] = {1,  1,  0, -1, -1, -1, 0, 1, 1};
static const int dinf_dy[9] = {0, -1, -1, -1,  0,  1, 1, 1, 0};

/*
We must convert the Dinf angle system to cells within the D8 system
I use the following grid for the D8 system
//234
//105
//876
To convert Dinf to this, take
(int)(flowdir/45)      D8
      0                4,5
      1                3,4
      2                2,3
      3                1,2
      4                1,8
      5                7,8
      6                6,7
      7                5,6
*/
//These arrays have a 9th element which repeats the 8th element because floating
//point rounding errors occassionally result in the 9th element being accessed.




//321
//4 0
//567
static void where_do_i_flow(float flowdir, int &nhigh, int &nlow){
  //If it is very close to being directed into just one cell
  //then we direct it into just one cell. If we mistakenly direct
  //it into 2 cells, then we may create unresolvable loops in the
  //flow accumulation algorithm, whereas nothing is lost if we
  //leave out one of the two cells (provided a negligible amount
  //of flow is directed to the one we leave out).
  assert(flowdir>=0 && flowdir<=2*M_PI+1e-6);

  flowdir /= (M_PI/4.);

  if(std::abs(flowdir-(int)flowdir)<1e-6){
    nlow  = -1;
    nhigh = (int)std::round(flowdir);
  } else {
    nlow  = (int)flowdir;
    nhigh = nlow+1;
  }

  //8 is not technically a direction, but, since things move in a circle,
  //it overlaps with 0. It should _never_ be greater than 8.
  assert(nhigh>=0 && nhigh<=8);
}

//This reacts correctly if the flow direction wedge number exceeds 7.
static void area_proportion(float flowdir, int nhigh, int nlow, float &phigh, float &plow){
  if(nlow==-1){
    phigh = 1;
    plow  = 0;
  } else {
    phigh = (nhigh*(M_PI/4.0)-flowdir)/(M_PI/4.0);
    plow  = 1-phigh;
  }

  assert(phigh+plow==1);  //TODO: This isn't necessarily so in floating-point... or is it?
}

/*//TODO: Debugging code used for checking for loops. Since loops should not occur in the output of the production code, this is not needed.
bool is_loop(const float_2d &flowdirs, int n, int x, int y, int c2x, int c2y){
  int nh,nl;
  if(! flowdirs.inGrid(c2x, c2y) || flowdirs(c2x,c2y)==flowdirs.noData() || flowdirs(c2x,c2y)==NO_FLOW)
    return false;
  where_do_i_flow(flowdirs(c2x,c2y),nh,nl);
  if(n==dinf_d8_inverse[nh] || (nl!=-1 && n==dinf_d8_inverse[nl])){
    printf("Beware dir %d (%d and %d).\n",n,nh,nl);
    flowdirs.surroundings(x,y,8);
    return true;
  }
  return false;
}*/



/**
  @brief  Calculate each cell's D-infinity flow accumulation value
  @author Tarboton (1997), Richard Barnes (rbarnes@umn.edu)

    TODO

  @param[in]  flowdirs   A grid of D-infinite flow directions
  @param[out] &area      A grid of flow accumulation values
*/
template <class T, class U>
void dinf_upslope_area(
  const Array2D<T> &flowdirs,
  Array2D<U> &area
){
  Array2D<int8_t> dependency;
  std::queue<GridCell> sources;
  ProgressBar progress;

  RDLOG_ALG_NAME<<"D-infinity Upslope Area";
  RDLOG_CITATION<<"Tarboton, D.G. 1997. A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water Resources Research. Vol. 33. pp 309-319.";

  RDLOG_PROGRESS<<"Setting up the dependency matrix...";
  dependency.resize(flowdirs);
  dependency.setAll(0);

  RDLOG_PROGRESS<<"Setting up the area matrix...";
  area.resize(flowdirs);
  area.setAll(0);
  area.setNoData(dinf_NO_DATA);

  bool has_cells_without_flow_directions=false;
  RDLOG_PROGRESS<<"Calculating dependency matrix & setting noData() cells...";
  progress.start( flowdirs.size() );

  ///////////////////////
  //Calculate the number of "dependencies" each cell has. That is, count the
  //number of cells which flow into each cell.

  #pragma omp parallel for reduction(|:has_cells_without_flow_directions)
  for(int y=0;y<flowdirs.height();y++){
    progress.update( y*flowdirs.width() );
    for(int x=0;x<flowdirs.width();x++){
      //If the flow direction of the cell is NoData, mark its area as NoData
      if(flowdirs.isNoData(x,y)){
        area(x,y)       = area.noData();
        dependency(x,y) = 9;  //TODO: This is an unnecessary safety precaution. This prevents the cell from ever being enqueued (an unnecessary safe guard? TODO)
        continue;             //Only necessary if there are bugs below (TODO)
      }

      //If the cell has no flow direction, note that so we can warn the user
      if(flowdirs(x,y)==NO_FLOW){
        has_cells_without_flow_directions=true;
        continue;
      }

      //TODO: More explanation of what's going on here
      int n_high, n_low;
      int nhx,nhy,nlx,nly;
      where_do_i_flow(flowdirs(x,y),n_high,n_low);
      nhx=x+dinf_dx[n_high];
      nhy=y+dinf_dy[n_high];
      if(n_low!=-1){
        nlx = x+dinf_dx[n_low];
        nly = y+dinf_dy[n_low];
      }
      if( n_low!=-1 && flowdirs.inGrid(nlx,nly) && flowdirs(nlx,nly)!=flowdirs.noData() )
        dependency(nlx,nly)++;
      if( flowdirs.inGrid(nhx,nhy) && flowdirs(nhx,nhy)!=flowdirs.noData() )
        dependency(nhx,nhy)++;
    }
  }
  RDLOG_TIME_USE<<"Succeeded in = "<<progress.stop()<<" s";

  if(has_cells_without_flow_directions)
    RDLOG_WARN<<"\033[91mNot all cells had defined flow directions! This implies that there will be digital dams!\033[39m";

  ///////////////////////
  //Find those cells which have no dependencies. These are the places to start
  //the flow accumulation calculation.

  RDLOG_PROGRESS<<"Locating source cells...";
  progress.start( flowdirs.size() );
  for(int y=0;y<flowdirs.height();y++){
    progress.update( y*flowdirs.width() );
    for(int x=0;x<flowdirs.width();x++)
      if(flowdirs(x,y)==flowdirs.noData())
        continue;
      else if(flowdirs(x,y)==NO_FLOW)
        continue;
      else if(dependency(x,y)==0)
        sources.emplace(x,y);
  }
  RDLOG_TIME_USE<<"Source cells located in = "<<progress.stop()<<" s";





  ///////////////////////
  //Calculate the flow accumulation by "pouring" a cell's flow accumulation
  //value into the cells below it, as indicated by the D-infinite flow routing
  //method.

  RDLOG_PROGRESS<<"Calculating up-slope areas...";
  progress.start( flowdirs.numDataCells() );
  long int ccount=0;
  while(sources.size()>0){
    auto c = sources.front();
    sources.pop();

    progress.update(ccount++);

    if(flowdirs.isNoData(c.x,c.y))  //TODO: This line shouldn't be necessary since NoData's do not get added below
      continue;

    area(c.x,c.y)+=1;

    if(flowdirs(c.x,c.y)==NO_FLOW)
      continue;

    int n_high,n_low,nhx,nhy,nlx,nly;
    where_do_i_flow(flowdirs(c.x,c.y),n_high,n_low);
    nhx = c.x+dinf_dx[n_high];
    nhy = c.y+dinf_dy[n_high];

    float phigh,plow;
    area_proportion(flowdirs(c.x,c.y), n_high, n_low, phigh, plow);
    if(flowdirs.inGrid(nhx,nhy) && flowdirs(nhx,nhy)!=flowdirs.noData())
      area(nhx,nhy)+=area(c.x,c.y)*phigh;

    if(n_low!=-1){
      nlx = c.x+dinf_dx[n_low];
      nly = c.y+dinf_dy[n_low];
      if(flowdirs.inGrid(nlx,nly) && flowdirs(nlx,nly)!=flowdirs.noData()){
        area(nlx,nly)+=area(c.x,c.y)*plow;
        if((--dependency(nlx,nly))==0)
          sources.emplace(nlx,nly);
      }
    }

    if( flowdirs.inGrid(nhx,nhy) && flowdirs(nhx,nhy)!=flowdirs.noData() && (--dependency(nhx,nhy))==0)
      sources.emplace(nhx,nhy);
  }
  RDLOG_PROGRESS<<"Succeeded in = "<<progress.stop()<<" s";
}

}

#endif
