/**
  @file
  @brief Functions for calculating flow directions according to a variety of authors

  Richard Barnes (rbarnes@umn.edu), 2016
*/
#ifndef _richdem_dall_flowdirs_hpp_
#define _richdem_dall_flowdirs_hpp_

#include "richdem/common/Array2D.hpp"
#include "richdem/common/ProgressBar.hpp"
#include "richdem/common/grid_cell.hpp"
#include "richdem/common/random.hpp"

#include <functional>

enum FDMode {
  CALC_DEPENDENCIES,
  CALC_ACCUM
};

typedef struct {
  double x;
} Params;

template<FDMode fd>
void foo(){
  std::cerr<<fd<<std::endl;
}

typedef Array2D<uint8_t> dep_t;

template<class E, class A>
static void KernelTarboton(
  const FDMode mode,
  const Array2D<E> &elevations,
  Array2D<A> &accum,
  dep_t &dep,
  std::queue<GridCell> &q,
  const int x,
  const int y,
  Array2D< std::pair<float,int8_t> > &fd
){
  //TODO: Assumes that the width and height of grid cells are equal and scaled
  //to 1.
  constexpr double d1   = 1;
  constexpr double d2   = 1;
  constexpr float  dang = std::atan2(d2,d1);

  auto nwrap = [](int n){ return (n==9)?1:n; };

  if(mode==FDMode::CALC_DEPENDENCIES){
    //Table 1 of Tarboton (1997)
    //          Column #  =   0    1    2    3    4    5   6    7
    // const int    dy_e1[8] = { 0 , -1 , -1 ,  0 ,  0 ,  1 , 1 ,  0 };
    // const int    dx_e1[8] = { 1 ,  0 ,  0 , -1 , -1 ,  0 , 0 ,  1 };
    // const int    dy_e2[8] = {-1 , -1 , -1 , -1 ,  1 ,  1 , 1 ,  1 };
    // const int    dx_e2[8] = { 1 ,  1 , -1 , -1 , -1 , -1 , 1 ,  1 };
    // const double ac[8]    = { 0.,  1.,  1.,  2.,  2.,  3., 3.,  4.};
    // const double af[8]    = { 1., -1.,  1., -1.,  1., -1., 1., -1.};

    //I remapped the foregoing table for ease of use with RichDEM. The facets
    //are renumbered as follows:
    //    3->1    2->2    1->3    0->4    7->5    6->6    5->7    4->8
    //This gives the following table
    //  Remapped Facet #  =  -   1    2     3    4    5   6    7    8  
    //  Tarboton Facet #  =  -   3    2     1    0    7   6    5    4  
    const int    dy_e1[9] = {0,  0 , -1 ,  -1 ,  0 ,  0 , 1 ,  1 ,  0  };
    const int    dx_e1[9] = {0, -1 ,  0 ,   0 ,  1 ,  1 , 0 ,  0 , -1  };
    const int    dy_e2[9] = {0, -1 , -1 ,  -1 , -1 ,  1 , 1 ,  1 ,  1  };
    const int    dx_e2[9] = {0, -1 , -1 ,   1 ,  1 ,  1 , 1 , -1 , -1  };
    //const double ac[9]    = {0,  2.,  1.,   1.,  0.,  4., 3.,  3.,  2. };
    const double af[9]    = {0, -1.,  1.,  -1.,  1., -1., 1., -1.,  1. };

    int8_t nmax = -1;
    double smax = 0;
    float  rmax = 0;

    for(int n=1;n<=8;n++){
      if(!elevations.inGrid (x+dx_e1[n],y+dy_e1[n]))
        continue;
      if(elevations.isNoData(x+dx_e1[n],y+dy_e1[n]))
        continue;
      if(!elevations.inGrid (x+dx_e2[n],y+dy_e2[n]))
        continue;
      if(elevations.isNoData(x+dx_e2[n],y+dy_e2[n]))
        continue;

      //Is is assumed that cells with a value of NoData have very negative
      //elevations with the result that they draw flow off of the grid.

      //Choose elevations based on Table 1 of Tarboton (1997), Barnes TODO
      const double e0 = elevations(x,y);
      const double e1 = elevations(x+dx_e1[n],y+dy_e1[n]);
      const double e2 = elevations(x+dx_e2[n],y+dy_e2[n]);

      const double s1 = (e0-e1)/d1;
      const double s2 = (e1-e2)/d2;

      double r = std::atan2(s2,s1);
      double s;

      if(r<1e-7){
        r = 0;
        s = s1;
      } else if(r>dang-1e-7){
        r = dang;
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

    if(nmax!=-1){
      if(af[nmax]==1 && rmax==0)
        rmax = dang;
      else if(af[nmax]==1 && rmax==dang)
        rmax = 0;
      else if(af[nmax]==1)
        rmax = M_PI/4-rmax;

      if(rmax==0){
        dep(x+dx[nmax],y+dy[nmax])++;
      } else if(rmax==dang){
        dep(x+dx[nwrap(nmax+1)],y+dy[nwrap(nmax+1)])++;
      } else {
        dep(x+dx[nmax],y+dy[nmax])++;
        dep(x+dx[nwrap(nmax+1)],y+dy[nwrap(nmax+1)])++;
      }
    }

    fd(x,y).first  = rmax;
    fd(x,y).second = nmax;

    //Code used by Tarboton to calculate the angle Rg. This should give the same
    //result despite the rearranged table
    // double rg = NO_FLOW;
    // if(nmax!=-1)
    //   rg = (af[nmax]*rmax+ac[nmax]*M_PI/2);
  } else {
    const float r = fd(x,y).first;
    const int   n = fd(x,y).second;

    if(n==-1)
      return;

    auto DoCell = [&](int n, float fraction){
      const int nx = x+dx[n];
      const int ny = y+dy[n];
      accum(nx,ny) += fraction*accum(x,y);
      if(--dep(nx,ny)==0)
        q.emplace(nx,ny);
    };

    if(r==0){
      DoCell(n, 1);
    } else if(r==dang){
      DoCell(nwrap(n+1), 1);
    } else {
      DoCell(n,            r/(M_PI/4.));
      DoCell(nwrap(n+1), 1-r/(M_PI/4.));
    }
  }
}

template<class E, class A>
static void KernelHolmgren(const FDMode mode, const Array2D<E> &elevations, Array2D<A> &accum, dep_t &dep, std::queue<GridCell> &q, const int x, const int y, const double xparam){
  const E e = elevations(x,y);

  constexpr double L1   = 0.5;
  constexpr double L2   = 0.354; //TODO: More decimal places
  constexpr double L[9] = {0,L1,L2,L1,L2,L1,L2,L1,L2};

  std::array<double,9> portions = {{0,0,0,0,0,0,0,0,0}};

  double C = 0;
  for(int n=1;n<=8;n++){
    const int nx = x+dx[n];
    const int ny = y+dy[n];

    if(!elevations.inGrid(nx,ny))
      continue;
    if(elevations.isNoData(nx,ny)) //TODO: Don't I want water to drain this way?
      continue;

    const E ne = elevations(nx,ny);

    if(ne<e){
      const double rise = e-ne;
      const double run  = dr[n];
      const double grad = rise/run;
      portions[n]       = std::pow(grad * L[n],xparam);
      C                += portions[n];
    }
  }

  C = accum(x,y)/C;

  for(int n=1;n<=8;n++){
    if(portions[n]>0){
      const int nx = x+dx[n];
      const int ny = y+dy[n];
      if(mode==FDMode::CALC_DEPENDENCIES){
        dep(nx,ny)++;
      } else {
        accum(nx,ny) += portions[n]*C;
        dep(nx,ny)--;
        if(dep(nx,ny)==0)
          q.emplace(nx,ny);
      }
    }
  }
}

template<class E, class A>
void KernelFairfieldLeymarie(const FDMode mode, const Array2D<E> &elevations, Array2D<A> &accum, dep_t &dep, std::queue<GridCell> &q, const int x, const int y, Array2D<d8_flowdir_t> &fd) {
  if(mode==FDMode::CALC_DEPENDENCIES){
    const E e = elevations(x,y);

    int    greatest_n     = 0;
    double greatest_slope = 0;
    for(int n=1;n<=8;n++){
      const int nx = x+dx[n];
      const int ny = y+dy[n];

      if(!elevations.inGrid(nx,ny))
        continue;
      if(elevations.isNoData(nx,ny)) //TODO: Don't I want water to drain this way?
        continue;

      const E ne = elevations(nx,ny);

      if(ne>=e)
        continue;

      double rho_slope = (e-ne);
      if(n_diag[n])
        rho_slope *= 1/(2-uniform_rand_real(0,1));

      if(rho_slope>greatest_slope){
        greatest_n     = n;
        greatest_slope = rho_slope;
      }
    }

    fd(x,y) = greatest_n;

    const int nx = x+dx[greatest_n];
    const int ny = y+dy[greatest_n];
    dep(nx,ny)++;
  } else {
    if(fd(x,y)==0)
      return;

    const int nx = x+dx[fd(x,y)];
    const int ny = y+dy[fd(x,y)];

    accum(nx,ny) += accum(x,y);
    dep(nx,ny)--;
    if(dep(nx,ny)==0)
      q.emplace(nx,ny);
  }
}


template<class F, class E, class A, typename... Args>
void KernelFlowdir(
  F f,
  const Array2D<E> &elevations,
  Array2D<A> &accum,
  Args&&... args
){
    ProgressBar progress;
    std::queue<GridCell> q;
    dep_t dep;

    Timer overall;
    overall.start();

    accum.setAll(0);
    accum.setNoData(ACCUM_NO_DATA);
    for(typename Array2D<E>::i_t i=0;i<elevations.size();i++)
      if(elevations.isNoData(i))
        accum(i) = ACCUM_NO_DATA;

    dep.resize(elevations);
    dep.setAll(0);

    std::cerr<<"p Calculating dependencies..."<<std::endl;
    progress.start(elevations.size());
    for(int y=0;y<elevations.height();y++)
    for(int x=0;x<elevations.width();x++){
      ++progress;
      if(!elevations.isNoData(x,y))
        f(FDMode::CALC_DEPENDENCIES,elevations,accum,dep,q,x,y,std::forward<Args>(args)...);
    }
    progress.stop();

    for(int y=0;y<dep.height();y++)
    for(int x=0;x<dep.width();x++)
      if(dep(x,y)==0 && !elevations.isNoData(x,y))
        q.emplace(x,y);

    std::cerr<<"p Calculating accumulation..."<<std::endl;
    progress.start(accum.numDataCells());
    while(!q.empty()){
      ++progress;

      const auto c = q.front();
      q.pop();

      accum(c.x,c.y) += 1;
      f(FDMode::CALC_ACCUM,elevations,accum,dep,q,c.x,c.y,std::forward<Args>(args)...);
    }
    progress.stop();

    for(int y=0;y<elevations.height();y++)
    for(int x=0;x<elevations.width();x++){
      if(accum(x,y)==0){
        std::cerr<<"x,y: "<<x<<","<<y<<std::endl;
        std::cerr<<"deps: "<<(int)dep(x,y)<<std::endl;
      }
    }

    std::cerr<<"m Data cells      = "<<elevations.numDataCells()<<std::endl;
    std::cerr<<"m Cells processed = "<<progress.cellsProcessed()<<std::endl;
    std::cerr<<"m Max accum       = "<<accum.max()              <<std::endl;
    std::cerr<<"m Min accum       = "<<accum.min()              <<std::endl;
    std::cerr<<"t Wall-time       = "<<overall.stop()<<" s"     <<std::endl;
  }

template<class E, class A>
void FA_FairfieldLeymarie(const Array2D<E> &elevations, Array2D<A> &accum){
  std::cerr<<"\nA Fairfield (Rho8) Flow Accumulation (TODO)"<<std::endl;
  Array2D<d8_flowdir_t> fd(elevations);
  KernelFlowdir(KernelFairfieldLeymarie<E,A>,elevations,accum,fd);
}

template<class E, class A>
void FA_Rho8(const Array2D<E> &elevations, Array2D<A> &accum){
  //Algorithm headers are taken care of in FA_FairfieldLeymarie()
  FA_FairfieldLeymarie(elevations,accum);
}

template<class E, class A>
void FA_Quinn(const Array2D<E> &elevations, Array2D<A> &accum){
  std::cerr<<"\nA Quinn 1991 Flow Accumulation (TODO)"<<std::endl;
  KernelFlowdir(KernelHolmgren<E,A>,elevations,accum,(double)1.0);
}

template<class E, class A>
void FA_Holmgren(const Array2D<E> &elevations, Array2D<A> &accum, double x){
  std::cerr<<"\nA Holmgren Flow Accumulation (TODO)"<<std::endl;
  KernelFlowdir(KernelHolmgren<E,A>,elevations,accum,x);
}

template<class E, class A>
void FA_Tarboton(const Array2D<E> &elevations, Array2D<A> &accum){
  std::cerr<<"\nA Tarboton (D-Infinity) Flow Accumulation (TODO)"<<std::endl;
  Array2D< std::pair<float,int8_t> > fd(elevations);
  KernelFlowdir(KernelTarboton<E,A>,elevations,accum,fd);
}







/**
  @brief  Calculates the D8 flow direction of a cell
  @author Richard Barnes (rbarnes@umn.edu)

  This calculates the D8 flow direction of a cell using the D8
  neighbour system, as defined in utility.h. Cells on the edge
  of the grid flow off the nearest edge.

  Helper function for d8_flow_directions().

  @param[in]  &elevations  A DEM
  @param[in]  x            x coordinate of cell
  @param[in]  y            y coordinate of cell

  @returns The D8 flow direction of the cell
*/
/*
template<class T>
static int dall_EdgeFlow(const Array2D<T> &elevations, const int x, const int y){
  T minimum_elevation = elevations(x,y);
  int flowdir         = NO_FLOW;

  if (elevations.isEdgeCell(x,y)){
    if(elevations.isTopLeft(x,y))
      return 2;
    else if(elevations.isBottomLeft(x,y))
      return 8;
    else if(elevations.isTopRight(x,y))
      return 4;
    else if(elevations.isBottomRight(x,y))
      return 6;
    else if(elevations.isLeftCol(x,y))
      return 1;
    else if(elevations.isRightCol(x,y))
      return 5;
    else if(elevations.isTopRow(x,y))
      return 3;
    else if(elevations.isBottomRow(x,y))
      return 7;
  }

  /*NOTE: Since the very edges of the DEM are defined to always flow outwards,
  if they have defined elevations, it is not necessary to check if a neighbour
  is IN_GRID in the following
  NOTE: It is assumed that the no_data datum is an extremely negative
  number, such that all water which makes it to the edge of the DEM's region
  of defined elevations is sucked directly off the grid, rather than piling up
  on the edges.*/
/*
  for(int n=1;n<=8;n++)
    if(
      elevations(x+dx[n],y+dy[n])<minimum_elevation
      || (elevations(x+dx[n],y+dy[n])==minimum_elevation
            && flowdir>0 && flowdir%2==0 && n%2==1) //TODO: What is this modulus stuff for?
    ){
      minimum_elevation=elevations(x+dx[n],y+dy[n]);
      flowdir=n;
    }

  return flowdir;
}
*/


//d8_flow_directions
/**
  @brief  Calculates the D8 flow directions of a DEM
  @author Richard Barnes (rbarnes@umn.edu)

  This calculates the D8 flow directions of a DEM. Its argument
  'flowdirs' will return a grid with flow directions using the D8
  neighbour system, as defined in utility.h. The choice of data type
  for array2d must be able to hold exact values for all neighbour
  identifiers (usually [-1,7]).

  Uses d8_FlowDir() as a helper function.

  @todo                    Combine dinf and d8 neighbour systems

  @param[in]  &elevations  A DEM
  @param[out] &flowdirs    Returns the flow direction of each cell
*/
/*
template<class T, class U>
void d8_flow_directions(
  const Array2D<T> &elevations,
        Array2D<U> &flowdirs
){
  ProgressBar progress;

  std::cerr<<"A D8 Flow Directions"<<std::endl;
  std::cerr<<"C TODO"<<std::endl;

  std::cerr<<"p Setting up the flow directions matrix..."<<std::endl;
  flowdirs.resize(elevations);
  flowdirs.setAll(NO_FLOW);
  flowdirs.setNoData(FLOWDIR_NO_DATA);

  std::cerr<<"p Calculating D8 flow directions..."<<std::endl;
  progress.start( elevations.width()*elevations.height() );
  #pragma omp parallel for
  for(int y=0;y<elevations.height();y++){
    progress.update( y*elevations.width() );
    for(int x=0;x<elevations.width();x++)
      if(elevations(x,y)==elevations.noData())
        flowdirs(x,y) = flowdirs.noData();
      else
        flowdirs(x,y) = d8_FlowDir(elevations,x,y);
  }
  std::cerr<<"t Succeeded in = "<<progress.stop()<<" s"<<std::endl;
}
*/
#endif