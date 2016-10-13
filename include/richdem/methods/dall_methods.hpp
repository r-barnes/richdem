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
class KernelFlowdir {
 private:
  virtual void doCell(const FDMode mode, const Array2D<E> &elevations, Array2D<A> &accum, const int x, const int y) = 0;

 protected:
  std::queue<GridCell> q;
  dep_t dep;

 public:
  void run(const Array2D<E> &elevations, Array2D<A> &accum){
    ProgressBar progress;

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
        doCell(FDMode::CALC_DEPENDENCIES,elevations,accum,x,y);
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
      doCell(FDMode::CALC_ACCUM,elevations,accum,c.x,c.y);
    }
    progress.stop();

    std::cerr<<"m Data cells      = "<<elevations.numDataCells()<<std::endl;
    std::cerr<<"m Cells processed = "<<progress.cellsProcessed()<<std::endl;
    std::cerr<<"m Max accum       = "<<accum.max()              <<std::endl;
    std::cerr<<"m Min accum       = "<<accum.min()              <<std::endl;
    std::cerr<<"t Wall-time       = "<<overall.stop()<<" s"     <<std::endl;
  }
};

template<class E, class A>
class KernelHolmgren : public KernelFlowdir<E,A> {
 public:
  double x;

  KernelHolmgren(double x){
    this->x = x;
  }

  void doCell(const FDMode mode, const Array2D<E> &elevations, Array2D<A> &accum, const int x, const int y){
    const E e = elevations(x,y);

    constexpr double L1   = 0.5;
    constexpr double L2   = 0.354; //TODO: More decimal places
    constexpr double L[9] = {0,L1,L2,L1,L2,L1,L2,L1,L2};

    std::array<double,9> portions = {0,0,0,0,0,0,0,0,0};

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
        portions[n]       = std::pow(grad * L[n],x);
        C                += grad * L[n];
      }
    }

    C = accum(x,y)/std::pow(C,x);

    for(int n=1;n<=8;n++){
      if(portions[n]>0){
        const int nx = x+dx[n];
        const int ny = y+dy[n];
        if(mode==FDMode::CALC_DEPENDENCIES){
          this->dep(nx,ny)++;
        } else {
          accum(nx,ny) += portions[n]*C;
          this->dep(nx,ny)--;
          if(this->dep(nx,ny)==0)
            this->q.emplace(nx,ny);
        }
      }
    }
  }
};

template<class E, class A>
class KernelFairfieldLeymarie : public KernelFlowdir<E,A> {
 public:
  Array2D<d8_flowdir_t> fd;

  KernelFairfieldLeymarie(const Array2D<E> &elevations) {
    fd.resize(elevations);
  }

  void doCell(const FDMode mode, const Array2D<E> &elevations, Array2D<A> &accum, const int x, const int y) {
    const E e = elevations(x,y);

    int    greatest_n     = 0;
    double greatest_slope = 0;
    if(mode==FDMode::CALC_DEPENDENCIES){
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
      this->dep(nx,ny)++;
    } else {
      const int nx = x+dx[fd(x,y)];
      const int ny = y+dy[fd(x,y)];

      accum(nx,ny) += accum(x,y);
      this->dep(nx,ny)--;
      if(this->dep(nx,ny)==0)
        this->q.emplace(nx,ny);
    }
  }
};



template<class E, class A>
void FA_Quinn(const Array2D<E> &elevations, Array2D<A> &accum){
  std::cerr<<"A Quinn 1991 Flow Accumulation (TODO)"<<std::endl;
  KernelHolmgren<E,A> kh(1.0);
  kh.run(elevations, accum);
}

template<class E, class A>
void FA_Holmgren(const Array2D<E> &elevations, Array2D<A> &accum, double x){
  std::cerr<<"A Holmgren Flow Accumulation (TODO)"<<std::endl;
  KernelHolmgren<E,A> kh(x);
  kh.run(elevations, accum);
}


template<class E, class A>
void FA_FairfieldLeymarie(const Array2D<E> &elevations, Array2D<A> &accum){
  std::cerr<<"A Fairfield (Rho8) Flow Accumulation (TODO)"<<std::endl;
  KernelFairfieldLeymarie<E,A> kfl(elevations);
  kfl.run(elevations, accum);
}

template<class E, class A>
void FA_Rho8(const Array2D<E> &elevations, Array2D<A> &accum){
  //Algorithm headers are taken care of in FA_FairfieldLeymarie()
  FA_FairfieldLeymarie(elevations,accum);
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