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
#include "richdem/common/Array2D.hpp"
#include "richdem/common/constants.hpp"
#include "richdem/common/ProgressBar.hpp"


/**
  @brief  Calculate the surface of a digital elevation model
  @author Jenness (2004), Richard Barnes (rbarnes@umn.edu)

    Calculates the surface area of a digital elevation model by connecting the
    central points of cells with triangles and then calculating the area of the
    portion of each triangle which falls within the focal cell. The method is
    described in detail in Jenness (2004)
    <doi:10.2193/0091-7648(2004)032[0829:CLSAFD]2.0.CO;2>

  @param[in]  &elevations A grid of elevations

  @return The surface area of the digital elevation model
*/
template <class T>
double dem_surface_area(
  const Array2D<T> &elevations
){
  ProgressBar progress;

  std::cerr<<"\nA DEM Surface Elevation"<<std::endl;
  std::cerr<<"C Jenness, J.S., 2004. Calculating landscape surface area from digital elevation models. Wildlife Society Bulletin 32, 829--839. doi:10.2193/0091-7648(2004)032[0829:CLSAFD]2.0.CO;2"<<std::endl;

  const auto euc_dist = [](const T a, const T b){ return std::sqrt(std::pow((double)a,2.0)+std::pow((double)b,2.0)); };

  double area = 0;

  const double xdist     = elevations.getCellLengthX();
  const double ydist     = elevations.getCellLengthY();
  const double diag_dist = euc_dist(xdist,ydist);

  progress.start(elevations.size());
  #pragma omp parallel for reduction(+:area)
  for(int y=0;y<elevations.height();y++){
    progress.update( y*elevations.width() );
    for(int x=0;x<elevations.width();x++){
      if(elevations.isNoData(x,y))
        continue;

      //Loop through neighbours
      for(int n=1;n<=8;n++){
        //This is the next neighbour, which forms part of the triangle
        int nn = n+1;
        if(nn==9) //Wrap around
          nn=1;

        //In each triangle one neighbour is in the diagonal direction and one is
        //in a straight direction.
        int dn  = n;   //Diagonal Neighbour
        int ndn = nn;  //Not-diagonal neighbour

        //As we walk around half of the time we'll misidentify the diagonal
        //neighbour, but we can fix that here by swapping the labels we just
        //gave the neighbours
        if(!n_diag[dn])
          std::swap(dn,ndn);

        const double my_elev = elevations(x,y);

        //Deal with the possibility that the neighbouring cells do not exist. In
        //this case, we pretend that they do exist and are at the same height as
        //the focal cell.
        double dn_elev;
        if(elevations.inGrid(x+dx[dn],y+dy[dn]) && !elevations.isNoData(x+dx[dn],y+dy[dn]))
          dn_elev = elevations(x+dx[dn],y+dy[dn]);
        else
          dn_elev = my_elev;

        //Do the same for the other neighbour
        double ndn_elev;
        if(elevations.inGrid(x+dx[ndn],y+dy[ndn]) && !elevations.isNoData(x+dx[ndn],y+dy[ndn]))
          ndn_elev = elevations(x+dx[ndn],y+dy[ndn]);
        else
          ndn_elev = my_elev;
   
        const double planar_dist_dn   = diag_dist;                   //Distance focal cell to diagonal neighbour
        const double planar_dist_ndn  = (dy[ndn] == 0)?xdist:ydist;  //Distance focal cell to non-diagonal neighbour
        const double planar_dist_bn   = (dy[ndn] == 0)?ydist:xdist;  //Distance between the neighbour cells
        
        const double elev_diff_dn     = dn_elev -my_elev; //Elevation drop between focal and diagonal
        const double elev_diff_ndn    = ndn_elev-my_elev; //Elevation drop between focal and non-diagonal
        const double elev_diff_bn     = ndn_elev-dn_elev; //Elevation drop between neighbours

        //Divide these distances by two to form a similar triangle constrained
        //by the boundary of the focal cell
        const double surf_dist_dn     = euc_dist(planar_dist_dn,elev_diff_dn)/2;   //3-space distance between center of focal and diagonal neighbour
        const double surf_dist_ndn    = euc_dist(planar_dist_ndn,elev_diff_ndn)/2; //3-space distance between center of focal and non-diagonal neighbour
        const double surf_dist_bn     = euc_dist(planar_dist_bn,elev_diff_bn)/2;   //3-space distance between neighbours

        //Used to get area of triangle
        const double s = (surf_dist_dn+surf_dist_ndn+surf_dist_bn)/2;

        //Accumulate area of triangle to
        area += std::sqrt(s*(s-surf_dist_dn)*(s-surf_dist_ndn)*(s-surf_dist_bn));
      }      
    }
  }
  std::cerr<<"p Succeeded in = "<<progress.stop()<<" s"<<std::endl;

  return area;
}

#endif
