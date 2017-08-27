#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include "richdem/common/version.hpp"
#include "richdem/common/router.hpp"
#include "richdem/misc/misc_methods.hpp"
#include "richdem/common/Array2D.hpp"

template<class T>
int PerformAlgorithm(const double zscale, std::string analysis, Array2D<T> elevation){
  elevation.loadData();

  const double sa_without_topo = elevation.numDataCells()*elevation.getCellArea();
  const double sa_with_topo    = dem_surface_area(elevation, zscale);
  const double sa_diff         = sa_with_topo-sa_without_topo;
  const double sa_increase     = 100*(sa_with_topo-sa_without_topo)/sa_without_topo;
  const double perim_len       = Perimeter(elevation, PerimType::SQUARE_EDGE);


  std::cout<<"Surface area with topography    = "
           <<std::fixed<<std::setprecision(10)<<sa_with_topo
           <<std::endl;
  
  std::cout<<"Surface area without topography = "
           <<std::fixed<<std::setprecision(10)<<sa_without_topo
           <<std::endl;

  std::cout<<"Surface area from topography    = "
           <<std::fixed<<std::setprecision(10)<<sa_diff
           <<std::endl;

  std::cout<<"Surface area increase           = "
           <<std::fixed<<std::setprecision(10)<<sa_increase<<" %"
           <<std::endl;

  std::cout<<"Perimeter                       = "
           <<std::fixed<<std::setprecision(10)<<perim_len
           <<std::endl;


  std::cout<<"Units                           = "
           <<"Input units squared"
           <<std::endl;

  return 0;
}

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=3){
    std::cerr<<"Calculate the surface are of a DEM accounting for topography."<<std::endl;
    std::cerr<<argv[0]<<" <Input> <z-scale>"<<std::endl;
    std::cerr<<"\t<z-scale> - All elevations are multiplied by the z-scale prior to calculation."<<std::endl;
    return -1;
  }

  return PerformAlgorithm(argv[1],std::stod(argv[2]),analysis);
}