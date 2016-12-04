#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/common/router.hpp"
#include "richdem/methods/d8_methods.hpp"
#include "richdem/common/Array2D.hpp"

template<class T>
int PerformAlgorithm(std::string output, int algorithm, float z_scale, std::string analysis, Array2D<T> dem){
  dem.loadData();

  Array2D<float> result(dem);

  switch(algorithm){
    case 1:  //d8_slope_riserun
      d8_slope_riserun(dem,result,z_scale);      break;
    case 2:  //d8_slope_percentage
      d8_slope_percentage(dem,result,z_scale);   break;
    case 3:  //d8_slope_degrees
      d8_slope_degrees(dem,result,z_scale);      break;
    case 4:  //d8_slope_radians
      d8_slope_radians(dem,result,z_scale);      break;
    case 5:  //d8_aspect
      d8_aspect(dem,result,z_scale);             break;
    case 6:  //d8_curvature
      d8_curvature(dem,result,z_scale);          break;
    case 7:  //d8_planform_curvature
      d8_planform_curvature(dem,result,z_scale); break;
    case 8: //d8_profile_curvature
      d8_profile_curvature(dem,result,z_scale);  break;
  }

  result.saveGDAL(output,analysis);

  return 0;
}

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);
  
  int   algorithm = 0;
  float z_scale   = 1;

  if(argc!=5){
    std::cerr<<"Calculate terrain attributes. Ensure that vertical and horizontal axes have the same units!"<<std::endl;
    std::cerr<<argv[0]<<" <DEM file> <Output File> <Algorithm #> <Z scaling factor>"<<std::endl;
    std::cerr<<"Algorithms:"<<std::endl;
    std::cerr<<" 1:  Slope (Rise/Run)   - Horn (1981)"<<std::endl;
    std::cerr<<" 2:  Slope (Percentage) - Horn (1981)"<<std::endl;
    std::cerr<<" 3:  Slope (Degrees)    - Horn (1981)"<<std::endl;
    std::cerr<<" 4:  Slope (Radians)    - Horn (1981)"<<std::endl;
    std::cerr<<" 5:  Apsect             - Horn (1981)"<<std::endl;
    std::cerr<<" 6:  Curvature          - Zevenbergen and Thorne (1987)"<<std::endl;
    std::cerr<<" 7:  Planform Curvature - Zevenbergen and Thorne (1987)"<<std::endl;
    std::cerr<<" 8: Profile Curvature  - Zevenbergen and Thorne (1987)"<<std::endl;
    return -1;
  }

  algorithm = std::stoi(argv[3]);
  z_scale   = std::stof(argv[4]);

  return PerformAlgorithm(std::string(argv[1]),std::string(argv[2]),algorithm,z_scale,analysis);
}