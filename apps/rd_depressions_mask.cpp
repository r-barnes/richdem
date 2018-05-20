#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/depressions/Barnes2014.hpp"
#include "richdem/common/Array2D.hpp"
using namespace richdem;


template<class T>
int PerformAlgorithm(std::string output, std::string analysis, Array2D<T> elevation){
  elevation.loadData();

  Array2D<uint8_t> mask(elevation);

  pit_mask<Topology::D8>(elevation, mask);

  mask.saveGDAL(output, analysis);

  return 0;
}

#include "router.hpp"



int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=3){
    std::cerr<<"Return a raster in which 1 indicates depressions, 0 indicates non-depressions, and 3 indicates NoData."<<std::endl;
    std::cerr<<argv[0]<<" <Input> <Output>"<<std::endl;
    return -1;
  }

  return PerformAlgorithm(std::string(argv[1]),std::string(argv[2]),analysis);
}
