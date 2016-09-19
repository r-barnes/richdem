#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/depressions/priority_flood.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/flats/flat_resolution.hpp"
#include "richdem/methods/d8_methods.hpp"
#include "richdem/common/router.hpp"

template<class T>
int PerformAlgorithm(std::string output_prefix, std::string analysis, Array2D<T> elevations){
  elevations.loadData();

  improved_priority_flood(elevations);

  Array2D<uint8_t> flowdirs;

  barnes_flat_resolution_d8(elevations,flowdirs,false);

  for(int32_t y=0;y<flowdirs.height();y++){
    for(int32_t x=0;x<flowdirs.width();x++)
      std::cerr<<(int)flowdirs(x,y)<<" ";
    std::cerr<<std::endl;
  }


  flowdirs.saveGDAL(output_prefix+"-flowdirs.tif", analysis);

  return 0;
}

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);

  if(argc!=3){
    std::cerr<<argv[0]<<" <INPUT> <OUTPUT_PREFIX>"<<std::endl;
    return -1;
  }

  PerformAlgorithm(argv[1],argv[2],analysis);

  return 0;
}