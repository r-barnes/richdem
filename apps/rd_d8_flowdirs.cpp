#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/depressions/Barnes2014.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/flats/flat_resolution.hpp"
#include "richdem/methods/d8_methods.hpp"

using namespace richdem;

template<class T>
int PerformAlgorithm(std::string output_prefix, std::string analysis, Array2D<T> elevations){
  elevations.loadData();

  PriorityFlood_Barnes2014<Topology::D8>(elevations);

  Array2D<uint8_t> flowdirs;

  barnes_flat_resolution_d8(elevations,flowdirs,false);

  flowdirs.saveGDAL(output_prefix+"-flowdirs.tif", analysis);

  return 0;
}

#include "router.hpp"

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);

  if(argc!=3){
    std::cerr<<argv[0]<<" <INPUT> <OUTPUT_PREFIX>"<<std::endl;
    return -1;
  }

  PerformAlgorithm(argv[1],argv[2],analysis);

  return 0;
}