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
int PerformAlgorithm(std::string outputname, std::string analysis, Array2D<T> elevations){
  bool flipH = false; //TODO
  bool flipV = false;
  std::string flip = "noflip";

  if(flip=="fliph")
    flipH=true;
  else if(flip=="flipv")
    flipV=true;
  else if(flip=="fliphv")
    flipH=flipV=true;
  else if(flip!="noflip"){
    std::cerr<<"Unrecognised flip directive."<<std::endl;
    return -1;
  }

  elevations.loadData();

  //Flip tiles if the geotransform demands it
  if( (elevations.geotransform[1]<0) ^ flipH)
    elevations.flipHorz();
  if( (elevations.geotransform[5]>0) ^ flipV)
    elevations.flipVert();


  PriorityFlood_Barnes2014<Topology::D8>(elevations);
  
  elevations.printStamp(5);

  Array2D<uint8_t> flowdirs;

  barnes_flat_resolution_d8(elevations,flowdirs,false);

  flowdirs.printStamp(5);

  //Flip tiles if the geotransform demands it
  if( (flowdirs.geotransform[1]<0) ^ flipH)
    flowdirs.flipHorz();
  if( (flowdirs.geotransform[5]>0) ^ flipV)
    flowdirs.flipVert();

  flowdirs.printStamp(5);

  flowdirs.saveGDAL(outputname, analysis);

  return 0;
}

#include "router.hpp"



int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);
  
  if(argc!=3){
    std::cerr<<argv[0]<<" <Input filename> <Output filename>"<<std::endl;
    return -1;
  }

  PerformAlgorithm(argv[1],argv[2],analysis);

  return 0;
}