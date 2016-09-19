#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/common/router.hpp"
#include "richdem/depressions/priority_flood.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/flats/flat_resolution.hpp"
#include "richdem/methods/d8_methods.hpp"

template<class T>
int PerformAlgorithm(std::string output, std::string flip, std::string analysis, Array2D<T> flowdirs){
  bool flipH = false; //TODO
  bool flipV = false;

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

  flowdirs.loadData();

  std::cerr<<"Geotransform: ";
  if(!flowdirs.geotransform.empty()){
    for(auto const x: flowdirs.geotransform)
      std::cerr<<std::setw(5)<<std::setprecision(2)<<x<<" ";
    std::cerr<<std::endl;
  }

  flowdirs.printStamp(5,"Stamp before reorientation");

  //Flip tiles if the geotransform demands it
  if( !flowdirs.geotransform.empty() && ((flowdirs.geotransform[1]<0) ^ flipH)){
    std::cerr<<"Flipping horizontally."<<std::endl;
    flowdirs.flipHorz();
  }
  if( !flowdirs.geotransform.empty() && ((flowdirs.geotransform[5]>0) ^ flipV)){
    std::cerr<<"Flipping vertically."<<std::endl;
    flowdirs.flipVert();
  }

  flowdirs.printStamp(5,"Stamp after reorientation");

  Array2D<int> area;
  d8_flow_accum(flowdirs, area);

  area.printStamp(5,"Output stamp before reorientation");

  if( !flowdirs.geotransform.empty() && ((area.geotransform[1]<0) ^ flipH))
    area.flipHorz();
  if( !flowdirs.geotransform.empty() && ((area.geotransform[5]>0) ^ flipV))
    area.flipVert();

  area.printStamp(5,"Output stamp after reorientation");

  area.saveGDAL(output,analysis);

  return 0;
}

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);
  
  if(argc!=4){
    std::cerr<<argv[0]<<" <Flowdirs input file> <Output filename> <noflip/fliph/flipv/fliphv>"<<std::endl;
    return -1;
  }

  return PerformAlgorithm(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),analysis);
}