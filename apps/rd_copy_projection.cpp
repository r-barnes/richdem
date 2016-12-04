#include <iostream>
#include <iomanip>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/router.hpp"

template<class T>
int PerformAlgorithm(std::string templatefile, std::string outputfile, std::string flip, std::string analysis, Array2D<T> raster){
  raster.loadData();
  Array2D<int8_t> temp  (templatefile,false,0,0,0,0,false,false); //Data type doesn't matter since we're not loading it

  raster.projection   = temp.projection;

  if(flip=="fliph" || flip=="fliphv")
    raster.flipHorz();
  if(flip=="flipv" || flip=="fliphv")
    raster.flipVert();
  if(flip!="fliph" && flip!="flipv" && flip!="fliphv" && flip!="noflip"){
    std::cerr<<"Unrecognised flip directive!"<<std::endl;
    return -1;
  }

  raster.saveGDAL(outputfile,analysis);

  return 0;
}

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);

  if(argc!=5){
    std::cerr<<argv[0]<<" <Input file> <File To Copy From> <Output File> <fliph/flipv/fliphv/noflip>"<<std::endl;
    return -1;
  }

  PerformAlgorithm(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4]),analysis);

  return 0;
}