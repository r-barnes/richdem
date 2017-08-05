#include <iostream>
#include <iomanip>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/router.hpp"

template<class T>
int PerformAlgorithm(std::string flip, Array2D<T> raster){
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

  raster.loadData();

  std::cerr<<"Geotransform: ";
  for(auto const x: raster.geotransform)
    std::cerr<<std::setw(6)<<std::setprecision(2)<<x<<" ";
  std::cerr<<std::endl;

  //Flip tiles if the geotransform demands it
  if( (raster.geotransform[1]<0) ^ flipH)
    raster.flipHorz();
  if( (raster.geotransform[5]>0) ^ flipV)
    raster.flipVert();

  raster.printStamp(5);

  for(int y=0;y<raster.height();y++){
    for(int x=0;x<raster.width();x++)
      std::cout<<std::setw(6)<<std::setprecision(2)<<(int)raster(x,y)<<" ";
    std::cout<<"\n";
  }

  return 0;
}

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);
  
  if(argc!=3){
    std::cerr<<"Print a raster to the terminal."<<std::endl;
    std::cerr<<argv[0]<<" <Input file> <noflip/fliph/flipv/fliphv>"<<std::endl;
    return -1;
  }

  PerformAlgorithm(argv[1],argv[2]);

  return 0;
}