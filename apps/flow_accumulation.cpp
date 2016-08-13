#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/depressions/priority_flood.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/flats/flat_resolution.hpp"
#include "richdem/methods/d8_methods.hpp"

template<class T>
int PerformAlgorithm(std::string filename, std::string output, std::string flip){
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

  Array2D<T> flowdirs(filename,false);

  std::cerr<<"Geotransform: ";
  for(auto const x: flowdirs.geotransform)
    std::cerr<<std::setw(5)<<std::setprecision(2)<<x<<" ";
  std::cerr<<std::endl;

  flowdirs.printStamp(5,"Stamp before reorientation");

  //Flip tiles if the geotransform demands it
  if( (flowdirs.geotransform[1]<0) ^ flipH){
    std::cerr<<"Flipping horizontally."<<std::endl;
    flowdirs.flipHorz();
  }
  if( (flowdirs.geotransform[5]>0) ^ flipV){
    std::cerr<<"Flipping vertically."<<std::endl;
    flowdirs.flipVert();
  }

  flowdirs.printStamp(5,"Stamp after reorientation");

  Array2D<int> area;
  d8_upslope_area(flowdirs, area);

  area.printStamp(5,"Output stamp before reorientation");

  if( (area.geotransform[1]<0) ^ flipH)
    area.flipHorz();
  if( (area.geotransform[5]>0) ^ flipV)
    area.flipVert();

  area.printStamp(5,"Output stamp after reorientation");

  area.saveGDAL(output,0,0);

  return 0;
}

int main(int argc, char **argv){
  PrintRichdemHeader();
  
  if(argc!=4){
    std::cerr<<argv[0]<<" <Flowdirs input file> <Output filename> <noflip/fliph/flipv/fliphv>"<<std::endl;
    return -1;
  }

  switch(peekGDALType(argv[1])){
    case GDT_Unknown:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[1]))<<std::endl;
      return -1;
    case GDT_Byte:
      return PerformAlgorithm<uint8_t >(argv[1],argv[2],argv[3]);
    case GDT_UInt16:
      return PerformAlgorithm<uint16_t>(argv[1],argv[2],argv[3]);
    case GDT_Int16:
      return PerformAlgorithm<int16_t >(argv[1],argv[2],argv[3]);
    case GDT_UInt32:
      return PerformAlgorithm<uint32_t>(argv[1],argv[2],argv[3]);
    case GDT_Int32:
      return PerformAlgorithm<int32_t >(argv[1],argv[2],argv[3]);
    case GDT_Float32:
      return PerformAlgorithm<float   >(argv[1],argv[2],argv[3]);
    case GDT_Float64:
      return PerformAlgorithm<double  >(argv[1],argv[2],argv[3]);
    case GDT_CInt16:
    case GDT_CInt32:
    case GDT_CFloat32:
    case GDT_CFloat64:
      std::cerr<<"Complex types are unsupported. Sorry!"<<std::endl;
      return -1;
    default:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[1]))<<std::endl;
      return -1;
  }

  return 0;
}