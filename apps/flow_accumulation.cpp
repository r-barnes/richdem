#include <iostream>
#include <string>
#include <cstdlib>
#include "../libs/depressions/priority_flood.hpp"
#include "../libs/common/Array2D.hpp"
#include "../libs/flats/flat_resolution.hpp"
#include "../libs/methods/d8_methods.hpp"

inline bool XOR(bool a, bool b){
  return a!=b;
}

template<class T>
int PerformAlgorithm(std::string filename, std::string output){
  bool flipH = false; //TODO
  bool flipV = false;

  Array2D<T> flowdirs(filename,false);


  //Flip tiles if the geotransform demands it
  if( XOR(flowdirs.geotransform[0]<0, flipH) )
    flowdirs.flipHorz();
  if( XOR(flowdirs.geotransform[5]<0, flipV) )
    flowdirs.flipVert();


  Array2D<int> area;
  d8_upslope_area(flowdirs, area);


  if( XOR(area.geotransform[0]<0, flipH) )
    area.flipHorz();
  if( XOR(area.geotransform[5]<0, flipV) )
    area.flipVert();

  area.saveGDAL(output,0,0);

  return 0;
}

int main(int argc, char **argv){
  if(argc!=3){
    std::cerr<<argv[0]<<" <Flowdirs input file> <Output filename>"<<std::endl;
    return -1;
  }

  switch(peekGDALType(argv[1])){
    case GDT_Unknown:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[1]))<<std::endl;
      return -1;
    case GDT_Byte:
      return PerformAlgorithm<uint8_t >(argv[1],argv[2]);
    case GDT_UInt16:
      return PerformAlgorithm<uint16_t>(argv[1],argv[2]);
    case GDT_Int16:
      return PerformAlgorithm<int16_t >(argv[1],argv[2]);
    case GDT_UInt32:
      return PerformAlgorithm<uint32_t>(argv[1],argv[2]);
    case GDT_Int32:
      return PerformAlgorithm<int32_t >(argv[1],argv[2]);
    case GDT_Float32:
      return PerformAlgorithm<float   >(argv[1],argv[2]);
    case GDT_Float64:
      return PerformAlgorithm<double  >(argv[1],argv[2]);
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