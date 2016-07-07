#include <iostream>
#include <string>
#include <cstdlib>
#include "../libs/depressions/priority_flood.hpp"
#include "../libs/common/Array2D.hpp"
#include "../libs/flats/flat_resolution.hpp"
#include "../libs/methods/d8_methods.hpp"

template<class T>
int PerformAlgorithm(std::string filename, std::string output_prefix){
  
  Array2D<T> elevations(filename,false);

  improved_priority_flood(elevations);

  Array2D<uint8_t> flowdirs;

  barnes_flat_resolution_d8(elevations,flowdirs,false);

  for(int y=0;y<flowdirs.viewHeight();y++){
    for(int x=0;x<flowdirs.viewWidth();x++)
      std::cerr<<(int)flowdirs(x,y)<<" ";
    std::cerr<<std::endl;
  }


  flowdirs.saveGDAL(output_prefix+"-flowdirs.tif", 0, 0);

  return 0;
}

int main(int argc, char **argv){
  if(argc!=3){
    std::cerr<<argv[0]<<" <INPUT> <OUTPUT_PREFIX>"<<std::endl;
    return -1;
  }

  auto data_type = peekGDALType(argv[1]);
  std::cerr<<"Loading data as type '"<<GDALGetDataTypeName(data_type)<<"'"<<std::endl;

  switch(data_type){
    case GDT_Unknown:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(data_type)<<std::endl;
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