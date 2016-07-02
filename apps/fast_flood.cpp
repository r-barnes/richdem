#include <iostream>
#include <string>
#include <cstdlib>
#include "../libs/depressions/priority_flood.hpp"
#include "../libs/depressions/Zhou2015pf.hpp"
#include "../libs/common/Array2D.hpp"
//#include "../libs/flats/flat_resolution.hpp"
//#include "../libs/methods/d8_methods.hpp"

template<class T>
int PerformAlgorithm(std::string filename, std::string output_prefix){
  Array2D<T> elevation(filename,false);

  auto elevation2 = elevation;
  auto elevation3 = elevation;


  //improved_priority_flood(elevation);  

  improved_priority_flood_original(elevation2);

  Timer timer;
  timer.start();
  Zhou2015Labels(elevation3);
  timer.stop();

  std::cout<<"Zhou time: "<<timer.accumulated()<<std::endl;

  //std::cout<<"Pf matches: "<<((int)(elevation2==elevation))<<std::endl;
  std::cout<<"Zhou matches: "<<((int)(elevation3==elevation2))<<std::endl;

  return 0;
}

int main(int argc, char **argv){
  if(argc!=3){
    std::cerr<<argv[0]<<" <INPUT> <OUTPUT_PREFIX>"<<std::endl;
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