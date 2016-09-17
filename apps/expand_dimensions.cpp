#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"

template<class T>
int PerformAlgorithm(std::string filename, std::string outname, int new_width, int new_height, std::string analysis){
  Array2D<T> inp(filename,false);

  if(new_width<inp.width()){
    std::cerr<<"Desired width is smaller than DEM's current width!"<<std::endl;
    return -1;
  }

  if(new_height<inp.height()){
    std::cerr<<"Desired height is smaller than DEM's current height!"<<std::endl;
    return -1;
  }

  inp.expand(new_width,new_height,inp.noData());

  inp.saveGDAL(outname,analysis);
  return 0;
}

template< typename... Arguments >
int Router(std::string inputfile, Arguments ... args){
  switch(peekGDALType(inputfile)){
    case GDT_Byte:
      return PerformAlgorithm<uint8_t >(args...);
    case GDT_UInt16:
      return PerformAlgorithm<uint16_t>(args...);
    case GDT_Int16:
      return PerformAlgorithm<int16_t >(args...);
    case GDT_UInt32:
      return PerformAlgorithm<uint32_t>(args...);
    case GDT_Int32:
      return PerformAlgorithm<int32_t >(args...);
    case GDT_Float32:
      return PerformAlgorithm<float   >(args...);
    case GDT_Float64:
      return PerformAlgorithm<double  >(args...);
    case GDT_CInt16:
    case GDT_CInt32:
    case GDT_CFloat32:
    case GDT_CFloat64:
      std::cerr<<"Complex types are unsupported. Sorry!"<<std::endl;
      return -1;
    case GDT_Unknown:
    default:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(inputfile))<<std::endl;
      return -1;
  }
}


int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);
  
  if(argc!=5){
    std::cerr<<argv[0]<<" <Input> <Output name> <Width> <Height>"<<std::endl;
    return -1;
  }

  int new_width  = std::stoi(argv[3]);
  int new_height = std::stoi(argv[4]);

  Router(argv[1],argv[1],argv[2],new_width,new_height,analysis);

  return 0;
}