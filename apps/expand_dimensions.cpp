#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/Array2D.hpp"

template<class T>
int PerformAlgorithm(std::string filename, std::string outname, int new_width, int new_height){
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

  inp.saveGDAL(outname,0,0);
  return 0;
}

int main(int argc, char **argv){
  if(argc!=5){
    std::cerr<<argv[0]<<" <Input> <Output name> <Width> <Height>"<<std::endl;
    return -1;
  }

  int new_width  = std::stoi(argv[3]);
  int new_height = std::stoi(argv[4]);

  switch(peekGDALType(argv[1])){
    case GDT_Unknown:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[1]))<<std::endl;
      return -1;
    case GDT_Byte:
      return PerformAlgorithm<uint8_t >(argv[1],argv[2],new_width,new_height);
    case GDT_UInt16:
      return PerformAlgorithm<uint16_t>(argv[1],argv[2],new_width,new_height);
    case GDT_Int16:
      return PerformAlgorithm<int16_t >(argv[1],argv[2],new_width,new_height);
    case GDT_UInt32:
      return PerformAlgorithm<uint32_t>(argv[1],argv[2],new_width,new_height);
    case GDT_Int32:
      return PerformAlgorithm<int32_t >(argv[1],argv[2],new_width,new_height);
    case GDT_Float32:
      return PerformAlgorithm<float   >(argv[1],argv[2],new_width,new_height);
    case GDT_Float64:
      return PerformAlgorithm<double  >(argv[1],argv[2],new_width,new_height);
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