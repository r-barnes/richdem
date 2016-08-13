#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"

template<class T>
int PerformAlgorithm(std::string filename, std::string outname, char *nodata){
  Array2D<T> inp(filename,false);
  inp.setNoData((T)std::stoi(nodata));
  inp.saveGDAL(outname,0,0);
  return 0;
}

int main(int argc, char **argv){
  PrintRichdemHeader();
  
  if(argc!=4){
    std::cerr<<argv[0]<<" <Input> <Output name> <NoData>"<<std::endl;
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