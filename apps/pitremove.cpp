#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/depressions/priority_flood.hpp"
#include "richdem/depressions/Zhou2015pf.hpp"
#include "richdem/common/Array2D.hpp"
//#include "richdem/flats/flat_resolution.hpp"
//#include "richdem/methods/d8_methods.hpp"

template<class T>
int PerformAlgorithm(std::string filename, std::string outputname, std::string analysis){
  Array2D<T> elevation(filename,false);

  improved_priority_flood(elevation);  

  elevation.saveGDAL(outputname,analysis);

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
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=3){
    std::cerr<<argv[0]<<" <Input> <Output name>"<<std::endl;
    return -1;
  }

  return Router(argv[1],argv[1],argv[2],analysis);
}