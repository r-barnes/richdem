#include <iostream>
#include <iomanip>
#include <unordered_map>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"

template<class T>
int PerformAlgorithm(std::string filename){
  std::unordered_map<T,int> counts;
  Array2D<T> rast(filename,false);
  for(int i=0;i<rast.size();i++)
    counts[rast(i)]++;

  std::cout<<"Nodata: "<<(int)rast.noData()<<std::endl;

  for(const auto &x: counts)
    std::cout<<std::setw(20)<<x.first<<" "<<std::setw(20)<<x.second<<"\n";

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
  std::string PrintRichdemHeader(argc,argv);
  
  if(argc!=2){
    std::cerr<<argv[0]<<" <Flowdirs input file>"<<std::endl;
    return -1;
  }

  Router(argv[1],argv[1]);

  return 0;
}