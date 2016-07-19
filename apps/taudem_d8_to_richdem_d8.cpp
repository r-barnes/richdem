#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/Array2D.hpp"

template<class T>
void PerformAlgorithm(std::string filename, std::string outname, char *nodata){
  const int taudem_to_richdem[9] = {0,5,4,3,2,1,8,7,6};
  Array2D<T> inp(filename,false);
  for(int y=0;y<inp.height();y++)
  for(int x=0;x<inp.width();x++)
    if(!inp.isNoData(x,y)){
      if(!(0<=inp(x,y) && inp(x,y)<=8)){
        std::cerr<<"Invalid flow direction '"<<inp(x,y)<<"' found!"<<std::endl;
        return;
      }
      inp(x,y) = taudem_to_richdem[(int)inp(x,y)];
    }
  inp.saveGDAL(outname,0,0);
  return;
}

int main(int argc, char **argv){
  if(argc!=3){
    std::cerr<<argv[0]<<" <Input> <Output name>"<<std::endl;
    return -1;
  }

  switch(peekGDALType(argv[1])){
    case GDT_Unknown:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[1]))<<std::endl;
      return -1;
    case GDT_Byte:
      PerformAlgorithm<uint8_t >(argv[1],argv[2],argv[3]);break;
    case GDT_UInt16:
      PerformAlgorithm<uint16_t>(argv[1],argv[2],argv[3]);break;
    case GDT_Int16:
      PerformAlgorithm<int16_t >(argv[1],argv[2],argv[3]);break;
    case GDT_UInt32:
      PerformAlgorithm<uint32_t>(argv[1],argv[2],argv[3]);break;
    case GDT_Int32:
      PerformAlgorithm<int32_t >(argv[1],argv[2],argv[3]);break;
    case GDT_Float32:
      PerformAlgorithm<float   >(argv[1],argv[2],argv[3]);break;
    case GDT_Float64:
      PerformAlgorithm<double  >(argv[1],argv[2],argv[3]);break;
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