#include <iostream>
#include <cstdint>
#include <string>
#include "richdem/common/Array2D.hpp"
#include "richdem/flats/flat_resolution.hpp"
#include "richdem/flats/garbrecht.hpp"

template<class T>
int PerformAlgorithm(std::string alg, std::string filename, std::string output){
  Timer overall;
  overall.start();

  Array2D<uint8_t> flowdirs;
  Array2D<T> elevations(filename,false);
  d8_flow_directions(elevations,flowdirs);

  if(alg=="1"){
    resolve_flats_barnes(elevations,flowdirs,false);
  } else if(alg=="2"){
    GarbrechtAlg(elevations,flowdirs);
  } else if(alg=="3"){
    resolve_flats_barnes(elevations,flowdirs,true);
  } else {
    std::cout<<"Unknown algorithm!"<<std::endl;
    return -1;
  }

  flowdirs.saveGDAL(output, filename, 0, 0);

  std::cout<<"Algorithm took "<<overall.stop()<<" seconds overall."<<std::endl;

  return 0;
}


int main(int argc, char **argv){
  std::string output_prefix;

  if(argc!=3 && argc!=4){
    std::cout<<"Syntax: "<<argv[0]<<" <Algorithm> <Input> [Output]"<<std::endl;
    std::cout<<"Algorithm choices are:"<<std::endl;
    std::cout<<"\t1: Barnes Flat Resolution Algorithm"<<std::endl;
    std::cout<<"\t2: Barnes Flat Resolution Algorithm: Modify Elevations"<<std::endl;
    std::cout<<"\t3: Garbrecht & Martz Reference Algorithm (slow!)"<<std::endl;
    return -1;
  }

  if(argc==4)
    output_prefix = argv[3];
  else
    output_prefix = "out";

  switch(peekGDALType(argv[2])){
    case GDT_Unknown:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[2]))<<std::endl;
      return -1;
    case GDT_Byte:
      return PerformAlgorithm<uint8_t >(argv[1],argv[2],output_prefix);
    case GDT_UInt16:
      return PerformAlgorithm<uint16_t>(argv[1],argv[2],output_prefix);
    case GDT_Int16:
      return PerformAlgorithm<int16_t >(argv[1],argv[2],output_prefix);
    case GDT_UInt32:
      return PerformAlgorithm<uint32_t>(argv[1],argv[2],output_prefix);
    case GDT_Int32:
      return PerformAlgorithm<int32_t >(argv[1],argv[2],output_prefix);
    case GDT_Float32:
      return PerformAlgorithm<float   >(argv[1],argv[2],output_prefix);
    case GDT_Float64:
      return PerformAlgorithm<double  >(argv[1],argv[2],output_prefix);
    case GDT_CInt16:
    case GDT_CInt32:
    case GDT_CFloat32:
    case GDT_CFloat64:
      std::cerr<<"Complex types are unsupported. Sorry!"<<std::endl;
      return -1;
    default:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[2]))<<std::endl;
      return -1;
  }

  return 0;
}
