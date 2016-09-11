#include <iostream>
#include <iomanip>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"

template<class T>
int PerformAlgorithm(
  std::string inputfile,
  std::string outputfile,
  std::string geo1,
  std::string geo2,
  std::string geo3,
  std::string geo4,
  std::string geo5,
  std::string geo6,
  std::string analysis
){
  Array2D<T> raster(inputfile,false);

  if(geo1!="x") raster.geotransform[0] = std::stod(geo1);
  if(geo2!="x") raster.geotransform[1] = std::stod(geo2);
  if(geo3!="x") raster.geotransform[2] = std::stod(geo3);
  if(geo4!="x") raster.geotransform[3] = std::stod(geo4);
  if(geo5!="x") raster.geotransform[4] = std::stod(geo5);
  if(geo6!="x") raster.geotransform[5] = std::stod(geo6);

  raster.saveGDAL(outputfile,analysis);

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

  if(argc==2){
    Array2D<int8_t> temp(argv[1],false,0,0,0,0,false,false); //Data type doesn't matter since we're not loading it
    for(auto v: temp.geotransform)
      std::cout<<std::setw(10)<<v<<" ";
    std::cout<<std::endl;
    return 0;
  }

  if(argc!=9){
    std::cerr<<argv[0]<<" <Display geotransform of this file>"<<std::endl;
    std::cerr<<argv[0]<<" <Input file> <Output File> <Geo1> <Geo2> <Geo3> <Geo4> <Geo5> <Geo6>"<<std::endl;
    return -1;
  }

  Router(argv[1],argv[1],argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],analysis);

  return 0;
}