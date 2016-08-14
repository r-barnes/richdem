#include <iostream>
#include <iomanip>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"

template<class T>
int PerformAlgorithm(std::string filename, std::string flip){
  bool flipH = false; //TODO
  bool flipV = false;

  if(flip=="fliph")
    flipH=true;
  else if(flip=="flipv")
    flipV=true;
  else if(flip=="fliphv")
    flipH=flipV=true;
  else if(flip!="noflip"){
    std::cerr<<"Unrecognised flip directive."<<std::endl;
    return -1;
  }


  Array2D<T> raster(filename,false);

  std::cerr<<"Geotransform: ";
  for(auto const x: raster.geotransform)
    std::cerr<<std::setw(6)<<std::setprecision(2)<<x<<" ";
  std::cerr<<std::endl;

  //Flip tiles if the geotransform demands it
  if( (raster.geotransform[1]<0) ^ flipH)
    raster.flipHorz();
  if( (raster.geotransform[5]>0) ^ flipV)
    raster.flipVert();

  raster.printStamp(5)

  for(int y=0;y<raster.height();y++){
    for(int x=0;x<raster.width();x++)
      std::cout<<std::setw(6)<<std::setprecision(2)<<(int)raster(x,y)<<" ";
    std::cout<<"\n";
  }

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
  
  if(argc!=3){
    std::cerr<<argv[0]<<" <Input file> <noflip/fliph/flipv/fliphv>"<<std::endl;
    return -1;
  }

  Router(argv[1],argv[1],argv[2]);

  return 0;
}