#include <iostream>
#include <string>
#include <cstdlib>
#include "../libs/flats/flat_resolution_dinf.hpp"
#include "../libs/depressions/priority_flood.hpp"
#include "../libs/common/Array2D.hpp"
#include "../libs/methods/d8_methods.hpp"
#include "../libs/methods/dinf_methods.hpp"

template<class T>
int PerformAlgorithm(std::string filename, std::string output_prefix, float zscale){

  Array2D<T> elevations(filename,false);

  improved_priority_flood(elevations);

  //Save pit-filled DEM
  elevations.saveGDAL(output_prefix+"-filled.tif", 0, 0);

  Array2D<float> flowdirs;

  //Get Dinf flow directions
  dinf_flow_directions(elevations,flowdirs);

  //Resolve flats without altering elevations
  resolve_flats_barnes_dinf(elevations,flowdirs);

  flowdirs.saveGDAL(output_prefix+"-flowdirs.tif", 0, 0);

  Array2D<float> slope;
  d8_slope(elevations,slope,TATTRIB_SLOPE_RISERUN,zscale);

  slope.saveGDAL(output_prefix+"-slope.tif",0,0);

  Array2D<float> area;
  flowdirs.countDataCells();
  dinf_upslope_area(flowdirs,area);

  area.saveGDAL(output_prefix+"-area.tif",0,0);

  Array2D<float> cti;
  d8_CTI(area, slope, cti);

  cti.saveGDAL(output_prefix+"-cti.tif",0,0);

  return 0;
}

int main(int argc, char **argv){
  if(argc!=4){
    std::cerr<<argv[0]<<" <ZSCALE> <INPUT> <OUTPUT_PREFIX>"<<std::endl;
    return -1;
  }

  float zscale = std::strtof(argv[1],NULL);

  switch(peekGDALType(argv[2])){
    case GDT_Unknown:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[2]))<<std::endl;
      return -1;
    case GDT_Byte:
      return PerformAlgorithm<uint8_t >(argv[2],argv[3],zscale);
    case GDT_UInt16:
      return PerformAlgorithm<uint16_t>(argv[2],argv[3],zscale);
    case GDT_Int16:
      return PerformAlgorithm<int16_t >(argv[2],argv[3],zscale);
    case GDT_UInt32:
      return PerformAlgorithm<uint32_t>(argv[2],argv[3],zscale);
    case GDT_Int32:
      return PerformAlgorithm<int32_t >(argv[2],argv[3],zscale);
    case GDT_Float32:
      return PerformAlgorithm<float   >(argv[2],argv[3],zscale);
    case GDT_Float64:
      return PerformAlgorithm<double  >(argv[2],argv[3],zscale);
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