#ifndef _richdem_router_hpp_
#define _richdem_router_hpp_

#include "gdal_priv.h"
#include "richdem/common/Array2D.hpp"

namespace richdem {

template<typename... Args>
int PerformAlgorithm(std::string inputfile, Args... args){
  switch(peekGDALType(inputfile)){
    case GDT_Byte:    {
      Array2D<uint8_t> arr(inputfile,false,0,0,0,0,false);
      return PerformAlgorithm(args..., arr);
    }
    case GDT_UInt16:  {
      Array2D<uint16_t> arr(inputfile,false,0,0,0,0,false);
      return PerformAlgorithm(args..., arr);
    }
    case GDT_Int16:   {
      Array2D<int16_t> arr(inputfile,false,0,0,0,0,false);
      return PerformAlgorithm(args..., arr);
    }
    case GDT_UInt32:  {
      Array2D<uint32_t> arr(inputfile,false,0,0,0,0,false);
      return PerformAlgorithm(args..., arr);
    }
    case GDT_Int32:   {
      Array2D<int32_t> arr(inputfile,false,0,0,0,0,false);
      return PerformAlgorithm(args..., arr);
    }
    case GDT_Float32: {
      Array2D<float> arr(inputfile,false,0,0,0,0,false);
      return PerformAlgorithm(args..., arr);
    }
    case GDT_Float64: {
      Array2D<double> arr(inputfile,false,0,0,0,0,false);
      return PerformAlgorithm(args..., arr);
    }
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

}

#endif
