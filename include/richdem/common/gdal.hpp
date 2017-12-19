#ifndef _richdem_gdal_hpp_
#define _richdem_gdal_hpp_

#ifdef USEGDAL

#include "gdal_priv.h"
#include <iostream>
#include <typeinfo>

namespace richdem {

/**
  @brief  Determine data type of a GDAL file's first layer
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  filename   Filename of file whose type should be determined
*/
GDALDataType peekGDALType(const std::string &filename) {
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  if(fin==NULL)
    throw std::runtime_error("Unable to open file '"+filename+"'!");

  GDALRasterBand *band   = fin->GetRasterBand(1);
  GDALDataType data_type = band->GetRasterDataType();

  GDALClose(fin);

  return data_type;
}



/**
  @brief  Retrieve height, width, NoData, and geotransform from a GDAL file
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  filename   GDAL file to peek at
  @param[out] height     Height of the raster in cells
  @param[out] width      Width of the raster in cells
  @param[out] no_data    Value of the raster's no_data property
  @param[out] geo_trans  Returns the SIX elements of the raster's geotransform
*/
template<class T>
void getGDALHeader(
  const   std::string &filename,
  int32_t &height,
  int32_t &width,
  T       &no_data,
  double  geotransform[6]
){
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  if(fin==NULL)
    throw std::runtime_error("Could not get GDAL header: file '" + filename + "'' did not open!");

  GDALRasterBand *band   = fin->GetRasterBand(1);

  height  = band->GetYSize();
  no_data = band->GetNoDataValue();
  width   = band->GetXSize();

  fin->GetGeoTransform(geotransform);

  GDALClose(fin);
}



/**
  @brief  Retrieve height, width, data type, and geotransform from a GDAL file
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  filename   GDAL file to peek at
  @param[out] height     Height of the raster in cells
  @param[out] width      Width of the raster in cells
  @param[out] dtype      Data type of the file in question
  @param[out] geo_trans  Returns the SIX elements of the raster's geotransform
*/
void getGDALDimensions(
  const   std::string &filename,
  int32_t &height,
  int32_t &width,
  GDALDataType &dtype,
  double geotransform[6]
){
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  if(fin==NULL)
    throw std::runtime_error("Could not open file '"+filename+"' to get dimensions.");

  GDALRasterBand *band = fin->GetRasterBand(1);
  
  dtype = band->GetRasterDataType();

  if(geotransform!=NULL && fin->GetGeoTransform(geotransform)!=CE_None)
    throw std::runtime_error("Error getting geotransform from '"+filename+"'.");

  height  = band->GetYSize();
  width   = band->GetXSize();

  GDALClose(fin);
}



/**
  @brief  Convert Array2D or any other template to its GDAL data type
  @author Richard Barnes (rbarnes@umn.edu)

  @return The GDAL datatype of T
*/
template<class T>
GDALDataType NativeTypeToGDAL() {
  if(typeid(T)==typeid(uint8_t))
    return GDT_Byte;
  else if(typeid(T)==typeid(uint16_t))
    return GDT_UInt16;
  else if(typeid(T)==typeid(int16_t))
    return GDT_Int16;
  else if(typeid(T)==typeid(uint32_t))
    return GDT_UInt32;
  else if(typeid(T)==typeid(int32_t))
    return GDT_Int32;
  else if(typeid(T)==typeid(float))
    return GDT_Float32;
  else if(typeid(T)==typeid(double))
    return GDT_Float64;
  else
    throw std::runtime_error("Could not map native type '"+std::string(typeid(T).name())+"' to GDAL type! (Use `c++filt -t` to decode.)");
  return GDT_Unknown;
}

}

#endif

#endif
