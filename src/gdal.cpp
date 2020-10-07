#include <richdem/common/gdal.hpp>

#ifdef USEGDAL

namespace richdem {

GDALDataType peekGDALType(const std::string &filename){
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  if(fin==NULL)
    throw std::runtime_error("Unable to open file '"+filename+"'!");

  GDALRasterBand *band   = fin->GetRasterBand(1);
  GDALDataType data_type = band->GetRasterDataType();

  GDALClose(fin);

  return data_type;
}

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

}

#endif