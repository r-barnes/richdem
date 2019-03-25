#ifndef _richdem_loaders_hpp_
#define _richdem_loaders_hpp_

#include <richdem/common/Array2D.hpp>
#include <richdem/common/gdal.hpp>

namespace richdem {

#ifdef USEGDAL
template<class T>
void LoadGDAL(
  const std::string &filename,
  Array2D<T> &arr, 
  uint64_t padx        = 0,
  uint64_t pady        = 0,
  uint64_t xOffset     = 0, 
  uint64_t yOffset     = 0, 
  uint64_t part_width  = 0, 
  uint64_t part_height = 0
){
  arr.clear();
  arr.filename = filename;

  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  if(fin==NULL)
    throw std::runtime_error("Could not open file '"+filename+"' with GDAL!");

  arr.geotransform.resize(6);
  if(fin->GetGeoTransform(arr.geotransform.data())!=CE_None){
    arr.geotransform = {{1000., 1., 0., 1000., 0., -1.}};
  }

  arr.metadata = ProcessMetadata(fin->GetMetadata());

  const char* projection_string=fin->GetProjectionRef();
  arr.projection = std::string(projection_string);

  GDALRasterBand *band = fin->GetRasterBand(1);

  uint64_t total_width  = band->GetXSize();         //Returns an int
  uint64_t total_height = band->GetYSize();         //Returns an int
  arr.setNoData((T)band->GetNoDataValue());

  //Ensure that the block we're grabbing doesn't extend outside the data
  if(xOffset+part_width>=total_width)
    part_width  = total_width-xOffset;
  if(yOffset+part_height>=total_height)
    part_height = total_height-yOffset;

  int view_width, view_height;
  if(part_width==0)
    part_width = total_width;
  view_width = part_width;

  if(part_height==0)
    part_height = total_height;
  view_height = part_height;

  arr.resize(view_width+2*padx,view_height+2*pady);
  auto temp = band->RasterIO( GF_Read, xOffset, yOffset, view_width, view_height, arr.data()+pady*view_width, view_width, view_height, NativeTypeToGDAL<T>(), 0, 0 );
  if(temp!=CE_None)
    throw std::runtime_error("An error occured while trying to read '"+filename+"' into RAM with GDAL.");

  if(padx>0){
    for(int y=0;y<arr.height();y++){
      for(int x=arr.width()-2*padx-1;x>=0;x--) //Shift cells to the right;
        arr(x+padx,y) = arr(x,y);
      for(unsigned int x=0;x<padx;x++)                  //Set halo pixels to 0
        arr(x,y) = 0;
    }
  }

  GDALClose(fin);

  //TODO
  // if(exact && (total_width-xOffset!=part_width || total_height-yOffset!=part_height))
    // throw std::runtime_error("Tile dimensions did not match expectations!");
}
#endif

}

#endif