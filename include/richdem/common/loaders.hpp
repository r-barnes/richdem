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

  uint64_t view_width, view_height;
  if(part_width==0)
    part_width = total_width;
  view_width = part_width;

  if(part_height==0)
    part_height = total_height;
  view_height = part_height;

  T *buf;
  std::vector<T> tempbuf;
  if(padx==0 && pady==0){
    arr.resize(view_width,view_height);
    buf = arr.data();
  } else {
    tempbuf.resize(view_width*view_height);
    arr.resize(view_width+2*padx,view_height+2*pady);
    buf = tempbuf.data();
  }

  auto temp = band->RasterIO( GF_Read, xOffset, yOffset, view_width, view_height, buf, view_width, view_height, NativeTypeToGDAL<T>(), 0, 0 );
  if(temp!=CE_None)
    throw std::runtime_error("An error occured while trying to read '"+filename+"' into RAM with GDAL.");

  if(padx>0 || pady>0){
    for(int y=0;y<view_height;y++)
    for(int x=0;x<view_width;x++)
      arr(x+padx,y+pady) = buf[y*view_width+x];
  }

  // if(padx>0){
  //   for(int y=0;y<arr.height();y++){
  //     for(int x=arr.width()-2*padx-1;x>=0;x--) //Shift cells to the right;
  //       arr(x+padx,y) = arr(x,y);
  //     for(unsigned int x=padx+arr.width();x<arr.width();x++)         //Set halo pixels to 0
  //       arr(x,y) = 0;
  //   }
  // }

  GDALClose(fin);

  //TODO
  // if(exact && (total_width-xOffset!=part_width || total_height-yOffset!=part_height))
    // throw std::runtime_error("Tile dimensions did not match expectations!");
}



template<class T>
void SaveGDAL(const Array2D<T> &arr, const std::string &filename, const std::string &metadata_str="", uint64_t xoffset=0, uint64_t yoffset=0, const uint64_t padx=0, const uint64_t pady=0, bool compress=false){
  char **papszOptions = NULL;
  if(compress){
    papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "DEFLATE" );
    papszOptions = CSLSetNameValue( papszOptions, "ZLEVEL",   "6" );
  }

  GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
  if(poDriver==NULL)
    throw std::runtime_error("Could not open GDAL driver!");
  GDALDataset *fout    = poDriver->Create(filename.c_str(), arr.width(), arr.height(), 1, NativeTypeToGDAL<T>(), papszOptions);
  if(fout==NULL)
    throw std::runtime_error("Could not open file '"+filename+"' for GDAL save!");

  GDALRasterBand *oband = fout->GetRasterBand(1);
  oband->SetNoDataValue(arr.noData());

  //This could be used to copy metadata
  //poDstDS->SetMetadata( poSrcDS->GetMetadata() );

  //TIFFTAG_SOFTWARE
  //TIFFTAG_ARTIST
  {
    std::time_t the_time = std::time(nullptr);
    char time_str[64];
    std::strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S UTC", std::gmtime(&the_time));
    fout->SetMetadataItem("TIFFTAG_DATETIME",   time_str);
    fout->SetMetadataItem("TIFFTAG_SOFTWARE",   program_identifier.c_str());

    std::string proc_hist;
    if(arr.metadata.count("PROCESSING_HISTORY")>0)
      proc_hist = arr.metadata.at("PROCESSING_HISTORY");
    
    if(!proc_hist.empty())
      proc_hist+="\n";

    proc_hist += std::string(time_str) + " | " + program_identifier + " | ";
    if(!metadata_str.empty())
      proc_hist += metadata_str;
    else
      proc_hist += "Unspecified Operation";

    fout->SetMetadataItem("PROCESSING_HISTORY", proc_hist.c_str());
  }

  for(const auto &kv: arr.metadata){
    if(kv.first!="PROCESSING_HISTORY")
      fout->SetMetadataItem(kv.first.c_str(), kv.second.c_str());
  }

  //The geotransform maps each grid cell to a point in an affine-transformed
  //projection of the actual terrain. The geostransform is specified as follows:
  //    Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
  //    Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
  //In case of north up images, the GT(2) and GT(4) coefficients are zero, and
  //the GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3))
  //position is the top left corner of the top left pixel of the raster.

  if(!arr.geotransform.empty()){
    auto out_geotransform = arr.geotransform;

    if(out_geotransform.size()!=6)
      throw std::runtime_error("Geotransform of output is not the right size. Found "+std::to_string(out_geotransform.size())+" expected 6.");

    //We shift the top-left pixel of hte image eastward to the appropriate
    //coordinate
    out_geotransform[0] += xoffset*arr.geotransform[1];

    //We shift the top-left pixel of the image southward to the appropriate
    //coordinate
    out_geotransform[3] += yoffset*arr.geotransform[5];

    fout->SetGeoTransform(out_geotransform.data());
  }

  if(!arr.projection.empty())
    fout->SetProjection(arr.projection.c_str());

  #ifdef DEBUG
    RDLOG_DEBUG<<"Filename: "<<std::setw(20)<<filename<<" Xoffset: "<<std::setw(6)<<xoffset<<" Yoffset: "<<std::setw(6)<<yoffset<<" Geotrans0: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotransform[0]<<" Geotrans3: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotransform[3];
  #endif


  T *buf       = const_cast<T*>(arr.data());
  auto owidth  = arr.width();
  auto oheight = arr.height();
  std::vector<T> temp_buf;
  if(padx>0 || pady>0){
    const uint64_t pwidth = arr.width();
    owidth          -= 2*padx;
    oheight         -= 2*pady;
    temp_buf.resize(owidth*oheight);
    for(auto y=pady;y<arr.height()-pady;y++)
    for(auto x=padx;x<arr.width()-padx;x++)
      temp_buf[(y-pady)*owidth+(x-padx)] = buf[y*pwidth+x];
    buf = temp_buf.data();
  }

  auto temp = oband->RasterIO(GF_Write, 0, 0, owidth, oheight, buf, owidth, oheight, NativeTypeToGDAL<T>(), 0, 0);
  if(temp!=CE_None)
    throw std::runtime_error("Error writing file with saveGDAL()!");

  GDALClose(fout);
}

#endif

}

#endif