#include "richdem/common/Layoutfile.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/constants.hpp"
#include "richdem/common/timer.hpp"
#include "richdem/tiled/A2Array2D.hpp"
#include <queue>
#include <unordered_set>
#include <iomanip>
#include <iostream>
#include <csignal>
#include <stdexcept>

template<class T>
void Master(std::string layoutfile, int cachesize, std::string outputname, std::string flip_style){
  Timer total_time;
  total_time.start();

  std::cerr<<"Loading results..."<<std::endl;
  A2Array2D<T> raster(layoutfile,cachesize);

  if(flip_style=="fliph" || flip_style=="fliphv")
    raster.flipH = true;
  if(flip_style=="flipv" || flip_style=="fliphv")
    raster.flipV = true;

  std::cerr<<"Saving results..."<<std::endl;

  raster.saveUnifiedGDAL(outputname);

  total_time.stop();

  std::cerr<<"Total time: "<<total_time.accumulated()<<"s ("<<(total_time.accumulated()/3600)<<"hr)"<<std::endl;

  std::cerr<<"Evictions: "<<raster.getEvictions()<<std::endl;
}

int main(int argc, char **argv){
  if(argc!=5){
    std::cerr<<"Syntax: "<<argv[0]<<" <Layout File> <Cache size> <Output File> <noflip/fliph/flipv/fliphv>"<<std::endl;
    std::cerr<<"\tor use 'table' for cache size"<<std::endl;
    return -1;
  }
  auto file_type         = peekLayoutType(argv[1]);

  if(argv[2]==std::string("table")){
    int dtype_size = GDALGetDataTypeSizeBytes(file_type);
    int tile_size  = peekLayoutTileSize(argv[1]);
    for(int i=2;i<500;i++)
      std::cerr<<std::setw(2)<<i<<"  "<<(dtype_size*tile_size*i/1000000.0)<<" MB"<<std::endl;
    return -1;
  }

  int cachesize = std::stoi(argv[2]);
  std::string output_filename(argv[3]);
  std::string flip_style(argv[4]);

  switch(file_type){
    case GDT_Byte:
      Master<uint8_t >(argv[1],cachesize,output_filename,flip_style);break;
    case GDT_UInt16:
      Master<uint16_t>(argv[1],cachesize,output_filename,flip_style);break;
    case GDT_Int16:
      Master<int16_t >(argv[1],cachesize,output_filename,flip_style);break;
    case GDT_UInt32:
      Master<uint32_t>(argv[1],cachesize,output_filename,flip_style);break;
    case GDT_Int32:
      Master<int32_t >(argv[1],cachesize,output_filename,flip_style);break;
    case GDT_Float32:
      Master<float   >(argv[1],cachesize,output_filename,flip_style);break;
    case GDT_Float64:
      Master<double  >(argv[1],cachesize,output_filename,flip_style);break;
    default:
      std::cerr<<"Unrecognised data type!"<<std::endl;
      return -1;
  }

  return 0;
}