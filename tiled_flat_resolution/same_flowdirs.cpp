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

typedef uint8_t flowdirs_t;
typedef int8_t  visited_t;

const flowdirs_t UNVISITED = 13;

const int fwidth=10;

template<class T>
void Master(std::string layoutfile, int cachesize, std::string merged_name, std::string flip_style){
  Timer total_time;
  total_time.start();

  A2Array2D<T> dem (layoutfile,cachesize);
  Array2D<T>   mdem(merged_name,false);

  if( (mdem.geotransform[1]<0) ^ (flip_style=="fliph" || flip_style=="fliphv")){
    std::cerr<<"Flipping mdem horizontally."<<std::endl;
    mdem.flipHorz();
  }
  if( (mdem.geotransform[5]>0) ^ (flip_style=="flipv" || flip_style=="fliphv")){
    std::cerr<<"Flipping mdem vertically."<<std::endl;
    mdem.flipVert();
  }

  if(flip_style=="fliph" || flip_style=="fliphv")
    dem.flipH = true;
  if(flip_style=="flipv" || flip_style=="fliphv")
    dem.flipV = true;

  if(dem.width()!=mdem.width())
    std::cerr<<"Widths differ dem="<<dem.width()<<", mdem="<<mdem.width()<<std::endl;
  if(dem.height()!=mdem.height())
    std::cerr<<"Heights differ dem="<<dem.height()<<", mdem="<<mdem.height()<<std::endl;

  int processed_cells = 0;

  for(int32_t ty=0;ty<dem.heightInTiles();ty++)
  for(int32_t tx=0;tx<dem.widthInTiles(); tx++){
    if(dem.isNullTile(tx,ty))
      continue;

    //std::cerr<<"Tile ("<<tx<<","<<ty<<") has dimensions "<<dem.tileHeight(tx,ty)<<" "<<dem.tileWidth(tx,ty)<<std::endl;

    int total_tiles       = dem.heightInTiles() * dem.widthInTiles();
    int processed_tiles   = ty*dem.widthInTiles()+tx;
    double est_total_time = (total_time.lap()/(double)processed_tiles)*(double)total_tiles;
    double time_left      = est_total_time-total_time.lap();
    std::cerr<<"Processed: "<<processed_tiles<<" of "<<total_tiles<<" tiles "
             <<time_left<<"s/"<<est_total_time<<"s ("<<(time_left/3600)<<"hr/"<<(est_total_time/3600)<<"hr)"<<std::endl;
    for(int py=0;py<dem.tileHeight(tx,ty);py++)
    for(int px=0;px<dem.tileWidth(tx,ty); px++){

      const int y = ty*dem.stdTileHeight()+py;
      const int x = tx*dem.stdTileWidth() +px;

      if(dem(x,y)!=mdem(x,y))
        std::cerr<<"Difference at ("<<x<<","<<y<<")"<<std::endl;

      processed_cells++;
    }
  }

  total_time.stop();

  std::cerr<<"Processed cells: "<<processed_cells<<std::endl;

  std::cerr<<"Total time: "<<total_time.accumulated()<<"s ("<<(total_time.accumulated()/3600)<<"hr)"<<std::endl;

  std::cerr<<"dem evictions: "<<dem.getEvictions()<<std::endl;
}

int main(int argc, char **argv){
  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <Layout File> <Cache size> <Merged File>"<<std::endl;
    std::cerr<<"\tor use 'table' for cache size"<<std::endl;
    return -1;
  }
  auto file_type         = peekLayoutType(argv[1]);
  std::string flip_style = argv[5];

  if(argv[2]==std::string("table")){
    int dtype_size = GDALGetDataTypeSizeBytes(file_type);
    int tile_size  = peekLayoutTileSize(argv[1]);
    for(int i=2;i<500;i++)
      std::cerr<<std::setw(2)<<i<<"  "<<((dtype_size+sizeof(flowdirs_t)+sizeof(visited_t))*tile_size*i/1000000.0)<<" MB"<<std::endl;
    return -1;
  }

  int cachesize = std::stoi(argv[2]);
  std::string merged_filename(argv[3]);

  switch(file_type){
    case GDT_Byte:
      Master<uint8_t >(argv[1],cachesize,merged_filename,flip_style);break;
    case GDT_UInt16:
      Master<uint16_t>(argv[1],cachesize,merged_filename,flip_style);break;
    case GDT_Int16:
      Master<int16_t >(argv[1],cachesize,merged_filename,flip_style);break;
    case GDT_UInt32:
      Master<uint32_t>(argv[1],cachesize,merged_filename,flip_style);break;
    case GDT_Int32:
      Master<int32_t >(argv[1],cachesize,merged_filename,flip_style);break;
    case GDT_Float32:
      Master<float   >(argv[1],cachesize,merged_filename,flip_style);break;
    case GDT_Float64:
      Master<double  >(argv[1],cachesize,merged_filename,flip_style);break;
    default:
      std::cerr<<"Unrecognised data type!"<<std::endl;
      return -1;
  }

  return 0;
}