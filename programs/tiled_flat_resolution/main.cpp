#include "richdem/common/version.hpp"
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

template<class T>
void ProcessFlat(
  A2Array2D<T>          &dem,
  A2Array2D<flowdirs_t> &fds,
  const int x0,
  const int y0
){
  std::queue< std::pair<int, int> > q;

  const T flat_height = dem(x0,y0);

  q.emplace(x0,y0);
  while(!q.empty()){
    const auto c = q.front();
    q.pop();

    for(int n=1;n<=8;n++){
      const int nx = c.first +dx[n];
      const int ny = c.second+dy[n];

      if(!dem.in_grid(nx,ny))
        continue;
      if(dem.isEdgeCell(nx,ny))
        continue;
      if(fds(nx,ny)!=NO_FLOW)
        continue;
      if(dem(nx,ny)!=flat_height)
        continue;

      if(dem(nx,ny)<flat_height || dem.isNoData(nx,ny)){
        fds(c.first,c.second) = n;
        continue;
      }

      fds(nx,ny) = d8_inverse[n];

      if(fds(c.first,c.second)==d8_inverse[fds(nx,ny)])
        std::cerr<<"Loop formed in flat resolution at ("<<c.first<<","<<c.second<<")"<<std::endl;

      q.emplace(nx,ny);
    }
  }
}

template<class T>
static int d8EdgeFlow(const A2Array2D<T> &elevations, const int x, const int y){
  if(x==0 && y==0)
    return 2;
  else if(x==0 && y==elevations.width()-1)
    return 8;
  else if(x==elevations.width()-1 && y==0)
    return 4;
  else if(x==elevations.width()-1 && y==elevations.height()-1)
    return 6;
  else if(x==0)
    return 1;
  else if(x==elevations.width()-1)
    return 5;
  else if(y==0)
    return 3;
  else if(y==elevations.height()-1)
    return 7;
  else
    std::cerr<<"Should never reach this point!"<<std::endl;
  throw std::runtime_error("Requested edge direction not on an edge!");
}


template<class T>
void Master(std::string layoutfile, int cachesize, std::string tempfile_name, std::string output_filename, std::string flip_style){
  Timer total_time;
  total_time.start();

  std::string temp_fds_name = tempfile_name;
  temp_fds_name.replace(temp_fds_name.find("%f"), 2, "%f-fds");

  A2Array2D<T>          dem(layoutfile,cachesize);

  if(flip_style=="fliph" || flip_style=="fliphv")
    dem.flipH = true;
  if(flip_style=="flipv" || flip_style=="fliphv")
    dem.flipV = true;

  dem.printStamp(5);

  A2Array2D<flowdirs_t> fds(temp_fds_name,dem,cachesize);

  fds.setNoData(FLOWDIR_NO_DATA);
  fds.setAll(NO_FLOW);

  int processed_cells = 0;
  int processed_tiles = 0;

  for(int32_t ty=0;ty<dem.heightInTiles();ty++)
  for(int32_t tx=0;tx<dem.widthInTiles(); tx++){
    if(dem.isNullTile(tx,ty))
      continue;

    double est_total_time = (total_time.lap()/(double)processed_tiles)*(double)dem.notNullTiles();
    double time_left      = est_total_time-total_time.lap();
    std::cerr<<"p Processed: "<<processed_tiles<<" of "<<dem.notNullTiles()<<" tiles "
             <<time_left<<"s/"<<est_total_time<<"s ("<<(time_left/3600)<<"hr/"<<(est_total_time/3600)<<"hr)"<<std::endl;

    processed_tiles++;

    for(int py=0;py<dem.tileHeight(tx,ty);py++)
    for(int px=0;px<dem.tileWidth(tx,ty); px++){

      const int y = ty*dem.stdTileHeight()+py;
      const int x = tx*dem.stdTileWidth() +px;

      processed_cells++;

      if(fds(tx,ty,px,py)!=NO_FLOW)
        continue;

      if(dem.isNoData(tx,ty,px,py)){
        fds(tx,ty,px,py) = FLOWDIR_NO_DATA;
        continue;
      }

      const auto myelev = dem(tx,ty,px,py);

      bool    drains       = false;
      bool    has_flat     = false;
      uint8_t nlowest      = 0;
      T       nlowest_elev = std::numeric_limits<T>::max();
      for(int n=1;n<=8;n++){
        const int nx = x+dx[n];
        const int ny = y+dy[n];
        if(!dem.in_grid(nx,ny) || dem.isNoData(nx,ny)){
          drains       = true;
          nlowest_elev = std::numeric_limits<T>::lowest();
          nlowest      = n;
          continue;
        }

        const auto nelev = dem(nx,ny);

        if(nelev==myelev){
          has_flat = true;
        } else if(nelev<myelev && nelev<nlowest_elev){
          drains       = true;
          nlowest_elev = nelev;
          nlowest      = n;
        }
      }

      if(nlowest!=0){
        fds(tx,ty,px,py) = nlowest;
        int nx = x+dx[nlowest];
        int ny = y+dy[nlowest];
        if(fds.in_grid(nx,ny) && fds(nx,ny)==d8_inverse[nlowest]){
          std::cerr<<"Two cell loop detected!"<<std::endl;
        }
      }

      if(dem.isEdgeCell(x,y))
        fds(tx,ty,px,py) = d8EdgeFlow(dem,x,y);

      if(drains && has_flat)
        ProcessFlat(dem,fds,x,y);
    }
  }

  // int no_flows = 0;
  // int loops    = 0;
  // for(int32_t ty=0;ty<dem.heightInTiles();ty++)
  // for(int32_t tx=0;tx<dem.widthInTiles(); tx++)
  // for(int py=0;py<dem.tileHeight(tx,ty);py++)
  // for(int px=0;px<dem.tileWidth(tx,ty); px++){

  //   const int y = ty*dem.stdTileHeight()+py;
  //   const int x = tx*dem.stdTileWidth() +px;

  //   auto my_fd = fds(tx,ty,px,py);
  //   if(my_fd==NO_FLOW){
  //     no_flows++;
  //     continue;
  //   }

  //   int nx = x+dx[my_fd];
  //   int ny = y+dy[my_fd];
  //   if(fds.in_grid(nx,ny) && my_fd==d8_inverse[fds(nx,ny)])
  //     loops++;
  // }
  // std::cerr<<"Found "<<loops<<" loops."<<std::endl;
  // std::cerr<<"Found "<<no_flows<<" cells with no flow."<<std::endl;

  fds.printStamp(5); //TODO

  std::cerr<<"p Saving results..."<<std::endl;
  fds.saveGDAL(output_filename);

  total_time.stop();

  std::cerr<<"m Processed cells = "<<processed_cells<<std::endl;

  std::cerr<<"t Total time = "<<total_time.accumulated()<<" s ("<<(total_time.accumulated()/3600)<<"hr)"<<std::endl;

  std::cerr<<"m dem evictions = "<<dem.getEvictions()<<std::endl;
  std::cerr<<"m fds evictions = "<<fds.getEvictions()<<std::endl;
}

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  std::cerr<<"A Tiled Flat Resolution (Greedy Resolution)"<<std::endl;
  std::cerr<<"C Barnes, R. 2016. RichDEM: Terrain Analysis Software. http://github.com/r-barnes/richdem"<<std::endl;
  if(argc!=6){
    std::cerr<<"Syntax: "<<argv[0]<<" <Layout File> <Cache size> <Temp Files> <Output Files> <noflip/fliph/flipv/fliphv>"<<std::endl;
    std::cerr<<"\tor use 'table' for cache size"<<std::endl;
    return -1;
  }
  auto file_type         = peekLayoutType(argv[1]);
  std::string flip_style = argv[5];

  if(argv[2]==std::string("table")){
    long dtype_size = GDALGetDataTypeSizeBytes(file_type);
    long tile_size  = peekLayoutTileSize(argv[1]);
    for(int i=2;i<500;i++)
      std::cerr<<std::setw(2)<<i<<"  "<<((dtype_size+sizeof(flowdirs_t))*tile_size*i/1000000.0)<<" MB"<<std::endl;
    return -1;
  }

  int cachesize = std::stoi(argv[2]);
  std::string tempfile_name(argv[3]);
  std::string output_filename(argv[4]);

  switch(file_type){
    case GDT_Byte:
      Master<uint8_t >(argv[1],cachesize,tempfile_name,output_filename,flip_style);break;
    case GDT_UInt16:
      Master<uint16_t>(argv[1],cachesize,tempfile_name,output_filename,flip_style);break;
    case GDT_Int16:
      Master<int16_t >(argv[1],cachesize,tempfile_name,output_filename,flip_style);break;
    case GDT_UInt32:
      Master<uint32_t>(argv[1],cachesize,tempfile_name,output_filename,flip_style);break;
    case GDT_Int32:
      Master<int32_t >(argv[1],cachesize,tempfile_name,output_filename,flip_style);break;
    case GDT_Float32:
      Master<float   >(argv[1],cachesize,tempfile_name,output_filename,flip_style);break;
    case GDT_Float64:
      Master<double  >(argv[1],cachesize,tempfile_name,output_filename,flip_style);break;
    default:
      std::cerr<<"Unrecognised data type!"<<std::endl;
      return -1;
  }

  return 0;
}