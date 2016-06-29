#include "Layoutfile.hpp"
#include "Array2D.hpp"
#include "common.hpp"
#include <queue>
#include <unordered_set>
#include <iomanip>
#include <iostream>
#include <csignal>

typedef uint8_t flowdirs_t;
typedef int8_t  visited_t;

#define UNVISITED    13
#define FLOW_NO_DATA 255

const int fwidth=10;

template<class T>
void print2dradius(A2Array2D<T> &arr, int xcen, int ycen, int radius){
  return;
  int minx  = std::max(xcen-radius,0);
  int maxx  = std::min(xcen+radius,(int)(arr.width()-1));
  int miny  = std::max(ycen-radius,0);
  int maxy  = std::min(ycen+radius,(int)(arr.height()-1));
  std::cout<<std::setw(fwidth)<<" ";
  for(int x=minx;x<=maxx;x++)
    std::cout<<std::setw(fwidth)<<x;
  std::cout<<std::endl;

  for(int y=miny;y<=maxy;y++){
    std::cout<<std::setw(fwidth)<<y;
    for(int x=minx;x<=maxx;x++){
      if(xcen==x && ycen==y)
        std::cout<<"\033[93m";

      if(NativeTypeToGDAL<T>()==GDT_Byte)
        std::cout<<std::setw(10)<<(int)arr(x,y);
      else
        std::cout<<std::setw(10)<<arr(x,y);

      if(xcen==x && ycen==y)
        std::cout<<"\033[39m";
    }
    std::cout<<std::endl;
  }
}

template<class T>
void InspectPoint(std::string header, A2Array2D<T> &dem, A2Array2D<flowdirs_t> &fds, int x, int y){
  return;
  if(x!=dem.width()-494 || y!=dem.height()-373)
    return;
  std::cerr<<header<<std::endl;
  print2dradius(dem, x, y, 5);
  print2dradius(fds, x, y, 5);
}

template<class T>
void ProcessFlat(
  A2Array2D<T>          &dem,
  A2Array2D<flowdirs_t> &fds,
  const int x0,
  const int y0
){
  std::queue< std::pair<int, int> > q;

  const auto flat_height = dem(x0,y0);

  InspectPoint("ProcessFlatTop",dem,fds,x0,y0);

  q.emplace(x0,y0);
  while(!q.empty()){
    const auto c = q.front();
    q.pop();

    for(int n=1;n<=8;n++){
      int nx = c.first +dx[n];
      int ny = c.second+dy[n];

      InspectPoint("ProcessFlatN",dem,fds,nx,ny);

      if(!dem.in_grid(nx,ny) || dem(nx,ny)<flat_height || dem.isNoData(nx,ny)){
        fds(c.first,c.second) = n;
        continue;
      }
      if(fds(nx,ny)!=UNVISITED || dem(nx,ny)!=flat_height)
        continue;
      fds(nx,ny) = d8_inverse[n];
      q.emplace(nx,ny);
    }
  }
}

template<class T>
void Master(std::string layoutfile, int cachesize, std::string tempfile_name, std::string output_filename){
  Timer total_time;
  total_time.start();

  std::string temp_fds_name = tempfile_name;
  temp_fds_name.replace(temp_fds_name.find("%f"), 2, "%f-fds");

  A2Array2D<T>          dem(layoutfile,cachesize);
  A2Array2D<flowdirs_t> fds(temp_fds_name,dem,cachesize);

  fds.setNoData(FLOW_NO_DATA);
  fds.setAll(UNVISITED);

  int processed_cells = 0;

  for(int ty=0;ty<dem.heightInTiles();ty++)
  for(int tx=0;tx<dem.widthInTiles(); tx++){
    if(dem.isNullTile(tx,ty))
      continue;

    std::cerr<<"Tile ("<<tx<<","<<ty<<") has dimensions "<<dem.tileHeight()<<" "<<dem.tileWidth()<<std::endl;

    int total_tiles       = dem.heightInTiles() * dem.widthInTiles();
    int processed_tiles   = ty*dem.widthInTiles()+tx;
    double est_total_time = (total_time.lap()/(double)processed_tiles)*(double)total_tiles;
    double time_left      = est_total_time-total_time.lap();
    std::cerr<<"Processed: "<<processed_tiles<<" of "<<total_tiles<<" tiles "
             <<time_left<<"s/"<<est_total_time<<"s ("<<(time_left/3600)<<"hr/"<<(est_total_time/3600)<<"hr)"<<std::endl;
    for(int py=0;py<dem.tileHeight();py++)
    for(int px=0;px<dem.tileWidth(); px++){

      const int y = ty*dem.tileHeight()+py;
      const int x = tx*dem.tileWidth() +px;

      processed_cells++;

      if(fds(tx,ty,px,py)!=UNVISITED)
        continue;

      if(dem.isNoData(tx,ty,px,py)){
        fds(tx,ty,px,py) = FLOW_NO_DATA; //TODO: global var
        continue;
      }

      const auto myelev = dem(tx,ty,px,py);

      bool    drains       = false;
      bool    has_flat     = false;
      uint8_t nlowest      = 0;
      T       nlowest_elev = std::numeric_limits<T>::max(); //TODO
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

      if(nlowest!=0)
        fds(tx,ty,px,py) = nlowest;

      InspectPoint("MainBottom",dem,fds,x,y);

      if(drains && has_flat)
        ProcessFlat(dem,fds,x,y);
    }
  }

  std::cerr<<"Saving results..."<<std::endl;
  fds.saveGDAL(output_filename, 0);

  total_time.stop();

  std::cerr<<"Processed cells: "<<processed_cells<<std::endl;

  std::cerr<<"Total time: "<<total_time.accumulated()<<"s ("<<(total_time.accumulated()/3600)<<"hr)"<<std::endl;

  std::cerr<<"dem evictions: "<<dem.getEvictions()<<std::endl;
  std::cerr<<"fds evictions: "<<fds.getEvictions()<<std::endl;
}

//TODO: Flip tiles where appropriate
int main(int argc, char **argv){
  if(argc!=5){
    std::cerr<<"Syntax: "<<argv[0]<<" <Layout File> <Cache size> <Temp Files> <Output Files>"<<std::endl;
    std::cerr<<"\tor use 'table' for cache size"<<std::endl;
    return -1;
  }
  auto file_type = peekLayoutType(argv[1]);

  if(argv[2]==std::string("table")){
    int dtype_size = GDALGetDataTypeSizeBytes(file_type);
    int tile_size  = peekLayoutTileSize(argv[1]);
    for(int i=2;i<500;i++)
      std::cerr<<std::setw(2)<<i<<"  "<<((dtype_size+sizeof(flowdirs_t)+sizeof(visited_t))*tile_size*i/1000000.0)<<" MB"<<std::endl;
    return -1;
  }

  int cachesize = std::stoi(argv[2]);
  std::string tempfile_name(argv[3]);
  std::string output_filename(argv[4]);

  switch(file_type){
    case GDT_Byte:
      Master<uint8_t >(argv[1],cachesize,tempfile_name,output_filename);break;
    case GDT_UInt16:
      Master<uint16_t>(argv[1],cachesize,tempfile_name,output_filename);break;
    case GDT_Int16:
      Master<int16_t >(argv[1],cachesize,tempfile_name,output_filename);break;
    case GDT_UInt32:
      Master<uint32_t>(argv[1],cachesize,tempfile_name,output_filename);break;
    case GDT_Int32:
      Master<int32_t >(argv[1],cachesize,tempfile_name,output_filename);break;
    case GDT_Float32:
      Master<float   >(argv[1],cachesize,tempfile_name,output_filename);break;
    case GDT_Float64:
      Master<double  >(argv[1],cachesize,tempfile_name,output_filename);break;
    default:
      std::cerr<<"Unrecognised data type!"<<std::endl;
      return -1;
  }

  return 0;
}