#include "Layoutfile.hpp"
#include "Array2D.hpp"
#include "common.hpp"
#include <queue>
#include <unordered_set>
#include <iomanip>
#include <iostream>

typedef uint8_t flowdirs_t;
typedef int8_t  visited_t;

template<class T>
void ProcessFlat(
  A2Array2D<T> &dem,
  A2Array2D<flowdirs_t> &fds,
  A2Array2D<visited_t> &visited, int x, int y
){
  std::queue< std::pair<int, int> > q;

  auto flat_height = dem(x,y);

  q.emplace(x,y);
  visited(x,y) = true;
  while(!q.empty()){
    auto c = q.front();
    q.pop();

    for(int n=1;n<=8;n++){
      int nx = c.first +dx[n];
      int ny = c.second+dy[n];
      if(!dem.in_grid(nx,ny))
        continue;
      if(visited(nx,ny))
        continue;
      if(dem(nx,ny)!=flat_height)
        continue;
      fds(nx,ny) = d8_inverse[n];
      q.emplace(nx,ny);
      visited(nx,ny) = true;
    }
  }
}

template<class T>
void Master(std::string layoutfile, int cachesize){
  A2Array2D<T>          dem(layoutfile,cachesize);
  A2Array2D<flowdirs_t> fds("/z/fdstemp",dem,cachesize);
  A2Array2D<visited_t>  visited("/z/vistemp",dem,cachesize);

  fds.setNoData(255);

  long total_cells     = dem.heightInTiles()*dem.widthInTiles()*dem.tileHeight()*dem.tileWidth();
  long processed_cells = 0;

  for(int ty=0;ty<dem.heightInTiles();ty++)
  for(int tx=0;tx<dem.widthInTiles(); tx++){
    std::cerr<<"Processed: "<<(processed_cells/1000)<<"k/"<<(total_cells/1000)<<"k"<<std::endl;
    for(int py=0;py<dem.tileHeight();py++)
    for(int px=0;px<dem.tileWidth(); px++){
      processed_cells++;

      const int y = ty*dem.tileHeight()+py;
      const int x = tx*dem.tileWidth() +px;

      if(visited(x,y))
        continue;

      auto myelev = dem(x,y);
      if(myelev==dem.noData()){
        fds(x,y) = 255;
        continue;
      }

      bool    drains       = false;
      bool    has_flat     = false;
      uint8_t nlowest      = 0;
      T       nlowest_elev = std::numeric_limits<T>::max(); //TODO
      for(int n=1;n<=8;n++){
        int nx = x+dx[n];
        int ny = y+dy[n];
        if(!dem.in_grid(nx,ny)){
          drains       = true;
          nlowest_elev = std::numeric_limits<T>::min();
          nlowest      = n;
          continue;
        }

        auto nelev = dem(nx,ny);
        if(nelev==dem.noData()){
          drains       = true;
          nlowest_elev = std::numeric_limits<T>::min();
          nlowest      = n;
        } else if(nelev==myelev){
          has_flat = true;
        } else if(nelev<myelev){
          drains = true;
          if(nelev<nlowest_elev){
            nlowest_elev = nelev;
            nlowest      = n;
          }
        }
      }

      fds(x,y) = nlowest;

      if(drains && has_flat)
        ProcessFlat(dem,fds,visited,x,y);
    }
  }

  std::cerr<<"Saving results..."<<std::endl;
  fds.saveGDAL("/z/fds-out");
}

//TODO: Flip tiles where appropriate
int main(int argc, char **argv){
  if(argc!=3){
    std::cerr<<"Syntax: "<<argv[0]<<" <Layout File> <Cache size>"<<std::endl;
    std::cerr<<"\tor use 'table' for cache size"<<std::endl;
    return -1;
  }
  auto file_type = peekLayoutType(argv[1]);

  if(argv[2]==std::string("table")){
    int dtype_size = GDALGetDataTypeSizeBytes(file_type);
    int tile_size  = peekLayoutTileSize(argv[1]);
    for(int i=2;i<30;i++)
      std::cerr<<std::setw(2)<<i<<"  "<<((dtype_size+sizeof(flowdirs_t)+sizeof(visited_t))*tile_size*i/1000000.0)<<" MB"<<std::endl;
    return -1;
  }


  int  cachesize = std::stoi(argv[2]);

  switch(file_type){
    case GDT_Byte:
      Master<uint8_t >(argv[1],cachesize);break;
    case GDT_UInt16:
      Master<uint16_t>(argv[1],cachesize);break;
    case GDT_Int16:
      Master<int16_t >(argv[1],cachesize);break;
    case GDT_UInt32:
      Master<uint32_t>(argv[1],cachesize);break;
    case GDT_Int32:
      Master<int32_t >(argv[1],cachesize);break;
    case GDT_Float32:
      Master<float   >(argv[1],cachesize);break;
    case GDT_Float64:
      Master<double  >(argv[1],cachesize);break;
    default:
      std::cerr<<"Unrecognised data type!"<<std::endl;
      return -1;
  }

  return 0;
}