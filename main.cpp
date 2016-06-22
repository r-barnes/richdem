#include "Layoutfile.hpp"
#include "Array2D.hpp"
#include "common.hpp"
#include <queue>
#include <unordered_set>

template<class T>
void ProcessFlat(A2Array2D<T> &dem, A2Array2D<uint8_t> &fds, A2Array2D<int8_t> &visited, int x, int y){
  std::queue< std::pair<int, int> > q;

  auto flat_height = dem(x,y);

  q.emplace(x,y);
  visited(x,y) = true;
  while(!q.empty()){
    auto c = q.front();
    q.pop();

    for(int n=1;n<=8;n++){
      int  nx = c.first+dx[n];
      int  ny = c.second+dy[n];
      if(!dem.in_grid(nx,ny))
        continue;
      if(visited(x,y))
        continue;
      if(dem(nx,ny)!=flat_height)
        continue;
      fds(nx,ny) = d8_inverse[n];
      q.emplace(nx,ny);
      visited(nx,ny);
    }
  }
}

template<class T>
void Master(std::string layoutfile){
  A2Array2D<T>       dem(layoutfile,30);
  A2Array2D<uint8_t> fds("/z/fdstemp",dem,30);
  A2Array2D<int8_t>  visited("/z/vistemp",dem,30);

  long total_cells     = dem.heightInTiles()*dem.widthInTiles()*dem.tileHeight()*dem.tileWidth();
  long processed_cells = 0;

  for(int ty=0;ty<dem.heightInTiles();ty++)
  for(int tx=0;tx<dem.widthInTiles(); tx++){
    std::cerr<<"Processed: "<<(processed_cells/1000)<<"k/"<<(total_cells/1000)<<"k"<<std::endl;
    for(int py=0;py<dem.tileHeight();py++)
    for(int px=0;px<dem.tileWidth(); px++){
      processed_cells++;

      int x = tx*dem.tileWidth()+px;
      int y = ty*dem.tileHeight()+py;

      if(visited(tx,ty,px,py))
        continue;

      auto myelev = dem(tx,ty,px,py);
      if(myelev==dem.noData()){
        fds(tx,ty,px,py) = 255;
        continue;
      }

      bool drains   = false;
      bool has_flat = false;
      for(int n=1;n<=8;n++){
        int nx = x+dx[n];
        int ny = y+dy[n];
        if(!dem.in_grid(nx,ny))
          continue;

        auto nelev = dem(nx,ny);
        if(nelev==dem.noData())
          continue;
        else if(nelev<myelev)
          drains = true;
        else if(nelev==myelev)
          has_flat = true;
      }

      if(drains && has_flat)
        ProcessFlat(dem,fds,visited,x,y);
    }
  }

  std::cerr<<"Saving results..."<<std::endl;
  fds.saveGDAL("/z/fds-out");
}

int main(int argc, char **argv){
  auto file_type = peekLayoutType(argv[1]);

  switch(file_type){
    //case GDT_Byte:
    //  Master<uint8_t >(argv[1]);break;
    // case GDT_UInt16:
    //   Master<uint16_t>(argv[1]);break;
    //case GDT_Int16:
    //  Master<int16_t >(argv[1]);break;
    // case GDT_UInt32:
    //   Master<uint32_t>(argv[1]);break;
    // case GDT_Int32:
    //   Master<int32_t >(argv[1]);break;
    case GDT_Float32:
      Master<float   >(argv[1]);break;
    // case GDT_Float64:
    //   Master<double  >(argv[1]);break;
    default:
      std::cerr<<"Unrecognised data type!"<<std::endl;
      return -1;
  }

  return 0;
}