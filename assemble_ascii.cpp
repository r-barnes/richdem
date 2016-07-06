#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "Array2D.hpp"

template<class T>
void Process(std::string path, int gridwidth, int gridheight, int tile_width, int tile_height){
  gridheight++;
  gridwidth++;
  
  std::vector< std::vector< Array2D<T> > > grid;
  grid.resize(gridheight, std::vector< Array2D<T> >(gridwidth));
  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++){
    std::string temp = path;
    temp.replace(temp.find("%n"), 2, std::to_string(x)+"_"+std::to_string(y));
    grid[y][x] = Array2D<T>(temp, false, 0, 0, 0, 0, false);
  }

  for(int y=0;y<gridheight*tile_height;y++){
    for(int x=0;x<gridwidth*tile_width;x++){
      int tx = x/tile_width;
      int ty = y/tile_height;
      int px = x%tile_width;
      int py = y%tile_height;
      std::cout<<std::setw(4)<<(int)grid[ty][tx](px,py);
    }
    std::cout<<std::endl;
  }
}

int main(int argc, char **argv){
  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <PATH> <WIDTH> <HEIGHT>"<<std::endl;
    return -1;
  }

  std::string path     = argv[1];
  const int gridwidth  = std::stoi(argv[2]);
  const int gridheight = std::stoi(argv[3]);

  int tile_width;
  int tile_height;
  GDALDataType data_type;
  std::vector<double> geotransform(6);

  std::string temp = path;
  temp.replace(temp.find("%n"), 2, "0_0");

  getGDALDimensions(temp, tile_height, tile_width, data_type, geotransform.data());

  switch(data_type){
    case GDT_Byte:
      Process<uint8_t >(path, gridwidth, gridheight, tile_width, tile_height);break;
    case GDT_UInt16:
      Process<uint16_t>(path, gridwidth, gridheight, tile_width, tile_height);break;
    case GDT_Int16:
      Process<int16_t >(path, gridwidth, gridheight, tile_width, tile_height);break;
    case GDT_UInt32:
      Process<uint32_t>(path, gridwidth, gridheight, tile_width, tile_height);break;
    case GDT_Int32:
      Process<int32_t >(path, gridwidth, gridheight, tile_width, tile_height);break;
    case GDT_Float32:
      Process<float   >(path, gridwidth, gridheight, tile_width, tile_height);break;
    case GDT_Float64:
      Process<double  >(path, gridwidth, gridheight, tile_width, tile_height);break;
    default:
      return -1;
  }
}