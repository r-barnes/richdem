#include <iostream>
#include <iomanip>
#include <unordered_set>
#include <richdem/common/version.hpp>
#include <richdem/common/Array2D.hpp>
using namespace richdem;

template<class flowdir_t>
void FollowPath(
  const int x0,                      //x-coordinate of initial cell
  const int y0,                      //y-coordinate of initial cell
  const Array2D<flowdir_t> &flowdirs //Flow directions matrix
){
  int x = x0;
  int y = y0;

  size_t path_len        = 0;
  size_t max_path_length = flowdirs.size(); //TODO: Should this have +1?
  max_path_length        = flowdirs.width();

  //Follow the flow path until it terminates
  while(path_len++<max_path_length){ //Follow the flow path until we reach its end
    const int n = flowdirs(x,y);     //Neighbour the current cell flows towards

    //Show the final part of the loop path (TODO)
    if(path_len>max_path_length-10)
      std::cout<<x<<","<<y<<" with flowdir "<<n<<"\n";

    if(flowdirs.isNoData(x,y) || n==NO_FLOW)
      return;

    const int nx = x+d8x[n]; //Get neighbour's x-coordinate.
    const int ny = y+d8y[n]; //Get neighbour's y-coordinate.

    if(!flowdirs.inGrid(nx,ny))
      return;

    x = nx;
    y = ny;
  }

  if(path_len>max_path_length-10)
    std::cout<<"\n\n";
}



template<class T>
int PerformAlgorithm(Array2D<T> inp){
  inp.loadData();

  for(int32_t y=0;y<inp.height();y++)
  for(int32_t x=0;x<inp.width();x++)
    FollowPath(x,y,inp);

  return 0;
}

#include "router.hpp"



int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);

  int32_t total_height;
  int32_t total_width;
  GDALDataType file_type;

  std::string inputfile;

  if(argc!=2){
    std::cerr<<argv[0]<<" <FILE>"<<std::endl;
    return -1;
  }

  inputfile = argv[1];

  //Get the total dimensions of the input file
  getGDALDimensions(inputfile, total_height, total_width, file_type, NULL);

  return PerformAlgorithm(inputfile);
}