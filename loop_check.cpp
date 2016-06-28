#include <iostream>
#include <iomanip>
#include <unordered_set>
#include "Array2D.hpp"

#define NO_FLOW          0 //TODO: Explain and ensure it fits in flowdir_t

template<class flowdir_t>
void FollowPath(
  const int x0,                      //x-coordinate of initial cell
  const int y0,                      //y-coordinate of initial cell
  const Array2D<flowdir_t> &flowdirs //Flow directions matrix
){
  int x = x0;
  int y = y0;

  std::vector< std::pair<int, int> >   path;
  std::unordered_set< int > visited;

  if(flowdirs.isNoData(x,y))
    return;

  visited.emplace(y*flowdirs.viewWidth()+x);
  path.emplace_back(x,y);

  //Follow the flow path until it terminates
  while(true){ //Follow the flow path until we reach its end
    const int n = flowdirs(x,y);     //Neighbour the current cell flows towards

    //If the neighbour this cell flows into is a no_data cell or this cell does
    //not flow into any neighbour, then mark the initial cell from which we
    //began this flow path as terminating somewhere unimportant: its flow cannot
    //pass to neighbouring segments/nodes for further processing.
    if(flowdirs.isNoData(x,y) || n==NO_FLOW)
      return;

    //Flow direction was valid. Where does it lead?
    const int nx = x+dx[n]; //Get neighbour's x-coordinate.
    const int ny = y+dy[n]; //Get neighbour's y-coordinate.

    //The neighbour cell is off one of the sides of the tile. Therefore, its
    //flow may be passed on to a neighbouring tile. Thus, we need to link this
    //flow path's initial cell to this terminal cell.
    if(!flowdirs.in_grid(nx,ny))
      return;

    //The flow path has not yet terminated. Continue following it.
    x = nx;
    y = ny;

    int serial = y*flowdirs.viewWidth()+x;

    if(visited.count(serial)){
      bool print = false;
      for(auto &p: path){
        if(p.second*flowdirs.viewWidth()+p.first==serial)
          print = true;
        if(print)
          std::cerr<<p.first<<","<<p.second<<" "<<flowdirs(p.first,p.second)<<std::endl;
      }
      std::cerr<<x<<","<<y<<" "<<flowdirs(x,y)<<std::endl;
      std::cerr<<"\n"<<std::endl;
      return;
    }

    visited.emplace(y*flowdirs.viewWidth()+x);
    path.emplace_back(x,y);
  }
}

template<class T>
void Master(std::string input){
  Array2D<T> inp(input, 0, 0, 0, 0, false);

  for(int y=0;y<inp.viewHeight();y++)
  for(int x=0;x<inp.viewWidth();x++)
    FollowPath(x,y,inp);
}

int main(int argc, char **argv){
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
  if(getGDALDimensions(inputfile, total_height, total_width, file_type, NULL)!=0){
    std::cerr<<"Error getting file information from '"<<inputfile<<"'!"<<std::endl;
    return -1;
  }

  switch(file_type){
    case GDT_Byte:
      Master<uint8_t >(inputfile);break;
    case GDT_UInt16:
      Master<uint16_t>(inputfile);break;
    case GDT_Int16:
      Master<int16_t >(inputfile);break;
    case GDT_UInt32:
      Master<uint32_t>(inputfile);break;
    case GDT_Int32:
      Master<int32_t >(inputfile);break;
    case GDT_Float32:
      Master<float   >(inputfile);break;
    case GDT_Float64:
      Master<double  >(inputfile);break;
    default:
      std::cerr<<"Unrecognised data type!"<<std::endl;
      return -1;
  }
}