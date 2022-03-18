#include "perimeters.hpp"

#include <richdem/common/Array2D.hpp>

#include <cassert>
#include <iomanip>
#include <iostream>

using namespace richdem;

//g++ -O3 -o test.exe `gdal-config --cflags` `gdal-config --libs` test.cpp -lgdal --std=c++17

//Are all serials unique?
//    ./test.exe | awk '{print $1}' | sort -n | uniq -c | sed 's/^\s*//' | sort -n -k 1
//Do all serials have a good bidirectional mapping?
//    ./test.exe | sort -n -k 1 | awk '{print $5}' | sort -r


int main(){
  Array2D<int> dem(5,7);
  int val=0;
  for(int y=0;y<dem.height();y++)
  for(int x=0;x<dem.width();x++)
    dem(x,y)=val;

  std::cerr<<"This is what happens when the x,y coordinate at the indicated spot is converted to serial. These should match expecations:"<<std::endl;
  for(int y=0;y<dem.height();y++){
    for(int x=0;x<dem.width();x++){
      dem(x,y) = -1;
      if(dem.isEdgeCell(x,y)){
        dem(x,y) = xyToSerial(x,y,dem.width(),dem.height());
        std::cout<<std::setw(3)<<dem(x,y);
      } else
        std::cout<<std::setw(3)<<" ";
    }
    std::cout<<std::endl;
  }

  std::cout<<"\nThese numbers should be the same on both rows:"<<std::endl;
  std::vector<int> edges;
  GridPerimToArray(dem,edges);
  for(size_t s=0;s<edges.size();s++)
    std::cout<<std::setw(3)<<s;
  std::cout<<std::endl;
  for(const auto &e: edges)
    std::cout<<std::setw(3)<<e;
  std::cout<<std::endl;

  for(int s=0;s<(int)edges.size();s++)
    if(s!=edges[s])
      std::cout<<"GridPerim to Order issue s="<<s<<" does not map to edge "<<edges[s]<<std::endl;

  std::cerr<<"\nTesting that serialToXY can get back to the xyToSerial input."<<std::endl;
  for(int y=0;y<dem.height();y++)
  for(int x=0;x<dem.width();x++){
    if(!dem.isEdgeCell(x,y))
      continue;
    int xr,yr;
    int s_good = xyToSerial(x,y,dem.width(),dem.height());
    serialToXY(s_good,xr,yr,dem.width(),dem.height());
    if(xr!=x || yr!=y)
      std::cout<<"("<<x<<","<<y<<") with serial "<<s_good<<" sadly maps to ("<<xr<<","<<yr<<")"<<std::endl;
  }

  std::cerr<<"Testing that xyToSerial can get back to the serialToXY input."<<std::endl;
  for(int s_good=0;s_good<(int)edges.size();s_good++){
    int xr,yr;
    serialToXY(s_good,xr,yr,dem.width(),dem.height());
    int s_to_test = xyToSerial(xr,yr,dem.width(),dem.height());
    if(s_to_test!=s_good)
      std::cout<<"Stest failed."<<std::endl;
  }

  std::cout<<"Unless there is sadness above, all inverse mappings checked out."<<std::endl;

}