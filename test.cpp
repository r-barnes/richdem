#include "Array2D.hpp"
#include <iostream>
#include <cassert>
#include <iomanip>

//g++ -O3 -o test.exe `gdal-config --cflags` `gdal-config --libs` test.cpp -lgdal --std=c++11

//Are all serials unique?
//    ./test.exe | awk '{print $1}' | sort -n | uniq -c | sed 's/^\s*//' | sort -n -k 1
//Do all serials have a good bidirectional mapping?
//    ./test.exe | sort -n -k 1 | awk '{print $5}' | sort -r

int xyToSerial(const int x, const int y, const int width, const int height){
  //Ensure cell is on the perimeter
  assert( (x==0 || x==width-1 || y==0 || y==height-1) && x>=0 && y>=0 && x<=width-1 && y<=height-1);

  if(y==0)                         //Top row
    return x;

  if(x==width-1)                   //Right hand side
    return (width-1)+y;

  if(y==height-1)                  //Bottom-row
    return (width-1)+(height)+x;   

  return 2*(width-1)+(height-1)+y; //Left-hand side
}

void serialToXY(const int serial, int &x, int &y, const int width, const int height){
  if(serial<width){                        //Top row
    x = serial;
    y = 0;
  } else if(serial<(width-1)+height){     //Right-hand side
    x = width-1;
    y = serial-(width-1);
  } else if(serial<2*(width-1)+(height)){ //Bottom row
    x = serial-(width-1)-(height-1)-1;
    y = height-1;
  } else {                                //Left-hand side
    x = 0;
    y = serial-2*(width-1)-(height-1);
  }

  //Ensure cell is on the perimeter
  assert( (x==0 || x==width-1 || y==0 || y==height-1) && x>=0 && y>=0 && x<=width-1 && y<=height-1);
}

template<class T>
void GridPerimToArray(const Array2D<T> &grid, std::vector<T> &vec){
  assert(vec.size()==0); //Ensure receiving array is empty

  std::vector<T> vec2copy;

  vec2copy = grid.getRowData(0);                         //Top
  vec.insert(vec.end(),vec2copy.begin(),vec2copy.end());

  vec2copy = grid.getColData(grid.viewWidth()-1);        //Right
  vec.insert(vec.end(),vec2copy.begin()+1,vec2copy.end());
  
  vec2copy = grid.getRowData(grid.viewHeight()-1);       //Bottom
  vec.insert(vec.end(),vec2copy.begin(),vec2copy.end()-1);
  
  vec2copy = grid.getColData(0);                         //Left
  vec.insert(vec.end(),vec2copy.begin()+1,vec2copy.end()-1);
}


int main(){
  Array2D<int> dem(13,17);
  int val=0;
  for(int y=0;y<dem.viewHeight();y++)
  for(int x=0;x<dem.viewWidth();x++)
    dem(x,y)=val;

  for(int y=0;y<dem.viewHeight();y++){
    for(int x=0;x<dem.viewWidth();x++){
      dem(x,y) = -1;
      if(dem.edge_grid(x,y)){
        dem(x,y) = xyToSerial(x,y,dem.viewWidth(),dem.viewHeight());
        std::cout<<std::setw(3)<<dem(x,y);
      } else
        std::cout<<std::setw(3)<<" ";
    }
    std::cout<<std::endl;
  }

  std::vector<int> edges;
  GridPerimToArray(dem,edges);
  for(int s=0;s<edges.size();s++)
    std::cout<<std::setw(3)<<s;
  std::cout<<std::endl;
  for(const auto &e: edges)
    std::cout<<std::setw(3)<<e;
  std::cout<<std::endl;

  for(int s=0;s<edges.size();s++)
    if(s!=edges[s])
      std::cout<<"GridPerim to Order issue s="<<s<<" does not map to edge "<<edges[s]<<std::endl;

  std::cerr<<"Testing that serialToXY can get back to the xyToSerial input."<<std::endl;
  for(int y=0;y<dem.viewHeight();y++)
  for(int x=0;x<dem.viewWidth();x++){
    if(!dem.edge_grid(x,y))
      continue;
    int xr,yr;
    int s_good = xyToSerial(x,y,dem.viewWidth(),dem.viewHeight());
    serialToXY(s_good,xr,yr,dem.viewWidth(),dem.viewHeight());
    if(xr!=x || yr!=y)
      std::cout<<"("<<x<<","<<y<<") with serial "<<s_good<<" sadly maps to ("<<xr<<","<<yr<<")"<<std::endl;
  }

  std::cerr<<"Testing that xyToSerial can get back to the serialToXY input."<<std::endl;
  for(int s_good=0;s_good<edges.size();s_good++){
    int xr,yr;
    serialToXY(s_good,xr,yr,dem.viewWidth(),dem.viewHeight());
    int s_to_test = xyToSerial(xr,yr,dem.viewWidth(),dem.viewHeight());
    if(s_to_test!=s_good)
      std::cout<<"Stest failed."<<std::endl;
  }

  std::cout<<"Unless there is sadness above, all inverse mappings checked out."<<std::endl;

}