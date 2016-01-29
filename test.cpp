#include "Array2D.hpp"
#include <iostream>
#include <cassert>

//Are all serials unique?
//    ./test.exe | awk '{print $1}' | sort -n | uniq -c | sed 's/^\s*//' | sort -n -k 1
//Do all serials have a good bidirectional mapping?
//    ./test.exe | sort -n -k 1 | awk '{print $5}' | sort -r

template<class T>
int xyToSerial(const int x, const int y, const Array2D<T> &grid){
  if(y==0)                                         //Top row
    return x;

  if(x==grid.viewWidth()-1)                        //Right hand side
    return grid.viewWidth()+y;

  if(y==grid.viewHeight()-1)      
    return grid.viewWidth()+grid.viewHeight()+x;   //Bottom-row

  return 2*grid.viewWidth()+grid.viewHeight()+y;   //Left-hand side
}

template<class T>
void serialToXY(const int serial, int &x, int &y, const Array2D<T> &grid){
  if(serial<grid.viewWidth()){                            //Top row
    x = serial;
    y = 0;
  } else if(serial<grid.viewWidth()+grid.viewHeight()){   //Right-hand side
    x = grid.viewWidth()-1;
    y = serial-grid.viewWidth();
  } else if(serial<2*grid.viewWidth()+grid.viewHeight()){ //Bottom row
    x = serial-grid.viewWidth()-grid.viewHeight();
    y = grid.viewHeight()-1;
  } else {                                                //Left-hand side
    x = 0;
    y = serial-2*grid.viewWidth()-grid.viewHeight(); 
  }
}

template<class T>
void testCombo(const int x, const int y, const Array2D<T> &grid){
  int xr, yr;
  int serial = xyToSerial(x,y,grid);
  std::cout<<serial<<" "<<x<<","<<y<<" = ";
  serialToXY(serial,xr,yr,grid);
  std::cout<<xr<<","<<yr<<" ";
  std::cout<<(x==xr && y==yr)<<std::endl;
}

template<class T>
void GridPerimToArray(const Array2D<T> &grid, std::vector<T> &vec){
  assert(vec.size()==0); //Ensure receiving array is empty

  std::vector<T> vec2copy;

  vec2copy = grid.getRowData(0);                         //Top
  vec.insert(vec.end(),vec2copy.begin(),vec2copy.end());

  vec2copy = grid.getColData(grid.viewWidth()-1);        //Right
  vec.insert(vec.end(),vec2copy.begin(),vec2copy.end());
  
  vec2copy = grid.getRowData(grid.viewHeight()-1);       //Bottom
  vec.insert(vec.end(),vec2copy.begin(),vec2copy.end());
  
  vec2copy = grid.getColData(0);                         //Left
  vec.insert(vec.end(),vec2copy.begin(),vec2copy.end());
}

template<class T>
void testXYserialization(const Array2D<T> &grid){
  for(int x=0;x<grid.viewWidth();x++){
    testCombo(x,0,grid);
    testCombo(x,grid.viewHeight()-1,grid);
  }

  for(int y=1;y<grid.viewHeight()-1;y++){
    testCombo(0,y,grid);
    testCombo(grid.viewWidth()-1,y,grid);
  }
}

int main(){
  Array2D<char> dem(3617,3301); //Prime numbers
  //testXYserialization(dem);

  std::vector<int> link_arr;
  for(int x=0;x<dem.viewWidth();x++)
    link_arr.push_back(xyToSerial(x,0,dem));
  for(int y=0;y<dem.viewHeight();y++)
    link_arr.push_back(xyToSerial(dem.viewWidth()-1,y,dem));
  for(int x=0;x<dem.viewWidth();x++)
    link_arr.push_back(xyToSerial(x,dem.viewHeight()-1,dem));
  for(int y=0;y<dem.viewHeight();y++)
    link_arr.push_back(xyToSerial(0,y,dem));

  for(int i=0;i<link_arr.size();i++){
    int xr,yr;
    serialToXY(link_arr[i],xr,yr,dem);
    std::cout<<link_arr[i]
             <<" "<<xr<<","<<yr<<" "
             <<xyToSerial(xr,yr,dem)<<" "
             <<((link_arr[i]==xyToSerial(xr,yr,dem))?"MATCHES":"NOPE")<<std::endl;
  }
}