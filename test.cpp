#include "Array2D.hpp"
#include <iostream>

//Are all serials unique?
//    ./test.exe | awk '{print $1}' | sort -n | uniq -c | sed 's/^\s*//' | sort -n -k 1

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
  testXYserialization(dem);
}