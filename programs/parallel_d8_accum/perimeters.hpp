#pragma once

#include <richdem/common/Array2D.hpp>

#include <cassert>

//TODO: Explain
int xyToSerial(const int x, const int y, const int width, const int height){
  //Ensure cell is on the perimeter
  assert( (x==0 || x==width-1 || y==0 || y==height-1) && x>=0 && y>=0 && x<width && y<height);

  if(y==0)                         //Top row
    return x;

  if(x==width-1)                   //Right hand side
    return (width-1)+y;

  if(y==height-1)                  //Bottom-row
    return (width-1)+(height)+x;

  return 2*(width-1)+(height-1)+y; //Left-hand side
}

//TODO: Explain
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
  assert( (x==0 || x==width-1 || y==0 || y==height-1) && x>=0 && y>=0 && x<width && y<height);
}

template<class U>
void GridPerimToArray(const richdem::Array2D<U> &grid, std::vector<U> &vec){
  assert(vec.size()==0); //Ensure receiving array is empty

  std::vector<U> vec2copy;

  vec2copy = grid.getRowData(0);                         //Top
  vec.insert(vec.end(),vec2copy.begin(),vec2copy.end());

  vec2copy = grid.getColData(grid.width()-1);        //Right
  vec.insert(vec.end(),vec2copy.begin()+1,vec2copy.end());

  vec2copy = grid.getRowData(grid.height()-1);       //Bottom
  vec.insert(vec.end(),vec2copy.begin(),vec2copy.end()-1);

  vec2copy = grid.getColData(0);                         //Left
  vec.insert(vec.end(),vec2copy.begin()+1,vec2copy.end()-1);
}