#include <iostream>
#include <string>
#include <cstdlib>
#include <richdem/common/version.hpp>
#include <richdem/common/Array2D.hpp>
using namespace richdem;

template<class T>
int PerformAlgorithm(std::string outname, int new_width, int new_height, std::string analysis, Array2D<T> inp){
  inp.loadData();

  if(new_width<inp.width()){
    std::cerr<<"Desired width is smaller than DEM's current width!"<<std::endl;
    return -1;
  }

  if(new_height<inp.height()){
    std::cerr<<"Desired height is smaller than DEM's current height!"<<std::endl;
    return -1;
  }

  inp.expand(new_width,new_height,inp.noData());

  inp.saveGDAL(outname,analysis);
  return 0;
}

#include "router.hpp"



int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);
  
  if(argc!=5){
    std::cerr<<argv[0]<<" <Input> <Output name> <Width> <Height>"<<std::endl;
    return -1;
  }

  int new_width  = std::stoi(argv[3]);
  int new_height = std::stoi(argv[4]);

  PerformAlgorithm(argv[1],argv[2],new_width,new_height,analysis);

  return 0;
}