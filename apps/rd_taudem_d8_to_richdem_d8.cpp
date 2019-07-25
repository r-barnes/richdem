#include <iostream>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <richdem/common/version.hpp>
#include <richdem/common/Array2D.hpp>
using namespace richdem;

template<class T>
int PerformAlgorithm(std::string outname, std::string analysis, Array2D<T> inp){
  const int taudem_to_richdem[9] = {0,5,4,3,2,1,8,7,6};
  inp.loadData();
  for(int y=0;y<inp.height();y++)
  for(int x=0;x<inp.width();x++)
    if(!inp.isNoData(x,y)){
      if(!(0<=inp(x,y) && inp(x,y)<=8))
        throw std::runtime_error("Invalid flow direction '"+std::to_string(inp(x,y))+"' found!");
      inp(x,y) = taudem_to_richdem[(int)inp(x,y)];
    }
  inp.saveGDAL(outname,analysis);
  
  return 0;
}

#include "router.hpp"



int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=3){
    std::cerr<<argv[0]<<" <Input> <Output name>"<<std::endl;
    return -1;
  }

  return PerformAlgorithm(std::string(argv[1]),std::string(argv[2]),analysis);
}
