#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/common/router.hpp"
#include "richdem/depressions/Zhou2016pf.hpp"
#include "richdem/common/Array2D.hpp"

template<class T>
int PerformAlgorithm(std::string outputname, std::string analysis, Array2D<T> elevation){
  elevation.loadData();

  Zhou2016(elevation);  

  elevation.saveGDAL(outputname,analysis);

  return 0;
}

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=3){
    std::cerr<<"Eliminate all depressions via flooding."<<std::endl;
    std::cerr<<argv[0]<<" <Input> <Output name>"<<std::endl;
    return -1;
  }

  return PerformAlgorithm(argv[1],argv[2],analysis);
}