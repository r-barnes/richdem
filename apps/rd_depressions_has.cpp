#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/depressions/Barnes2014.hpp"
#include "richdem/common/Array2D.hpp"

using namespace richdem;

template<class T>
int PerformAlgorithm(std::string analysis, Array2D<T> elevation){
  elevation.loadData();

  if(HasDepressions<Topology::D8>(elevation))
    std::cout<<"m Depressions found."<<std::endl;
  else
    std::cout<<"m No depressions found."<<std::endl;

  return 0;
}

#include "router.hpp"

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=2){
    std::cerr<<argv[0]<<" <Input>"<<std::endl;
    return -1;
  }

  return PerformAlgorithm(argv[1],analysis);
}
