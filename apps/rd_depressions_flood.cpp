#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/common/router.hpp"
#include "richdem/depressions/Zhou2016pf.hpp"
#include "richdem/depressions/priority_flood.hpp"
#include "richdem/common/Array2D.hpp"

template<class T>
int PerformAlgorithm(std::string outputname, uint32_t max_dep_size, std::string analysis, Array2D<T> elevation){
  elevation.loadData();

  if(max_dep_size==0)
    Zhou2016(elevation);
  else
    improved_priority_flood_max_dep(elevation,max_dep_size);

  elevation.saveGDAL(outputname,analysis);

  return 0;
}

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=4){
    std::cerr<<"Eliminate all depressions via flooding."<<std::endl;
    std::cerr<<argv[0]<<" <Input> <Output name> <Maximum Depression Size>"<<std::endl;
    std::cerr<<"\t<Maximum Depression Size> - Depressions larger than this are not flooded."<<std::endl;
    std::cerr<<"                              Use `0` to flood all depressions.            "<<std::endl;
    return -1;
  }

  uint32_t max_dep_size = std::stoul(argv[3]);

  return PerformAlgorithm(argv[1],argv[2],max_dep_size,analysis);
}