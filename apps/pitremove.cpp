#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/common/router.hpp"
#include "richdem/depressions/priority_flood.hpp"
#include "richdem/depressions/Zhou2016pf.hpp"
#include "richdem/common/Array2D.hpp"
//#include "richdem/flats/flat_resolution.hpp"
//#include "richdem/methods/d8_methods.hpp"

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
    std::cerr<<argv[0]<<" <Input> <Output name>"<<std::endl;
    return -1;
  }

  return PerformAlgorithm(argv[1],argv[2],analysis);
}