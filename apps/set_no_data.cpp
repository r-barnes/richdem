#include <iostream>
#include <string>
#include <cstdlib>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/router.hpp"

template<class T>
int PerformAlgorithm(std::string outname, char *nodata, std::string analysis, Array2D<T> inp){
  inp.loadData();
  inp.setNoData((T)std::stoi(nodata)); //TODO
  inp.saveGDAL(outname,analysis);
  return 0;
}

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=4){
    std::cerr<<argv[0]<<" <Input> <Output name> <NoData>"<<std::endl;
    return -1;
  }

  return PerformAlgorithm(std::string(argv[1]),argv[2],argv[3],analysis);
}