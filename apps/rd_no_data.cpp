#include <iostream>
#include <string>
#include <cstdlib>
#include <richdem/common/version.hpp>
#include <richdem/common/Array2D.hpp>
using namespace richdem;

template<class T>
int PerformAlgorithm(std::string outname, char *nodata, std::string analysis, Array2D<T> inp){
  inp.loadData();
  inp.setNoData((T)std::stoi(nodata)); //TODO
  inp.saveGDAL(outname,analysis);
  return 0;
}

#include "router.hpp"



int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=2 && argc!=4){
    std::cerr<<"Get or set a raster's NoData value."<<std::endl;
    std::cerr<<argv[0]<<" <Input>"<<std::endl;
    std::cerr<<argv[0]<<" <Input> <Output name> <NoData>"<<std::endl;
    return -1;
  }

  if(argc==2){
    Array2D<int8_t> temp(argv[1],false,0,0,0,0,false,false); //Data type doesn't matter since we're not loading it
    std::cerr<<temp.noData()<<std::endl;
  } else {
    PerformAlgorithm(std::string(argv[1]),argv[2],argv[3],analysis);
  }

  return 0;
}
