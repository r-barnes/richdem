#include <iostream>
#include <iomanip>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"
using namespace richdem;

template<class T>
int PerformAlgorithm(std::string outputfile, std::string templatefile, std::string analysis, Array2D<T> raster){
  raster.loadData();

  if(templatefile[0]==':'){
    raster.projection = templatefile.substr(1); //Leave off the colon
  } else {
    Array2D<int8_t> temp(templatefile,false,0,0,0,0,false,false); //Data type doesn't matter since we're not loading it
    raster.projection = temp.projection;
  }

  raster.saveGDAL(outputfile,analysis);

  return 0;
}

#include "router.hpp"



int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);

  if(argc!=4){
    std::cerr<<"Alter a given raster's projection."<<std::endl;
    std::cerr<<"<New Projection> may be either"<<std::endl;
    std::cerr<<"  *A filename, in which case the projection is copied from the file."<<std::endl;
    std::cerr<<"  *A projection string, with a colon preceding it."<<std::endl;
    std::cerr<<argv[0]<<" <Input file> <Output File> <New Projection>"<<std::endl;
    return -1;
  }

  PerformAlgorithm(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),analysis);

  return 0;
}
