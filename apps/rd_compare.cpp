#include <iostream>
#include <iomanip>
#include <richdem/common/version.hpp>
#include <richdem/common/Array2D.hpp>
using namespace richdem;

template<typename T, typename U>
int PerformAlgorithm(Array2D<T> a, Array2D<U> b){
  a.loadData();
  b.loadData();
  if(a.geotransform!=b.geotransform)
    std::cout<<"Geotransforms differ."<<std::endl;
  if(a.noData()!=b.noData())
    std::cout<<"NoData values differ."<<std::endl;
  if(a.projection!=b.projection)
    std::cout<<"Projections differ."<<std::endl;
  if(a.width()!=b.width())
    std::cout<<"Widths differ."<<std::endl;
  if(a.height()!=b.height())
    std::cout<<"Heights differ."<<std::endl;

  uint64_t differs = 0;
  for(uint64_t i=0;i<a.size();i++)
    if(a(i)!=b(i))
      differs++;

  std::cout<<"Number of values which differ = "<<differs<<std::endl;

  return 0;
}

#include "router.hpp"

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc, argv);

  if(argc!=3){
    std::cerr<<"Determine whether, and in what ways, two rasters differ. Geotransform, NoData, Projection, Width, Height, and all data values are checked."<<std::endl;
    std::cerr<<argv[0]<<" <Input File A> <Input File B>"<<std::endl;
    return -1;
  }

  PerformAlgorithm(argv[1],argv[2]);

  return 0;
}