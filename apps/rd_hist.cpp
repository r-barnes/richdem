#include <iostream>
#include <iomanip>
#include <unordered_map>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"
using namespace richdem;

template<class T>
int PerformAlgorithm(Array2D<T> rast){
  std::unordered_map<T,int> counts;

  rast.loadData();

  for(unsigned int i=0;i<rast.size();i++)
    counts[rast(i)]++;

  std::cout<<"Nodata: "<<(int)rast.noData()<<std::endl;

  for(const auto &x: counts)
    std::cout<<std::setw(20)<<x.first<<" "<<std::setw(20)<<x.second<<"\n";

  return 0;
}

#include "router.hpp"



int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=2){
    std::cerr<<argv[0]<<" <Flowdirs input file>"<<std::endl;
    return -1;
  }

  PerformAlgorithm(std::string(argv[1]));

  return 0;
}