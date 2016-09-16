#include <iostream>
#include <vector>
#include <string>
#include "richdem/common/version.hpp"
#include "richdem/common/Array2D.hpp"

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);
  
  if(argc!=2){
    std::cerr<<argv[0]<<" <Input>"<<std::endl;
    return -1;
  }

  int width;
  int height;
  GDALDataType dtype;
  std::vector<double> geotransform(6);

  getGDALDimensions(argv[1],height,width,dtype,geotransform.data());

  std::cout<<"Geotransform for '"<<argv[1]<<"': ";
  for(const auto x: geotransform)
    std::cout<<x<<" ";
  std::cout<<width<<" "<<height<<std::endl;

  return 0;
}