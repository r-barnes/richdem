#include "PerlinNoise.h"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/version.hpp"
#include <iostream>

int main(int argc, char **argv){
  std::string analysis = PrintRichdemHeader(argc,argv);

  if(argc!=3){
    std::cerr<<"Syntax: "<<argv[0]<<" <Output Name> <Size>"<<std::endl;
    return -1;
  }

  const int tsize = std::stoi(argv[2]);

  PerlinNoise pn;

  Array2D<float> terrain(tsize,tsize);
  for(int y=0;y<tsize;y++)
  for(int x=0;x<tsize;x++)
    terrain(x,y) = pn.noise(10*x/(double)tsize,10*y/(double)tsize,0.8);

  terrain.saveGDAL(argv[1],analysis);

  return 0;
}