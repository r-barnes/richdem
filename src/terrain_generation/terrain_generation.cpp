#include <richdem/common/Array2D.hpp>
#include <richdem/common/version.hpp>
#include <richdem/terrain_generation.hpp>

#include "PerlinNoise.h"

namespace richdem {

Array2D<double> perlin(const int size, const uint32_t seed){
  PerlinNoise pn(seed);

  Array2D<double> terrain(size,size);
  for(int y=0;y<size;y++)
  for(int x=0;x<size;x++)
    terrain(x,y) = pn.noise(10*x/(double)size,10*y/(double)size,0.8);

  return terrain;
}

}