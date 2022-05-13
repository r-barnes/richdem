#include "PerlinNoise.h"

#include <richdem/common/Array2D.hpp>
#include <richdem/common/version.hpp>
#include <richdem/terrain_generation.hpp>

#include <stdexcept>

namespace richdem {

void generate_perlin_terrain(Array2D<double>& arr, const uint32_t seed){
  PerlinNoise pn(seed);

  if(arr.width() != arr.height()){
    throw std::runtime_error("Perlin noise array must be square!");
  }

  const auto size = arr.width();

  for(int y=0;y<size;y++)
  for(int x=0;x<size;x++){
    arr(x,y) = pn.noise(10*x/static_cast<double>(size),10*y/static_cast<double>(size),0.8);
  }
}

}