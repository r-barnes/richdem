#pragma once

#include <richdem/common/Array2D.hpp>

namespace richdem {

void generate_perlin_terrain(Array2D<double>& arr, const uint32_t seed);

}