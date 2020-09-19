#pragma once

#include <richdem/common/Array2D.hpp>

namespace richdem {

Array2D<double> perlin(const int size, const uint32_t seed);

}