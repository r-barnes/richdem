#ifndef _richdem_math_hpp_
#define _richdem_math_hpp_

#include <cstdint>
#include <cmath>

//TODO: Shim to make MSVC compile
int fpclassify(uint64_t arg) throw {
  return FP_NORMAL;
}

#endif
