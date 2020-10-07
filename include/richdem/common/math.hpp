#pragma once

#include <cstdint>
#include <cmath>

//TODO: Shim to make MSVC compile
#ifdef _MSC_VER
int fpclassify(uint64_t arg) {
  return FP_NORMAL;
}
#endif

namespace richdem {

static const auto RICHDEM_FP_COMPARISON_ERROR = 1e-6;

///////////////////////////////////
//Floating-Point Comparisons
///////////////////////////////////

///Is a<=b?
inline bool fp_le(const double &a, const double &b){
  return a<b || std::abs(a-b)<RICHDEM_FP_COMPARISON_ERROR;
}

///Is a>=b?
inline bool fp_ge(const double &a, const double &b){
  return a>b || std::abs(a-b)<RICHDEM_FP_COMPARISON_ERROR;
}

///Does a==b?
inline bool fp_eq(const double &a, const double &b){
  return std::abs(a-b)<RICHDEM_FP_COMPARISON_ERROR;
}

}