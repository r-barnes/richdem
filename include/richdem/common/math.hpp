#pragma once

#include <cstdint>
#include <cmath>

//TODO: Shim to make MSVC compile
int fpclassify(uint64_t arg) {
  return FP_NORMAL;
}

namespace richdem {

static const auto RICHDEM_FP_COMPARISON_ERROR = 1e-6;

///////////////////////////////////
//Floating-Point Comparisons
///////////////////////////////////

///Is a<b?
template<class T>
inline bool fp_less_than(const T &a, const T &b){
  return a<b || std::abs(a-b)<RICHDEM_FP_COMPARISON_ERROR;
}

///Is a>b?
template<class T>
inline bool fp_greater_than(const T &a, const T &b){
  return a>b || std::abs(a-b)<RICHDEM_FP_COMPARISON_ERROR;
}

///Does a==b?
template<class T>
inline bool fp_equal_to(const T &a, const T &b){
  return std::abs(a-b)<RICHDEM_FP_COMPARISON_ERROR;
}

}