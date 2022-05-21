#pragma once

#include <richdem/common/Array2D.hpp>
#include <functional>

namespace richdem {
  template<class T>
  void iterate_flat(Array2D<T> &arr, std::function<void(Array2D<T>&, typename Array2D<T>::i_t)> func){
    for(auto i=arr.i0();i<arr.size();i++){
      func(arr, i);
    }
  }

  template<class T>
  void iterate_2d(Array2D<T> &arr, std::function<void(typename Array2D<T>::xy_t, typename Array2D<T>::xy_t)> func){
    for(int y=0;y<arr.height();y++)
    for(int x=0;x<arr.width();x++)
      func(x, y);
  }
}