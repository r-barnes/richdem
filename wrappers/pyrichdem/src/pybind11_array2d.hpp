#ifndef _richdem_pybind11_array2d_
#define _richdem_pybind11_array2d_

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


//Dark magic that wraps an Array2D structure of type T around a Python array
//type (such as NumPy) in such a way that the original data can be modified in-
//place. The following also takes Array2D outputs and converts them into NumPy
//arrays.
namespace pybind11 { namespace detail {
  namespace py = pybind11;
  namespace rd = richdem;

  template <typename T> struct type_caster<rd::Array2D<T>>
  {
    public:

      PYBIND11_TYPE_CASTER(rd::Array2D<T>, _("rd::Array2D<T>"));

      // Conversion part 1 (Python -> C++)
      bool load(py::handle src, bool convert) 
      {
        if (!convert && !py::array_t<T>::check_(src))
          return false;

        auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
          return false;

        auto dims = buf.ndim();
        if (dims != 2 )
          return false;

        value = rd::Array2D<T>((T*) buf.data(),buf.shape()[1],buf.shape()[0]);

        return true;
      }

      //Conversion part 2 (C++ -> Python)
      static py::handle cast(const rd::Array2D<T>& src, py::return_value_policy policy, py::handle parent) 
      {

        std::vector<size_t> shape  (2);
        std::vector<size_t> strides(2);

        shape[0] = src.height();
        shape[1] = src.width();

        strides[0] = src.width();
        strides[1] = strides[0]*sizeof(T);

        py::array a(std::move(shape), std::move(strides), src.data() );

        return a.release();
      }
  };
}}

#endif
