#include <richdem/depressions/Zhou2016pf.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;



//Dark magic that wraps an Array2D structure of type T around a Python array
//type (such as NumPy) in such a way that the original data can be modified in-
//place. The following also takes Array2D outputs and converts them into NumPy
//arrays.
namespace pybind11 { namespace detail {
  template <typename T> struct type_caster<Array2D<T>>
  {
    public:

      PYBIND11_TYPE_CASTER(Array2D<T>, _("Array2D<T>"));

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

        value = Array2D<T>((T*) buf.data(),buf.shape()[1],buf.shape()[0]);

        return true;
      }

      //Conversion part 2 (C++ -> Python)
      static py::handle cast(const Array2D<T>& src, py::return_value_policy policy, py::handle parent) 
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



//Tutorials
//http://www.benjack.io/2017/06/12/python-cpp-tests.html
//https://pybind11.readthedocs.io/en/stable/classes.html
//http://people.duke.edu/~ccc14/sta-663-2016/18G_C++_Python_pybind11.html

//Passing data around
//https://github.com/pybind/pybind11/issues/27

//Note:
//py::array_t<double, py::array::c_style | py::array::forcecast>
//forcecast forces a conversion. We don't use it here in order to ensure that memory is not unnecessarily copied


PYBIND11_MODULE(_richdem, m) {
  m.doc() = "Internal library used by pyRichDEM for calculations";

  m.def("rdFillDepressions",&Zhou2016<float>,"Fill all depressions.");
  m.def("rdFillDepressions",&Zhou2016<double>,"Fill all depressions.");

  // m.def(
  //   "getBoundedScoresForGeoJSON",
  //   &getBoundedScoresForGeoJSON,
  //   "Takes a GeoJSON string as input, calculates all scores, return a JSON dictionary of the scores keyed to values identified 'id' which are expected to be properties of the GeoJSON objects. If id='', then the object's 0-indexed order is used.",
  //   py::arg("gj_subunit"),
  //   py::arg("gj_superunit"),
  //   py::arg("join_on")="",
  //   py::arg("join_id")="",
  //   py::arg("score_list")=""
  // );
}