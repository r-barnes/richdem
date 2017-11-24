#include <pybind11/pybind11.h>
#include <richdem/richdem.hpp>
#include <pybind11/numpy.h>
// #include "pybind11_array2d.hpp"
#include <string>

namespace py = pybind11;

using namespace richdem;

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

  m.def("rdFillDepressions",&Zhou2016<float>, "@@depressions/Zhou2016pf.hpp:Zhou2016@@");
  m.def("rdFillDepressions",&Zhou2016<double>,"Fill all depressions.");
  m.def("rdPFepsilon",      &priority_flood_epsilon<float>,"Fill all depressions with epsilon.");

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

  py::class_<Array2D<float>>(m, "Array2Dfloat", py::buffer_protocol(), py::dynamic_attr())
      .def(py::init<>())
      //.def(py::init<const std::string&>())
      .def(py::init<Array2D<float>::xy_t,Array2D<float>::xy_t,float>())
      .def("size",      &Array2D<float>::size)
      .def("width",     &Array2D<float>::width)
      .def("height",    &Array2D<float>::height)
      .def("empty",     &Array2D<float>::empty)
      .def("noData",    &Array2D<float>::noData)
      .def("min",       &Array2D<float>::min)
      .def("max",       &Array2D<float>::max)
      .def("setNoData", &Array2D<float>::setNoData)
      .def_readwrite("projection", &Array2D<float>::projection)
      .def_readwrite("processing_history", &Array2D<float>::processing_history)
      .def("copy", [](const Array2D<float> a){
        return a;
      })
      .def("fromArray", [](Array2D<float> &a, py::handle src){
        // if(!py::array_t<float>::check_(src))
          // return false;

        auto buf = py::array_t<float, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
          throw std::runtime_error("Unable to convert array to RichDEM object!");

        auto dims = buf.ndim();
        if (dims != 2 )
          throw std::runtime_error("Array must have two dimensions!");

        a.clear();
        a.resize(buf.shape()[1], buf.shape()[0]);
        float* dat = (float*)buf.data();
        for(Array2D<float>::i_t i=0;i<a.size();i++)
          a(i) = dat[i];
      })
      .def_buffer([](Array2D<float> &arr) -> py::buffer_info {
        return py::buffer_info(
          arr.getData(),
          sizeof(float),
          py::format_descriptor<float>::format(),
          2,                                           //Dimensions
          {arr.height(), arr.width()},                 //Shape
          {sizeof(float) * arr.width(), sizeof(float)} //Stride (in bytes)
        );
      })
      .def("__repr__",
        [](const Array2D<float> &a) {
            return "<RichDEM array: type=float, width="+std::to_string(a.width())+", height="+std::to_string(a.height())+">";
        }
      );
}
 