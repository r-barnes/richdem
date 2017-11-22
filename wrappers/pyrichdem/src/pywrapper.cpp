#include <pybind11/pybind11.h>
#include <richdem/richdem.hpp>
#include "pybind11_array2d.hpp"

namespace py = pybind11;
namespace rd = richdem;

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

  m.def("rdFillDepressions",&rd::Zhou2016<float>,"Fill all depressions.");
  m.def("rdFillDepressions",&rd::Zhou2016<double>,"Fill all depressions.");

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
