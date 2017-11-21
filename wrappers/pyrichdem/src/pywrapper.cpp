#include <richdem/depressions/Zhou2016pf.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

//Tutorials
//http://www.benjack.io/2017/06/12/python-cpp-tests.html
//https://pybind11.readthedocs.io/en/stable/classes.html
//http://people.duke.edu/~ccc14/sta-663-2016/18G_C++_Python_pybind11.html

//Passing data around
//https://github.com/pybind/pybind11/issues/27

//Note:
//py::array_t<double, py::array::c_style | py::array::forcecast>
//forcecast forces a conversion. We don't use it here in order to ensure that memory is not unnecessarily copied



/* Bind MatrixXd (or some other Eigen type) to Python */
// typedef Eigen::MatrixXd Matrix;

// typedef Matrix::Scalar Scalar;
// constexpr bool rowMajor = Matrix::Flags & Eigen::RowMajorBit;

// py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
//     .def("__init__", [](Matrix &m, py::buffer b) {
//         typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

//         /* Request a buffer descriptor from Python */
//         py::buffer_info info = b.request();

//         /* Some sanity checks ... */
//         if (info.format != py::format_descriptor<Scalar>::format())
//             throw std::runtime_error("Incompatible format: expected a double array!");

//         if (info.ndim != 2)
//             throw std::runtime_error("Incompatible buffer dimension!");

//         auto strides = Strides(
//             info.strides[rowMajor ? 0 : 1] / (py::ssize_t)sizeof(Scalar),
//             info.strides[rowMajor ? 1 : 0] / (py::ssize_t)sizeof(Scalar));

//         auto map = Eigen::Map<Matrix, 0, Strides>(
//             static_cast<Scalar *>(info.ptr), info.shape[0], info.shape[1], strides);

//         new (&m) Matrix(map);
//     });


template<class T>
void rdFillDepressions(py::array_t<T, py::array::c_style> elev){
  auto elevbuf = elev.request();
  Array2D<T> a_elev((T*) elevbuf.ptr, elevbuf.shape[1], elevbuf.shape[0]);
  Zhou2016(a_elev);
}




PYBIND11_PLUGIN(_richdem) {
  py::module m("_richdem", "Internal library used by pyRichDEM for calculations");

  m.def(
    "rdFillDepressions",
    &rdFillDepressions<double>,
    "Fill all depressions."
  );

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

  return m.ptr();
}