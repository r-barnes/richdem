#include "pywrapper.hpp"

#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include <string>

namespace py = pybind11;

using namespace richdem;

PYBIND11_MODULE(_richdem, m) {
  m.doc() = "Internal library used by pyRichDEM for calculations";

  //py::bind_vector<std::vector<double>>(m, "VecDouble");
  py::bind_map<std::map<std::string, std::string>>(m, "MapStringString");

  TemplatedFunctionsWrapper<float   >(m, "float"   );
  TemplatedFunctionsWrapper<double  >(m, "double"  );
  TemplatedFunctionsWrapper<int8_t  >(m, "int8_t"  );
  TemplatedFunctionsWrapper<int16_t >(m, "int16_t" );
  TemplatedFunctionsWrapper<int32_t >(m, "int32_t" );
  TemplatedFunctionsWrapper<int64_t >(m, "int64_t" );
  TemplatedFunctionsWrapper<uint8_t >(m, "uint8_t" );
  TemplatedFunctionsWrapper<uint16_t>(m, "uint16_t");
  TemplatedFunctionsWrapper<uint32_t>(m, "uint32_t");
  TemplatedFunctionsWrapper<uint64_t>(m, "uint64_t");

  TemplatedArrayWrapper<float   >(m, "float"   );
  TemplatedArrayWrapper<double  >(m, "double"  );
  TemplatedArrayWrapper<int8_t  >(m, "int8_t"  );
  TemplatedArrayWrapper<int16_t >(m, "int16_t" );
  TemplatedArrayWrapper<int32_t >(m, "int32_t" );
  TemplatedArrayWrapper<int64_t >(m, "int64_t" );
  TemplatedArrayWrapper<uint8_t >(m, "uint8_t" );
  TemplatedArrayWrapper<uint16_t>(m, "uint16_t");
  TemplatedArrayWrapper<uint32_t>(m, "uint32_t");
  TemplatedArrayWrapper<uint64_t>(m, "uint64_t");

  m.def("rdHash",        &rdHash,        "Git hash of previous commit");
  m.def("rdCompileTime", &rdCompileTime, "Commit time of previous commit");

  m.def("FlowAccumulation", &FlowAccumulation<double>, "TODO");

  py::class_<Array3D<float>>(m, "Array3D_float", py::buffer_protocol(), py::dynamic_attr())
      .def(py::init<>())
      .def(py::init<typename Array3D<float>::xy_t, typename Array3D<float>::xy_t,float>())
      
      // .def(py::init<const Array2D<float   >&, T>())
      // .def(py::init<const Array2D<double  >&, T>())
      // .def(py::init<const Array2D<int8_t  >&, T>())
      // .def(py::init<const Array2D<int16_t >&, T>())
      // .def(py::init<const Array2D<int32_t >&, T>())
      // .def(py::init<const Array2D<int64_t >&, T>())
      // .def(py::init<const Array2D<uint8_t >&, T>())
      // .def(py::init<const Array2D<uint16_t>&, T>())
      // .def(py::init<const Array2D<uint32_t>&, T>())
      // .def(py::init<const Array2D<uint64_t>&, T>())

      //NOTE: This does not do reference counting. For that we would want
      //py::object and a wrapped derived class of Array2D
      .def(py::init([](py::handle src){
        // if(!py::array_t<T>::check_(src)) //TODO: What's this about?
          // return false;

        auto buf = py::array_t<float, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
          throw std::runtime_error("Unable to convert array to RichDEM object!");

        //TODO: CHeck stride
        auto dims = buf.ndim();
        if (dims != 3 )
          throw std::runtime_error("Array must have three dimensions!");

        //Array comes to us in (y,x,z) form
        return new Array3D<float>((float*)buf.data(), buf.shape()[1], buf.shape()[0]);
      }))

      .def("size",      &Array3D<float>::size)
      .def("width",     &Array3D<float>::width)
      .def("height",    &Array3D<float>::height)
      .def("empty",     &Array3D<float>::empty)
      .def("noData",    &Array3D<float>::noData)
      
      //TODO: Simplify by casting to double in Python
      .def("setNoData", [](Array3D<float> &a, const float    ndval){ a.setNoData((float)ndval); })
      .def("setNoData", [](Array3D<float> &a, const double   ndval){ a.setNoData((float)ndval); })
      .def("setNoData", [](Array3D<float> &a, const int8_t   ndval){ a.setNoData((float)ndval); })
      .def("setNoData", [](Array3D<float> &a, const int16_t  ndval){ a.setNoData((float)ndval); })
      .def("setNoData", [](Array3D<float> &a, const int32_t  ndval){ a.setNoData((float)ndval); })
      .def("setNoData", [](Array3D<float> &a, const int64_t  ndval){ a.setNoData((float)ndval); })
      .def("setNoData", [](Array3D<float> &a, const uint8_t  ndval){ a.setNoData((float)ndval); })
      .def("setNoData", [](Array3D<float> &a, const uint16_t ndval){ a.setNoData((float)ndval); })
      .def("setNoData", [](Array3D<float> &a, const uint32_t ndval){ a.setNoData((float)ndval); })
      .def("setNoData", [](Array3D<float> &a, const uint64_t ndval){ a.setNoData((float)ndval); })
            
      .def_readwrite("geotransform", &Array3D<float>::geotransform)
      .def_readwrite("projection",   &Array3D<float>::projection)
      .def_readwrite("metadata",     &Array3D<float>::metadata)
      .def("copy", [](const Array3D<float> a){
        return a;
      })

      // .def_buffer([](Array3D<float> &arr) -> py::buffer_info {
      //   return py::buffer_info(
      //     arr.getData(),
      //     sizeof(T),
      //     py::format_descriptor<T>::format(),
      //     2,                                           //Dimensions
      //     {arr.height(), arr.width()},                 //Shape
      //     {sizeof(T) * arr.width(), sizeof(T)} //Stride (in bytes)
      //   );
      // })
      .def("__repr__",
        [=](const Array3D<float> &a) {
            return "<RichDEM 3D array: type=float, width="+std::to_string(a.width())+", height="+std::to_string(a.height())+", owned="+std::to_string(a.owned())+">";
        }
      )
      .def("__call__",
        [](Array3D<float> &a, const int x, const int y, const int n) -> float& {
          return a(x,y,n);
        }
      )
      .def("getIN",
        [](Array3D<float> &a, const int i, const int n) -> float& {
          return a.getIN(i,n);
        }
      );
}
