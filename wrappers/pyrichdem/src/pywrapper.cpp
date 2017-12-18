#include <pybind11/pybind11.h>
#include <richdem/richdem.hpp>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>
#include <richdem/methods/dall_methods.hpp>
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

// #include "pybind11_array2d.hpp"

template<class T>
class pyArray2D : public Array2D<T> {
 private:
  py::object myobj;

 public:
  pyArray2D() = default;

  pyArray2D(typename Array2D<T>::xy_t width, typename Array2D<T>::xy_t height, T val) : Array2D<T>(width,height,val){}

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

  pyArray2D(py::object obj) : myobj(obj) {
    // if(!py::array_t<T>::check_(src)) //TODO: What's this about?
      // return false;

    auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(obj);
    if (!buf)
      throw std::runtime_error("Unable to convert array to RichDEM object!");

    //TODO: CHeck stride
    auto dims = buf.ndim();
    if (dims != 2 )
      throw std::runtime_error("Array must have two dimensions!");

    Array2D<T>::data        = ManagedVector<T>((T*)buf.data(), buf.shape()[1]*buf.shape()[0]);
    Array2D<T>::view_width  = buf.shape()[1];
    Array2D<T>::view_height = buf.shape()[0];    
  }


};


template<class T>
void TemplatedWrapper(py::module &m, std::string tname){
  m.def("rdFillDepressions",     &Zhou2016<T>,                     "@@depressions/Zhou2016pf.hpp:Zhou2016@@"); //TODO
  m.def("rdPFepsilon",           &priority_flood_epsilon<T>,       "Fill all depressions with epsilon."); //TODO

  m.def("rdBreach",              &Lindsay2016<T>,                  "TODO");
  
  //m.def("rdBreach",              [](Array2D<T> &dem, const int mode, bool fill_depressions){&Lindsay2016<T>(dem,mode,fill_depressions);}, "TODO");

  m.def("TA_SPI",                &TA_SPI<T, float, double>,       "TODO");         
  m.def("TA_CTI",                &TA_CTI<T, float, double>,       "TODO");         
  m.def("TA_slope_riserun",      &TA_slope_riserun     <T>,       "TODO");                   
  m.def("TA_slope_percentage",   &TA_slope_percentage  <T>,       "TODO");                      
  m.def("TA_slope_degrees",      &TA_slope_degrees     <T>,       "TODO");                   
  m.def("TA_slope_radians",      &TA_slope_radians     <T>,       "TODO");                   
  m.def("TA_aspect",             &TA_aspect            <T>,       "TODO");            
  m.def("TA_curvature",          &TA_curvature         <T>,       "TODO");               
  m.def("TA_planform_curvature", &TA_planform_curvature<T>,       "TODO");                        
  m.def("TA_profile_curvature",  &TA_profile_curvature <T>,       "TODO");                       

  m.def("FA_Tarboton",           &FA_Tarboton         <T,double>, "TODO");
  m.def("FA_Holmgren",           &FA_Holmgren         <T,double>, "TODO");
  m.def("FA_Quinn",              &FA_Quinn            <T,double>, "TODO");
  m.def("FA_Freeman",            &FA_Freeman          <T,double>, "TODO");
  m.def("FA_FairfieldLeymarie",  &FA_FairfieldLeymarie<T,double>, "TODO");
  m.def("FA_Rho8",               &FA_Rho8             <T,double>, "TODO");
  m.def("FA_D8",                 &FA_D8               <T,double>, "TODO");
  m.def("FA_OCallaghan",         &FA_OCallaghan       <T,double>, "TODO");

  //TODO
  // m.def(
  //   "WrapNumPy",
  //   [](py::handle src){
  //     if (!py::array_t<T>::check_(src))
  //       throw std::runtime_error("WrapNumPy: Array check failed!");

  //     auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
  //     if (!buf)
  //       throw std::runtime_error("WrapNumPy: Could not generate a buffer!")

  //     auto dims = buf.ndim();
  //     if (dims != 2 )
  //       throw std::runtime_error("WrapNumPy: NumPy array must have 2 dimensions.")

  //     value = rd::Array2D<T>((T*) buf.data(),buf.shape()[1],buf.shape()[0]);

  //     return true;
  //   },
  //   "TODO"
  // );

  py::class_<pyArray2D<T>>(m, ("Array2D_" + tname).c_str(), py::buffer_protocol(), py::dynamic_attr())
      .def(py::init<>())
      .def(py::init<typename pyArray2D<T>::xy_t, typename pyArray2D<T>::xy_t,T>())
      
      // .def(py::init<const pyArray2D<float   >&, T>())
      // .def(py::init<const pyArray2D<double  >&, T>())
      // .def(py::init<const pyArray2D<int8_t  >&, T>())
      // .def(py::init<const pyArray2D<int16_t >&, T>())
      // .def(py::init<const pyArray2D<int32_t >&, T>())
      // .def(py::init<const pyArray2D<int64_t >&, T>())
      // .def(py::init<const pyArray2D<uint8_t >&, T>())
      // .def(py::init<const pyArray2D<uint16_t>&, T>())
      // .def(py::init<const pyArray2D<uint32_t>&, T>())
      // .def(py::init<const pyArray2D<uint64_t>&, T>())

      .def(py::init<py::object>())


      // //TODO: Turn this into an appropriate wrapping function so we can remove WrapNumPy in Python?
      // .def(py::init([](py::handle src){
      //   // if(!py::array_t<T>::check_(src)) //TODO: What's this about?
      //     // return false;

      //   auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
      //   if (!buf)
      //     throw std::runtime_error("Unable to convert array to RichDEM object!");

      //   //TODO: CHeck stride
      //   auto dims = buf.ndim();
      //   if (dims != 2 )
      //     throw std::runtime_error("Array must have two dimensions!");

      //   return new pyArray2D<T>((T*)buf.data(), buf.shape()[1], buf.shape()[0]);
      // }))

      .def("size",      &pyArray2D<T>::size)
      .def("width",     &pyArray2D<T>::width)
      .def("height",    &pyArray2D<T>::height)
      .def("empty",     &pyArray2D<T>::empty)
      .def("noData",    &pyArray2D<T>::noData)
      .def("min",       &pyArray2D<T>::min)
      .def("max",       &pyArray2D<T>::max)
      
      //TODO: Simplify by casting to double in Python
      .def("setNoData", [](pyArray2D<T> &a, const float    ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](pyArray2D<T> &a, const double   ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](pyArray2D<T> &a, const int8_t   ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](pyArray2D<T> &a, const int16_t  ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](pyArray2D<T> &a, const int32_t  ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](pyArray2D<T> &a, const int64_t  ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](pyArray2D<T> &a, const uint8_t  ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](pyArray2D<T> &a, const uint16_t ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](pyArray2D<T> &a, const uint32_t ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](pyArray2D<T> &a, const uint64_t ndval){ a.setNoData((T)ndval); })
            
      .def_readwrite("geotransform", &pyArray2D<T>::geotransform)
      .def_readwrite("projection",   &pyArray2D<T>::projection)
      .def_readwrite("metadata",     &pyArray2D<T>::metadata)
      .def("copy", [](const pyArray2D<T> a){
        return a;
      })

      .def("copyFromArray", [](pyArray2D<T> &a, py::handle src){
        // if(!py::array_t<T>::check_(src)) //TODO: What's this about?
          // return false;

        auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
          throw std::runtime_error("Unable to convert array to RichDEM object!");

        //TODO: CHeck stride
        auto dims = buf.ndim();
        if (dims != 2 )
          throw std::runtime_error("Array must have two dimensions!");

        a.clear();
        a.resize(buf.shape()[1], buf.shape()[0]);
        T* dat = (T*)buf.data();
        for(typename pyArray2D<T>::i_t i=0;i<a.size();i++)
          a(i) = dat[i];
      })
      .def_buffer([](pyArray2D<T> &arr) -> py::buffer_info {
        return py::buffer_info(
          arr.getData(),
          sizeof(T),
          py::format_descriptor<T>::format(),
          2,                                           //Dimensions
          {arr.height(), arr.width()},                 //Shape
          {sizeof(T) * arr.width(), sizeof(T)} //Stride (in bytes)
        );
      })
      .def("__repr__",
        [](const pyArray2D<T> &a) {
            return "<RichDEM array: type=T, width="+std::to_string(a.width())+", height="+std::to_string(a.height())+", owned="+std::to_string(a.owned())+">";
        }
      )
      .def("__call__",
        [](pyArray2D<T> &a, const int x, const int y) -> T& {
          return a(x,y);
        }
      )
      .def("__call__",
        [](pyArray2D<T> &a, const int i) -> T& {
          return a(i);
        }
      );      
}



PYBIND11_MODULE(_richdem, m) {
  m.doc() = "Internal library used by pyRichDEM for calculations";

  py::bind_vector<std::vector<double>>(m, "VecDouble");
  py::bind_map<std::map<std::string, std::string>>(m, "MapStringString");

  TemplatedWrapper<float   >(m, "float"   );
  TemplatedWrapper<double  >(m, "double"  );
  TemplatedWrapper<int8_t  >(m, "int8_t"  );
  TemplatedWrapper<int16_t >(m, "int16_t" );
  TemplatedWrapper<int32_t >(m, "int32_t" );
  TemplatedWrapper<int64_t >(m, "int64_t" );
  TemplatedWrapper<uint8_t >(m, "uint8_t" );
  TemplatedWrapper<uint16_t>(m, "uint16_t");
  TemplatedWrapper<uint32_t>(m, "uint32_t");
  TemplatedWrapper<uint64_t>(m, "uint64_t");

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
 