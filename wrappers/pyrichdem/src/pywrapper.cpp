#include <pybind11/pybind11.h>
#include <richdem/depressions/depressions.hpp>
#include <richdem/methods/terrain_attributes.hpp>
#include <richdem/methods/flow_accumulation.hpp>
#include <richdem/flats/flats.hpp>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
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
void TemplatedWrapper(py::module &m, std::string tname){
  m.def("rdFillDepressionsD8",   &PriorityFlood_Zhou2016<T>,                "@@depressions/Zhou2016pf.hpp:Zhou2016@@"); //TODO
  m.def("rdFillDepressionsD4",   &PriorityFlood_Barnes2014<Topology::D4,T>, "@@depressions/Zhou2016pf.hpp:Zhou2016@@"); //TODO
  m.def("rdPFepsilonD8",         &PriorityFloodEpsilon_Barnes2014<Topology::D8,T>, "Fill all depressions with epsilon."); //TODO
  m.def("rdPFepsilonD4",         &PriorityFloodEpsilon_Barnes2014<Topology::D4,T>, "Fill all depressions with epsilon."); //TODO

  m.def("rdResolveFlatsEpsilon", &ResolveFlatsEpsilon<T>,         "TODO");

  m.def("rdBreachDepressionsD8",   &BreachDepressions<Topology::D8,T>,               "@@depressions/Lindsay2016.hpp:Lindsay2016@@"); //TODO  
  m.def("rdBreachDepressionsD4",   &BreachDepressions<Topology::D4,T>,               "@@depressions/Lindsay2016.hpp:Lindsay2016@@"); //TODO  

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

  m.def("FA_Tarboton",            &FA_Tarboton           <T,double>, "TODO");
  m.def("FA_Dinfinity",           &FA_Dinfinity          <T,double>, "TODO");
  m.def("FA_Holmgren",            &FA_Holmgren           <T,double>, "TODO");
  m.def("FA_Quinn",               &FA_Quinn              <T,double>, "TODO");
  m.def("FA_Freeman",             &FA_Freeman            <T,double>, "TODO");
  m.def("FA_FairfieldLeymarieD8", &FA_FairfieldLeymarieD8<T,double>, "TODO");
  m.def("FA_FairfieldLeymarieD4", &FA_FairfieldLeymarieD4<T,double>, "TODO");
  m.def("FA_Rho8",                &FA_Rho8               <T,double>, "TODO");
  m.def("FA_Rho4",                &FA_Rho4               <T,double>, "TODO");
  m.def("FA_D8",                  &FA_D8                 <T,double>, "TODO");
  m.def("FA_D4",                  &FA_D4                 <T,double>, "TODO");
  m.def("FA_OCallaghanD8",        &FA_OCallaghanD8       <T,double>, "TODO");
  m.def("FA_OCallaghanD4",        &FA_OCallaghanD4       <T,double>, "TODO");

  m.def("FM_Tarboton",            &FM_Tarboton          <T>,              "TODO");
  m.def("FM_Dinfinity",           &FM_Dinfinity         <T>,              "TODO");
  m.def("FM_Holmgren",            &FM_Holmgren          <T>,              "TODO");
  m.def("FM_Quinn",               &FM_Quinn             <T>,              "TODO");
  m.def("FM_Freeman",             &FM_Freeman           <T>,              "TODO");
  m.def("FM_FairfieldLeymarieD8", &FM_FairfieldLeymarie <Topology::D8,T>, "TODO");
  m.def("FM_FairfieldLeymarieD4", &FM_FairfieldLeymarie <Topology::D4,T>, "TODO");
  m.def("FM_Rho8",                &FM_Rho8              <T>,              "TODO");
  m.def("FM_Rho4",                &FM_Rho4              <T>,              "TODO");
  m.def("FM_OCallaghanD8",        &FM_OCallaghan        <Topology::D8,T>, "TODO");
  m.def("FM_OCallaghanD4",        &FM_OCallaghan        <Topology::D4,T>, "TODO");
  m.def("FM_D8",                  &FM_D8                <T>,              "TODO");
  m.def("FM_D4",                  &FM_D4                <T>,              "TODO");

  py::class_<Array2D<T>>(m, ("Array2D_" + tname).c_str(), py::buffer_protocol(), py::dynamic_attr())
      .def(py::init<>())
      .def(py::init<typename Array2D<T>::xy_t, typename Array2D<T>::xy_t,T>())
      
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

        auto buf = py::array_t<T, py::array::c_style | py::array::forcecast>::ensure(src);
        if (!buf)
          throw std::runtime_error("Unable to convert array to RichDEM object!");

        //TODO: CHeck stride
        auto dims = buf.ndim();
        if (dims != 2 )
          throw std::runtime_error("Array must have two dimensions!");

        return new Array2D<T>((T*)buf.data(), buf.shape()[1], buf.shape()[0]);
      }))

      .def("size",      &Array2D<T>::size)
      .def("width",     &Array2D<T>::width)
      .def("height",    &Array2D<T>::height)
      .def("empty",     &Array2D<T>::empty)
      .def("noData",    &Array2D<T>::noData)
      .def("min",       &Array2D<T>::min)
      .def("max",       &Array2D<T>::max)
      
      //TODO: Simplify by casting to double in Python
      .def("setNoData", [](Array2D<T> &a, const float    ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](Array2D<T> &a, const double   ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](Array2D<T> &a, const int8_t   ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](Array2D<T> &a, const int16_t  ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](Array2D<T> &a, const int32_t  ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](Array2D<T> &a, const int64_t  ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](Array2D<T> &a, const uint8_t  ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](Array2D<T> &a, const uint16_t ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](Array2D<T> &a, const uint32_t ndval){ a.setNoData((T)ndval); })
      .def("setNoData", [](Array2D<T> &a, const uint64_t ndval){ a.setNoData((T)ndval); })
            
      .def_readwrite("geotransform", &Array2D<T>::geotransform)
      .def_readwrite("projection",   &Array2D<T>::projection)
      .def_readwrite("metadata",     &Array2D<T>::metadata)
      .def("copy", [](const Array2D<T> a){
        return a;
      })

      // .def_buffer([](Array2D<T> &arr) -> py::buffer_info {
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
        [=](const Array2D<T> &a) {
            return "<RichDEM array: type="+tname+", width="+std::to_string(a.width())+", height="+std::to_string(a.height())+", owned="+std::to_string(a.owned())+">";
        }
      )
      .def("__call__",
        [](Array2D<T> &a, const int x, const int y) -> T& {
          return a(x,y);
        }
      )
      .def("__call__",
        [](Array2D<T> &a, const int i) -> T& {
          return a(i);
        }
      );
}



PYBIND11_MODULE(_richdem, m) {
  m.doc() = "Internal library used by pyRichDEM for calculations";

  //py::bind_vector<std::vector<double>>(m, "VecDouble");
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
