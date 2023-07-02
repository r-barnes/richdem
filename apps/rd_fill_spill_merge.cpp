#include <richdem/common/Array2D.hpp>
#include <richdem/common/gdal.hpp>
#include <richdem/depressions/fill_spill_merge.hpp>
#include <richdem/misc/misc_methods.hpp>
#include <richdem/ui/cli_options.hpp>

#include <iostream>
#include <string>
#include <stdexcept>

namespace rd = richdem;
namespace dh = richdem::dephier;

int main(int argc, char **argv){
  CLI::App app("Fill-Spill-Merge Example Program");

  std::string topography_filename;
  std::string output_prefix;
  double      surface_water_level = std::numeric_limits<double>::quiet_NaN();
  std::string surface_water_filename;
  double      ocean_level;

  size_t num;
  size_t len;
  app.add_option("topography", topography_filename, "Topography to run FSM on")->required();
  app.add_option("output",     output_prefix,       "Path of GeoTiff output file")->required();
  const auto swl_ptr = app.add_option("--swl", surface_water_level, "Surface water level as a numeric constant");
  app.add_option("--swf", surface_water_filename, "File containing surface water levels")->excludes(swl_ptr);
  app.add_option("ocean_level",         ocean_level, "Elevation of the ocean")->required();

  CLI11_PARSE(app, argc, argv);

  std::cout<<"m Input DEM           = "<<topography_filename   <<std::endl;
  std::cout<<"m Output prefix       = "<<output_prefix         <<std::endl;
  std::cout<<"m Surface water level = "<<surface_water_level   <<std::endl;
  std::cout<<"m Surface water file  = "<<surface_water_filename<<std::endl;
  std::cout<<"m Ocean level         = "<<ocean_level           <<std::endl;

  std::cout<<"p Reading topography..."<<std::endl;
  rd::Timer timer_io;
  timer_io.start();
  rd::Array2D<double> topo(topography_filename);
  timer_io.stop();

  std::cout<<"m Data width  = "<<topo.width ()<<std::endl;
  std::cout<<"m Data height = "<<topo.height()<<std::endl;
  std::cout<<"m Data cells  = "<<topo.numDataCells()<<std::endl;

  rd::Array2D<double> wtd;
  if(surface_water_filename.empty()){
    //All cells have the same amount of water
    wtd.resize(topo.width(), topo.height(), surface_water_level);
  } else {
    wtd = rd::Array2D<double>(surface_water_filename);
  }

  rd::Array2D<dh::dh_label_t> label   (topo.width(), topo.height(), dh::NO_DEP );      //No cells are part of a depression
  rd::Array2D<rd::flowdir_t>  flowdirs(topo.width(), topo.height(), rd::NO_FLOW);      //No cells flow anywhere

  wtd.setNoData(topo.noData());

  //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  std::cout<<"p Performing bucket fill..."<<std::endl;
  rd::BucketFillFromEdges<rd::Topology::D8>(topo, label, ocean_level, dh::OCEAN);

  //Make NoData cells also ocean cells. Ocean has no water on it to begin with.
  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++){
    if(topo.isNoData(i) || label(i)==dh::OCEAN){
      label(i) = dh::OCEAN;
      wtd  (i) = 0;
    }
  }

  rd::Timer timer_calc;
  timer_calc.start();

  //Generate flow directions, label all the depressions, and get the hierarchy
  //connecting them
  std::cout<<"p Getting depression hierarchy..."<<std::endl;
  auto deps = dh::GetDepressionHierarchy<double,rd::Topology::D8>(topo, label, flowdirs);

  std::cout<<"p Performing FillSpillMerge..."<<std::endl;
  dh::FillSpillMerge(topo, label, flowdirs, deps, wtd);

  timer_calc.stop();

  timer_io.start();

  //Output the water table depth
  wtd.saveGDAL(output_prefix+"-wtd.tif");

  for(unsigned int i=0;i<topo.size();i++)
    if(!topo.isNoData(i))
      wtd(i) += topo(i);

  //Output the new height of the hydraulic surface
  wtd.saveGDAL(output_prefix+"-hydrologic-surface-height.tif");

  timer_io.stop();

  std::cout<<"Finished."<<std::endl;
  std::cout<<"IO time   = "<<timer_io.accumulated()  <<" s"<<std::endl;
  std::cout<<"Calc time = "<<timer_calc.accumulated()<<" s"<<std::endl;

  return 0;
}
