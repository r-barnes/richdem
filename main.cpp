#define  RICHDEM_VERSION  "1.0.0"
#include "utility.hpp"
#include "data_structures.hpp"
#include "data_io.hpp"
#include "d8_methods.hpp"
#include "dinf_methods.hpp"
#include "pit_fill.hpp"
#include "interface.hpp"
#include "flat_resolution.hpp"
#include "debug.hpp"
#include <string>

#include "tclap/CmdLine.h"

int main(int argc, char **argv){
  setbuf ( stderr , NULL );
#ifdef ARCGIS
  setbuf ( stdout , NULL );
#endif

  Timer running_calc_time,running_io_time,total_time;

  total_time.start();

  TCLAP::CmdLine cmd("RichDEM is a suite of DEM analysis functions for determining hydrologic properties. It has been developed by Richard Barnes (rbarnes@umn.edu). Find RichDEM on the web at \"http://www.richdem.com\".", ' ', RICHDEM_VERSION);
  TCLAP::SwitchArg cl_d8("8","d8","Use the D8 flow metric (Dinf is default)", cmd, false);
  TCLAP::SwitchArg cl_fill_pits("p","pits","Perform pit-filling prior to other operations", cmd, false);
  TCLAP::SwitchArg cl_fill_pfdirs("r","pfdirs","Use Priority-Flood+Flow Directions algorithm.", cmd, false);
  TCLAP::ValueArg<std::string> cl_output_pit_filled("l","pitfilled","Output pit-filled DEM - only applicable if -p is specified",false,"","file",cmd);
  TCLAP::ValueArg<std::string> cl_output_unresolved_flowdirs("u","uflowdirs","Output flow directions before flat resolution",false,"","file",cmd);
  TCLAP::ValueArg<std::string> cl_output_resolved_flowdirs("f","flowdirs","Output flow directions after flat resolution",false,"","file",cmd);
  TCLAP::ValueArg<std::string> cl_output_flow_acculm("a","acculm","Output flow accumulation (aka: contributing area, upslope area)",false,"","file",cmd); //TODO: Are these really all equivalent?
  TCLAP::ValueArg<std::string> cl_output_spi("","spi","Output SPI (stream power index)",false,"","file",cmd);
  TCLAP::ValueArg<std::string> cl_output_cti("","cti","Output CTI (compound topographic index, wetness index)",false,"","file",cmd);
//  TCLAP::SwitchArg cl_are_there_dams("d","dams","Determines if the input file has digital dams. All other options ignored.", cmd, false); //TODO
  TCLAP::ValueArg<std::string> cl_output_tikz("","tikz","Output TikZ flow directions",false,"","file",cmd);
  TCLAP::ValueArg<float> cl_output_tikzx("","tikzx","TikZ X scaling",false,1.,"file",cmd);
  TCLAP::ValueArg<float> cl_output_tikzy("","tikzy","TikZ Y scaling",false,1.,"file",cmd);
  TCLAP::ValueArg<float> cl_output_tikzxo("","tikzxo","TikZ X offset",false,0.,"file",cmd);
  TCLAP::ValueArg<float> cl_output_tikzyo("","tikzyo","TikZ X offset",false,0.,"file",cmd);
  TCLAP::SwitchArg cl_tikz_omit_edges("","tikzne","Omits edges of the grid in TikZ output", cmd, true);
  TCLAP::UnlabeledValueArg<std::string> cl_inputDEM("inputDEM", "The DEM to be processed (must be in ArcGrid ASCII format)", true, "", "Input DEM", cmd);
  try {
    cmd.parse( argc, argv );
  } catch (TCLAP::ArgException &e) {
    std::cerr << "Command line argument (" << e.argId() << ") error: " << e.error() << std::endl;
  }

  diagnostic_arg("Running RichDEMv%s.\n",RICHDEM_VERSION);

  float_2d elevations;
  running_io_time.start();
  load_ascii_data(cl_inputDEM.getValue().c_str(),elevations);
  running_io_time.stop();

//  if(cl_are_there_dams.getValue()){
//    diagnostic_arg("Found %d digital dams.\n",digital_dams(elevations));
//    return 0;
//  }

  if(cl_fill_pfdirs.getValue()){
    char_2d flowdirs;

    running_calc_time.start();
    priority_flood_flowdirs(elevations,flowdirs);
    running_calc_time.stop();

    if(!cl_output_resolved_flowdirs.getValue().empty()){
      running_io_time.start();
      output_ascii_data(cl_output_resolved_flowdirs.getValue().c_str(), flowdirs);
      running_io_time.stop();
    }

    if(!cl_output_tikz.getValue().empty()){
      running_io_time.start();
      tikz_flowdir_print(flowdirs, cl_output_tikz.getValue(), cl_output_tikzx.getValue(), cl_output_tikzy.getValue(), cl_output_tikzxo.getValue(), cl_output_tikzyo.getValue(), cl_tikz_omit_edges.getValue());
      running_io_time.stop();
    }

    if(!cl_output_flow_acculm.getValue().empty()){
      if(cl_d8.getValue()){
        int_2d area(elevations);

        running_calc_time.start();
        d8_upslope_area(flowdirs, area);
        running_calc_time.stop();

        running_io_time.start();
        output_ascii_data(cl_output_flow_acculm.getValue().c_str(),area);
        running_io_time.stop();
      } else
        diagnostic("Dinf has not been implemented yet for Priority-Flood+FlowDirections.\n"); //TODO
    }

    goto end_main;
  }

  if(cl_fill_pits.getValue()){
    running_calc_time.start();
    barnes_flood(elevations);
    running_calc_time.stop();
    if(!cl_output_pit_filled.getValue().empty()){
      running_io_time.start();
      output_ascii_data(cl_output_pit_filled.getValue().c_str(),elevations);
      running_io_time.stop();
    }
  }


  if(cl_d8.getValue()){
    char_2d flowdirs(elevations);

    running_calc_time.start();
    d8_flow_directions(elevations,flowdirs);
    running_calc_time.stop();

    if(!cl_output_unresolved_flowdirs.getValue().empty()){
      running_io_time.start();
      output_ascii_data(cl_output_unresolved_flowdirs.getValue().c_str(),flowdirs);
      running_io_time.stop();
    }

    {  //Used to try to free memory after the process
      int_2d flat_resolution_mask(elevations), groups;
      resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);
      d8_flow_flats(flat_resolution_mask,groups,flowdirs);

      flat_resolution_mask.clear();
      groups.clear();

      if(!cl_output_resolved_flowdirs.getValue().empty())
        output_ascii_data(cl_output_resolved_flowdirs.getValue().c_str(),flowdirs);
    }

    if(!cl_output_flow_acculm.getValue().empty()){
      int_2d area(elevations);

      running_calc_time.start();
      d8_upslope_area(flowdirs, area);
      running_calc_time.stop();

      running_io_time.start();
      output_ascii_data(cl_output_flow_acculm.getValue().c_str(),area);
      running_io_time.stop();
    }
  } else {
    float_2d flowdirs(elevations);
    dinf_flow_directions(elevations,flowdirs);
    if(!cl_output_unresolved_flowdirs.getValue().empty())
      output_ascii_data(cl_output_unresolved_flowdirs.getValue().c_str(),flowdirs);

    {  //Used to try to free memory after the process
      int_2d flat_resolution_mask(elevations), groups;
      resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);
      dinf_flow_flats(flat_resolution_mask,groups,flowdirs);

      flat_resolution_mask.clear();
      groups.clear();

      if(!cl_output_resolved_flowdirs.getValue().empty())
        output_ascii_data(cl_output_resolved_flowdirs.getValue().c_str(),flowdirs);
    }

    if(!cl_output_flow_acculm.getValue().empty()){
      float_2d area(elevations);

      running_calc_time.start();
      dinf_upslope_area(flowdirs, area);
      running_calc_time.stop();

      running_io_time.start();
      output_ascii_data(cl_output_flow_acculm.getValue().c_str(),area);
      running_io_time.stop();
    }
  }


end_main:
  total_time.stop();
  diagnostic_arg("Total time was: %lfs\n", total_time.accumulated());
  diagnostic_arg("Time spent in computation: %lfs\n",running_calc_time.accumulated());
  diagnostic_arg("Time spent in I/O: %lfs\n",running_io_time.accumulated());
  diagnostic_arg("Ran RichDEMv%s.\n",RICHDEM_VERSION);

  return 0;
}
