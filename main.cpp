#define  RICHDEM_VERSION  "0r~219"
#include "utility.h"
#include "data_structures.h"
#include "data_io.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"
#include "flat_resolution.h"
#include "watershed.h"
//#include "debug.h"
#include <string>

#include "tclap/CmdLine.h"

int main(int argc, char **argv){
	setbuf ( stderr , NULL );
#ifdef ARCGIS
	setbuf ( stdout , NULL );
#endif
	TCLAP::CmdLine cmd("RichDEM is a suite of DEM analysis functions for determining hydrologic properties. It has been developed by Richard Barnes (rbarnes@umn.edu). Find RichDEM on the web at \"http://www.richdem.com\".", ' ', RICHDEM_VERSION);
	TCLAP::SwitchArg cl_d8("8","d8","Use the D8 flow metric (Dinf is default)", cmd, false);
	TCLAP::SwitchArg cl_fill_pits("p","pits","Perform pit-filling prior to other operations", cmd, false);
	TCLAP::ValueArg<std::string> cl_output_pit_filled("l","pitfilled","Output pit-filled DEM - only applicable if -p is specified",false,"","file",cmd);
	TCLAP::ValueArg<std::string> cl_output_unresolved_flowdirs("u","uflowdirs","Output flow directions before flat resolution",false,"","file",cmd);
	TCLAP::ValueArg<std::string> cl_output_resolved_flowdirs("f","flowdirs","Output flow directions after flat resolution",false,"","file",cmd);
	TCLAP::ValueArg<std::string> cl_output_flow_acculm("a","acculm","Output flow accumulation (aka: contributing area, upslope area)",false,"","file",cmd); //TODO: Are these really all equivalent?
	TCLAP::UnlabeledValueArg<std::string> cl_inputDEM("inputDEM", "The DEM to be processed (must be in ArcGrid ASCII format)", true, "", "Input DEM", cmd);
	try {
		cmd.parse( argc, argv );
	//	int bob = nameArg.getValue();
	//	bool reverseName = reverseSwitch.getValue();
	} catch (TCLAP::ArgException &e) {
		std::cerr << "Command line argument (" << e.argId() << ") error: " << e.error() << std::endl;
	}

	diagnostic_arg("Running RichDEMv%s.\n",RICHDEM_VERSION);

	float_2d elevations;
	load_ascii_data(cl_inputDEM.getValue().c_str(),elevations);

	if(cl_fill_pits.getValue()){
		barnes_flood(elevations);
		if(!cl_output_pit_filled.getValue().empty())
			output_ascii_data(cl_output_pit_filled.getValue().c_str(),elevations);
	}


	if(cl_d8.getValue()){	//TODO: D8 functions haven't been tested lately
		char_2d flowdirs(elevations);
		d8_flow_directions(elevations,flowdirs);
		if(!cl_output_unresolved_flowdirs.getValue().empty())
			output_ascii_data(cl_output_unresolved_flowdirs.getValue().c_str(),flowdirs);

		{
			int_2d flat_resolution_mask(elevations), groups;
			resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);

			d8_flow_flats(flat_resolution_mask,groups,flowdirs);
			if(!cl_output_resolved_flowdirs.getValue().empty())
				output_ascii_data(cl_output_resolved_flowdirs.getValue().c_str(),flowdirs);
		}

		int_2d area(elevations);
		d8_upslope_area(flowdirs, area);
		if(!cl_output_flow_acculm.getValue().empty())
			output_ascii_data(cl_output_flow_acculm.getValue().c_str(),area);
	} else {
		float_2d flowdirs(elevations);
		dinf_flow_directions(elevations,flowdirs);
		if(!cl_output_unresolved_flowdirs.getValue().empty())
			output_ascii_data(cl_output_unresolved_flowdirs.getValue().c_str(),flowdirs);

		{
			int_2d flat_resolution_mask(elevations), groups;
			resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);

			dinf_flow_flats(flat_resolution_mask,groups,flowdirs);
			if(!cl_output_resolved_flowdirs.getValue().empty())
				output_ascii_data(cl_output_resolved_flowdirs.getValue().c_str(),flowdirs);
		}

		float_2d area(elevations);
		dinf_upslope_area(flowdirs, area);
		if(!cl_output_flow_acculm.getValue().empty())
			output_ascii_data(cl_output_flow_acculm.getValue().c_str(),area);
	}

	return 0;
}
