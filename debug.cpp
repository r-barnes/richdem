#include <deque>
#include <fstream>
#include <string>
#include "data_structures.h"
#include "interface.h"
#include "flat_resolution.h"
#include "dinf_methods.h"

void print_edges(const float_2d &elevations, const std::deque<grid_cell> &low_edges, const std::deque<grid_cell> &high_edges){
	for(int y=0;y<elevations.height();y++){
		for(int x=0;x<elevations.width();x++){
			for(std::deque<grid_cell>::const_iterator i=low_edges.begin();i!=low_edges.end();i++)
				if(x==i->x && y==i->y){
					diagnostic("\033[36m");
					goto printedges_done;
				}
			for(std::deque<grid_cell>::const_iterator i=high_edges.begin();i!=high_edges.end();i++)
				if(x==i->x && y==i->y){
					diagnostic("\033[31m");
					goto printedges_done;
				}
			printedges_done: diagnostic_arg("%2.0f\033[39m ",elevations(x,y));
		}
		diagnostic("\n");
	}
}



void dinf_pit_flows(const float_2d &elevations, float_2d &flowdirs){
	flowdirs.copyprops(elevations);
	dinf_flow_directions(elevations,flowdirs);

	int_2d flat_resolution_mask(elevations), groups;
	resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);

	dinf_flow_flats(flat_resolution_mask,groups,flowdirs);
}

void tikz_flowdir_print(const char_2d &flowdirs, std::string filename, float x_scale, float y_scale, float x_offset, float y_offset, bool omit_edges){
	std::ofstream fout;

	diagnostic_arg("%f %f %f %f\n",x_scale,y_scale,x_offset,y_offset);

	diagnostic_arg("Opening TikZ flowdirs output file \"%s\"...",filename.c_str());
	fout.open(filename.c_str());
	if(!fout.is_open()){
		diagnostic("failed!\n");
		exit(-1);
	}
	diagnostic("succeeded.\n");

	for(int y=(omit_edges?1:0);y<flowdirs.height()-(omit_edges?1:0);y++)
	for(int x=(omit_edges?1:0);x<flowdirs.width()-(omit_edges?1:0);x++){
		fout<<"\\node at ("<<((x)*x_scale+x_offset)<<","<<(((flowdirs.height()-1)-y)*y_scale+y_offset)<<") {$\\";
		switch(flowdirs(x,y)){
			case 1:
				fout<<"leftarrow";break;
			case 2:
				fout<<"nwarrow";break;
			case 3:
				fout<<"uparrow";break;
			case 4:
				fout<<"nearrow";break;
			case 5:
				fout<<"eastarrow";break;
			case 6:
				fout<<"searrow";break;
			case 7:
				fout<<"downarrow";break;
			case 8:
				fout<<"swarrow";break;
		}
		fout<<"$};"<<std::endl;
	}
	fout.close();
}
