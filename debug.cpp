#include <deque>
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
