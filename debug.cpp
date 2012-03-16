#include <deque>
#include "data_structures.h"
#include "interface.h"

void print_edges(float_2d &elevations, std::deque<grid_cell> &low_edges, std::deque<grid_cell> &high_edges){
	for(int y=0;y<elevations.height();y++){
		for(int x=0;x<elevations.width();x++){
			for(std::deque<grid_cell>::iterator i=low_edges.begin();i!=low_edges.end();i++)
				if(x==i->x && y==i->y){
					diagnostic("\033[36m");
					goto printedges_done;
				}
			for(std::deque<grid_cell>::iterator i=high_edges.begin();i!=high_edges.end();i++)
				if(x==i->x && y==i->y){
					diagnostic("\033[31m");
					goto printedges_done;
				}
			printedges_done: diagnostic_arg("%2.0f\033[39m ",elevations(x,y));
		}
		diagnostic("\n");
	}
}
