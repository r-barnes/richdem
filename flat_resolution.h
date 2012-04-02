#ifndef _flat_resolution_included
#define _flat_resolution_included

#include "flat_resolution.h"
#include "utility.h"
#include "interface.h"
#include "data_structures.h"
#include <omp.h>
#include <deque>
#include <vector>
#include <queue>
#include "debug.h"

template <class T, class U>
int BarnesStep(const array2d<T> &elevations, const array2d<U> &flowdirs, int_2d &incrementations, std::deque<grid_cell> edges, std::vector<int> &flat_height, const int_2d &groups){
	int loops=1;
	int x,y,nx,ny;
	bool past_edge=false;
	bool do_increment;
	grid_cell iteration_marker(-1,-1);

	//Incrementation
	edges.push_back(iteration_marker);
	while(edges.size()!=1){	//Only iteration marker is left in the end
		x=edges.front().x;
		y=edges.front().y;
		edges.pop_front();

		if(x==-1){	//I'm an iteration marker
			loops++;
			edges.push_back(iteration_marker);
			past_edge=true;
			continue;
		}

		if(incrementations(x,y)>0) continue;	//I've already been incremented!

		//Should I increment?
		if(past_edge){
			do_increment=false;
			for(int n=1;n<=8;n++){
				nx=x+dx[n];
				ny=y+dy[n];
				if(incrementations(nx,ny)>0){
					do_increment=true;
					break;
				}
			}
		} else
			do_increment=true;

		//If I incremented, maybe my neighbours should too
		if(do_increment){
			incrementations(x,y)=loops;
			flat_height[groups(x,y)]=loops;
			for(int n=1;n<=8;n++){
				nx=x+dx[n];	
				ny=y+dy[n];
				if(!(IN_GRID(nx,ny,elevations.width(),elevations.height()) 
						&& elevations(nx,ny)==elevations(x,y) 
						&& flowdirs(nx,ny)==NO_FLOW)) continue; //TODO: Shouldn't need to check if elevation is equal.
				edges.push_back(grid_cell(nx,ny));
			}
		}
	}

	return loops;
}

template <class T>
void BarnesStep3(const array2d<T> &elevations, int_2d &inc1, int_2d &inc2, int_2d &flat_resolution_mask, std::deque<grid_cell> &edge, const std::vector<int> &flat_height, const int_2d &groups){
	int x,y,nx,ny;

	while(edge.size()!=0){
		x=edge.front().x;
		y=edge.front().y;
		edge.pop_front();

		if(inc1(x,y)==-1) continue;

		for(int n=1;n<=8;n++){
			nx=x+dx[n];
			ny=y+dy[n];
			if(!(IN_GRID(nx,ny,elevations.width(),elevations.height()) && elevations(nx,ny)==elevations(x,y))) continue;
			edge.push_back(grid_cell(nx,ny));
		}
		if(inc2(x,y)>0)
			flat_resolution_mask(x,y)=2*(inc1(x,y)-1)+(flat_height[groups(x,y)]-inc2(x,y)+1);
		else
			flat_resolution_mask(x,y)=0;
		inc1(x,y)=-1;
	}
}

template<class T>
void label_this(int x, int y, const T target_elevation, const int label, int_2d &groups, const array2d<T> &elevations){
	std::queue<grid_cell> to_fill; //TODO: Queue versus stack
	to_fill.push(grid_cell(x,y));

	while(to_fill.size()>0){
		x=to_fill.front().x;
		y=to_fill.front().y;
		to_fill.pop();
		if(elevations(x,y)!=target_elevation || groups(x,y)>-1) continue;
		groups(x,y)=label;
		for(int n=1;n<=8;n++)
			if(IN_GRID(x+dx[n],y+dy[n],groups.width(),groups.height()))
				to_fill.push(grid_cell(x+dx[n],y+dy[n]));
	}
}

template <class T, class U>
void find_flat_edges(std::deque<grid_cell> &low_edges, std::deque<grid_cell> &high_edges, const array2d<U> &flowdirs, const array2d<T> &elevations){
	int nx,ny;
	diagnostic("\r\033[2KSearching for flats...\n");
	progress_bar(-1);
	for(int x=0;x<flowdirs.width();x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.height()*100/(flowdirs.width()*flowdirs.height()));
		for(int y=0;y<flowdirs.height();y++){
			if(flowdirs(x,y)==flowdirs.no_data)
				continue;
			for(int n=1;n<=8;n++){
				nx=x+dx[n];
				ny=y+dy[n];

				if(!IN_GRID(nx,ny,flowdirs.width(),flowdirs.height())) continue;
				if(flowdirs(nx,ny)==flowdirs.no_data) continue; //TODO: This isn't really necessary, but it is safe.

				if(flowdirs(x,y)!=NO_FLOW && flowdirs(nx,ny)==NO_FLOW && elevations(nx,ny)==elevations(x,y))
					low_edges.push_back(grid_cell(x,y));
				else if(flowdirs(x,y)==NO_FLOW && elevations(x,y)<elevations(nx,ny))
					high_edges.push_back(grid_cell(x,y));
			}
		}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}

template <class T, class U>
int resolve_flats(const array2d<T> &elevations, const array2d<U> &flowdirs, int_2d &flat_resolution_mask, int_2d &groups){
	std::deque<grid_cell> low_edges,high_edges;	//TODO: Need estimate of size

	diagnostic_arg("The groups matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(int)/1024/1024);
	diagnostic("Resizing groups matrix...");
	groups.resize(flowdirs.width(),flowdirs.height(),false);
	diagnostic("succeeded.\n");
	diagnostic("Initializing groups matrix...");
	groups.init(-1);
	diagnostic("succeeded.\n");

	diagnostic_arg("The flat resolution mask will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(int)/1024/1024);
	diagnostic("Resizing flat resolution mask...");
	flat_resolution_mask.resize(flowdirs.width(),flowdirs.height(),false);
	diagnostic("succeeded.\n");
	diagnostic("Initializing flat resolution mask...");
	flat_resolution_mask.init(2147483000);
	diagnostic("succeeded.\n");
	diagnostic("Setting flat resolution mask's no_data...");
	flat_resolution_mask.no_data=2147483000;
	diagnostic("succeeded!\n");

	find_flat_edges(low_edges, high_edges, flowdirs, elevations);

	if(low_edges.size()==0){
		if(high_edges.size()>0)
			diagnostic("There were flats, but none of them had outlets!\n");
		else
			diagnostic("There were no flats!\n");
		return 0;
	}

	diagnostic("Labeling flats...");
	int group_number=0;
	for(std::deque<grid_cell>::iterator i=low_edges.begin();i!=low_edges.end();i++)
		if(groups(i->x,i->y)==-1)
			label_this(i->x, i->y, elevations(i->x,i->y), group_number++, groups, elevations);
	diagnostic("succeeded!\n");

	diagnostic_arg("Found %d unique flats.\n",group_number);

	diagnostic("Removing flats without outlets from the queue...");
	bool flat_without_outlet=false;
	std::deque<grid_cell> temp;
	for(std::deque<grid_cell>::iterator i=high_edges.begin();i!=high_edges.end();i++)	//i++ in 'else'
		if(groups(i->x,i->y)!=-1)
			temp.push_back(*i);
		else
			flat_without_outlet=true;
	high_edges=temp;
	temp.clear();
	diagnostic("succeeded.\n");

	if(flat_without_outlet)	//TODO: Prompt user for intervention?
		diagnostic("\033[91mNot all flats have outlets; the DEM contains sinks/pits/depressions!\033[39m\n");

	diagnostic_arg("The incrementation matricies will require approximately %ldMB of RAM.\n",3*flowdirs.width()*flowdirs.height()*sizeof(int)/1024/1024);
	diagnostic("Creating incrementation matricies...");
	int_2d inc1(elevations,true);
	int_2d inc2(elevations,true);
	diagnostic("succeeded!\n");

	diagnostic("Initializing incrementation matricies...");
	inc1.init(0);
	inc2.init(0);
	diagnostic("succeeded!\n");

	diagnostic_arg("The flat height vector will require approximately %ldMB of RAM.\n",group_number*sizeof(int)/1024/1024);
	diagnostic("Creating flat height vector...");
	std::vector<int> flat_height(group_number);
	diagnostic("succeeded!\n");

	diagnostic("Performing Barnes flat resolution...\n");
	//Flat_height is used with low_edges, but the results will be overwritten by high_edges
	BarnesStep(elevations, flowdirs, inc1, low_edges, flat_height, groups);
	BarnesStep(elevations, flowdirs, inc2, high_edges, flat_height, groups);

	diagnostic("Combining Barnes flat resolution steps...\n");
	BarnesStep3(elevations, inc1, inc2, flat_resolution_mask, low_edges, flat_height, groups);
	diagnostic("succeeded!\n");

	return 0;
}

#endif
