#ifndef _flat_resolution_included
#define _flat_resolution_included

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
	int x,y,nx,ny;
	int loops=1;
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
			continue;
		}

		if(incrementations(x,y)>0) continue;	//I've already been incremented!

		//Should I increment?
		if(loops>1){
			int n;
			for(n=1;n<=8;n++)
				if(incrementations(x+dx[n],y+dy[n])>0)
					break;
			if(n>8) continue;
		}

		//If I incremented, maybe my neighbours should too
		incrementations(x,y)=loops;
		flat_height[groups(x,y)]=loops;
		for(int n=1;n<=8;n++){
			nx=x+dx[n];	
			ny=y+dy[n];
			if(IN_GRID(nx,ny,elevations.width(),elevations.height()) 
					&& elevations(nx,ny)==elevations(x,y) 
					&& flowdirs(nx,ny)==NO_FLOW) //TODO: Need to check if elevation is equal because this could be a low edge abutting another flat's high edge
				edges.push_back(grid_cell(nx,ny));
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
			if(IN_GRID(nx,ny,elevations.width(),elevations.height()) && elevations(nx,ny)==elevations(x,y))
				edge.push_back(grid_cell(nx,ny));
		}
		if(inc1(x,y)>0){
			flat_resolution_mask(x,y)=2*(inc1(x,y)-1);
			if(inc2(x,y)>0)
				flat_resolution_mask(x,y)+=flat_height[groups(x,y)]-inc2(x,y)+1;
		}
			
		inc1(x,y)=-1;
	}
}

template<class T>
void label_this(int x, int y, const T target_elevation, const int label, int_2d &groups, const array2d<T> &elevations){
	std::queue<grid_cell> to_fill;
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
void resolve_flats(const array2d<T> &elevations, const array2d<U> &flowdirs, int_2d &flat_resolution_mask, int_2d &groups){
	std::deque<grid_cell> low_edges,high_edges;	//TODO: Need estimate of size

	diagnostic_arg("The groups matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(int)/1024/1024);
	diagnostic("Setting up groups matrix...");
	groups.resize(flowdirs.width(),flowdirs.height(),false);
	groups.init(-1);
	diagnostic("succeeded.\n");

	diagnostic_arg("The flat resolution mask will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(int)/1024/1024);
	diagnostic("Setting up flat resolution mask...");
	flat_resolution_mask.resize(flowdirs.width(),flowdirs.height(),false);
	flat_resolution_mask.init(-1);
	flat_resolution_mask.no_data=-1;
	diagnostic("succeeded!\n");

	find_flat_edges(low_edges, high_edges, flowdirs, elevations);

	if(low_edges.size()==0){
		if(high_edges.size()>0)
			diagnostic("There were flats, but none of them had outlets!\n");
		else
			diagnostic("There were no flats!\n");
		return;
	}

	diagnostic("Labeling flats...");
	int group_number=0;
	for(std::deque<grid_cell>::iterator i=low_edges.begin();i!=low_edges.end();i++)
		if(groups(i->x,i->y)==-1)
			label_this(i->x, i->y, elevations(i->x,i->y), group_number++, groups, elevations);
	diagnostic("succeeded!\n");

	diagnostic_arg("Found %d unique flats.\n",group_number);

	diagnostic("Removing flats without outlets from the queue...");
	std::deque<grid_cell> temp;
	for(std::deque<grid_cell>::iterator i=high_edges.begin();i!=high_edges.end();i++)
		if(groups(i->x,i->y)!=-1)
			temp.push_back(*i);
	diagnostic("succeeded.\n");

	if(temp.size()<high_edges.size())	//TODO: Prompt for intervention?
		diagnostic("\033[91mNot all flats have outlets; the DEM contains sinks/pits/depressions!\033[39m\n");
	high_edges=temp;
	temp.clear();

	diagnostic_arg("The incrementation matricies will require approximately %ldMB of RAM.\n",3*flowdirs.width()*flowdirs.height()*sizeof(int)/1024/1024);
	diagnostic("Setting up incrementation matricies...");
	int_2d inc1(elevations,true);
	int_2d inc2(elevations,true);
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
	diagnostic("succeeded!\n");

	diagnostic("Combining Barnes flat resolution steps...\n");
	BarnesStep3(elevations, inc1, inc2, flat_resolution_mask, low_edges, flat_height, groups);
	diagnostic("succeeded!\n");
}

#endif
