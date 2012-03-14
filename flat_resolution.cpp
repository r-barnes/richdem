#include "flat_resolution.h"
#include "utility.h"
#include "interface.h"
#include <omp.h>
#include <string>
#include <deque>
#include <vector>

/*void print_edges(float_2d &elevations, std::deque<grid_cell> &low_edges, std::deque<grid_cell> &high_edges){
	for(int y=0;y<elevations.height();y++){
		for(int x=0;x<elevations.width();x++){
			for(std::deque<grid_cell>::iterator i=low_edges.begin();i!=low_edges.end();i++)
				if(x==i->x && y==i->y){
					printf("\033[36m");
					goto printedges_done;
				}
			for(std::deque<grid_cell>::iterator i=high_edges.begin();i!=high_edges.end();i++)
				if(x==i->x && y==i->y){
					printf("\033[31m");
					goto printedges_done;
				}
			printedges_done: printf("%2.0f\033[39m ",elevations(x,y));
		}
		printf("\n");
	}
}*/

int BarnesStep(const float_2d &elevations, const char_2d &flowdirs, int_2d &incrementations, std::deque<grid_cell> edges, std::vector<int> &flat_height, const int_2d &groups){
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

		if(incrementations(x,y)!=0) continue;	//I've already been incremented!

		//Should I increment?
		if(past_edge){
			do_increment=false;
			for(int n=1;n<=8;n++){
				nx=x+dx[n];
				ny=y+dy[n];
				if(incrementations(nx,ny)>0 || elevations(nx,ny)>elevations(x,y)){
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
				if(!(IN_GRID(nx,ny,elevations.width(),elevations.height()) && elevations(nx,ny)==elevations(x,y) && flowdirs(nx,ny)==NO_FLOW)) continue; //TODO: Shouldn't need to check if elevation is equal.
				edges.push_back(grid_cell(nx,ny));
			}
		}
	}

	return loops;
}

void BarnesStep3(float_2d &elevations, int_2d &inc1, int_2d &inc2, std::deque<grid_cell> &edge, const std::vector<int> &flat_height, const int_2d &groups, float epsilon){
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
		if(inc2(x,y)>0){
			elevations(x,y)+=(epsilon+epsilon/2)*((inc1(x,y)-1))+epsilon*(flat_height[groups(x,y)]-inc2(x,y)+1);
			inc2(x,y)=(inc1(x,y)-1)+(flat_height[groups(x,y)]-inc2(x,y)+1);
		}
		inc1(x,y)=-1;
	}
}

int label_this(const int x, const int y, const float target_elevation, const int label, int_2d &groups, const float_2d &elevations){
	if(!INTERIOR_GRID(x,y,elevations.width(),elevations.height())) return 0;
	if(elevations(x,y)!=target_elevation || groups(x,y)!=-1) return 0;
	groups(x,y)=label;
	label_this(x+1,y,target_elevation,label,groups,elevations);
	label_this(x-1,y,target_elevation,label,groups,elevations);
	label_this(x,y+1,target_elevation,label,groups,elevations);
	label_this(x,y-1,target_elevation,label,groups,elevations);
	label_this(x+1,y+1,target_elevation,label,groups,elevations);
	label_this(x-1,y+1,target_elevation,label,groups,elevations);
	label_this(x+1,y-1,target_elevation,label,groups,elevations);
	label_this(x-1,y-1,target_elevation,label,groups,elevations);
	return 1;
}

int find_flat_edges(std::deque<grid_cell> &low_edges, std::deque<grid_cell> &high_edges, const char_2d &flowdirs, const float_2d &elevations, int_2d &groups){
	int nx,ny;
	int group_number=0;
	diagnostic("\r\033[2KSearching for flats...\n");
	progress_bar(-1);
	for(int x=1;x<flowdirs.width()-1;x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.height()*100/(flowdirs.width()*flowdirs.height()));
		for(int y=1;y<flowdirs.height()-1;y++)
			for(int n=1;n<=8;n++){
				nx=x+dx[n];
				ny=y+dy[n];

				if(!IN_GRID(nx,ny,flowdirs.width(),flowdirs.height())) continue;
				if(flowdirs(nx,ny)==d8_NO_DATA) continue;

				if(flowdirs(x,y)!=NO_FLOW && flowdirs(nx,ny)==NO_FLOW && elevations(nx,ny)==elevations(x,y) && INTERIOR_GRID(nx,ny,flowdirs.width(),flowdirs.height())){
					low_edges.push_back(grid_cell(x,y));
					group_number+=label_this(x, y, elevations(x,y), group_number, groups, elevations);
					break;
				} else if(flowdirs(x,y)==NO_FLOW && elevations(x,y)<elevations(nx,ny)){
					high_edges.push_back(grid_cell(x,y));
					group_number+=label_this(x, y, elevations(x,y), group_number, groups, elevations);
					break;
				}	//TODO: You could think about using a GOTO in the code above to move the label_this statements to the same line somewhere
			}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	return group_number;
}

int resolve_flats(float_2d &elevations, const char_2d &flowdirs){
	std::deque<grid_cell> low_edges,high_edges;	//TODO: Need estimate of size
	int_2d groups(flowdirs);
	int group_max=-1;

	diagnostic_arg("The groups matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(int)/1024/1024);
	diagnostic("Resizing groups matrix...");
	groups.resize(flowdirs.width(),flowdirs.height(),false);
	diagnostic("succeeded.\n");
	diagnostic("Initializing groups matrix...");
	groups.init(-1);
	diagnostic("succeeded.\n");

	diagnostic("Entering find_flat_edges function...");
	group_max=find_flat_edges(low_edges, high_edges, flowdirs, elevations, groups);

	diagnostic_arg("Found %d unique flats.\n",group_max);
//	print2d("%2.0f ",elevations);	//TODO
//	print2d("%2d ",flowdirs);
//	print2d("%2d ",groups);

	if(group_max==0){
		diagnostic("No flats to resolve!\n");
		return 0;
	}

	diagnostic_arg("The boolean outlets vector will require approximately %ldMB of RAM.\n",group_max*sizeof(bool)/1024/1024);
	diagnostic("Creating boolean outlets vector...");
	std::vector<int> has_outlet(group_max,false);
	diagnostic("succeeded!\n");

	diagnostic("Determining which flats have outlets...");
	int outlet_count=0;
	for(std::deque<grid_cell>::iterator i=low_edges.begin();i!=low_edges.end();i++){
		if(!has_outlet[groups(i->x,i->y)]){
			has_outlet[groups(i->x,i->y)]=true;
			outlet_count++;
		}
	}
	diagnostic("succeeded!\n");
	diagnostic_arg("%d out of %d flats had outlets.\n",outlet_count,group_max);
	if(outlet_count!=group_max)	//TODO: Prompt user for intervention?
		diagnostic("\033[91mNot all flats had outlets; the DEM contains sinks/pits/depressions!\033[39m\n");


	diagnostic("Removing flats without outlets from the queue...");
	for(std::deque<grid_cell>::iterator i=low_edges.begin();i!=low_edges.end();)	//i++ in 'else'
		if(!has_outlet[groups(i->x,i->y)])
			i=low_edges.erase(i);
		else
			i++;
	for(std::deque<grid_cell>::iterator i=high_edges.begin();i!=high_edges.end();)	//i++ in 'else'
		if(!has_outlet[groups(i->x,i->y)])
			i=high_edges.erase(i);
		else
			i++;
	diagnostic("succeeded.\n");

//	printedges(elevations,low_edges,high_edges);	//TODO

	diagnostic_arg("The incrementation matricies will require approximately %ldMB of RAM.\n",3*flowdirs.width()*flowdirs.height()*sizeof(int)/1024/1024);
	diagnostic("Creating incrementation matricies...");
	int_2d inc1(elevations,true);
	int_2d inc2(elevations,true);
	int_2d inc3(elevations,true);
	diagnostic("succeeded!\n");

	diagnostic("Initializing incrementation matricies...");
	inc1.init(0);
	inc2.init(0);
	inc3.init(0);
	diagnostic("succeeded!\n");

	diagnostic_arg("The flat height vector will require approximately %ldMB of RAM.\n",group_max*sizeof(int)/1024/1024);
	diagnostic("Creating flat height vector...");
	std::vector<int> flat_height(group_max);
	diagnostic("succeeded!\n");

	diagnostic("Performing Barnes flat resolution...\n");
	BarnesStep(elevations, flowdirs, inc1, low_edges, flat_height, groups);		//Flat_height is used here, but the results will be overwritten by the next line
	BarnesStep(elevations, flowdirs, inc2, high_edges, flat_height, groups);
//	print2d("%d ", inc1);	//TODO
//	print2d("%d ", inc2);	//TODO

	diagnostic("Combining Barnes flat resolution steps...\n");
	BarnesStep3(elevations, inc1, inc2, low_edges, flat_height, groups, 1e-1);
	//print2d("%d ", inc2);	//TODO
	diagnostic("succeeded!\n");

	return 0;
}
