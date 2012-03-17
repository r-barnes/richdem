#include "flat_resolution.h"
#include "utility.h"
#include "interface.h"
#include "data_structures.h"
#include <omp.h>
#include <string>
#include <deque>
#include <vector>
#include <queue>
#include "debug.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_queue.h"

int BarnesStep(const float_2d &elevations, const char_2d &flowdirs, 
				int_2d &incrementations, std::vector<grid_cell> edges, 
				std::vector<int> &flat_height, const int_2d &groups){
	int loops=0;
	int x,y,nx,ny;
	bool past_edge=false;
	bool do_increment;

	std::vector<grid_cell> new_edges[omp_get_max_threads()];

	while(!edges.empty()){
		loops++;

		#pragma omp parallel num_threads(omp_get_max_threads()) private(x, y, nx, ny, do_increment)
	for(unsigned int i=omp_get_thread_num()*edges.size()/omp_get_num_threads();
		i<edges.size() && i<(omp_get_thread_num()+1)*edges.size()/omp_get_num_threads();
		i++){
			x=edges[i].x;
			y=edges[i].y;

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
					if(!(IN_GRID(nx,ny,elevations.width(),elevations.height()) 
							&& elevations(nx,ny)==elevations(x,y) 
							&& flowdirs(nx,ny)==NO_FLOW)) continue; //TODO: Shouldn't need to check if elevation is equal.
					new_edges[omp_get_thread_num()].push_back(grid_cell(nx,ny));
				}
			}
		}
		past_edge=true;
		edges.clear();
		for(int i=0;i<omp_get_max_threads();i++){
			edges.insert(edges.end(),new_edges[i].begin(),new_edges[i].end());
			new_edges[i].clear();
		}
	}

	return loops;
}

void BarnesStep3(float_2d &elevations, int_2d &inc1, int_2d &inc2, 
					tbb::concurrent_queue<grid_cell> &edge,
					const std::vector<int> &flat_height, 
					const int_2d &groups, float epsilon){
	int x,y,nx,ny;
	grid_cell c(0,0);

	#pragma omp parallel private(x, y, c, nx,ny)
	while(edge.try_pop(c)){
		x=c.x;
		y=c.y;

		if(inc1(x,y)==-1) continue;

		for(int n=1;n<=8;n++){
			nx=x+dx[n];
			ny=y+dy[n];
			if(!(IN_GRID(nx,ny,elevations.width(),elevations.height())
					&& elevations(nx,ny)==elevations(x,y))) continue;
			edge.push(grid_cell(nx,ny));
		}
		if(inc2(x,y)>0){
			elevations(x,y)+= 	(epsilon+epsilon/2)*((inc1(x,y)-1))
								+epsilon*(flat_height[groups(x,y)]-inc2(x,y)+1);
			inc2(x,y)=(inc1(x,y)-1)+(flat_height[groups(x,y)]-inc2(x,y)+1);
		}
		inc1(x,y)=-1;
	}
}

int label_this(int x, int y, const float target_elevation,
					const int label, int_2d &groups, const float_2d &elevations){
	std::queue<grid_cell> to_fill; //TODO: Queue versus stack
	if(elevations(x,y)!=target_elevation || groups(x,y)!=-1) return 0;
	to_fill.push(grid_cell(x,y));

	#pragma omp critical
	while(to_fill.size()>0){
		x=to_fill.front().x;
		y=to_fill.front().y;
		to_fill.pop();
		if(elevations(x,y)!=target_elevation || groups(x,y)!=-1) continue;
		groups(x,y)=label;
		for(int n=1;n<=8;n++)
			if(IN_GRID(x+dx[n],y+dy[n],groups.width(),groups.height()))
				to_fill.push(grid_cell(x+dx[n],y+dy[n]));
	}

	return 1;
}

int find_flat_edges(std::vector<grid_cell> &low_edges, std::vector<grid_cell> &high_edges, 
						const char_2d &flowdirs, const float_2d &elevations, int_2d &groups){
	int nx,ny;
	int group_number=0;

	tbb::concurrent_vector<grid_cell> plow_edges,phigh_edges;

	diagnostic("\r\033[2KSearching for flats...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=1;x<flowdirs.width()-1;x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.height()*100/(flowdirs.width()*flowdirs.height()));
		for(int y=1;y<flowdirs.height()-1;y++)
			for(int n=1;n<=8;n++){
				nx=x+dx[n];
				ny=y+dy[n];

				if(!IN_GRID(nx,ny,flowdirs.width(),flowdirs.height())) continue;
				if(flowdirs(nx,ny)==flowdirs.no_data) continue;

				if(flowdirs(x,y)!=NO_FLOW && flowdirs(nx,ny)==NO_FLOW
							&& elevations(nx,ny)==elevations(x,y) 
							&& INTERIOR_GRID(nx,ny,flowdirs.width(),flowdirs.height())){
					plow_edges.push_back(grid_cell(x,y));
					group_number+=label_this(x, y, elevations(x,y), group_number, groups, elevations);
					break;
				} else if(flowdirs(x,y)==NO_FLOW && elevations(x,y)<elevations(nx,ny)){
					phigh_edges.push_back(grid_cell(x,y));
					group_number+=label_this(x, y, elevations(x,y), group_number, groups, elevations);
					break;
				}	//TODO: You could think about using a GOTO in the code above to move the label_this statements to the same line somewhere
			}
	}
	for(unsigned int i=0;i<plow_edges.size();i++)
		low_edges.push_back(plow_edges[i]);
	for(unsigned int i=0;i<phigh_edges.size();i++)
		high_edges.push_back(phigh_edges[i]);
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	return group_number;
}

int resolve_flats(float_2d &elevations, const char_2d &flowdirs){
	std::vector<grid_cell> llow_edges,lhigh_edges,temp;	//TODO: Need estimate of size
	tbb::concurrent_queue<grid_cell> qlow_edges,qhigh_edges;
	int_2d groups(flowdirs);
	int group_max=-1;

	diagnostic_arg("The groups matrix will require approximately %ldMB of RAM.\n",
						flowdirs.width()*flowdirs.height()*sizeof(int)/1024/1024);
	diagnostic("Resizing groups matrix...");
	groups.resize(flowdirs.width(),flowdirs.height(),false);
	diagnostic("succeeded.\n");
	diagnostic("Initializing groups matrix...");
	groups.init(-1);
	diagnostic("succeeded.\n");

	diagnostic("Entering find_flat_edges function...");
	group_max=find_flat_edges(llow_edges, lhigh_edges, flowdirs, elevations, groups);

	diagnostic_arg("Found %d unique flats.\n",group_max);

	if(group_max==0){
		diagnostic("No flats to resolve!\n");
		return 0;
	}

	diagnostic_arg("The boolean outlets vector will require approximately %ldMB of RAM.\n",
							group_max*sizeof(bool)/1024/1024);
	diagnostic("Creating boolean outlets vector...");
	std::vector<int> has_outlet(group_max,false);
	diagnostic("succeeded!\n");

	diagnostic("Determining which flats have outlets...");
	int outlet_count=0;
	for(std::vector<grid_cell>::iterator i=llow_edges.begin();i!=llow_edges.end();i++){
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
	//The next two for blocks could be run in parallel, but the majority of the time in this algorithm is not spent here.
	for(std::vector<grid_cell>::iterator i=llow_edges.begin();i!=llow_edges.end();i++)	//i++ in 'else'
		if(has_outlet[groups(i->x,i->y)])
			temp.push_back(*i);
	llow_edges.clear();
	llow_edges.insert(llow_edges.end(),temp.begin(),temp.end());
	temp.clear();
	for(std::vector<grid_cell>::iterator i=lhigh_edges.begin();i!=lhigh_edges.end();i++)	//i++ in 'else'
		if(has_outlet[groups(i->x,i->y)])
			temp.push_back(*i);
	llow_edges.clear();
	llow_edges.insert(llow_edges.end(),temp.begin(),temp.end());
	diagnostic("succeeded.\n");

	diagnostic_arg("The incrementation matricies will require approximately %ldMB of RAM.\n",
							3*flowdirs.width()*flowdirs.height()*sizeof(int)/1024/1024);
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

	diagnostic_arg("The flat height vector will require approximately %ldMB of RAM.\n",
							group_max*sizeof(int)/1024/1024);
	diagnostic("Creating flat height vector...");
	std::vector<int> flat_height(group_max);
	diagnostic("succeeded!\n");

	diagnostic("Performing Barnes flat resolution...\n");
	//Flat_height is used with low_edges, but the results will be overwritten by high_edges
	BarnesStep(elevations, flowdirs, inc1, llow_edges, flat_height, groups);
	BarnesStep(elevations, flowdirs, inc2, lhigh_edges, flat_height, groups);

	diagnostic("Combining Barnes flat resolution steps...\n");
	BarnesStep3(elevations, inc1, inc2, qlow_edges, flat_height, groups, 1e-1);
	diagnostic("succeeded!\n");

	return 0;
}
