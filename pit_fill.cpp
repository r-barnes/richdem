#include "utility.h"
#include "data_structures.h"
#include "interface.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <omp.h>
#include <limits>
#include <vector>
#include <queue>

//D8 Directions
static int const dx[9]={0,-1,-1,0,1,1,1,0,-1};
static int const dy[9]={0,0,-1,-1,-1,0,1,1,1};
//Neighbour directions
//234
//105
//876






/*
@article{planchon2002fast,
  title={A fast, simple and versatile algorithm to fill the depressions of digital elevation models},
  author={Planchon, O. and Darboux, F.},
  journal={Catena},
  volume={46},
  number={2-3},
  pages={159--176},
  year={2002}
}
*/
int pit_fill_planchon_direct(float_2d &elevations, float epsilon_straight, float epsilon_diagonal){
	float_2d w;

	float epsilon[9]={0,epsilon_straight, epsilon_diagonal, epsilon_straight, epsilon_diagonal, epsilon_straight, epsilon_straight, epsilon_diagonal, epsilon_straight};

	diagnostic_arg("The intermediate elevation surface 'W' will require approximately %ldMB of RAM.\n",elevations.size1()*elevations.size2()*sizeof(float)/1024/1024);
	diagnostic("Resizing intermediate elevation surface 'W'...");
	try{
		w.resize(elevations.size1(),elevations.size2());
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Initializing intermediate elevation surface 'W'...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.size1();x++){
		progress_bar(x*omp_get_num_threads()*elevations.size2()*100/(elevations.size1()*elevations.size2()));
		for(int y=0;y<elevations.size2();y++)
			if(EDGE_GRID(x,y,elevations.size1(),elevations.size2()))
				w(x,y)=elevations(x,y);
			else
				w(x,y)=std::numeric_limits<float>::infinity();
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	diagnostic("Performing Planchon pit fill Stage 2, direct implementation...\n");
	progress_bar(-1);
	bool something_done=true;
	while(something_done){
		something_done=false;

		#pragma omp parallel for reduction(|:something_done)
		for(int x=1;x<elevations.size1()-1;x++)
		for(int y=1;y<elevations.size2()-1;y++)
			if(w(x,y)>elevations(x,y))
				for(int n=1;n<=8;n++){
					if(elevations(x,y)>=w(x+dx[n],y+dy[n])+epsilon[n]){
						w(x,y)=elevations(x,y);
						something_done=true;
					} else if(w(x,y)>w(x+dx[n],y+dy[n])+epsilon[n]){
						w(x,y)=w(x+dx[n],y+dy[n])+epsilon[n];
						something_done=true;
					}
				}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}











/*
@article{planchon2002fast,
  title={A fast, simple and versatile algorithm to fill the depressions of digital elevation models},
  author={Planchon, O. and Darboux, F.},
  journal={Catena},
  volume={46},
  number={2-3},
  pages={159--176},
  year={2002}
}
*//*
int Dry_upward_cell(int x, int y, int recursive_depth, int max_recursive_depth, float_2d &elevations, float_2d &w, float epsilon[]){ //TODO: This could be implemented as an iterative routine with a queue, but then it could use potentially O(N) memory. This way, the memory usage is limited.
	if(++recursive_depth>max_recursive_depth) return;

	for(int n=1;n<=8;n++){
		if(!IN_GRID(x+dx[n],y+dy[n],elevations.size1(),elevations.size2())) continue;
		if(elevations(x+dx[n],y+dy[n])!=std::numeric_limits<float>::infinity()) continue;
		if(elevations(x+dx[n],y+dy[n])>=w(x,y)+epsilon[n]){
			w(x+dx[n],y+dy[n])=elevations(x+dx[n],y+dy[n]);
			Dry_upward_cell(x+dx[n],y+dy[n],recursive_depth,max_recursive_depth,elevations,w,epsilon);
		}
	}
}

int pit_fill_planchon_optimized(float_2d &elevations, float epsilon_straight, float epsilon_diagonal){
	float_2d w;

	float epsilon[9]={0,epsilon_straight, epsilon_diagonal, epsilon_straight, epsilon_diagonal, epsilon_straight, epsilon_straight, epsilon_diagonal, epsilon_straight};

	diagnostic_arg("The intermediate elevation surface 'W' will require approximately %ldMB of RAM.\n",elevations.size1()*elevations.size2()*sizeof(float)/1024/1024);
	diagnostic("Resizing intermediate elevation surface 'W'...");
	try{
		w.resize(elevations.size1(),elevations.size2());
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Initializing intermediate elevation surface 'W'...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.size1();x++){
		progress_bar(x*omp_get_num_threads()*elevations.size2()*100/(elevations.size1()*elevations.size2()));
		for(int y=0;y<elevations.size2();y++)
			if(EDGE_GRID(x,y,elevations.size1(),elevations.size2()))
				w(x,y)=elevations(x,y);
			else
				w(x,y)=std::numeric_limits<float>::infinity();
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	diagnostic("Performing Planchon pit fill Stage 2, direct implementation...\n");
	progress_bar(-1);
	bool something_done=true;
	while(something_done){
		something_done=false;

		#pragma omp parallel for reduction(|:something_done)
		for(int x=1;x<elevations.size1()-1;x++)
		for(int y=1;y<elevations.size2()-1;y++)
			if(w(x,y)>elevations(x,y))
				for(int n=1;n<=8;n++){
					if(elevations(x,y)>=w(x+dx[n],y+dy[n])+epsilon[n]){
						w(x,y)=elevations(x,y);
						something_done=true;
					} else if(w(x,y)>w(x+dx[n],y+dy[n])+epsilon[n]){
						w(x,y)=w(x+dx[n],y+dy[n])+epsilon[n];
						something_done=true;
					}
				}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}
*/











/*
@article{wang2006efficient,
  title={An efficient method for identifying and filling surface depressions in digital elevation models for hydrologic analysis and modelling},
  author={Wang, L. and Liu, H.},
  journal={International Journal of Geographical Information Science},
  volume={20},
  number={2},
  pages={193--213},
  year={2006},
  publisher={Taylor \& Francis}
}

This algorithm is essentially reproduced (but in poorer quality) by
@article{liu2009another,
  title={Another fast and simple DEM depression-filling algo-rithm based on priority queue structure},
  author={LIU, Y.H. and ZHANG, W.C. and XU, J.W.},
  journal={Atmos. Oceanic Sci. Lett},
  volume={2},
  pages={214--219},
  year={2009}
}

This algorithm is first, but also poorly-described, though the pseudo-code appears to work
@article{soille1994efficient,
  title={An efficient algorithm for drainage network extraction on DEMs},
  author={Soille, P. and Gratin, C.},
  journal={Journal of Visual Communication and Image Representation},
  volume={5},
  number={2},
  pages={181--189},
  year={1994},
  publisher={Elsevier}
}
*/

typedef struct grid_cell_type {
	int x;
	int y;
	float z;
	grid_cell_type(int x0, int y0, float z0){
		x=x0;
		y=y0;
		z=z0;
	}
} grid_cell;

class grid_cell_compare{
	bool reverse;

	public:
		grid_cell_compare(const bool& revparam=false){reverse=revparam;}

		bool operator() (const grid_cell *lhs, const grid_cell *rhs) const{
			if (reverse) return (lhs->z<rhs->z);
			else return (lhs->z>rhs->z);
		}
};

int pit_fill_wang(float_2d &elevations){
	std::priority_queue<grid_cell*, std::vector<grid_cell*>, grid_cell_compare> open;
	grid_cell *c;
	char_2d closed;	//TODO: This could probably be made into a boolean
	unsigned long processed_cells=0;

	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.size1()*elevations.size2()*sizeof(char)/1024/1024);
	diagnostic("Resizing boolean flood array matrix...");
	try{
		closed.resize(elevations.size1(),elevations.size2());
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");
	diagnostic("Initializing closed matrix...");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<closed.size1();x++){
		progress_bar(x*omp_get_num_threads()*elevations.size2()*100/(elevations.size1()*elevations.size2()));
		for(int y=0;y<closed.size2();y++)
			closed(x,y)=0;
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.size1()*2+elevations.size2()*2)*sizeof(grid_cell)/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.size1();x++){
		open.push(new grid_cell(x,0,elevations(x,0) ));
		open.push(new grid_cell(x,elevations.size2()-1,elevations(x,elevations.size2()-1) ));
		closed(x,0)=1;
		closed(x,elevations.size2()-1)=1;
	}
	for(int y=1;y<elevations.size2()-1;y++){
		open.push(new grid_cell(0,y,elevations(0,y)	));
		open.push(new grid_cell(elevations.size1()-1,y,elevations(elevations.size1()-1,y) ));
		closed(0,y)=1;
		closed(elevations.size1()-1,y)=1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Performing the Wang fill...\n");
	progress_bar(-1);
	while(open.size()>0){
		c=open.top();
		open.pop();
		closed(c->x,c->y)=2;

		for(int n=1;n<=8;n++){
			if(!IN_GRID(c->x+dx[n],c->y+dy[n],elevations.size1(),elevations.size2())) continue;
			if(closed(c->x+dx[n],c->y+dy[n])>0)
				continue;
			else{
				elevations(c->x+dx[n],c->y + dy[n])=MAX(elevations(c->x+dx[n],c->y+dy[n]),c->z);
				open.push(new grid_cell(c->x+dx[n],c->y+dy[n],elevations(c->x+dx[n],c->y+dy[n]) ));
				closed(c->x+dx[n],c->y+dy[n])=1;
			}
		}
		progress_bar(processed_cells*100/(elevations.size1()*elevations.size2()));
		processed_cells++;
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

//	print_dem(elevations);
}
