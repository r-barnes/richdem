#include "utility.h"
#include "data_structures.h"
#include "interface.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <omp.h>
#include <limits>
#include <vector>
#include <queue>


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
void pit_fill_planchon_direct(float_2d &elevations, float epsilon_straight, float epsilon_diagonal){
	float_2d w;

	float epsilon[9]={0,epsilon_straight, epsilon_diagonal, epsilon_straight, epsilon_diagonal, epsilon_straight, epsilon_straight, epsilon_diagonal, epsilon_straight};

	diagnostic_arg("The intermediate elevation surface 'W' will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(float)/1024/1024);
	diagnostic("Resizing intermediate elevation surface 'W'...");
	w.resize(elevations.width(),elevations.height());
	diagnostic("succeeded.\n");

	diagnostic("Initializing intermediate elevation surface 'W'...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++){
		progress_bar(x*omp_get_num_threads()*elevations.height()*100/(elevations.width()*elevations.height()));
		for(int y=0;y<elevations.height();y++)
			if(EDGE_GRID(x,y,elevations.width(),elevations.height()))
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
		for(int x=1;x<elevations.width()-1;x++)
		for(int y=1;y<elevations.height()-1;y++)
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
		if(!IN_GRID(x+dx[n],y+dy[n],elevations.width(),elevations.height())) continue;
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

	diagnostic_arg("The intermediate elevation surface 'W' will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(float)/1024/1024);
	diagnostic("Resizing intermediate elevation surface 'W'...");
	w.resize(elevations.width(),elevations.height());
	diagnostic("succeeded.\n");

	diagnostic("Initializing intermediate elevation surface 'W'...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++){
		progress_bar(x*omp_get_num_threads()*elevations.height()*100/(elevations.width()*elevations.height()));
		for(int y=0;y<elevations.height();y++)
			if(EDGE_GRID(x,y,elevations.width(),elevations.height()))
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
		for(int x=1;x<elevations.width()-1;x++)
		for(int y=1;y<elevations.height()-1;y++)
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
31 citations. None appear to be improvements.
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
0 citations listed on Google Scholar.
@article{liu2009another,
  title={Another fast and simple DEM depression-filling algo-rithm based on priority queue structure},
  author={LIU, Y.H. and ZHANG, W.C. and XU, J.W.},
  journal={Atmos. Oceanic Sci. Lett},
  volume={2},
  pages={214--219},
  year={2009}
}

This algorithm is first, but also poorly-described, though the pseudo-code appears to work
44 citations. None appear to be improvements.
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

class grid_cell_compare{
	bool reverse;

	public:
		grid_cell_compare(const bool& revparam=false){reverse=revparam;}

		bool operator() (const grid_cellz *lhs, const grid_cellz *rhs) const{
			if (reverse) return (lhs->z<rhs->z);
			else return (lhs->z>rhs->z);
		}
};

class grid_cell_compare2{
	bool reverse;

	public:
		grid_cell_compare2(const bool& revparam=false){reverse=revparam;}

		bool operator() (const grid_cellz &lhs, const grid_cellz &rhs) const{
			if (reverse) return (lhs.z<rhs.z);
			else return (lhs.z>rhs.z);
		}
};





void pit_fill_wang(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare2> open;
	char_2d closed;	//TODO: This could probably be made into a boolean
	unsigned long processed_cells=0;

	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Resizing boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	diagnostic("succeeded.\n");
	diagnostic("Initializing closed matrix...\n");
	closed.init(0);
	diagnostic("\tsucceeded.\n");

	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*sizeof(grid_cell)/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.width();x++){
		open.push(grid_cellz(x,0,elevations(x,0) ));
		open.push(grid_cellz(x,elevations.height()-1,elevations(x,elevations.height()-1) ));
		closed(x,0)=1;
		closed(x,elevations.height()-1)=1;
	}
	for(int y=1;y<elevations.height()-1;y++){
		open.push(grid_cellz(0,y,elevations(0,y)	));
		open.push(grid_cellz(elevations.width()-1,y,elevations(elevations.width()-1,y) ));
		closed(0,y)=1;
		closed(elevations.width()-1,y)=1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Performing the Wang fill...\n");
	progress_bar(-1);
	while(open.size()>0){
		grid_cellz c=open.top();
		closed(c.x,c.y)=2;
		processed_cells++;

		for(int n=1;n<=8;n++){
			int nx=c.x+dx[n];
			int ny=c.y+dy[n];
			if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
			if(closed(nx,ny)>0)
				continue;
			else{
				elevations(nx,ny)=MAX(elevations(nx,ny),c.z);
				open.push(grid_cellz(nx,ny,elevations(nx,ny) ));
				closed(nx,ny)=1;
			}
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
		open.pop();
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}









void pit_fill_barnes1(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare2> open;
	std::queue<grid_cell> meander;
	char_2d closed;	//TODO: This could probably be made into a boolean
	unsigned long processed_cells=0;

	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Resizing boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	diagnostic("succeeded.\n");
	diagnostic("Initializing closed matrix...\n");
	closed.init(0);
	diagnostic("\tsucceeded.\n");

	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*sizeof(grid_cell)/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.width();x++){
		open.push(grid_cellz(x,0,elevations(x,0) ));
		open.push(grid_cellz(x,elevations.height()-1,elevations(x,elevations.height()-1) ));
		closed(x,0)=1;
		closed(x,elevations.height()-1)=1;
	}
	for(int y=1;y<elevations.height()-1;y++){
		open.push(grid_cellz(0,y,elevations(0,y)	));
		open.push(grid_cellz(elevations.width()-1,y,elevations(elevations.width()-1,y) ));
		closed(0,y)=1;
		closed(elevations.width()-1,y)=1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Performing the Barnes fill...\n");
	progress_bar(-1);
	while(open.size()>0 || meander.size()>0){
		if(meander.size()>0){
			grid_cell c=meander.front();
			closed(c.x,c.y)=2;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)>0)
					continue;
				else if(elevations(nx,ny)<=elevations(c.x,c.y)){
					elevations(nx,ny)=elevations(c.x,c.y);
					meander.push(grid_cell(nx,ny));
					closed(nx,ny)=1;
					continue;
				} else {
					open.push(grid_cellz(nx,ny,elevations(nx,ny)));
					closed(nx,ny)=1;
					continue;
				}
			}
			progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
			meander.pop();
		} else {
			grid_cellz c=open.top();
			closed(c.x,c.y)=2;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)>0)
					continue;
				else if(elevations(nx,ny)>elevations(c.x,c.y)){
					open.push(grid_cellz(nx,ny,elevations(nx,ny)));
					closed(nx,ny)=1;
					continue;
				} else if(elevations(nx,ny)<elevations(c.x,c.y))
					elevations(nx,ny)=elevations(c.x,c.y);
				meander.push(grid_cell(nx,ny));
				closed(nx,ny)=1;
			}
			progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
			open.pop();
		}
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}









void pit_fill_barnes2(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare2> open;
	std::queue<grid_cell> meander;
	std::queue<grid_cell> climb;
	char_2d closed;	//TODO: This could probably be made into a boolean
	unsigned long processed_cells=0;

	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Resizing boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	diagnostic("succeeded.\n");
	diagnostic("Initializing closed matrix...\n");
	closed.init(0);
	diagnostic("\tsucceeded.\n");

	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*sizeof(grid_cell)/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.width();x++){
		open.push(grid_cellz(x,0,elevations(x,0) ));
		open.push(grid_cellz(x,elevations.height()-1,elevations(x,elevations.height()-1) ));
		closed(x,0)=1;
		closed(x,elevations.height()-1)=1;
	}
	for(int y=1;y<elevations.height()-1;y++){
		open.push(grid_cellz(0,y,elevations(0,y)	));
		open.push(grid_cellz(elevations.width()-1,y,elevations(elevations.width()-1,y) ));
		closed(0,y)=1;
		closed(elevations.width()-1,y)=1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Performing the Barnes fill...\n");
	progress_bar(-1);
	while(open.size()>0 || meander.size()>0 || climb.size()>0){
		if(meander.size()>0){
			grid_cell c=meander.front();
			closed(c.x,c.y)=2;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)>0)
					continue;
				else if(elevations(nx,ny)<=elevations(c.x,c.y)){
					elevations(nx,ny)=elevations(c.x,c.y);
					meander.push(grid_cell(nx,ny));
					closed(nx,ny)=1;
					continue;
				} else {
					climb.push(grid_cell(nx,ny));
					closed(nx,ny)=1;
					continue;
				}
			}
			meander.pop();
		} else if (climb.size()>0){
			grid_cell c=climb.front();
			closed(c.x,c.y)=2;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)>0)
					continue;
				else if(elevations(nx,ny)<elevations(c.x,c.y)){
					open.push(grid_cellz(nx,ny,elevations(c.x,c.y)));
					continue;
				} else {
					climb.push(grid_cell(nx,ny));
					closed(nx,ny)=1;
					continue;
				}
			}
			climb.pop();
		} else {
			grid_cellz c=open.top();

			if(closed(c.x,c.y)==2){
				open.pop();
				continue;
			}

			closed(c.x,c.y)=2;
			processed_cells++;

			elevations(c.x,c.y)=c.z;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)>0)
					continue;
				else if(elevations(nx,ny)>elevations(c.x,c.y)){
					climb.push(grid_cell(nx,ny));
					closed(nx,ny)=1;
					continue;
				} else if(elevations(nx,ny)<=elevations(c.x,c.y)){
					elevations(nx,ny)=elevations(c.x,c.y);
					meander.push(grid_cell(nx,ny));
					closed(nx,ny)=1;
					continue;
				}
			}
			open.pop();
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}










void pit_fill_barnes3(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare2> open;
	std::queue<grid_cell> meander;
	std::queue<grid_cell> climb;
	char_2d closed;	//TODO: This could probably be made into a boolean
	float_2d info;
	unsigned long processed_cells=0;
	const int CLOSED=1, QUEUED=2, OPEN=3;

	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Resizing boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	info.resize(elevations.width(),elevations.height());
	diagnostic("succeeded.\n");
	diagnostic("Initializing closed matrix...\n");
	closed.init(OPEN);
	info.init(999999);
	diagnostic("\tsucceeded.\n");

	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*sizeof(grid_cell)/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.width();x++){
		open.push(grid_cellz(x,0,elevations(x,0) ));
		open.push(grid_cellz(x,elevations.height()-1,elevations(x,elevations.height()-1) ));
		closed(x,0)=QUEUED;
		closed(x,elevations.height()-1)=QUEUED;
	}
	for(int y=1;y<elevations.height()-1;y++){
		open.push(grid_cellz(0,y,elevations(0,y)	));
		open.push(grid_cellz(elevations.width()-1,y,elevations(elevations.width()-1,y) ));
		closed(0,y)=QUEUED;
		closed(elevations.width()-1,y)=QUEUED;
	}
	diagnostic("succeeded.\n");

	diagnostic("Performing the Barnes fill...\n");
	progress_bar(-1);
	while(open.size()>0 || meander.size()>0 || climb.size()>0){
		if(meander.size()>0){
			grid_cell c=meander.front();
			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)<=QUEUED)
					continue;
				else if(elevations(nx,ny)<=elevations(c.x,c.y)){
					elevations(nx,ny)=elevations(c.x,c.y);
					meander.push(grid_cell(nx,ny));
					closed(nx,ny)=QUEUED;
					continue;
				} else {
					climb.push(grid_cell(nx,ny));
					closed(nx,ny)=QUEUED;
					continue;
				}
			}
			meander.pop();
		} else if (climb.size()>0){
			grid_cell c=climb.front();
			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)<=QUEUED || elevations(c.x,c.y)>=info(nx,ny))
					continue;
				else if(elevations(nx,ny)<elevations(c.x,c.y)){
					open.push(grid_cellz(nx,ny,elevations(c.x,c.y)));
					info(nx,ny)=elevations(c.x,c.y);
					continue;
				} else {
					climb.push(grid_cell(nx,ny));
					closed(nx,ny)=QUEUED;
					continue;
				}
			}
			climb.pop();
		} else {
			grid_cellz c=open.top();

			if(closed(c.x,c.y)==CLOSED){
				open.pop();
				continue;
			}

			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			elevations(c.x,c.y)=c.z;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)<=QUEUED)
					continue;
				else if(elevations(nx,ny)>elevations(c.x,c.y)){
					climb.push(grid_cell(nx,ny));
					closed(nx,ny)=QUEUED;
					continue;
				} else if(elevations(nx,ny)<=elevations(c.x,c.y)){
					elevations(nx,ny)=elevations(c.x,c.y);
					meander.push(grid_cell(nx,ny));
					closed(nx,ny)=QUEUED;
					continue;
				}
			}
			open.pop();
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}
