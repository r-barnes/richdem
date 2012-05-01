#include "pit_fill.h"
#include "utility.h"
#include "data_structures.h"
#include "interface.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <omp.h>
#include <limits>
#include <vector>
#include <queue>
#include <limits>
#include <stack>
#include <algorithm>







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

		bool operator() (const grid_cellz &lhs, const grid_cellz &rhs) const{
			if (reverse) return (lhs.z<rhs.z);
			else return (lhs.z>rhs.z);
		}
};

void push_edges(const float_2d &elevations, std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> &open, char_2d &closed, int closed_val=1){
	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*sizeof(grid_cellz)/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.width();x++){
		open.push(grid_cellz(x,0,elevations(x,0) ));
		open.push(grid_cellz(x,elevations.height()-1,elevations(x,elevations.height()-1) ));
		closed(x,0)=closed_val;
		closed(x,elevations.height()-1)=closed_val;
	}
	for(int y=1;y<elevations.height()-1;y++){
		open.push(grid_cellz(0,y,elevations(0,y)	));
		open.push(grid_cellz(elevations.width()-1,y,elevations(elevations.width()-1,y) ));
		closed(0,y)=closed_val;
		closed(elevations.width()-1,y)=closed_val;
	}
	diagnostic("succeeded.\n");
}


//Wang runs everything through a priority queue. This the standard against which to test everything. We will make improvements by reducing the number of items going through the PQ.
void pit_fill_wang(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	char_2d closed;	//TODO: This could probably be made into a boolean
	unsigned long processed_cells=0;

	diagnostic("\n###Wang (2006) Pit Fill\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Resizing boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	diagnostic("succeeded.\n");
	diagnostic("Initializing closed matrix...");
	closed.init(0);
	diagnostic("succeeded.\n");

	push_edges(elevations, open, closed);

	diagnostic("Performing the Wang fill...\n");
	progress_bar(-1);
	while(open.size()>0){
		grid_cellz c=open.top();open.pop();
		closed(c.x,c.y)=2;
		processed_cells++;

		for(int n=1;n<=8;n++){
			int nx=c.x+dx[n];
			int ny=c.y+dy[n];
			if(!elevations.in_grid(nx,ny)) continue;
			if( !(closed(nx,ny)>0) ){
				elevations(nx,ny)=MAX(elevations(nx,ny),c.z);
				open.push(grid_cellz(nx,ny,elevations(nx,ny) ));
				closed(nx,ny)=1;
			}
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);	//TODO
}




/*
//WangND explores the possibility that the majority of the speed gains the Barnes variants make over Wang stem from efficient processing of the NoData cells. Therefore, it runs the NoData cells through a standard queue and everything else through a PQ.
void pit_fill_wangND(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	char_2d closed;	//TODO: This could probably be made into a boolean
	unsigned long processed_cells=0;

	diagnostic("\n###Wang (2006) Pit Fill (Richard's ND version)\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(0);
	diagnostic("succeeded.\n");

	push_edges(elevations, open, closed);

	std::queue<grid_cell> nd_cells;

	diagnostic("Performing the WangND fill...\n");
	progress_bar(-1);
	while(nd_cells.size()>0 || open.size()>0){
		if(nd_cells.size()>0){
			grid_cell c=nd_cells.front();nd_cells.pop();
			closed(c.x,c.y)=2;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if( !(closed(nx,ny)>0) ){
					if(elevations(nx,ny)==elevations.no_data)
						nd_cells.push(grid_cell(nx,ny));
					else
						open.push(grid_cellz(nx,ny,elevations(nx,ny) ));
					closed(nx,ny)=1;
				}
			}

		} else {
			grid_cellz c=open.top();open.pop();
			closed(c.x,c.y)=2;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if( !(closed(nx,ny)>0) ){
					if(elevations(nx,ny)==elevations.no_data)
						nd_cells.push(grid_cell(nx,ny));
					else {
						elevations(nx,ny)=MAX(elevations(nx,ny),c.z);
						open.push(grid_cellz(nx,ny,elevations(nx,ny) ));
					}
					closed(nx,ny)=1;
				}
			}
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);	//TODO
}

*/




bool sort_cells (const grid_cellz &lhs,const grid_cellz &rhs) { return (lhs.z>rhs.z); }

void pit_fill_barneslehman(float_2d &elevations){
	float max,min;
	bool_2d closed;
	max=elevations.max();
	min=elevations.min();

	std::vector<std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> > open((int)(10*max), std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare>() );
	std::queue<grid_cellz> pit;
	int lvl=((int)(10*min));
	int processed_cells=0;
	const int tcells=elevations.width()*elevations.height();

	diagnostic("\n###Barnes-Lehman Pit Fill v1\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(false);
	diagnostic("succeeded.\n");

	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*sizeof(grid_cellz)/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.width();x++){
		if(elevations(x,0)==elevations.no_data)
			pit.push(grid_cellz(x,0,elevations.no_data));
		else{
			open[(int)(elevations(x,0)*10)].push(grid_cellz(x,0,elevations(x,0) ));
		}
		if(elevations(x,elevations.height()-1)==elevations.no_data)
			pit.push(grid_cellz(x,elevations.height()-1,elevations.no_data));
		else{
			open[(int)(elevations(x,elevations.height()-1)*10)].push(grid_cellz(x,elevations.height()-1,elevations(x,elevations.height()-1) ));
		}
		closed(x,0)=true;
		closed(x,elevations.height()-1)=true;
	}
	for(int y=1;y<elevations.height()-1;y++){
		if(elevations(0,y)==elevations.no_data)
			pit.push(grid_cellz(0,y,elevations.no_data));
		else{
			open[(int)(elevations(0,y)*10)].push(grid_cellz(0,y,elevations(0,y) ));
		}
		if(elevations(elevations.width()-1,y)==elevations.no_data)
			pit.push(grid_cellz(elevations.width()-1,y,elevations.no_data));
		else{
			open[(int)(elevations(elevations.width()-1,y)*10)].push(grid_cellz(elevations.width()-1,y,elevations(elevations.width()-1,y) ));
		}
		closed(0,y)=true;
		closed(elevations.width()-1,y)=true;
	}
	diagnostic("succeeded.\n");

	diagnostic("Performing the Barnes1 fill...\n");
	progress_bar(-1);
	while(pit.size()>0 || processed_cells<tcells){
		if(pit.size()>0){
			grid_cellz c=pit.front();
			pit.pop();
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)) continue;
				closed(nx,ny)=true;
				if(elevations(nx,ny)<=c.z){
					elevations(nx,ny)=c.z;
					pit.push(grid_cellz(nx,ny,c.z));
				} else
					open[(int)(10*elevations(nx,ny))].push(grid_cellz(nx,ny,elevations(nx,ny)));
			}
		} else {
			while(open[lvl].size()==0)	//Find next bin
				lvl++;

			grid_cellz c=open[lvl].top();
			open[lvl].pop();
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)) continue;
				closed(nx,ny)=true;
				if(elevations(nx,ny)<=c.z){
					elevations(nx,ny)=c.z;
					pit.push(grid_cellz(nx,ny,c.z));
				} else
					open[(int)(10*elevations(nx,ny))].push(grid_cellz(nx,ny,elevations(nx,ny)));
			}
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs\033[39m\n",progress_bar(-1));
}










//Barnes1 explores depressions and flats by pushing them onto the meander queue. When meander encounters cells at <= elevation, it meanders over them. Otherwise, it pushes them onto the open queue.
void pit_fill_barnes1(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	std::queue<grid_cellz> meander;
//	std::stack<grid_cellz, std::vector<grid_cellz> > meander;
	bool_2d closed;	//TODO: This could probably be made into a boolean
	unsigned long processed_cells=0;
	unsigned long pitc=0,openc=0;

	diagnostic("\n###Barnes Pit Fill v1\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(false);
	diagnostic("succeeded.\n");

	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*sizeof(grid_cellz)/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.width();x++){
		open.push(grid_cellz(x,0,elevations(x,0) ));
		open.push(grid_cellz(x,elevations.height()-1,elevations(x,elevations.height()-1) ));
		closed(x,0)=true;
		closed(x,elevations.height()-1)=true;
	}
	for(int y=1;y<elevations.height()-1;y++){
		open.push(grid_cellz(0,y,elevations(0,y)	));
		open.push(grid_cellz(elevations.width()-1,y,elevations(elevations.width()-1,y) ));
		closed(0,y)=true;
		closed(elevations.width()-1,y)=true;
	}
	diagnostic("succeeded.\n");

	diagnostic("Performing the Barnes1 fill...\n");
	progress_bar(-1);
	while(open.size()>0 || meander.size()>0){
		grid_cellz c;
		if(meander.size()>0){
			c=meander.front();
			meander.pop();
			pitc++;
		} else {
			c=open.top();
			open.pop();
			openc++;
		}
		processed_cells++;

		for(int n=1;n<=8;n++){
			int nx=c.x+dx[n];
			int ny=c.y+dy[n];
			if(!elevations.in_grid(nx,ny)) continue;
			if(closed(nx,ny)) 
				continue;

			closed(nx,ny)=true;
			if(elevations(nx,ny)<=c.z){
				elevations(nx,ny)=c.z;
				meander.push(grid_cellz(nx,ny,c.z));
			} else
				open.push(grid_cellz(nx,ny,elevations(nx,ny)));
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
	printf("Pit=%ld, Open=%ld\n",pitc,openc);
}




/*


//Barnes1b works in the same way as Barnes1, but does a flood fill over cells which are <= the elevation at which a pit was entered. Barnes1 does a flood fill as well, but with slightly different checks.
//Tests on Steele County showed that Barnes1 took 155.463474s, while Barnes1b took 152.906847s.
void flood_barnes1b(int x, int y, const float_2d &elevations, char_2d &closed, std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> &open, unsigned long &processed_cells){
	std::queue<grid_cell> to_flood;
	const float elevation=elevations(x,y);
	to_flood.push(grid_cell(x,y));

	while(to_flood.size()>0){
		grid_cell c=to_flood.front();to_flood.pop();
		closed(c.x,c.y)=2;
		processed_cells++;
		for(int n=1;n<=8;n++){
			int nx=c.x+dx[n];
			int ny=c.y+dy[n];
			if(!elevations.in_grid(nx,ny)) continue;
			if(closed(nx,ny)>0) continue;

			closed(nx,ny)=1;
			if(elevations(nx,ny)>elevation)
				open.push(grid_cellz(nx,ny,elevations(nx,ny)));
			else
				to_flood.push(grid_cell(nx,ny));
		}
	}
}

void pit_fill_barnes1b(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	std::queue<grid_cell> meander;
	char_2d closed;	//TODO: This could probably be made into a boolean
	unsigned long processed_cells=0;

	diagnostic("\n###Barnes Pit Fill v1b\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(0);
	diagnostic("succeeded.\n");

	push_edges(elevations, open, closed);

	diagnostic("Performing the Barnes1b fill...\n");
	progress_bar(-1);
	while(open.size()>0){
		grid_cellz c=open.top();open.pop();
		if(closed(c.x,c.y)==2) continue;
		closed(c.x,c.y)=2;
		processed_cells++;

		for(int n=1;n<=8;n++){
			int nx=c.x+dx[n];
			int ny=c.y+dy[n];
			if(!elevations.in_grid(nx,ny)) continue;
			if(closed(nx,ny)>0)
				continue;
			else if(elevations(nx,ny)>elevations(c.x,c.y)){
				open.push(grid_cellz(nx,ny,elevations(nx,ny)));
				closed(nx,ny)=1;
				continue;
			} else
				flood_barnes1b(nx,ny,elevations,closed,open,processed_cells);
				continue;
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}








//Barnes2 uses a meander queue in the same way as in Barnes1, but introduces a climb queue. The climb queue can only ascend from cells which have a known outlet (i.e. those popped off the PQ or Meander queue). If a cell is lower than a cell popped off the climb queue, then that cell must be added to the PQ. This means that a cell may be placed onto the PQ multiple times.
void pit_fill_barnes2(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	std::queue<grid_cellz> meander;
	std::queue<grid_cellz> climb;
	char_2d closed;	//TODO: This could probably be made into a boolean
	unsigned long processed_cells=0;

	diagnostic("\n###Barnes Pit Fill v2\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(0);
	diagnostic("succeeded.\n");

	push_edges(elevations, open, closed);

	diagnostic("Performing the Barnes2 fill...\n");
	progress_bar(-1);
	while(open.size()>0 || meander.size()>0 || climb.size()>0){
		if(meander.size()>0){
			grid_cellz c=meander.front();meander.pop();
			closed(c.x,c.y)=2;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)>0) continue;

				closed(nx,ny)=1;
				if(elevations(nx,ny)<=c.z){
					elevations(nx,ny)=c.z;
					meander.push(grid_cellz(nx,ny,c.z));
				} else
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
			}
		} else if (climb.size()>0){
			grid_cellz c=climb.front();climb.pop();
			closed(c.x,c.y)=2;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)>0) continue;

				if(elevations(nx,ny)<c.z)
					open.push(grid_cellz(nx,ny,c.z));
				else {
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
					closed(nx,ny)=1;
				}
			}
		} else {
			grid_cellz c=open.top();open.pop();

			if(closed(c.x,c.y)==2) continue;

			closed(c.x,c.y)=2;
			processed_cells++;

			elevations(c.x,c.y)=c.z;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)>0) continue;
				closed(nx,ny)=1;
				if(elevations(nx,ny)>elevations(c.x,c.y))
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
				else {
					elevations(nx,ny)=c.z;
					meander.push(grid_cellz(nx,ny,c.z));
				}
			}
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}









//Barnes3 uses a meander queue in the same way as in Barnes1, but introduces a climb queue. The climb queue can only ascend from cells which have a known outlet (i.e. those popped off the PQ or Meander queue). If a cell is lower than a cell popped off the climb queue, then that cell must be added to the PQ. In Barnes2, this meant that a cell could be placed onto the PQ multiple times. In Barnes3, this is reduced through the use of the info array which has the minimum encountered elevation a cell can be. A popped climb cell will only place a cell into the PQ if its estimated elevation is lower than that already in the info array.
void pit_fill_barnes3(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	std::queue<grid_cellz> meander;
	std::queue<grid_cellz> climb;
	char_2d closed;	//TODO: This could probably be made into a boolean
	float_2d info;
	unsigned long processed_cells=0;
	const int CLOSED=1, QUEUED=2, OPEN=3;

	diagnostic("\n###Barnes Pit Fill v3\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(OPEN);
	diagnostic("succeeded.\n");
	info.resize(elevations.width(),elevations.height());
	info.init(9e12);

	push_edges(elevations, open, closed, QUEUED);

	diagnostic("Performing the Barnes3 fill...\n");
	progress_bar(-1);
	while(open.size()>0 || meander.size()>0 || climb.size()>0){
		if(meander.size()>0){
			grid_cellz c=meander.front();meander.pop();
			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)!=OPEN) continue;

				closed(nx,ny)=QUEUED;
				if(elevations(nx,ny)<=c.z){
					elevations(nx,ny)=c.z;
					meander.push(grid_cellz(nx,ny,c.z));
				} else
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
			}
		} else if (climb.size()>0){
			grid_cellz c=climb.front();climb.pop();
			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)!=OPEN || c.z>=info(nx,ny))
					continue;
				else if(elevations(nx,ny)<c.z){
					open.push(grid_cellz(nx,ny,c.z));
					info(nx,ny)=elevations(c.x,c.y);
				} else {
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
					closed(nx,ny)=QUEUED;
				}
			}
		} else {
			grid_cellz c=open.top();open.pop();

			if(closed(c.x,c.y)==CLOSED) continue;

			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			elevations(c.x,c.y)=c.z; //TODO?

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)!=OPEN) continue;

				closed(nx,ny)=QUEUED;
				if(elevations(nx,ny)>c.z)
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
				else if(elevations(nx,ny)<=c.z){
					elevations(nx,ny)=c.z;
					meander.push(grid_cellz(nx,ny,c.z));
				}
			}
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}
















//Barnes4 works in the same way as Barnes3, but knows which QUEUE a cell has been added to - the meander queue takes precedence.
void pit_fill_barnes4(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	std::queue<grid_cellz> meander;
	std::queue<grid_cellz> climb;
	char_2d closed;	//TODO: This could probably be made into a boolean
	float_2d info;
	unsigned long processed_cells=0;
	const int CLOSED=1, PITQUEUED=2, CLIMBQUEUED=3, OPENQUEUED=4, OPEN=5;

	diagnostic("\n###Barnes Pit Fill v4\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(OPEN);
	diagnostic("succeeded.\n");
	info.resize(elevations.width(),elevations.height());
	info.init(9e12);

	push_edges(elevations, open, closed, OPENQUEUED);

	diagnostic("Performing the Barnes4 fill...\n");
	progress_bar(-1);
	while(open.size()>0 || meander.size()>0 || climb.size()>0){
		if(meander.size()>0){
			grid_cellz c=meander.front();meander.pop();
			if(closed(c.x,c.y)==CLOSED) continue;
			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)<PITQUEUED) continue;

				if(elevations(nx,ny)<=c.z){
					elevations(nx,ny)=c.z;
					meander.push(grid_cellz(nx,ny,c.z));
					closed(nx,ny)=PITQUEUED;
				} else{
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
					closed(nx,ny)=CLIMBQUEUED;
				}
			}
		} else if (climb.size()>0){
			grid_cellz c=climb.front();climb.pop();
			if(closed(c.x,c.y)==CLOSED) continue;
			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)<CLIMBQUEUED || c.z>=info(nx,ny))
					continue;
				else if(elevations(nx,ny)<c.z){
					open.push(grid_cellz(nx,ny,c.z));
					info(nx,ny)=elevations(c.x,c.y);
					closed(nx,ny)=OPENQUEUED;
				} else {
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
					closed(nx,ny)=CLIMBQUEUED;
				}
			}
		} else {
			grid_cellz c=open.top();open.pop();

			if(closed(c.x,c.y)==CLOSED) continue;

			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			elevations(c.x,c.y)=c.z; //TODO?

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)<OPENQUEUED) continue;

				if(elevations(nx,ny)>c.z){
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
					closed(nx,ny)=CLIMBQUEUED;
				} else if(elevations(nx,ny)<=c.z){
					elevations(nx,ny)=c.z;
					meander.push(grid_cellz(nx,ny,c.z));
					closed(nx,ny)=PITQUEUED;
				}
			}
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}










//Barnes5 works in the same way as Barnes3, but attempts to climb better
void pit_fill_barnes5(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	std::queue<grid_cellz> meander;
	std::queue<grid_cellz> climb;
	char_2d closed;	//TODO: This could probably be made into a boolean
	float_2d info;
	unsigned long processed_cells=0;
	const int CLOSED=1, QUEUED=2, OPEN=3;

	diagnostic("\n###Barnes Pit Fill v5\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(OPEN);
	diagnostic("succeeded.\n");
	info.resize(elevations.width(),elevations.height());
	info=elevations;

	push_edges(elevations, open, closed, QUEUED);

	diagnostic("Performing the Barnes5 fill...\n");
	progress_bar(-1);
	while(open.size()>0 || meander.size()>0 || climb.size()>0){
		if(meander.size()>0){
			grid_cellz c=meander.front();meander.pop();
			if(closed(c.x,c.y)==CLOSED) continue;
			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)!=OPEN) continue;

				closed(nx,ny)=QUEUED;
				if(elevations(nx,ny)<=c.z){
					elevations(nx,ny)=c.z;
					meander.push(grid_cellz(nx,ny,c.z));
				} else
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
			}
		} else if (climb.size()>0){
			grid_cellz c=climb.front();climb.pop();
			if(closed(c.x,c.y)==CLOSED) continue;
			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)==CLOSED) continue;
				if(closed(nx,ny)==QUEUED && c.z<=info(nx,ny)) continue;

				if(elevations(nx,ny)<c.z){
					open.push(grid_cellz(nx,ny,c.z));
					info(nx,ny)=elevations(c.x,c.y);
				} else {
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
					closed(nx,ny)=QUEUED;
				}
			}
		} else {
			grid_cellz c=open.top();open.pop();

			if(closed(c.x,c.y)==CLOSED) continue;

			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			elevations(c.x,c.y)=c.z; //TODO?

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!elevations.in_grid(nx,ny)) continue;
				if(closed(nx,ny)!=OPEN) continue;

				closed(nx,ny)=QUEUED;
				if(elevations(nx,ny)>c.z)
					climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
				else if(elevations(nx,ny)<=c.z){
					elevations(nx,ny)=c.z;
					meander.push(grid_cellz(nx,ny,c.z));
				}
			}
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}




















//Barnes6 works the same way as Barnes3, but doesn't immediately added cells found in climbing to the Open queue.
void pit_fill_barnes6(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	std::queue<grid_cellz> meander;
	std::queue<grid_cellz> climb;
	std::vector<grid_cellz> potential_open;
	char_2d closed;	//TODO: This could probably be made into a boolean
	float_2d info;
	unsigned long processed_cells=0;
	const int CLOSED=1, QUEUED=2, OPEN=3;
	unsigned long pitc=0,openc=0,climbc=0,climbpc=0;

	diagnostic("\n###Barnes Pit Fill v6\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(OPEN);
	diagnostic("succeeded.\n");
	info.resize(elevations.width(),elevations.height());
	info.init(9e12);

	push_edges(elevations, open, closed, QUEUED);

	diagnostic("Performing the Barnes6 fill...\n");
	progress_bar(-1);
	do{
		while(open.size()>0 || meander.size()>0 || climb.size()>0){
			if(meander.size()>0){
				grid_cellz c=meander.front();meander.pop();
				closed(c.x,c.y)=CLOSED;
				processed_cells++;
				pitc++;

				for(int n=1;n<=8;n++){
					int nx=c.x+dx[n];
					int ny=c.y+dy[n];
					if(!elevations.in_grid(nx,ny)) continue;
					if(closed(nx,ny)!=OPEN) continue;

					closed(nx,ny)=QUEUED;
					if(elevations(nx,ny)<=c.z){
						elevations(nx,ny)=c.z;
						meander.push(grid_cellz(nx,ny,c.z));
					} else
						climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
				}
			} else if (climb.size()>0){
				grid_cellz c=climb.front();climb.pop();
				closed(c.x,c.y)=CLOSED;
				processed_cells++;
				climbc++;

				for(int n=1;n<=8;n++){
					int nx=c.x+dx[n];
					int ny=c.y+dy[n];
					if(!elevations.in_grid(nx,ny)) continue;
					if(closed(nx,ny)!=OPEN || c.z>=info(nx,ny))
						continue;
					else if(elevations(nx,ny)<c.z){
						potential_open.push_back(grid_cellz(nx,ny,c.z));
						info(nx,ny)=elevations(c.x,c.y);
						climbpc++;
					} else {
						climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
						closed(nx,ny)=QUEUED;
					}
				}
			} else {
				for(std::vector<grid_cellz>::iterator i=potential_open.begin();i!=potential_open.end();i++)
					if(closed(i->x,i->y)!=2)
						open.push(*i);
				potential_open.clear();

				grid_cellz c=open.top();open.pop();

				if(closed(c.x,c.y)==CLOSED) continue;

				openc++;

				closed(c.x,c.y)=CLOSED;
				processed_cells++;

				elevations(c.x,c.y)=c.z; //TODO?

				for(int n=1;n<=8;n++){
					int nx=c.x+dx[n];
					int ny=c.y+dy[n];
					if(!elevations.in_grid(nx,ny)) continue;
					if(closed(nx,ny)!=OPEN) continue;

					closed(nx,ny)=QUEUED;
					if(elevations(nx,ny)>c.z)
						climb.push(grid_cellz(nx,ny,elevations(nx,ny)));
					else if(elevations(nx,ny)<=c.z){
						elevations(nx,ny)=c.z;
						meander.push(grid_cellz(nx,ny,c.z));
					}
				}
			}
			progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
		}
	} while (open.size()>0);
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
	diagnostic_arg("Pit=%ld, Climb=%ld, Open=%ld, Climb Push=%ld\n",pitc,climbc,openc,climbpc);
}*/
