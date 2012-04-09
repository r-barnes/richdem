#include "utility.h"
#include "data_structures.h"
#include "interface.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <omp.h>
#include <limits>
#include <vector>
#include <queue>
#include <limits>


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
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));

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
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
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
*/
/*
const int dR[8]={0,0,1,-1,0,0,1,-1};
const int dC[8]={1,-1,0,0,-1,1,0,0};
bool Next_Cell(int &C, int &R,int i,const float_2d &e){
	const int fR[8]={1,-1,-e.height()+1,e.height()-1,1,-1,-e.height()+1,e.height()-1};
	const int fC[8]={-e.width()+1,e.width()-1,-1,1,e.width()-1,-e.width()+1,1,-1};
	R+=dr[i];
	C+=dC[i];
	if(R<0 || C<0 || R>=e.height() || C>=e.width()){
		R+=fR[i];
		C+=fC[i];
		if(R<0 || C<0 || R>=e.height() || C>=e.width())
			return false;
	}
	return true;
}

void PlanchonStage1(float_2d &w, float_2d &elevations){
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++)
	for(int y=0;y<elevations.height();y++)
		if(EDGE_GRID(x,y,elevations.width(),elevations.height())
			w(x,y)=elevations(x,y);
		else
			w(x,y)=std::numeric_limits<float>::infinity();
}


int Dry_upward_cell(int x, int y, const float_2d &elevations, float_2d &w, float epsilon){
	std::queue<grid_cell> to_dry;
	to_dry.push(grid_cell(x,y));

	while(to_dry.size()>0){
		grid_cell c=to_dry.front();to_dry.pop();
		for(int n=1;n<=8;n++){
			int nx=x+dx[n];
			int ny=y+dy[n];
			if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
			if(elevations(nx,ny)!=std::numeric_limits<float>::infinity()) continue;
			if(elevations(nx,ny)>=w(x,y)+epsilon[n]){
				w(nx,ny)=elevations(nx,ny);
				Dry_upward_cell(nx,ny,elevations,w,epsilon);
			}
		}
	}
}
/*
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
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));

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
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
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
			if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
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
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
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
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
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








//Barnes1 explores depressions and flats by pushing them onto the meander queue. When meander encounters cells at <= elevation, it meanders over them. Otherwise, it pushes them onto the open queue.
void pit_fill_barnes1(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	std::queue<grid_cellz> meander;
	char_2d closed;	//TODO: This could probably be made into a boolean
	unsigned long processed_cells=0;

	diagnostic("\n###Barnes Pit Fill v1\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(char)/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(0);
	diagnostic("succeeded.\n");

	push_edges(elevations, open, closed);

	diagnostic("Performing the Barnes1 fill...\n");
	progress_bar(-1);
	while(open.size()>0 || meander.size()>0){
		grid_cellz c;
		if(meander.size()>0){
			c=meander.front();
			meander.pop();
		} else {
			c=open.top();
			open.pop();
		}
		closed(c.x,c.y)=2;
		processed_cells++;

		for(int n=1;n<=8;n++){
			int nx=c.x+dx[n];
			int ny=c.y+dy[n];
			if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
			if(closed(nx,ny)>0) 
				continue;

			closed(nx,ny)=1;
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
}







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
			if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
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
			if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
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
	std::queue<grid_cell> meander;
	std::queue<grid_cell> climb;
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









//Barnes3 uses a meander queue in the same way as in Barnes1, but introduces a climb queue. The climb queue can only ascend from cells which have a known outlet (i.e. those popped off the PQ or Meander queue). If a cell is lower than a cell popped off the climb queue, then that cell must be added to the PQ. In Barnes2, this meant that a cell could be placed onto the PQ multiple times. In Barnes3, this is reduced through the use of the info array which has the minimum encountered elevation a cell can be. A popped climb cell will only place a cell into the PQ if its estimated elevation is lower than that already in the info array.
void pit_fill_barnes3(float_2d &elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cell_compare> open;
	std::queue<grid_cell> meander;
	std::queue<grid_cell> climb;
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
			grid_cell c=meander.front();meander.pop();
			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)!=OPEN)
					continue;
				else if(elevations(nx,ny)<=elevations(c.x,c.y)){
					elevations(nx,ny)=elevations(c.x,c.y);
					meander.push(grid_cell(nx,ny));
					closed(nx,ny)=QUEUED;
				} else {
					climb.push(grid_cell(nx,ny));
					closed(nx,ny)=QUEUED;
				}
			}
		} else if (climb.size()>0){
			grid_cell c=climb.front();climb.pop();
			closed(c.x,c.y)=CLOSED;
			processed_cells++;

			for(int n=1;n<=8;n++){
				int nx=c.x+dx[n];
				int ny=c.y+dy[n];
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)!=OPEN || elevations(c.x,c.y)>=info(nx,ny))
					continue;
				else if(elevations(nx,ny)<elevations(c.x,c.y)){
					open.push(grid_cellz(nx,ny,elevations(c.x,c.y)));
					info(nx,ny)=elevations(c.x,c.y);
				} else {
					climb.push(grid_cell(nx,ny));
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
				if(!IN_GRID(nx,ny,elevations.width(),elevations.height())) continue;
				if(closed(nx,ny)!=OPEN)
					continue;
				else if(elevations(nx,ny)>elevations(c.x,c.y)){
					climb.push(grid_cell(nx,ny));
					closed(nx,ny)=QUEUED;
				} else if(elevations(nx,ny)<=elevations(c.x,c.y)){
					elevations(nx,ny)=elevations(c.x,c.y);
					meander.push(grid_cell(nx,ny));
					closed(nx,ny)=QUEUED;
				}
			}
		}
		progress_bar(processed_cells*100/(elevations.width()*elevations.height()));
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
	diagnostic_arg("%ld cells processed.\n",processed_cells);
}
