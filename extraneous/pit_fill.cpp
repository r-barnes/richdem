#include "utility.h"
#include "data_structures.h"
#include "interface.h"
#include <queue>
#include <vector>
using namespace std;

//D8 Directions
static int const dx[9]={0,-1,-1,0,1,1,1,0,-1};
static int const dy[9]={0,0,-1,-1,-1,0,1,1,1};

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
			if (reverse) return (lhs->z>rhs->z);
			else return (lhs->z<rhs->z);
		}
};
/*
int pit_fill_yonghe2009(float_2d &elevations){
	priority_queue<grid_cell*, vector<grid_cell*>, grid_cell_compare> pq;
	grid_cell *temp;
	bool_2d marks;

	diagnostic_arg("The boolean flood array will require approximately %ldMB of RAM.\n",elevations.size1()*elevations.size2()*((long)sizeof(bool))/1024/1024);
	diagnostic("Resizing boolean flood array matrix...");
	try{
		marks.resize(elevations.size1(),elevations.size2());
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");
	diagnostic("Initializing boolean flood array matrix...");
	#pragma omp parallel for
	for(int x=0;x<marks.size1();x++)
	for(int y=0;y<marks.size2();y++)
		marks(x,y)=1;
	diagnostic("succeeded.\n");

	diagnostic_arg("The flood queue will require approximately %ldMB of RAM.\n",(elevations.size1()*2+elevations.size2()*2)*((long)sizeof(grid_cell))/1024/1024);
	diagnostic("Adding cells to flood queue...");
	for(int x=0;x<elevations.size1();x++){
		pq.push(new grid_cell( x,0,elevations(x,0) ));
		pq.push(new grid_cell( x,elevations.size2()-1,elevations(x,elevations.size2()-1) ));
	}
	for(int y=0;y<elevations.size2();y++){
		pq.push(new grid_cell( 0,y,elevations(0,y) ));
		pq.push(new grid_cell( elevations.size1()-1,y,elevations(elevations.size1()-1,y) ));
	}
	diagnostic("succeeded.\n");

	double waterlevel=0;
	int ccount=0;
	progress_bar(-1);
	while(pq.size()>0){
		ccount++;
		progress_bar(ccount*100/(elevations.size1()*elevations.size2()));

		temp=pq.top();
		pq.pop();

		waterlevel=temp->z;
		for(int n=1;n<=8;n++){
			if(!IN_GRID(temp->x+dx[n],temp->y+dy[n],elevations.size1(),elevations.size2())) continue;
			if(elevations(temp->x+dx[n],temp->y+dy[n])==elev_NO_DATA) continue;
			if(marks(temp->x+dx[n],temp->y+dy[n])==0)
				continue;
			else
				marks(temp->x+dx[n],temp->y+dy[n])=0;
			else if (elevations(temp->x+dx[n],temp->y+dy[n])<=waterlevel)
				elevations(temp->x+dx[n],temp->y+dy[n])=waterlevel;
			pq.push(new grid_cell( temp->x+dx[n],temp->y+dy[n],elevations(temp->x+dx[n],temp->y+dy[n]) ));
		}
	}
	progress_bar(-1);
}*/
