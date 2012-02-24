#ifndef _utility_included
#define _utility_included

#include <string>
#include "data_structures.h"
//Neighbour directions
//234
//105
//876

//Inverse neighbour directions
//678
//501
//432

#define IN_GRID(x,y,x_max,y_max) (x>=0 && y>=0 && x<x_max && y<y_max)
#define EDGE_GRID(x,y,x_max,y_max) (x==0 || y==0 || x==x_max-1 || y==y_max-1)
#define d8_NO_DATA		-100
#define dinf_NO_DATA	-9999
#define elev_NO_DATA	-9999
#define NO_FLOW			-1

#define MAX(A,B)	(((A)>(B))?(A):(B))
//int dx[9]={0,-1,-1,0,1,1,1,0,-1};
//int dy[9]={0,0,-1,-1,-1,0,1,1,1};
//std::string fd[9]={"·","←","↖","↑","↗","→","↘","↓","↙"};

void print_dem(float_2d &elevations, int mark_x, int mark_y);
void print_flow(float_2d &flowdirs);
void print_char_2d(const char_2d &arr);

typedef struct grid_cell_typez {
	int x;
	int y;
	float z;
	grid_cell_typez(int x0, int y0, float z0){
		x=x0;
		y=y0;
		z=z0;
	}
} grid_cellz;

typedef struct grid_cell_type {
	int x;
	int y;
	grid_cell_type(int x0, int y0){
		x=x0;
		y=y0;
	}
} grid_cell;

#endif
