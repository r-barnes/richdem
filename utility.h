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

//D8 Directions
const int dx[9]={0,-1,-1,0,1,1,1,0,-1};
const int dy[9]={0,0,-1,-1,-1,0,1,1,1};
const int inverse_flow[9]={0,5,6,7,8,1,2,3,4}; //Inverse of a given n from chart below
//derived from the following neighbour directions
//234
//105
//876


//std::string fd[9]={"·","←","↖","↑","↗","→","↘","↓","↙"};

void print_dem(const float_2d &elevations, int mark_x=-1, int mark_y=-1, int colour=91);
void print_flow(const float_2d &flowdirs);
void print_char_2d(const char_2d &arr);

#endif
