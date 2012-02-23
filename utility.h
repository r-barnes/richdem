#ifndef _utility_included
#define _utility_included

#include <string>
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

#endif
