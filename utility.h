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

#define d8_NO_DATA		-100
#define dinf_NO_DATA	-100
#define NO_FLOW			-1
#define SQRT2			1.414213562373095048801688724209698078569671875376948

#define MIN(A,B) 	(((A)<(B))?(A):(B))
#define MAX(A,B) 	(((A)>(B))?(A):(B))
#define DEG(A)		((A)*180.0/M_PI)

#define ROUND(A)	floor((A) + 0.5)

//D8 Directions
const int dx[9]={0,-1,-1,0,1,1,1,0,-1};	//TODO: These should be merged with my new dinf_d8 to reflect a uniform and intelligent directional system
const int dy[9]={0,0,-1,-1,-1,0,1,1,1};
const int dr[9]={0,1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2};
const int inverse_flow[9]={0,5,6,7,8,1,2,3,4}; //Inverse of a given n from chart below
const std::string fd[9]={"·","←","↖","↑","↗","→","↘","↓","↙"};
//derived from the following neighbour directions
//234
//105
//876

void print_dem(const float_2d &elevations, int mark_x=-1, int mark_y=-1, int colour=91);

//TODO: Don't think I'm using the must_be functionality any more
//Used for checking input from files for structured tags
struct must_be{
	std::string val;
    must_be(const std::string val):val(val){}
};

std::istream& operator>>(std::istream& inputstream, const must_be &a);

#endif
