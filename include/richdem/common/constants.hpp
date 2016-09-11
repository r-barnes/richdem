#ifndef _richdem_constants_
#define _richdem_constants_

//D8 Neighbour Directions
//234
//105
//876

///x and y offsets of D8 neighbours, from a central cell
const int dx[9]={0, -1, -1,  0,  1, 1, 1, 0, -1};
const int dy[9]={0,  0, -1, -1, -1, 0, 1, 1,  1};

///dx[] and dy[] offsets are labeled 0-8. This maps the inverse path from each of those cells.
const int d8_inverse[9] = {0,5,6,7,8,1,2,3,4};

///Arrows indicating flow directions
const wchar_t fdarrows[9]={L'·',L'←',L'↖',L'↑',L'↗',L'→',L'↘',L'↓',L'↙'};

///sqrt(2), used to generate distances from a central cell to its neighbours
const double SQRT2 = 1.414213562373095048801688724209698078569671875376948;

///Distances from a central cell to each of its 8 neighbours
const double dr[9]={0,1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2}; //TODO: Check these for new D8 directions

//ArcGIS uses:
//32 64 128
//16  0   1
// 8  4   2
const uint8_t d8_arcgis[9] = {0,16,32,64,128,1,2,4,8};

const uint8_t FLOWDIR_NO_DATA = 255;

///Value used to indicate that a cell does not have a defined flow direction
//(i.e. that it has no local gradient)
const uint8_t NO_FLOW = 0;

const int32_t ACCUM_NO_DATA = -1;

//Used for indicating whether a block is on the edge of the larger DEM and which
//edges it is adjacent to
const uint8_t GRID_LEFT   = 1;
const uint8_t GRID_TOP    = 2;
const uint8_t GRID_RIGHT  = 4;
const uint8_t GRID_BOTTOM = 8;

#endif