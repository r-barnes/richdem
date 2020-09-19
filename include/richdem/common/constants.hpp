/**
  @file
  @brief Defines a number of constants used by many of the algorithms.

  RichDEM uses the following D8 neighbourhood. This is used by the dx[] and dy[]
  variables, among many others.

      234
      105
      876

  ArcGIS uses the following bits to indicate flow toward a given neighbour:

      32 64 128
      16  0   1
       8  4   2

  D4 directions
    2
   103
    4
*/
#ifndef _richdem_constants_hpp_
#define _richdem_constants_hpp_

#include <cstdint>
#include <stdexcept>
#include <string>

namespace richdem {

///sqrt(2), used to generate distances from a central cell to its neighbours
const double SQRT2 = 1.414213562373095048801688724209698078569671875376948;


//Constant used to hold D8 flow directions
typedef uint8_t d8_flowdir_t;
typedef int8_t  flowdir_t;

//D8 Neighbour Directions

//Facet                 0   1   2   3   4  5  6  7   8
const int  dx[9]     = {0, -1, -1,  0,  1, 1, 1, 0, -1}; ///< x offsets of D8 neighbours, from a central cell
const int  dy[9]     = {0,  0, -1, -1, -1, 0, 1, 1,  1}; ///< y offsets of D8 neighbours, from a central cell
const double d8r[9]  = {0,1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2};
const bool n_diag[9] = {0,  0,  1,  0,  1, 0, 1, 0,  1}; ///< True along diagonal directions, false along north, south, east, west
const int D8_WEST    = 1;
const int D8_NORTH   = 3;
const int D8_EAST    = 5;
const int D8_SOUTH   = 7;

const int *const d8x = dx;
const int *const d8y = dy;
const int d4x[5]   = {0, -1,  0, 1, 0}; ///< x offsets of D4 neighbours, from a central cell
const int d4y[5]   = {0,  0, -1, 0, 1}; ///< y offsets of D4 neighbours, from a central cell
const double d4r[5] = {0, 1, 1, 1, 1};
const int D4_WEST  = 1;
const int D4_NORTH = 2;
const int D4_EAST  = 3;
const int D4_SOUTH = 4;

///@brief Directions from neighbours to the central cell.
///Neighbours are labeled 0-8. This is the inverse direction leading from a
///neighbour to the central cell.
const int d8_inverse[9] = {0,5,6,7,8,1,2,3,4};

const int d4_inverse[5] = {0, 3, 4, 1, 2};

///Distances from a central cell to each of its 8 neighbours
const double dr[9] = {0,1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2}; //TODO: Check these for new D8 directions

///Convert from RichDEM flowdirs to ArcGIS flowdirs
const uint8_t d8_arcgis[9] = {0,16,32,64,128,1,2,4,8};

///Used to indicate that a flowdir cell is NoData
const uint8_t FLOWDIR_NO_DATA = 255;

///Value used to indicate that a cell does not have a defined flow direction
//(i.e. that it has no local gradient)
const flowdir_t NO_FLOW = 0;

///Value used to indicate NoFlow in generic flow metric outputs
const float NO_FLOW_GEN  = -1;
const float HAS_FLOW_GEN =  0;
const float NO_DATA_GEN  = -2;

///Value used to indicate that a flow accumulation cell is NoData
const int32_t ACCUM_NO_DATA = -1;

//These are used predominantly by the parallel algorithms/programs for working
//on tiled datasets.
const uint8_t GRID_LEFT   = 1; ///< Indicates a tile is on the LHS of a DEM
const uint8_t GRID_TOP    = 2; ///< Indicates a tile is on the top of a DEM
const uint8_t GRID_RIGHT  = 4; ///< Indicates a tile is on the RHS of a DEM
const uint8_t GRID_BOTTOM = 8; ///< Indicates a tile is on the bottom of a DEM

enum class Topology {
  D8,
  D4
};

inline std::string TopologyName(Topology topo){
  switch(topo){
    case Topology::D8: return "D8";
    case Topology::D4: return "D4";
    default:
      throw std::runtime_error("Unrecognised topology!");
  }
}

template<Topology topo>
void TopologicalResolver(
  const int*&    dx,
  const int*&    dy,
  const double*& dr,
  const int*&    dinverse,
  int&       neighbours
){
  if(topo==Topology::D4){
    dx         = d4x;
    dy         = d4y;
    dr         = d4r;
    dinverse   = d4_inverse;
    neighbours = 4;
  } else if(topo==Topology::D8){
    dx         = d8x;
    dy         = d8y;
    dr         = d8r;
    dinverse   = d8_inverse;
    neighbours = 8;
  } else {
    throw std::runtime_error("Unrecognised topology!");
  }
}

}

#endif
