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

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace richdem {

///sqrt(2), used to generate distances from a central cell to its neighbours
constexpr double SQRT2 = 1.414213562373095048801688724209698078569671875376948;


//Constant used to hold D8 flow directions
typedef uint8_t d8_flowdir_t;
typedef int8_t  flowdir_t;

//D8 Neighbour Directions

//Facet                 0   1   2   3   4  5  6  7   8
constexpr std::array<int, 9> d8x = {0, -1, -1,  0,  1, 1, 1, 0, -1}; ///< x offsets of D8 neighbours, from a central cell
constexpr std::array<int, 9> d8y = {0,  0, -1, -1, -1, 0, 1, 1,  1}; ///< y offsets of D8 neighbours, from a central cell
constexpr double d8r[9]  = {0,1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2};
constexpr std::array<bool, 9> n8_diag = {false,  false,  true,  false,  true, false, true, false,  true}; ///< True along diagonal directions, false along north, south, east, west
constexpr int D8_WEST    = 1;
constexpr int D8_NORTH   = 3;
constexpr int D8_EAST    = 5;
constexpr int D8_SOUTH   = 7;

constexpr std::array<int, 5> d4x = {0, -1,  0, 1, 0}; ///< x offsets of D4 neighbours, from a central cell
constexpr std::array<int, 5> d4y = {0,  0, -1, 0, 1}; ///< y offsets of D4 neighbours, from a central cell
constexpr std::array<double, 5> d4r = {0, 1, 1, 1, 1};
constexpr std::array<bool, 5> n4_diag = {false, false, false, false, false};
constexpr int D4_WEST  = 1;
constexpr int D4_NORTH = 2;
constexpr int D4_EAST  = 3;
constexpr int D4_SOUTH = 4;

///@brief Directions from neighbours to the central cell.
///Neighbours are labeled 0-8. This is the inverse direction leading from a
///neighbour to the central cell.
constexpr std::array<int, 9> d8_inverse = {0,5,6,7,8,1,2,3,4};

constexpr std::array<int, 5> d4_inverse = {0, 3, 4, 1, 2};

///Distances from a central cell to each of its 8 neighbours
constexpr std::array<double, 9> dr = {0,1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2}; //TODO: Check these for new D8 directions

///Convert from RichDEM flowdirs to ArcGIS flowdirs
constexpr uint8_t d8_arcgis[9] = {0,16,32,64,128,1,2,4,8};

///Used to indicate that a flowdir cell is NoData
constexpr uint8_t FLOWDIR_NO_DATA = 255;

///Value used to indicate that a cell does not have a defined flow direction
//(i.e. that it has no local gradient)
constexpr flowdir_t NO_FLOW = 0;

///Value used to indicate NoFlow in generic flow metric outputs
constexpr float NO_FLOW_GEN  = -1;
constexpr float HAS_FLOW_GEN =  0;
constexpr float NO_DATA_GEN  = -2;

///Value used to indicate that a flow accumulation cell is NoData
constexpr int32_t ACCUM_NO_DATA = -1;

//These are used predominantly by the parallel algorithms/programs for working
//on tiled datasets.
constexpr uint8_t GRID_LEFT   = 1; ///< Indicates a tile is on the LHS of a DEM
constexpr uint8_t GRID_TOP    = 2; ///< Indicates a tile is on the top of a DEM
constexpr uint8_t GRID_RIGHT  = 4; ///< Indicates a tile is on the RHS of a DEM
constexpr uint8_t GRID_BOTTOM = 8; ///< Indicates a tile is on the bottom of a DEM

enum class Topology {
  D8,
  D4
};

template<Topology topo>
constexpr auto get_dx_for_topology() {
  if constexpr (topo==Topology::D8){
    return d8x;
  } else if constexpr (topo==Topology::D4){
    return d4x;
  } else {
    //static_assert(false, "Unknown topology!");
  }
}

template<Topology topo>
constexpr auto get_dy_for_topology() {
  if constexpr (topo==Topology::D8){
    return d8y;
  } else if constexpr (topo==Topology::D4){
    return d4y;
  } else {
    //static_assert(false, "Unknown topology!");
  }
}

template<Topology topo>
constexpr auto get_nmax_for_topology() {
  if constexpr (topo==Topology::D8){
    return 8;
  } else if constexpr (topo==Topology::D4){
    return 4;
  } else {
    //static_assert(false, "Unknown topology!");
  }
}

template<Topology topo>
constexpr auto get_dinverse_for_topology() {
  if constexpr (topo==Topology::D8){
    return d8_inverse;
  } else if constexpr (topo==Topology::D4){
    return d4_inverse;
  } else {
    //static_assert(false, "Unknown topology!");
  }
}

template<Topology topo>
constexpr auto get_n_diag_for_topology() {
  if constexpr (topo==Topology::D8){
    return n8_diag;
  } else if constexpr (topo==Topology::D4){
    return n4_diag;
  } else {
    //static_assert(false, "Unknown topology!");
  }
}

inline std::string TopologyName(Topology topo){
  switch(topo){
    case Topology::D8: return "D8";
    case Topology::D4: return "D4";
    default:
      throw std::runtime_error("Unrecognised topology!");
  }
}

}

#endif
