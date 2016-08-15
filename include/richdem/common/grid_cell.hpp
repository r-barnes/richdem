/**
  @file
  Defines grid cells data structures and the queues that hold them.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_grid_hpp_
#define _richdem_grid_hpp_

#include <vector>
#include <queue>
#include <cmath>

/// Stores the (x,y) coordinates of a grid cell
class GridCell {
 public:
  int x; ///< Grid cell's x-coordinate
  int y; ///< Grid cell's y-coordinate
  /// Initiate the grid cell without coordinates; should generally be avoided
  GridCell(){}
  /// Initiate the grid cell to the coordinates (x0,y0)
  GridCell(int x, int y):x(x),y(y){}
};


/// Stores the (x,y,z) coordinates of a grid cell; useful for priority sorting
/// with \ref grid_cellz_compare.
/// TODO: Specialize for integers
template<class elev_t>
class GridCellZ : public GridCell {
 public:
  elev_t z;         ///< Grid cell's z-coordinate
  GridCellZ(int x, int y, elev_t z): GridCell(x,y), z(z) {}
  GridCellZ(){}
  //TODO: Add tests of isnan
  bool operator> (const GridCellZ<elev_t>& a) const { return z>a.z; }
};





template<>
class GridCellZ<double> : public GridCell {
 public:
  double z;         ///< Grid cell's z-coordinate
  GridCellZ(int x, int y, double z): GridCell(x,y), z(z) {}
  GridCellZ(){}
  //TODO: Add tests of isnan
  bool operator< (const GridCellZ<double>& a) const { return ( std::isnan(z) && !std::isnan(a.z)) || z< a.z; }
  bool operator> (const GridCellZ<double>& a) const { return (!std::isnan(z) &&  std::isnan(a.z)) || z> a.z; }
  bool operator>=(const GridCellZ<double>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || (!std::isnan(z) &&  std::isnan(a.z)) || z>=a.z; }
  bool operator<=(const GridCellZ<double>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || ( std::isnan(z) && !std::isnan(a.z)) || z<=a.z; }
  bool operator==(const GridCellZ<double>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || z==a.z; }
  bool operator!=(const GridCellZ<double>& a) const { return !(std::isnan(z) &&  std::isnan(a.z)) && z!=a.z; }
};


template<>
class GridCellZ<float>: public GridCell {
 public:
  float z;         ///< Grid cell's z-coordinate
  GridCellZ(int x, int y, float z): GridCell(x,y), z(z) {}
  GridCellZ(){}
  //TODO: Add tests of isnan
  bool operator< (const GridCellZ<float>& a) const { return ( std::isnan(z) && !std::isnan(a.z)) || z< a.z; }
  bool operator> (const GridCellZ<float>& a) const { return (!std::isnan(z) &&  std::isnan(a.z)) || z> a.z; }
  bool operator>=(const GridCellZ<float>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || (!std::isnan(z) &&  std::isnan(a.z)) || z>=a.z; }
  bool operator<=(const GridCellZ<float>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || ( std::isnan(z) && !std::isnan(a.z)) || z<=a.z; }
  bool operator==(const GridCellZ<float>& a) const { return ( std::isnan(z) &&  std::isnan(a.z)) || z==a.z; }
  bool operator!=(const GridCellZ<float>& a) const { return !(std::isnan(z) &&  std::isnan(a.z)) && z!=a.z; }
};













/// Stores the (x,y,z) coordinates of a grid cell and a priority indicator k;
/// used by grid_cellz_pq.
template<class elev_t>
class GridCellZk : public GridCellZ<elev_t> {
  public:
    int k;           ///< Used to store an integer to make sorting stable
    GridCellZk(int x, int y, elev_t z, int k): GridCellZ<elev_t>(x,y,z), k(k) {}
    GridCellZk(){}
    bool operator< (const GridCellZk<elev_t>& a) const { return GridCellZ<elev_t>::z< a.z || ( std::isnan(GridCellZ<elev_t>::z) && !std::isnan(a.z)) || (GridCellZ<elev_t>::z==a.z && k<a.k) || (std::isnan(GridCellZ<elev_t>::z) && std::isnan(a.z) && k<a.k); }
    bool operator> (const GridCellZk<elev_t>& a) const { return GridCellZ<elev_t>::z> a.z || (!std::isnan(GridCellZ<elev_t>::z) &&  std::isnan(a.z)) || (GridCellZ<elev_t>::z==a.z && k>a.k) || (std::isnan(GridCellZ<elev_t>::z) && std::isnan(a.z) && k>a.k); }
};





///A priority queue of grid_cells, sorted by ascending height
template<typename elev_t>
using GridCellZ_pq = std::priority_queue<GridCellZ<elev_t>, std::vector<GridCellZ<elev_t> >, std::greater<GridCellZ<elev_t> > >;


///A priority queue of grid_cells, sorted by ascending height or, if heights
///are equal, by the order of insertion
template<typename elev_t>
using GridCellZk_pq = std::priority_queue<GridCellZk<elev_t>, std::vector<GridCellZk<elev_t> >, std::greater<GridCellZk<elev_t> > >;

#endif
