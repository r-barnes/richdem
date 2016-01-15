/**
  @file
  Template code defining data structures used throughout the package.
  Defines grid cells.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_grid_hpp_
#define _richdem_grid_hpp_

#include <vector>
#include <queue>
#include <cmath>

/// Stores the (x,y) coordinates of a grid cell
class grid_cell {
  public:
    int x; ///< Grid cell's x-coordinate
    int y; ///< Grid cell's y-coordinate
    /// Initiate the grid cell without coordinates; should generally be avoided
    grid_cell(){}
    /// Initiate the grid cell to the coordinates (x0,y0)
    grid_cell(int x, int y):x(x),y(y){}
};


/// Stores the (x,y,z) coordinates of a grid cell; useful for priority sorting
/// with \ref grid_cellz_compare.
template<class elev_t>
class grid_cellz : public grid_cell {
  public:
    elev_t z;         ///< Grid cell's z-coordinate
    grid_cellz(int x, int y, elev_t z): grid_cell(x,y), z(z) {}
    grid_cellz(){}
    bool operator< (const grid_cellz& a) const { return ( isnan(z) && !isnan(a.z)) || z< a.z; }
    bool operator> (const grid_cellz& a) const { return (!isnan(z) &&  isnan(a.z)) || z> a.z; }
    bool operator>=(const grid_cellz& a) const { return ( isnan(z) &&  isnan(a.z)) || (!isnan(z) &&  isnan(a.z)) || z>=a.z; }
    bool operator<=(const grid_cellz& a) const { return ( isnan(z) &&  isnan(a.z)) || ( isnan(z) && !isnan(a.z)) || z<=a.z; }
    bool operator==(const grid_cellz& a) const { return ( isnan(z) &&  isnan(a.z)) || z==a.z; }
    bool operator!=(const grid_cellz& a) const { return !(isnan(z) &&  isnan(a.z)) && z!=a.z; }
};



/// Stores the (x,y,z) coordinates of a grid cell and a priority indicator k;
/// used by grid_cellz_pq.
template<class elev_t>
class grid_cellzk : public grid_cellz<elev_t> {
  public:
    int k;           ///< Used to store an integer to make sorting stable
    grid_cellzk(int x, int y, elev_t z, int k): grid_cellz<elev_t>(x,y,z), k(k) {}
    grid_cellzk(){}
    bool operator< (const grid_cellzk<elev_t>& a) const { return grid_cellz<elev_t>::z< a.z || ( isnan(grid_cellz<elev_t>::z) && !isnan(a.z)) || (grid_cellz<elev_t>::z==a.z && k<a.k) || (isnan(grid_cellz<elev_t>::z) && isnan(a.z) && k<a.k); }
    bool operator> (const grid_cellzk<elev_t>& a) const { return grid_cellz<elev_t>::z> a.z || (!isnan(grid_cellz<elev_t>::z) &&  isnan(a.z)) || (grid_cellz<elev_t>::z==a.z && k>a.k) || (isnan(grid_cellz<elev_t>::z) && isnan(a.z) && k>a.k); }
};

///A priority queue of grid_cells, sorted by ascending height
template<class elev_t>
class grid_cellz_pq : public std::priority_queue<grid_cellz<elev_t>, std::vector<grid_cellz<elev_t> >, std::greater<grid_cellz<elev_t> > > {
  public:
    void push_cell(int x, int y, elev_t z){
      std::priority_queue<grid_cellz<elev_t>, std::vector<grid_cellz<elev_t> >, std::greater<grid_cellz<elev_t> > >::push(grid_cellz<elev_t>(x,y,z));
    }
};

///A priority queue of grid_cells, sorted by ascending height or, if heights
///are equal, by the order of insertion
template<class elev_t>
class grid_cellzk_pq : public std::priority_queue<grid_cellzk<elev_t>, std::vector<grid_cellzk<elev_t> >, std::greater<grid_cellzk<elev_t> > > {
  private:
    int count;
  public:
    grid_cellzk_pq() : count(0) {}
    void push(const grid_cellz<elev_t> &a){
      std::priority_queue<grid_cellzk<elev_t>, std::vector<grid_cellzk<elev_t> >, std::greater<grid_cellzk<elev_t> > >::push(grid_cellzk<elev_t>(a.x,a.y,a.z,count++));
    }
    void push_cell(int x, int y, elev_t z){
      std::priority_queue<grid_cellzk<elev_t>, std::vector<grid_cellzk<elev_t> >, std::greater<grid_cellzk<elev_t> > >::push(grid_cellzk<elev_t>(x,y,z,count++));
    }
};

#endif
