/**
  @file
  @brief Defines a 3D array object with convenient methods for working raster 
         data where information about neighbours needs to be stored and 
         processed

  Richard Barnes (rbarnes@umn.edu), 2018
*/
#ifndef _richdem_array_3d_hpp_
#define _richdem_array_3d_hpp_

#include "gdal.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <typeinfo>
#include <stdexcept>
#include <limits>
#include <ctime>         //Used for timestamping output files
#include <unordered_set> //For printStamp
#include <stdexcept>
#include <map>
#include "richdem/common/Array2D.hpp"
#include "richdem/common/logger.hpp"
#include "richdem/common/version.hpp"
#include "richdem/common/constants.hpp"
#include "richdem/common/ManagedVector.hpp"


namespace richdem {

template<typename> class Array2D;

/**
  @brief  Class to hold and 2D rasters with neighbour information
  @author Richard Barnes (rbarnes@umn.edu)

  Array3D manages a two-dimensional raster dataset with information about
  neighbours.

  Array3D implements two addressing schemes: "xyn" and "i". All methods are
  available in each scheme; users may use whichever is convenient. The xyn-
  scheme accesses raster cells by their xyn-coordinates. The i-scheme accesses
  cells by their address in a flat array. Internally, xyn-addresses are
  converted to i-addresses. i-addressing is frequently faster because it reduces
  the space needed to store coordinates and requires no addressing mathematics;
  however, xyn-addressing may be more intuitive. It is suggested to develop
  algorithms using xyn-addressing and then convert them to i-addressing if
  additional speed is desired. The results of the two versions can then be
  compared against each other to verify that using i-addressing has not
  introduced any errors.
*/
template<class T>
class Array3D {
 public:
  std::string filename;             ///< File, if any, from which the data was loaded
  std::string basename;             ///< Filename without path or extension
  std::vector<double> geotransform; ///< Geotransform of the raster
  std::string projection;           ///< Projection of the raster
  std::map<std::string, std::string> metadata; ///< Raster's metadata in key-value pairs

  //Using uint32_t for i-addressing allows for rasters of ~65535^2. These 
  //dimensions fit easily within an int32_t xy-address.
  typedef int32_t     xy_t;         ///< xy-addressing data type
  typedef std::size_t i_t;          ///< i-addressing data type
  typedef uint8_t     n_t;          ///< neighbour addressing data type

  static const i_t NO_I = std::numeric_limits<i_t>::max(); //TODO: What is this?

 private:
  template<typename> friend class Array2D;
  template<typename> friend class Array3D;

  ManagedVector<T> data;            ///< Holds the raster data in a 1D array
                                    ///< this improves caching versus a 2D array

  T   no_data;                       ///< NoData value of the raster
  mutable i_t num_data_cells = NO_I; ///< Number of cells which are not NoData

  xy_t view_width  = 0;              ///< Height of raster in cells
  xy_t view_height = 0;              ///< Width of raster in cells

  ///@{ A rectangular subregion of a larger raster can be extracted. These
  ///   variables store the offsets of this subregion in case the subregion
  ///   needs to be saved into a raster with other subregions
  xy_t view_xoff = 0;
  xy_t view_yoff = 0;
  ///@}

 public:
  Array3D() = default;

  /**
    @brief Creates a raster of the specified dimensions

    @param[in] width   Width of the raster
    @param[in] height  Height of the raster
    @param[in] val     Initial value of all the raster's cells. Defaults to the
                       Array3D template type's default value
  */
  Array3D(xy_t width, xy_t height, const T& val = T()) : Array3D() {
    resize(width,height,val);
  }

  /**
    @brief Wraps a flat array in an Array3D object.

    Wraps a flat array in an Array3D object. The Array3D does not take ownership
    of the data.

    @param[in] data0   Pointer to data to wrap
    @param[in] width   Width of the data
    @param[in] height  Height of the data
  */
  Array3D(T *data0, const xy_t width, const xy_t height) : Array3D() {
    data        = ManagedVector<T>(data0, 9*width*height);
    view_width  = width;
    view_height = height;
  }

  /**
    @brief Create a raster with the same properties and dimensions as another
           raster. No data is copied between the two.

    @param[in] other   Raster whose properties and dimensions should be copied
    @param[in] val     Initial value of all the raster's cells.
  */
  template<class U>
  Array3D(const Array3D<U> &other, const T& val=T()) : Array3D() {
    view_width         = other.view_width;
    view_height        = other.view_height;
    view_xoff          = other.view_xoff;
    view_yoff          = other.view_yoff;
    geotransform       = other.geotransform;
    metadata           = other.metadata;
    projection         = other.projection;
    basename           = other.basename;    
    resize(other.width(), other.height(), val);
  }

  /**
    @brief Create a raster with the same properties and dimensions as another
           raster. No data is copied between the two.

    @param[in] other   Raster whose properties and dimensions should be copied
    @param[in] val     Initial value of all the raster's cells.
  */
  template<class U>
  Array3D(const Array2D<U> &other, const T& val=T()) : Array3D() {
    view_width         = other.view_width;
    view_height        = other.view_height;
    view_xoff          = other.view_xoff;
    view_yoff          = other.view_yoff;
    geotransform       = other.geotransform;
    metadata           = other.metadata;
    projection         = other.projection;
    basename           = other.basename;    
    resize(other.width(), other.height(), val);
  }

  ///Returns a pointer to the internal data array
  T* getData() { return data.data(); }

  ///@brief Number of cells in the DEM
  i_t size() const { return view_width*view_height; }

  ///Width of the raster
  xy_t width() const { return view_width; }

  ///Height of the raster
  xy_t height() const { return view_height; }

  ///X-Offset of this subregion of whatever raster we loaded from
  xy_t viewXoff() const { return view_xoff; }

  ///Y-Offset of this subregion of whatever raster we loaded from
  xy_t viewYoff() const { return view_yoff; }

  ///Returns TRUE if no data is present in RAM
  bool empty() const { return data.empty(); }

  ///Returns the NoData value of the raster. Cells equal to this value sould
  ///generally not be used in calculations. But note that the isNoData() method
  ///is a much better choice for testing whether a cell is NoData or not.
  T noData() const { return no_data; }

  //TODO
  inline i_t i0() const {
    return (i_t)0;
  }

  /**
    @brief Convert from x,y coordinates to index coordinates

    @param[in]  x   X-coordinate to convert
    @param[in]  y   Y-coordinate to convert

    @return Returns the index coordinate i of (x,y)
  */
  inline i_t xyToI(xy_t x, xy_t y, n_t n) const {
    assert(0<=n && n<=9);
    return (i_t)y*9*(i_t)view_width+9*(i_t)x+n;
  }

  /**
    @brief Determine if two rasters are equivalent based on dimensions,
           NoData value, and their data
  */
  bool operator==(const Array3D<T> &o) const {
    if(width()!=o.width() || height()!=o.height())
      return false;
    for(unsigned int i=0;i<o.data.size();i++)
      if(data[i]!=o.data[i])
        return false;
    return true;
  }

  /**
    @brief Whether or not a cell is NoData using x,y coordinates

    @param[in]  x   X-coordinate of cell to test
    @param[in]  y   Y-coordinate of cell to test

    @return Returns TRUE if the cell is NoData
  */
  inline bool isNoData(xy_t x, xy_t y) const {
    assert(0<=x && x<view_width);
    assert(0<=y && y<view_height);
    return operator()(x,y,0)==no_data;
  }

  /**
    @brief Whether or not a cell is NoData using i coordinates

    @param[in]  i   i-coordinate of cell to test

    @return Returns TRUE if the cell is NoData
  */
  inline bool isNoData(i_t i) const {
    assert(0<=i && i<size());
    return getIN(i,0)==no_data;
  }

  /**
    @brief Test whether a cell lies within the boundaries of the raster

    @param[in]  x   X-coordinate of cell to test
    @param[in]  y   Y-coordinate of cell to test

    @return TRUE if cell lies within the raster
  */
  inline bool inGrid(xy_t x, xy_t y) const {
    return 0<=x && x<view_width && 0<=y && y<view_height;
  }

  /*
    @brief Test whether a cell lies within the boundaries of the raster.

    Obviously this bears some difference from `inGrid(x,y)`.

    @param[in]  i   i-coordinate of cell to test

    @return TRUE if cell lies within the raster
  */
  // bool inGrid(i_t i) const {
  //   return 0<=i && i<size();
  // }

  /**
    @brief Sets the NoData value of the raster

    @param[in]   ndval    Value to change NoData to
  */
  void setNoData(const T &ndval){
    no_data = ndval;
  }

  /**
    @brief Sets all of the raster's cells to 'val'

    @param[in]   val      Value to change the cells to
  */
  void setAll(const T val){
    for(i_t i=0;i<data.size();i++)
      data[i] = val;
  }

  /**
    @brief Resize the raster. Note: this clears all the raster's data.

    @param[in]   width0    New width of the raster
    @param[in]   height0   New height of the raster
    @param[in]   val0      Value to set all the cells to. Defaults to the 
                          raster's template type default value
  */
  void resize(const xy_t width0, const xy_t height0, const T& val0 = T()){
    data.resize((i_t)9*(i_t)width0*(i_t)height0);

    view_width  = width0;
    view_height = height0;

    setAll(val0);
  }

  /*
    @brief Resize a raster to copy another raster's dimensions. Copy properies.

    @param[in]   other    Raster to match sizes with
    @param[in]   val      Value to set all the cells to. Defaults to the
                          raster's template type default value
  */
  template<class U>
  void resize(const Array3D<U> &other, const T& val = T()){
    resize(other.width(), other.height(), val);
  }

  /**
    @brief Counts the number of cells which are not NoData.
  */
  void countDataCells() const {
    num_data_cells = 0;
    for(unsigned int i=0;i<size();i++)
      if(data[i]!=no_data)
        num_data_cells++;
  }

  /**
    @brief Returns the number of cells which are not NoData. May count them.

    @return Returns the number of cells which are not NoData.
  */
  i_t numDataCells() const {
    if(num_data_cells==NO_I)
      countDataCells();
    return num_data_cells;
  }


  /**
    @brief Return cell value based on x,y coordinates

    @param[in]   x    X-coordinate of cell whose data should be fetched.
    @param[in]   y    Y-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by x,y
  */
  T& operator()(xy_t x, xy_t y, n_t n){
    assert(x>=0);
    assert(y>=0);
    assert(n>=0);
    assert(x<width());
    assert(y<height());
    assert(n<9);
    return data[xyToI(x,y,n)];
  }

  /**
    @brief Return cell value based on x,y coordinates

    @param[in]   x    X-coordinate of cell whose data should be fetched.
    @param[in]   y    Y-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by x,y
  */
  T operator()(xy_t x, xy_t y, n_t n) const {
    assert(x>=0);
    assert(y>=0);
    assert(n>=0);
    assert(x<width());
    assert(y<height());
    assert(n<9);
    return data[xyToI(x,y,n)];
  }

  T getIN(i_t i, n_t n) const {
    assert(0<=i);
    assert(i<size());
    assert(0<=n);
    assert(n<9);
    return data[9*i+n];
  }

  T& getIN(i_t i, n_t n){
    assert(0<=i);
    assert(i<size());
    assert(0<=n);
    assert(n<9);
    return data[9*i+n];
  }

  ///Clears all raster data from RAM
  void clear(){
    data = ManagedVector<T>();
  }

  //TODO
  inline bool owned() const {
    return data.owned();
  }
};

}

#endif
