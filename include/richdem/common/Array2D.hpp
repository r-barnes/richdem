/**
  @file
  @brief Defines a 2D array object with many convenient methods for working with raster data, along with several functions for checking file data types.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_array_2d_hpp_
#define _richdem_array_2d_hpp_

#include "gdal.hpp"
#include <array>
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
#include <cmath>
#include <unordered_set> //For printStamp
#include <stdexcept>
#include <map>
#include "richdem/common/Array3D.hpp"
#include "richdem/common/logger.hpp"
#include "richdem/common/version.hpp"
#include "richdem/common/constants.hpp"
#include "richdem/common/ManagedVector.hpp"

//These enable compression in the loadNative() and saveNative() methods
#ifdef WITH_COMPRESSION
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#endif

namespace richdem {

template<typename> class Array3D;

std::map<std::string, std::string> ProcessMetadata(char **metadata){
  std::map<std::string, std::string> ret;
  if(metadata==NULL)
    return ret;

  for(int metstri=0;metadata[metstri]==NULL;metstri++){
    std::string metstr = metadata[metstri];
    const auto equals  = metstr.find("=");
    if(equals==std::string::npos){
      RDLOG_WARN<<"Skipping improper metadata string: '"<<metstr<<"'";
      continue;
    }
    std::string keystr = metstr.substr(0,equals);
    std::string valstr = metstr.substr(equals+1);
    if(ret.count(keystr)>0){
      RDLOG_WARN<<"Duplicate key '"<<keystr<<"' found in metadata. Only latter value will be kept.";
    }
    ret[keystr] = valstr;
  }

  return ret;
}

/**
  @brief  Class to hold and manipulate GDAL and native rasters
  @author Richard Barnes (rbarnes@umn.edu)

  Array2D manages a two-dimensional raster dataset. Passed a request to load
  such data, it peeks at the file header and can either load data on
  construction or wait until a later point. It can also offload data to disk.

  Array2D permits simple copy construction as well as templated copies, which
  transfer projections and geotransforms, but not the actual data. This is
  useful for say, create a flow directions raster which is homologous to a DEM.

  Array2D implements two addressing schemes: "xy" and "i". All methods are
  available in each scheme; users may use whichever is convenient. The xy-scheme
  accesses raster cells by their xy-coordinates. The i-scheme accesses cells by
  their address in a flat array. Internally, xy-addresses are converted to
  i-addresses. i-addressing is frequently faster because it reduces the space
  needed to store coordinates and requires no addressing mathematics; however,
  xy-addressing may be more intuitive. It is suggested to develop algorithms
  using xy-addressing and then convert them to i-addressing if additional speed
  is desired. The results of the two versions can then be compared against each
  other to verify that using i-addressing has not introduced any errors.
*/
template<class T>
class Array2D {
 public:
  std::string filename;             ///< File, if any, from which the data was loaded
  std::string basename;             ///< Filename without path or extension
  std::vector<double> geotransform; ///< Geotransform of the raster
  std::string projection;           ///< Projection of the raster
  std::map<std::string, std::string> metadata; ///< Raster's metadata in key-value pairs

  //Using uint32_t for i-addressing allows for rasters of ~65535^2. These 
  //dimensions fit easily within an int32_t xy-address.
  typedef int32_t  xy_t;            ///< xy-addressing data type
  typedef uint32_t i_t;             ///< i-addressing data type

  static const i_t NO_I = std::numeric_limits<i_t>::max(); //TODO: What is this?

 private:
  template<typename> friend class Array2D;
  template<typename> friend class Array3D;

  std::array<int, 9> _nshift;       ///< Offset to neighbouring cells;

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
  
  ///If TRUE, loadData() loads data from the cache assuming  the Native format.
  ///Otherwise, it assumes it is loading from a GDAL file.
  bool from_cache;

  #ifdef USEGDAL
  ///TODO
  void loadGDAL(const std::string &filename, xy_t xOffset=0, xy_t yOffset=0, xy_t part_width=0, xy_t part_height=0, bool exact=false, bool load_data=true){
    assert(empty());

    from_cache = false;

    this->filename = filename;

    RDLOG_PROGRESS<<"Trying to open file '"<<filename<<"'...";

    GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
    if(fin==NULL)
      throw std::runtime_error("Could not open file '"+filename+"' with GDAL!");

    geotransform.resize(6);
    if(fin->GetGeoTransform(geotransform.data())!=CE_None){
      RDLOG_WARN<<"Could not get a geotransform from '"<<filename<<"'! Setting to an arbitrary standard geotransform.";
      geotransform = {{1000., 1., 0., 1000., 0., -1.}};
    }

    metadata = ProcessMetadata(fin->GetMetadata());

    const char* projection_string=fin->GetProjectionRef();
    projection = std::string(projection_string);

    GDALRasterBand *band = fin->GetRasterBand(1);

    xy_t total_width  = band->GetXSize();         //Returns an int
    xy_t total_height = band->GetYSize();         //Returns an int
    no_data           = band->GetNoDataValue();

    if(exact && (total_width-xOffset!=part_width || total_height-yOffset!=part_height))
      throw std::runtime_error("Tile dimensions did not match expectations!");

    //TODO: What's going on here?

    if(xOffset+part_width>=total_width)
      part_width  = total_width-xOffset;
    if(yOffset+part_height>=total_height)
      part_height = total_height-yOffset;

    if(part_width==0)
      part_width = total_width;
    view_width = part_width;

    if(part_height==0)
      part_height = total_height;
    view_height = part_height;

    view_xoff = xOffset;
    view_yoff = yOffset;

    GDALClose(fin);

    if(load_data)
      loadData();
  }
  #endif


  #ifdef USEGDAL
  ///Returns the GDAL data type of the Array2D template type
  GDALDataType myGDALType() const {
    return NativeTypeToGDAL<T>();
  }
  #endif


  /**
    @brief Saves raster to a simply-structure file on disk, possibly using
           compression.

    @post  Using loadData() after running this function will result in data
           being loaded from the cache, rather than the original file (if any).
  */
  //TODO: Should save metadata
  void saveToCache(const std::string &filename){
    std::fstream fout;

    from_cache     = true;
    this->filename = filename;

    fout.open(filename, std::ios_base::binary | std::ios_base::out | std::ios::trunc);
    if(!fout.good())
      throw std::logic_error("Failed to open cache file '"+filename+"'.");

    #ifdef WITH_COMPRESSION
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::zlib_compressor());
      out.push(fout);
    #else
      auto &out = fout;
    #endif

    out.write(reinterpret_cast<const char*>(&view_height),    sizeof(xy_t));
    out.write(reinterpret_cast<const char*>(&view_width),     sizeof(xy_t));
    out.write(reinterpret_cast<const char*>(&view_xoff),      sizeof(xy_t));
    out.write(reinterpret_cast<const char*>(&view_yoff),      sizeof(xy_t));
    out.write(reinterpret_cast<const char*>(&num_data_cells), sizeof(i_t));
    out.write(reinterpret_cast<const char*>(&no_data),        sizeof(T  ));

    out.write(reinterpret_cast<const char*>(geotransform.data()), 6*sizeof(double));
    std::string::size_type projection_size = projection.size();
    out.write(reinterpret_cast<const char*>(&projection_size), sizeof(std::string::size_type));
    out.write(reinterpret_cast<const char*>(projection.data()), projection.size()*sizeof(const char));

    out.write(reinterpret_cast<const char*>(data.data()), size()*sizeof(T));
  }

  ///TODO
  void loadNative(const std::string &filename, bool load_data=true){
    std::ifstream fin(filename, std::ios::in | std::ios::binary);
    assert(fin.good());

    this->filename = filename;
    from_cache    = true;

    #ifdef WITH_COMPRESSION
      boost::iostreams::filtering_istream in;
      in.push(boost::iostreams::zlib_decompressor());
      in.push(fin);
    #else
      auto &in = fin;
    #endif

    in.read(reinterpret_cast<char*>(&view_height),    sizeof(xy_t));
    in.read(reinterpret_cast<char*>(&view_width),     sizeof(xy_t));
    in.read(reinterpret_cast<char*>(&view_xoff),      sizeof(xy_t));
    in.read(reinterpret_cast<char*>(&view_yoff),      sizeof(xy_t));
    in.read(reinterpret_cast<char*>(&num_data_cells), sizeof(i_t));
    in.read(reinterpret_cast<char*>(&no_data),        sizeof(T  ));
    geotransform.resize(6);
    in.read(reinterpret_cast<char*>(geotransform.data()), 6*sizeof(double));

    std::string::size_type projection_size;
    in.read(reinterpret_cast<char*>(&projection_size), sizeof(std::string::size_type));
    projection.resize(projection_size,' ');
    in.read(reinterpret_cast<char*>(&projection[0]), projection.size()*sizeof(char));

    if(load_data){
      resize(view_width,view_height);
      in.read(reinterpret_cast<char*>(data.data()), size()*sizeof(T));
    }
  }

 public:
  Array2D(){
    #ifdef USEGDAL
    GDALAllRegister();
    #endif
  }

  /**
    @brief Creates a raster of the specified dimensions

    @param[in] width   Width of the raster
    @param[in] height  Height of the raster
    @param[in] val     Initial value of all the raster's cells. Defaults to the
                       Array2D template type's default value
  */
  Array2D(xy_t width, xy_t height, const T& val = T()) : Array2D() {
    resize(width,height,val);
  }

  /**
    @brief Wraps a flat array in an Array2D object.

    Wraps a flat array in an Array2D object. The Array2D does not take ownership
    of the data.

    @param[in] data0   Pointer to data to wrap
    @param[in] width   Width of the data
    @param[in] height  Height of the data
  */
  Array2D(T *data0, const xy_t width, const xy_t height) : Array2D() {
    data        = ManagedVector<T>(data0, width*height);
    view_width  = width;
    view_height = height;
    _nshift     = {{0,-1,-width-1,-width,-width+1,1,width+1,width,width-1}};
  }

  /**
    @brief Create a raster with the same properties and dimensions as another
           raster. No data is copied between the two.

    @param[in] other   Raster whose properties and dimensions should be copied
    @param[in] val     Initial value of all the raster's cells.
  */
  template<class U>
  Array2D(const Array2D<U> &other, const T& val=T()) : Array2D() {
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
  Array2D(const Array3D<U> &other, const T& val=T()) : Array2D() {
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

  Array2D(const std::string &filename) : Array2D(filename, false, 0,0,0,0, false, true) {}

  ///TODO
  Array2D(const std::string &filename, bool native, xy_t xOffset=0, xy_t yOffset=0, xy_t part_width=0, xy_t part_height=0, bool exact=false, bool load_data=true) : Array2D() {
    if(native){
      loadNative(filename, load_data);
    } else {
      #ifdef USEGDAL
      loadGDAL(filename, xOffset, yOffset, part_width, part_height, exact, load_data);
      #else
        throw std::runtime_error("RichDEM was not compiled with GDAL!");
      #endif
    }
  }

  void setCacheFilename(const std::string &filename){
    this->filename = filename;
  }

  /**
    @brief Caches the raster data and all its properties to disk. Data is then
           purged from RAM.

    @post  Calls to loadData() after this will result in data being loaded from
           the cache.
  */
  void dumpData(){
    saveToCache(filename);
    clear();
  }

  /**
    @brief Loads data from disk into RAM. 

    If dumpData() has been previously called, data is loaded from the cache; 
    otherwise, it is loaded from a GDAL file. No data is loaded if data is
    already present in RAM.
  */
  void loadData() {
    if(!data.empty())
      throw std::runtime_error("Data already loaded!");

    if(from_cache){
      loadNative(filename, true);
    } else {
      #ifdef USEGDAL
      GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
        if(fin==NULL)
          throw std::runtime_error("Failed to loadData() into tile from '"+filename+"'");

      GDALRasterBand *band = fin->GetRasterBand(1);

      resize(view_width,view_height);
      auto temp = band->RasterIO( GF_Read, view_xoff, view_yoff, view_width, view_height, data.data(), view_width, view_height, myGDALType(), 0, 0 );
        if(temp!=CE_None)
          throw std::runtime_error("An error occured while trying to read '"+filename+"' into RAM with GDAL.");

      GDALClose(fin);
      #else
        throw std::runtime_error("RichDEM was not compiled with GDAL!");
      #endif
    }
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

  ///Finds the minimum value of the raster, ignoring NoData cells
  T min() const {
    T minval = std::numeric_limits<T>::max();
    for(unsigned int i=0;i<size();i++){
      if(data[i]==no_data)
        continue;
      minval = std::min(minval,data[i]);
    }
    return minval;
  }

  ///Finds the maximum value of the raster, ignoring NoData cells
  T max() const {
    T maxval = std::numeric_limits<T>::min();
    for(unsigned int i=0;i<size();i++){
      if(data[i]==no_data)
        continue;
      maxval = std::max(maxval,data[i]);
    }
    return maxval;
  }

  /**
    @brief Replace one cell value with another throughout the raster. Can
           operate on NoData cells.

    @param[in] oldval   Value to be replaced
    @param[in] newval   Value to replace 'oldval' with
  */
  void replace(const T oldval, const T newval){
    for(unsigned int i=0;i<size();i++)
      if(data[i]==oldval)
        data[i] = newval;
  }

  /**
    @brief Counts the number of occurrences of a particular value in the raster.
           Can operate on NoData cells.

    @param[in] val   Value to be be counted

    @return The number of times 'val' appears in the raster. Will be 0 if raster
            is not loaded in RAM.
  */
  i_t countval(const T val) const {
    //TODO: Warn if raster is empty?
    i_t count=0;
    for(unsigned int i=0;i<size();i++)
      if(data[i]==val)
        count++;
    return count;
  }

  //TODO
  inline i_t i0() const {
    return (i_t)0;
  }

  /**
    @brief Convert from index coordinates to x,y coordinates

    @param[in]  i   Index coordinate
    @param[out] x   X-coordinate of i
    @param[out] y   Y-coordinate of i
  */
  void iToxy(const i_t i, xy_t &x, xy_t &y) const {
    x = i%view_width;
    y = i/view_width;
  }

  /**
    @brief Convert from x,y coordinates to index coordinates

    @param[in]  x   X-coordinate to convert
    @param[in]  y   Y-coordinate to convert

    @return Returns the index coordinate i of (x,y)
  */
  inline i_t xyToI(xy_t x, xy_t y) const {
    return (i_t)y*(i_t)view_width+(i_t)x;
  }

  /**
    @brief Given a cell identified by an i-coordinate, return the i-coordinate
           of the neighbour identified by dx,dy

    @param[in]  i   i-coordinate of cell whose neighbour needs to be identified
    @param[in] dx   x-displacement of the neighbour from i
    @param[in] dy   y-displacement of the neighbour from i

    @return i-coordinate of the neighbour. Usually referred to as 'ni'
  */
  i_t nToI(i_t i, xy_t dx, xy_t dy) const {
    int32_t x=i%view_width+dx;
    int32_t y=i/view_width+dy;
    if(x<0 || y<0 || x>=view_width || y>=view_height)
      return NO_I;
    return xyToI(x,y);
  }

  /**
    @brief Given a cell identified by an i-coordinate, return the i-coordinate
           of the neighbour identified by n

    @param[in]  i   i-coordinate of cell whose neighbour needs to be identified
    @param[in]  n   Neighbour to be identified

    @return i-coordinate of the neighbour. Usually referred to as 'ni'
  */
  i_t getN(i_t i, uint8_t n) const {
    assert(0<=n && n<=8);
    xy_t x = i%view_width+(xy_t)dx[n];
    xy_t y = i/view_width+(xy_t)dy[n];
    if(x<0 || y<0 || x>=view_width || y>=view_height)
      return NO_I;
    return xyToI(x,y);
  }

  /**
    @brief Return the offset of the neighbour cell identified by n

    @param[in]  n   Neighbour for which offset should be retrieved

    @return Offset of the neighbour n
  */
  inline int nshift(const uint8_t n) const {
    assert(0<=n && n<=8);
    return _nshift[n];
  }

  /**
    @brief Determine if two rasters are equivalent based on dimensions,
           NoData value, and their data
  */
  bool operator==(const Array2D<T> &o) const {
    if(width()!=o.width() || height()!=o.height())
      return false;
    if(noData()!=o.noData())
      return false;
    for(unsigned int i=0;i<o.size();i++)
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
    return data[xyToI(x,y)]==no_data;
  }

  /**
    @brief Whether or not a cell is NoData using i coordinates

    @param[in]  i   i-coordinate of cell to test

    @return Returns TRUE if the cell is NoData
  */
  inline bool isNoData(i_t i) const {
    assert(0<=i && i<size());
    return data[i]==no_data;
  }

  /**
    @brief Flips the raster from top to bottom
  */
  void flipVert(){
    for(xy_t y=0;y<view_height/2;y++)
    for(xy_t x=0;x<view_width;x++)
      std::swap(data[xyToI(x,y)], data[xyToI(x,view_height-1-y)]);
  }

  /**
    @brief Flips the raster from side-to-side
  */
  void flipHorz(){
    for(xy_t y=0;y<view_height;y++){
      T* start = &data[xyToI(0,y)];
      T* end   = &data[xyToI(view_width,y)];
      while(start<end){
        std::swap(*start,*end);
        start++;
        end--;
      }
    }
  }

  /**
    @brief Flips the raster about its diagonal axis, like a matrix tranpose.
  */
  void transpose(){
    RDLOG_WARN<<"transpose() is an experimental feature.";
    std::vector<T> new_data(view_width*view_height);
    for(xy_t y=0;y<view_height;y++)
    for(xy_t x=0;x<view_width;x++)
      std::swap(data[(i_t)x*(i_t)view_height+(i_t)y], data[xyToI(x,y)]);
    std::swap(view_width,view_height);
    //TODO: Offsets?
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
    @brief Test whether a cell lies on the boundary of the raster

    @param[in]  x   X-coordinate of cell to test
    @param[in]  y   X-coordinate of cell to test

    @return TRUE if cell lies on the raster's boundary
  */
  inline bool isEdgeCell(xy_t x, xy_t y) const {
    return x==0 || y==0 || x==view_width-1 || y==view_height-1;
  }

  ///@brief Determines whether an (x,y) pair is the top left of the DEM
  ///@return True, if the (x,y) pair is the top left of the DEM; otherwise, false
  bool isTopLeft    (xy_t x, xy_t y) const { return x==0         && y==0;          }
  ///@brief Determines whether an (x,y) pair is the top right of the DEM
  ///@return True, if the (x,y) pair is the top right of the DEM; otherwise, false
  bool isTopRight   (xy_t x, xy_t y) const { return x==width()-1 && y==0;          }
  ///@brief Determines whether an (x,y) pair is the bottom left of the DEM
  ///@return True, if the (x,y) pair is the bottom left of the DEM; otherwise, false
  bool isBottomLeft (xy_t x, xy_t y) const { return x==0         && y==height()-1; }
  ///@brief Determines whether an (x,y) pair is the bottom right of the DEM
  ///@return True, if the (x,y) pair is the bottom right of the DEM; otherwise, false
  bool isBottomRight(xy_t x, xy_t y) const { return x==width()-1 && y==height()-1; }

  ///@brief Determines whether an (x,y) pair is in the top row of the DEM
  ///@return True, if the (x,y) pair is in the top row of the DEM; otherwise, false
  bool isTopRow    (xy_t x, xy_t y) const { return y==0;          }
  ///@brief Determines whether an (x,y) pair is in the bottom row of the DEM
  ///@return True, if the (x,y) pair is in the bottom row of the DEM; otherwise, false
  bool isBottomRow (xy_t x, xy_t y) const { return y==height()-1; }
  ///@brief Determines whether an (x,y) pair is in the left column of the DEM
  ///@return True, if the (x,y) pair is in the left column of the DEM; otherwise, false
  bool isLeftCol   (xy_t x, xy_t y) const { return x==0;          }
  ///@brief Determines whether an (x,y) pair is in the right column of the DEM
  ///@return True, if the (x,y) pair is in the right column of the DEM; otherwise, false
  bool isRightCol  (xy_t x, xy_t y) const { return x==width()-1;  }

  /**
    @brief Test whether a cell lies on the boundary of the raster

    @param[in]  i   i-coordinate of cell to test

    @return TRUE if cell lies on the raster's boundary
  */
  bool isEdgeCell(i_t i) const {
    xy_t x,y;
    iToxy(i,x,y);
    return isEdgeCell(x,y);
  }

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
    for(i_t i=0;i<size();i++)
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
    data.resize(width0*height0);

    _nshift     = {{0,-1,-width0-1,-width0,-width0+1,1,width0+1,width0,width0-1}};

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
  void resize(const Array2D<U> &other, const T& val = T()){
    resize(other.width(), other.height(), val);
    geotransform       = other.geotransform;
    projection         = other.projection;
  }

  /**
    @brief Makes a raster larger and retains the raster's old data, similar to resize.

    Note: Using this command requires RAM equal to the sum of the old raster and
    the new raster. The old raster is placed in the upper-left of the new
    raster.

    @param[in] new_width  New width of the raster. Must be >= the old width.
    @param[in] new_height New height of the raster. Must be >= the old height.
    @param[in] val        Value to set the new cells to
  */
  void expand(xy_t new_width, xy_t new_height, const T val){
    RDLOG_DEBUG<<"Array2D::expand(width,height,val)";

    if(new_width==view_width && new_height==view_height)
      return;    
    if(!owned())
      throw std::runtime_error("RichDEM can only expand memory it owns!");

    if(new_width<view_width)
      throw std::runtime_error("expand(): new_width<view_width");
    if(new_height<view_height)
      throw std::runtime_error("expand(): new_height<view_height");
    
    xy_t old_width  = width();
    xy_t old_height = height();

    auto old_data = std::move(data);   //This gets the pointer to the old data before it is replaced

    resize(new_width,new_height,val);

    for(xy_t y=0;y<old_height;y++)
    for(xy_t x=0;x<old_width;x++)
      data[y*new_width+x] = old_data[y*old_width+x];
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
    @brief Return cell value based on i-coordinate

    @param[in]   i    i-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by 'i'
  */
  T& operator()(i_t i){
    assert(i>=0);
    assert(i<(i_t)view_width*view_height);
    return data[i];
  }

  /**
    @brief Return cell value based on i-coordinate

    @param[in]   i    i-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by 'i'
  */
  T operator()(i_t i) const {
    assert(i>=0);
    assert(i<(i_t)view_width*view_height);
    return data[i];
  }

  /**
    @brief Return cell value based on x,y coordinates

    @param[in]   x    X-coordinate of cell whose data should be fetched.
    @param[in]   y    Y-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by x,y
  */
  T& operator()(xy_t x, xy_t y){
    assert(x>=0);
    assert(y>=0);
    assert(x<width());
    assert(y<height());
    return data[xyToI(x,y)];
  }

  /**
    @brief Return cell value based on x,y coordinates

    @param[in]   x    X-coordinate of cell whose data should be fetched.
    @param[in]   y    Y-coordinate of cell whose data should be fetched.

    @return The value of the cell identified by x,y
  */
  T operator()(xy_t x, xy_t y) const {
    assert(x>=0);
    assert(y>=0);
    assert(x<width());
    assert(y<height());
    return data[xyToI(x,y)];
  }

  /**
    @brief Returns a copy of the top row of the raster

    @return A vector containing a copy of the top row of the raster
  */
  std::vector<T> topRow() const {  
    return getRowData(0);  
  }

  /**
    @brief Returns a copy of the bottom row of the raster

    @return A vector containing a copy of the bottom row of the raster
  */
  std::vector<T> bottomRow() const {
    return getRowData(view_height-1);
  }

  /**
    @brief Returns a copy of the left column of the raster

    Top to bottom is reoriented as left to right.

    @return A vector containing a copy of the left column of the raster
  */
  std::vector<T> leftColumn() const { 
    return getColData(0); 
  }

  /**
    @brief Returns a copy of the right column of the raster

    Top to bottom is reoriented as left to right.

    @return A vector containing a copy of the right column of the raster
  */
  std::vector<T> rightColumn() const { 
    return getColData(view_width-1);
  }

  /**
    @brief Sets an entire row of a raster to a given value.

    @param[in]   y    The row to be set
    @param[in] val    The value to set the row to
  */
  void setRow(xy_t y, const T &val){
    for(xy_t x=0;x<view_width;x++)
      data[xyToI(x,y)] = val;
  }

  /**
    @brief Sets an entire column of a raster to a given value.

    @param[in]   x    The column to be set
    @param[in] val    The value to set the column to
  */
  void setCol(xy_t x, const T &val){
    for(xy_t y=0;y<view_height;y++)
      data[xyToI(x,y)] = val;
  }

  /**
    @brief Returns a copy of an arbitrary row of the raster

    @param[in]   y    The row to retrieve

    @return A vector containing a copy of the selected row
  */
  std::vector<T> getRowData(xy_t y) const {
    return std::vector<T>(data.data()+xyToI(0,y),data.data()+xyToI(0,y)+view_width);
  }

  /**
    @brief Returns a copy of an arbitrary column of the raster

    @param[in]   x    The column to retrieve

    @return A vector containing a copy of the selected column
  */
  std::vector<T> getColData(xy_t x) const {
    std::vector<T> temp(view_height);
    for(xy_t y=0;y<view_height;y++)
      temp[y]=data[xyToI(x,y)];
    return temp;
  }

  ///Clears all raster data from RAM
  void clear(){
    data = ManagedVector<T>();
  }

  /**
    @brief Copies the geotransform, projection, and basename of another raster

    @param[in]    other    Raster to copy from
  */
  template<class U>
  void templateCopy(const Array2D<U> &other){
    geotransform       = other.geotransform;
    projection         = other.projection;
    basename           = other.basename;
    metadata     = other.metadata;
  }


  #ifdef USEGDAL
  void saveGDAL(const std::string &filename, const std::string &metadata_str="", xy_t xoffset=0, xy_t yoffset=0, bool compress=false){
    char **papszOptions = NULL;
    if(compress){
      papszOptions = CSLSetNameValue( papszOptions, "COMPRESS", "DEFLATE" );
      papszOptions = CSLSetNameValue( papszOptions, "ZLEVEL",   "6" );
    }

    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if(poDriver==NULL)
      throw std::runtime_error("Could not open GDAL driver!");
    GDALDataset *fout    = poDriver->Create(filename.c_str(), width(), height(), 1, myGDALType(), papszOptions);
    if(fout==NULL)
      throw std::runtime_error("Could not open file '"+filename+"' for GDAL save!");

    GDALRasterBand *oband = fout->GetRasterBand(1);
    oband->SetNoDataValue(no_data);

    //This could be used to copy metadata
    //poDstDS->SetMetadata( poSrcDS->GetMetadata() );

    //TIFFTAG_SOFTWARE
    //TIFFTAG_ARTIST
    {
      std::time_t the_time = std::time(nullptr);
      char time_str[64];
      std::strftime(time_str, sizeof(time_str), "%Y-%m-%d %H:%M:%S UTC", std::gmtime(&the_time));
      fout->SetMetadataItem("TIFFTAG_DATETIME",   time_str);
      fout->SetMetadataItem("TIFFTAG_SOFTWARE",   program_identifier.c_str());

      //TODO: `metadata_str` may need removing
      metadata["PROCESSING_HISTORY"] += "\n" + std::string(time_str) + " | " + program_identifier + " | ";
      if(!metadata.empty())
        metadata["PROCESSING_HISTORY"] += metadata_str;
      else
        metadata["PROCESSING_HISTORY"] += "Unspecified Operation";
    }

    for(const auto &kv: metadata)
      fout->SetMetadataItem(kv.first.c_str(), kv.second.c_str());

    //The geotransform maps each grid cell to a point in an affine-transformed
    //projection of the actual terrain. The geostransform is specified as follows:
    //    Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
    //    Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
    //In case of north up images, the GT(2) and GT(4) coefficients are zero, and
    //the GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3))
    //position is the top left corner of the top left pixel of the raster.

    if(!geotransform.empty()){
      auto out_geotransform = geotransform;

      if(out_geotransform.size()!=6)
        throw std::runtime_error("Geotransform of output is not the right size. Found "+std::to_string(out_geotransform.size())+" expected 6.");

      //We shift the top-left pixel of hte image eastward to the appropriate
      //coordinate
      out_geotransform[0] += xoffset*geotransform[1];

      //We shift the top-left pixel of the image southward to the appropriate
      //coordinate
      out_geotransform[3] += yoffset*geotransform[5];

      fout->SetGeoTransform(out_geotransform.data());
    }

    if(!projection.empty())
      fout->SetProjection(projection.c_str());

    #ifdef DEBUG
      RDLOG_DEBUG<<"Filename: "<<std::setw(20)<<filename<<" Xoffset: "<<std::setw(6)<<xoffset<<" Yoffset: "<<std::setw(6)<<yoffset<<" Geotrans0: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotransform[0]<<" Geotrans3: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotransform[3];
    #endif

    auto temp = oband->RasterIO(GF_Write, 0, 0, view_width, view_height, data.data(), view_width, view_height, myGDALType(), 0, 0);
    if(temp!=CE_None)
      throw std::runtime_error("Error writing file with saveGDAL()!");

    GDALClose(fout);
  }
  #endif



  /**
    @brief Output a square of cells useful for determining raster orientation.

    This method prints out a square block of cells whose upper-left corner is
    the (integer-division) center of the raster.

    Stamps are only shown if the SHOW_STAMPS preprocessor variable is set.

    Since algorithms may have to flip rasters horizontally or vertically before
    manipulating them, it is important that all algorithms work on data in the
    same orientation. This method, used in testing, helps a user ensure that 
    their algorithm is orientating data correctly.

    @param[in]  size   Output stamp will be size x size
    @param[in]  msg    Message to print prior to the stamp

  */
  void printStamp(int size, std::string msg="") const {
    #ifdef SHOW_STAMPS
      const xy_t sx = width()/2;
      const xy_t sy = height()/2;

      if(msg.size()>0)
        std::cout<<msg<<std::endl;
      std::cout<<"Stamp for basename='"<<basename
               <<"', filename='"<<filename
               #ifdef USEGDAL
               <<"', dtype="<<GDALGetDataTypeName(myGDALType())
               #endif
               <<" at "<<sx<<","<<sy<<"\n";

      const xy_t sxmax = std::min(width(), sx+size);
      const xy_t symax = std::min(height(),sy+size);

      for(xy_t y=sy;y<symax;y++){
        for(xy_t x=sx;x<sxmax;x++)
          std::cout<<std::setw(5)<<std::setprecision(3)<<(int)data[xyToI(x,y)]<<" ";
        std::cout<<"\n";
      }
    #endif
  }


  /**
    @brief Prints a square of cells centered at x,y. Useful for debugging.

    @param[in]  radius   Output stamp will be 2*radius x 2*radius
    @param[in]      x0   X-coordinate of block center
    @param[in]      y0   Y-coordinate of block center
    @param[in]   color   Print the (x,y) cell in colour?
    @param[in]     msg   Optional message to print above the block
  */
  void printBlock(const int radius, const xy_t x0, const xy_t y0, bool color=false, const std::string msg="") const {
    if(msg.size()!=0)
      std::cout<<msg<<std::endl;

    xy_t xmin = std::max(0,x0-radius);
    xy_t ymin = std::max(0,y0-radius);
    xy_t xmax = std::min(width(),x0+radius);
    xy_t ymax = std::min(height(),y0+radius);

    for(xy_t y=ymin;y<ymax;y++){
      for(xy_t x=xmin;x<xmax;x++){
        if(color && x==x0 && y==y0)
          std::cout<<"\033[92m";
        std::cout<<std::setw(5)<<(int)data[xyToI(x,y)]<<" ";
        if(color && x==x0 && y==y0)
          std::cout<<"\033[39m";
      }
      std::cout<<std::endl;
    }
  }

  /**
    @brief Prints the entire array

    @param[in]     msg   Optional message to print above the block
  */
  void printAll(const std::string msg="") const {
    if(!msg.empty())
      std::cout<<msg<<std::endl;

    for(xy_t y=0;y<height();y++){
      for(xy_t x=0;x<width();x++)
        std::cout<<std::setw(5)<<data[xyToI(x,y)]<<" ";
      std::cout<<std::endl;
    }
  }

  /**
    @brief Get the area of an individual cell in square projection units

    @return The area of the cell in square projection units
  */
  double getCellArea() const {
    assert(geotransform.size()>0);
    return std::abs(geotransform[1]*geotransform[5]);
  }

  /**
    @brief Get the length of a cell along the raster's horizontal axis
    @return The length of the cell along the raster's horizontal axis
  */
  double getCellLengthX() const {
    assert(geotransform.size()>0);
    return std::abs(geotransform[1]);
  }

  /**
    @brief Get the length of a cell along the raster's horizontal axis
    @return The length of the cell along the raster's horizontal axis
  */
  double getCellLengthY() const {
    assert(geotransform.size()>0);
    return std::abs(geotransform[5]);
  }

  /**
    @brief Multiplies the entire array by a scalar

    @param[in]     x     Value to multiply array by
  */
  void scale(const double x) {
    for(i_t i=0;i<size();i++)
      if(data[i]!=no_data)
        data[i] *= x;
  }

  //TODO
  inline bool owned() const {
    return data.owned();
  }
};

}

#endif
