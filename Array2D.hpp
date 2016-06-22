#ifndef _array_2d_hpp_
#define _array_2d_hpp_

#include "gdal_priv.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <typeinfo>
#include <stdexcept>
#include "Layoutfile.hpp"
#include "lru.hpp"

//These enable compression in the loadNative() and saveNative() methods
#ifdef WITH_COMPRESSION
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#endif

GDALDataType peekGDALType(const std::string &filename) {
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  assert(fin!=NULL);

  GDALRasterBand *band   = fin->GetRasterBand(1);
  GDALDataType data_type = band->GetRasterDataType();

  GDALClose(fin);

  return data_type;
}

GDALDataType peekLayoutType(const std::string &layout_filename) {
  LayoutfileReader lf(layout_filename);

  while(lf.next()){
    if(lf.getFilename().size()==0)
      continue;

    GDALAllRegister();
    std::string tile_path = lf.getPath()+lf.getFilename();
    GDALDataset *fin = (GDALDataset*)GDALOpen(tile_path.c_str(), GA_ReadOnly);
    if(fin==NULL)
      std::cerr<<"Could not open file'"<<(lf.getPath()+lf.getFilename())<<"' to determine layout type."<<std::endl;

    GDALRasterBand *band   = fin->GetRasterBand(1);
    GDALDataType data_type = band->GetRasterDataType();

    GDALClose(fin);

    return data_type;
  }

  throw std::runtime_error("Empty layout file!");
}

template<class T>
T getGDALnodata(const std::string filename){
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  assert(fin!=NULL);

  GDALRasterBand *band = fin->GetRasterBand(1);
  return band->GetNoDataValue();
}

//Get the dimensions of a GDAL file
int getGDALDimensions(
  const std::string &filename,
  int &height,
  int &width,
  GDALDataType &dtype,
  double geotransform[6]
){
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  assert(fin!=NULL);

  GDALRasterBand *band = fin->GetRasterBand(1);
  
  dtype   = band->GetRasterDataType();

  if(geotransform!=NULL && fin->GetGeoTransform(geotransform)!=CE_None){
    std::cerr<<"Error getting geotransform from '"<<filename<<"'!"<<std::endl;
    return -1; //TODO: Set error code
  }

  height  = band->GetYSize();
  width   = band->GetXSize();

  GDALClose(fin);

  return 0;
}


template<class T>
GDALDataType NativeTypeToGDAL() {
  if(typeid(T)==typeid(uint8_t))
    return GDT_Byte;
  else if(typeid(T)==typeid(uint16_t))
    return GDT_UInt16;
  else if(typeid(T)==typeid(int16_t))
    return GDT_Int16;
  else if(typeid(T)==typeid(uint32_t))
    return GDT_UInt32;
  else if(typeid(T)==typeid(int32_t))
    return GDT_Int32;
  else if(typeid(T)==typeid(float))
    return GDT_Float32;
  else if(typeid(T)==typeid(double))
    return GDT_Float64;
  else 
    std::cerr<<"Could not map native type "<<typeid(T).name()<<" to GDAL type! (Use `c++filt -t` to decode.)"<<std::endl;
    throw std::runtime_error("Could not map native data type to GDAL type!");
  return GDT_Unknown;
}


template<class T>
class Array2D {
 public:
  typedef std::vector<T>   Row;
  std::string filename;
  std::vector<double> geotrans;
  std::string projection;

 private:
  typedef std::vector<Row> InternalArray;
  InternalArray data;

  GDALDataType data_type;

  static const int HEADER_SIZE = 7*sizeof(int) + sizeof(T);

  int total_height;
  int total_width;
  int view_height;
  int view_width;
  int view_xoff;
  int view_yoff;
  int num_data_cells = -1;
  bool file_native;

  T   no_data;

  void loadGDAL(const std::string &filename, int xOffset=0, int yOffset=0, int part_width=0, int part_height=0, bool exact=false, bool load_data=true){
    assert(empty());
    assert(xOffset>=0);
    assert(yOffset>=0);

    this->filename = filename;

    GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
    assert(fin!=NULL);

    if(fin->GetGeoTransform(geotrans.data())!=CE_None){
      std::cerr<<"Error getting geotransform from '"<<filename<<"'!"<<std::endl;
      throw std::runtime_error("Error getting geotransform!");
    }

    const char* projection_string=fin->GetProjectionRef();
    projection = std::string(projection_string);

    GDALRasterBand *band = fin->GetRasterBand(1);
    data_type            = band->GetRasterDataType();

    total_width  = band->GetXSize();
    total_height = band->GetYSize();
    no_data      = band->GetNoDataValue();

    if(exact && xOffset==0 && yOffset==0 && (part_width!=total_width || part_height!=total_height))
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

    //std::cerr<<"Allocating: "<<view_height<<" rows by "<<view_width<<" columns"<<std::endl;
    if(load_data)
      loadData();
  }

  GDALDataType myGDALType() const {
    return NativeTypeToGDAL<T>();
  }

 public:
  Array2D(){
    GDALAllRegister();
    total_height = 0;
    total_width  = 0;
    view_width   = 0;
    view_height  = 0;
    view_xoff    = 0;
    view_yoff    = 0;
    geotrans.resize(6);
  }

  //Create an internal array
  Array2D(int width, int height, const T& val = T()) : Array2D() {
    resize(width,height,val);
  }

  //Create internal array from a file
  Array2D(const std::string &filename, bool native, int xOffset=0, int yOffset=0, int part_width=0, int part_height=0, bool exact=false, bool load_data=true) : Array2D() {
    if(native)
      loadNative(filename, load_data);
    else
      loadGDAL(filename, xOffset, yOffset, part_width, part_height, exact, load_data);
  }

  void saveNative(const std::string &filename){
    std::fstream fout;

    std::cerr<<"Saving native to: "<<filename<<std::endl;

    fout.open(filename, std::ios_base::binary | std::ios_base::out | std::ios::trunc);
    if(!fout.good()){
      std::cerr<<"Failed to open file '"<<filename<<"'."<<std::endl;
      throw std::logic_error("Failed to open a file!");
    }

    #ifdef WITH_COMPRESSION
      boost::iostreams::filtering_ostream out;
      out.push(boost::iostreams::zlib_compressor());
      out.push(fout);
    #else
      auto &out = fout;
    #endif

    out.write(reinterpret_cast<char*>(&total_height),   sizeof(int));
    out.write(reinterpret_cast<char*>(&total_width),    sizeof(int));
    out.write(reinterpret_cast<char*>(&view_height),    sizeof(int));
    out.write(reinterpret_cast<char*>(&view_width),     sizeof(int));
    out.write(reinterpret_cast<char*>(&view_xoff),      sizeof(int));
    out.write(reinterpret_cast<char*>(&view_yoff),      sizeof(int));
    out.write(reinterpret_cast<char*>(&num_data_cells), sizeof(int));
    out.write(reinterpret_cast<char*>(&no_data),        sizeof(T  ));

    for(int y=0;y<view_height;y++)
      out.write(reinterpret_cast<char*>(data[y].data()), view_width*sizeof(T));
  }

  void setFilename(const std::string &filename){
    this->filename = filename;
  }

  void dumpData() {
    saveNative(filename);
    clear();
  }

  void loadData() {
    if(!data.empty())
      return;

    if(file_native){
      std::ifstream fin(filename, std::ios::in | std::ios::binary);
      assert(fin.good());

      #ifdef WITH_COMPRESSION
        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::zlib_decompressor());
        in.push(fin);
      #else
        auto &in = fin;
      #endif

      data = InternalArray(view_height, Row(view_width));

      in.seekg(7*sizeof(int)+sizeof(T));

      for(int y=0;y<view_height;y++)
        in.read(reinterpret_cast<char*>(data[y].data()), view_width*sizeof(T));
    } else {
      GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
      if(fin==NULL){
        std::cerr<<"Failed to loadData() into tile from '"<<filename<<"'"<<std::endl;
        throw std::runtime_error("Failed to loadData() into tile.");
      }

      GDALRasterBand *band = fin->GetRasterBand(1);

      data = InternalArray(view_height, Row(view_width));
      for(int y=view_yoff;y<view_yoff+view_height;y++){
        auto temp = band->RasterIO( GF_Read, view_xoff, y, view_width, 1, data[y-view_yoff].data(), view_width, 1, data_type, 0, 0 ); //TODO: Check for success
        if(temp!=CE_None)
          throw std::runtime_error("Error reading file with GDAL!");
      }

      GDALClose(fin);
    }
  }

  void loadNative(const std::string &filename, bool load_data=true){
    std::ifstream fin(filename, std::ios::in | std::ios::binary);
    assert(fin.good());

    this->filename = filename;
    file_native    = true;

    #ifdef WITH_COMPRESSION
      boost::iostreams::filtering_istream in;
      in.push(boost::iostreams::zlib_decompressor());
      in.push(fin);
    #else
      auto &in = fin;
    #endif

    in.read(reinterpret_cast<char*>(&total_height),   sizeof(int));
    in.read(reinterpret_cast<char*>(&total_width),    sizeof(int));
    in.read(reinterpret_cast<char*>(&view_height),    sizeof(int));
    in.read(reinterpret_cast<char*>(&view_width),     sizeof(int));
    in.read(reinterpret_cast<char*>(&view_xoff),      sizeof(int));
    in.read(reinterpret_cast<char*>(&view_yoff),      sizeof(int));
    in.read(reinterpret_cast<char*>(&num_data_cells), sizeof(int));
    in.read(reinterpret_cast<char*>(&no_data),        sizeof(T  ));

    if(load_data)
      loadData();
  }

  //Note: The following functions return signed integers, which make them
  //generally easier to work with. If your DEM has a dimension which exceeds
  //2147483647, some other modifications to this program will probably be
  //necessary.
  int  totalWidth () const { return total_width;    }
  int  totalHeight() const { return total_height;   }
  int  viewWidth  () const { return view_width;     }
  int  viewHeight () const { return view_height;    }
  int  viewXoff   () const { return view_xoff;      }
  int  viewYoff   () const { return view_yoff;      }
  bool empty      () const { return data.empty();   }
  T    noData     () const { return no_data;        }

  bool operator==(const Array2D<T> &o){
    if(viewWidth()!=o.viewWidth() || viewHeight()!=o.viewHeight())
      return false;
    if(noData()!=o.noData())
      return false;
    for(int y=0;y<viewHeight();y++)
    for(int x=0;x<viewWidth();x++)
      if(data[y][x]!=o.data[y][x])
        return false;
    return true;
  }

  bool isNoData(int x, int y) const {
    return data[y][x]==no_data;
  }

  void flipVert(){
    std::reverse(data.begin(),data.end());
  }

  void flipHorz(){
    for(auto &row: data)
      std::reverse(row.begin(),row.end());
  }

  bool in_grid(int x, int y) const {
    return 0<=x && x<viewWidth() && 0<=y && y<viewHeight();
  }

  void setNoData(const T &ndval){
    no_data = ndval;
  }

  void setAll(const T &val){
    for(auto &row: data)
      std::fill(row.begin(),row.end(),val);
  }

  void init(T val){
    setAll(val);
  }

  //Destructively resizes the array. All data will die!
  void resize(int width, int height, const T& val = T()){
    data         = InternalArray(height, Row(width, val));
    total_height = view_height = height;
    total_width  = view_width  = width;
  }

  void countDataCells(){
    num_data_cells = 0;
    for(int y=0;y<viewHeight();y++)
    for(int x=0;x<viewWidth();x++)
      if(data[y][x]!=no_data)
        num_data_cells++;
  }

  int numDataCells(){
    if(num_data_cells==-1)
      countDataCells();
    return num_data_cells;
  }

  T& operator()(int x, int y){
    assert(x>=0);
    assert(y>=0);
    //std::cerr<<"Width: "<<viewWidth()<<" Height: "<<viewHeight()<<" x: "<<x<<" y: "<<y<<std::endl;
    assert(x<viewWidth());
    assert(y<viewHeight());
    return data[y][x];
  }

  const T& operator()(int x, int y) const {
    assert(x>=0);
    assert(y>=0);
    assert(x<viewWidth());
    assert(y<viewHeight());
    return data[y][x];
  }

  Row&       topRow   ()       { return data.front(); }
  Row&       bottomRow()       { return data.back (); }
  const Row& topRow   () const { return data.front(); }
  const Row& bottomRow() const { return data.back (); }

  Row leftColumn() const {
    Row temp(data.size());
    for(size_t y=0;y<data.size();y++)
      temp[y] = data[y][0];
    return temp;
  }

  Row rightColumn() const {
    Row temp(data.size());
    size_t right = data[0].size()-1;
    for(size_t y=0;y<data.size();y++)
      temp[y] = data[y][right];
    return temp;
  }

  void emplaceRow(Row &row){
    data.emplace_back(row);
  }

  Row& rowRef(int rownum){
    return data[rownum];
  }

  void setRow(int rownum, const T &val){
    std::fill(data[rownum].begin(),data[rownum].end(),val);
  }

  void setCol(int colnum, const T &val){
    for(int y=0;y<viewHeight();y++)
      data[y][colnum] = val;
  }

  const std::vector<T>& getRowData(int rownum){
    return data[rownum].data();
  }

  void clear(){
    data.clear();
    data.shrink_to_fit();
  }

  template<class U>
  void templateCopy(const Array2D<U> &other){
    geotrans   = other.geotrans;
    projection = other.projection;
  }

  void saveGDAL(const std::string &filename, int xoffset, int yoffset){
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if(poDriver==NULL){
      std::cerr<<"Could not load GTiff driver!"<<std::endl;
      throw std::runtime_error("Could not load GTiff driver!");
    }
    GDALDataset *fout    = poDriver->Create(filename.c_str(), viewWidth(), viewHeight(), 1, myGDALType(), NULL);
    if(fout==NULL){
      std::cerr<<"Could not open output file!"<<std::endl;
      throw std::runtime_error("Could not open output file!");
    }

    GDALRasterBand *oband = fout->GetRasterBand(1);
    oband->SetNoDataValue(no_data);

    //The geotransform maps each grid cell to a point in an affine-transformed
    //projection of the actual terrain. The geostransform is specified as follows:
    //    Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
    //    Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
    //In case of north up images, the GT(2) and GT(4) coefficients are zero, and
    //the GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3))
    //position is the top left corner of the top left pixel of the raster.

    //We shift the top-left pixel of hte image eastward to the appropriate
    //coordinate
    geotrans[0] += xoffset*geotrans[1];

    //We shift the top-left pixel of the image southward to the appropriate
    //coordinate
    geotrans[3] += yoffset*geotrans[5];

    fout->SetGeoTransform(geotrans.data());
    fout->SetProjection(projection.c_str());

    #ifdef DEBUG
      std::cerr<<"Filename: "<<std::setw(20)<<filename<<" Xoffset: "<<std::setw(6)<<xoffset<<" Yoffset: "<<std::setw(6)<<yoffset<<" Geotrans0: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotrans[0]<<" Geotrans3: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotrans[3]<< std::endl;
    #endif

    for(int y=0;y<view_height;y++){
      auto temp = oband->RasterIO(GF_Write, 0, y, viewWidth(), 1, data[y].data(), viewWidth(), 1, myGDALType(), 0, 0); //TODO: Check for success
      if(temp!=CE_None)
        std::cerr<<"Error writing file! Continuing in the hopes that some work can be salvaged."<<std::endl;
    }

    GDALClose(fout);
  }
};

















template<class T>
class A2Array2D {
 private:
  template<typename U> friend class A2Array2D;

  class WrappedArray2D : public Array2D<T> {
   public:
    using Array2D<T>::Array2D;
    bool null_tile = false;
    bool loaded    = false;
    bool created   = true;
  };
  std::vector< std::vector< WrappedArray2D > > data;

  LRU< WrappedArray2D* > lru;

  int32_t total_width_in_cells  = 0;
  int32_t total_height_in_cells = 0;
  int32_t per_tile_width        = -1;
  int32_t per_tile_height       = -1;
  T       no_data;

  bool readonly = true;

  void LoadTile(int tile_x, int tile_y){
    if(data[tile_y][tile_x].loaded){
      lru.insert(&data[tile_y][tile_x]);
      return;
    }

    if(lru.full()){
      auto tile_to_unload = lru.back();

      if(readonly)
        tile_to_unload->clear();
      else
        tile_to_unload->dumpData();

      tile_to_unload->loaded = false;
      std::cout<<"Loading "<<tile_x<<","<<tile_y<<std::endl;
      lru.pop_back();
    }

    if(data[tile_y][tile_x].created){
      data[tile_y][tile_x].loadData();
    } else {
      data[tile_y][tile_x].resize(per_tile_width,per_tile_height);
      data[tile_y][tile_x].created = true;
    }
    data[tile_y][tile_x].loaded = true;
    lru.insert(&data[tile_y][tile_x]);
  }

 public:

  A2Array2D(std::string layoutfile, int cachesize){
    lru.setCapacity(cachesize);
    readonly = true;

    long    cell_count     = 0;
    int     not_null_tiles = 0;
    std::vector<double> chunk_geotransform(6);

    LayoutfileReader lf(layoutfile);
    while(lf.next()){
      if(lf.newRow()) //Add a row to the grid of chunks
        data.emplace_back();

      if(lf.getFilename().size()==0){
        data.back().emplace_back();
        data.back().back().null_tile = true;
        continue;
      }

      not_null_tiles++;

      if(per_tile_height==-1){
        //Retrieve information about this chunk. All chunks must have the same
        //dimensions, which we could check here, but opening and closing
        //thousands of files is expensive. Therefore, we rely on the user to
        //check this beforehand if they want to. We will, however, verify that
        //things are correct in Consumer() as we open the files for reading.
        GDALDataType file_type;
        std::cerr<<(lf.getPath() + lf.getFilename())<<std::endl;
        if(getGDALDimensions(
            lf.getPath() + lf.getFilename(),
            per_tile_height,
            per_tile_width,
            file_type,
            chunk_geotransform.data()
        )!=0){
          std::cerr<<"Error getting file information from '"<<(lf.getPath() + lf.getFilename())<<"'!"<<std::endl;
          throw std::runtime_error("Error get information from a DEM file!");
        }

        no_data = getGDALnodata<T>(lf.getPath() + lf.getFilename());
        std::cerr<<"NoData: "<<no_data<<std::endl;
      }

      cell_count += per_tile_width*per_tile_height;

      data.back().emplace_back(
        lf.getPath()+lf.getFilename(),
        false,
        0,
        0,
        0,
        0,
        false,
        false
      );

      if(no_data!=data.back().back().noData()){
        std::cerr<<"NoData found: "<<data.back().back().noData()<<" Expected: "<<no_data<<std::endl;
        throw std::runtime_error("Tiles did not all have the same NoData value!");
      }
      if(per_tile_width!=data.back().back().viewWidth())
        throw std::runtime_error("Tiles were not all the same width!");
      if(per_tile_height!=data.back().back().viewHeight())
        throw std::runtime_error("Tiles were not all the same width!");

      if(lf.getY()==0)
        total_width_in_cells  += data.back().back().viewWidth();
      if(lf.newRow())
        total_height_in_cells += data.back().back().viewHeight();
    }

    std::cerr<<"TWIC: "<<total_width_in_cells<<std::endl;
  }

  A2Array2D(std::string prefix, int per_tile_width, int per_tile_height, int width, int height, int cachesize){
    lru.setCapacity(cachesize);

    readonly = false;

    this->per_tile_width  = per_tile_width;
    this->per_tile_height = per_tile_height;

    this->total_width_in_cells  = per_tile_width*width;
    this->total_height_in_cells = per_tile_height*height;

    int tile=0;
    for(int y=0;y<height;y++){
      data.emplace_back();
      for(int x=0;x<width;x++){
        tile++;
        data.back().emplace_back();
        data.back().back().setFilename(prefix+std::to_string(tile)+".native");
        data.back().back().created  = false;
      }
    }
  }

  template<class U>
  A2Array2D(std::string prefix, const A2Array2D<U> &other, int cachesize) : A2Array2D(prefix, other.tileWidth(), other.tileHeight(), other.widthInTiles(), other.heightInTiles(), cachesize) {
    for(int y=0;y<heightInTiles();y++)
    for(int x=0;x<widthInTiles();x++)
      data[y][x].templateCopy(other.data[y][x]);
  }

  T& getn(int tx, int ty, int x, int y, int dx, int dy){
    x += dx;
    y += dy;
    if(x<0){
      tx--;
      x = per_tile_width-1;
    } else if (x==per_tile_width) {
      tx++;
      x=0;
    }
    if(y<0){
      ty--;
      y = per_tile_height-1;
    } else if (y==per_tile_height){
      ty++;
      y=0;
    }
    assert(x>=0);
    assert(y>=0);
    assert(x<per_tile_width);
    assert(y<per_tile_height);
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<widthInTiles());
    assert(ty<heightInTiles());
    return data[ty][tx](x,y);
  }

  T& operator()(int tx, int ty, int x, int y){
    assert(x>=0);
    assert(y>=0);
    assert(x<per_tile_width);
    assert(y<per_tile_height);
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<widthInTiles());
    assert(ty<heightInTiles());
    LoadTile(tx, ty);
    return data[ty][tx](x,y);
  }

  const T& operator()(int tx, int ty, int x, int y) const {
    assert(x>=0);
    assert(y>=0);
    assert(x<per_tile_width);
    assert(y<per_tile_height);
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<widthInTiles());
    assert(ty<heightInTiles());
    LoadTile(tx, ty);
    return data[ty][tx](x,y);
  }

  T& operator()(int x, int y){
    assert(x>=0);
    assert(y>=0);
    assert(x<total_width_in_cells);
    assert(y<total_height_in_cells);
    int tile_x = x/per_tile_width;
    int tile_y = y/per_tile_height;
    x          = x%per_tile_width;
    y          = y%per_tile_height;
    LoadTile(tile_x, tile_y);
    return data[tile_y][tile_x](x,y);
  }

  const T& operator()(int x, int y) const {
    assert(x>=0);
    assert(y>=0);
    assert(x<total_width_in_cells);
    assert(y<total_height_in_cells);
    int tile_x = x/per_tile_width;
    int tile_y = y/per_tile_height;
    x          = x%per_tile_width;
    y          = y%per_tile_height;
    LoadTile(tile_x, tile_y);
    return data[tile_y][tile_x](x,y);
  }

  int width() const {
    return total_width_in_cells;
  }

  int height() const {
    return total_height_in_cells;
  }

  int widthInTiles() const {
    return data.back().size();
  }

  int heightInTiles() const {
    return data.size();
  }

  int tileWidth() const {
    return per_tile_width;
  }

  int tileHeight() const {
    return per_tile_height;
  }

  void setAll(const T &val){
    for(auto &row: data)
    for(auto &tile: row)
      tile.setAll(val);
  }

  T noData() const {
    return no_data;
  }

  bool in_grid(int x, int y) const {
    return (x>=0 && y>=0 && x<total_width_in_cells && y<total_height_in_cells);
  }

  bool isReadonly() const {
    return readonly;
  }

  void saveGDAL(std::string prefix) {
    int tile_i = 0;
    for(auto &row: data)
    for(auto &tile: row){
      std::cerr<<"Trying to save: "<<(prefix+std::to_string(tile_i)+".tif")<<std::endl;
      tile_i++;
      if(!tile.loaded)
        tile.loadData();
      tile.saveGDAL(prefix+std::to_string(tile_i)+".tif", 0, 0);
      tile.clear();
    }
  }
};

#endif