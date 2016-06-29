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
#include "common.hpp"

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
      std::cerr<<"Could not open '"<<(lf.getPath()+lf.getFilename())<<"' to determine layout type."<<std::endl;

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
  
  dtype = band->GetRasterDataType();

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
  else {
    std::cerr<<"Could not map native type '"<<typeid(T).name()<<"' to GDAL type! (Use `c++filt -t` to decode.)"<<std::endl;
    throw std::runtime_error("Could not map native data type to GDAL type!");
  }
  return GDT_Unknown;
}


int peekLayoutTileSize(const std::string &layout_filename) {
  LayoutfileReader lf(layout_filename);

  GDALAllRegister();
  while(lf.next()){
    if(lf.getFilename().size()==0)
      continue;

    int height;
    int width;
    GDALDataType dtype;
    getGDALDimensions(lf.getPath()+lf.getFilename(),height,width,dtype,NULL);
    return width*height;
  }

  throw std::runtime_error("Empty layout file!");
}


template<class T>
class Array2D {
 public:
  std::string filename;
  std::string basename;
  std::vector<double> geotransform;
  std::string projection;

 private:
  template<typename U> friend class Array2D;

  std::vector<T> data;

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

    file_native = false;

    this->filename = filename;

    GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
    if(fin==NULL){
      std::cerr<<"Could not open file '"<<filename<<"'!"<<std::endl;
      throw std::runtime_error("Could not open a GDAL file!");
    }

    geotransform.resize(6);
    if(fin->GetGeoTransform(geotransform.data())!=CE_None){
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
  }

  //Create an internal array
  Array2D(int width, int height, const T& val = T()) : Array2D() {
    resize(width,height,val);
  }

  template<class U>
  Array2D(const Array2D<U> &other, const T& val=T()) : Array2D() {
    total_height = other.total_height;
    total_width  = other.total_width;
    view_width   = other.view_width;
    view_height  = other.view_height;
    view_xoff    = other.view_xoff;
    view_yoff    = other.view_yoff;
    geotransform = other.geotransform;
    projection   = other.projection;
    basename     = other.basename;
    resize(other.viewWidth(), other.viewHeight(), val);
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

    file_native = true;

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

    out.write(reinterpret_cast<char*>(data.data()), data.size()*sizeof(T));
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

      data.resize(view_height*view_width);

      in.seekg(7*sizeof(int)+sizeof(T));

      in.read(reinterpret_cast<char*>(data.data()), view_width*view_height*sizeof(T));
    } else {
      GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
      if(fin==NULL){
        std::cerr<<"Failed to loadData() into tile from '"<<filename<<"'"<<std::endl;
        throw std::runtime_error("Failed to loadData() into tile.");
      }

      GDALRasterBand *band = fin->GetRasterBand(1);

      data.resize(view_width*view_height);
      auto temp = band->RasterIO( GF_Read, view_xoff, view_yoff, view_width, view_height, data.data(), view_width, view_height, data_type, 0, 0 );
      if(temp!=CE_None)
        throw std::runtime_error("Error reading file with GDAL!");

      std::cerr<<"Geotrans: ";
      for(auto x: geotransform)
        std::cerr<<x<<" ";
      std::cerr<<std::endl;

      if(geotransform[0]<0)
        flipHorz();
      if(geotransform[5]<0){
        std::cerr<<"Flip vert on load"<<std::endl;
        flipVert();
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

  T min() const {
    T minval = std::numeric_limits<T>::max();
    for(auto const x: data){
      if(x==no_data)
        continue;
      minval = std::min(minval,x);
    }
    return minval;
  }

  T max() const {
    T maxval = std::numeric_limits<T>::min();
    for(auto const x: data){
      if(x==no_data)
        continue;
      maxval = std::max(maxval,x);
    }
    return maxval;
  }

  void replace(const T oldval, const T newval){
    for(auto &x: data)
      if(x==oldval)
        x = newval;
  }

  int countval(const T val) const {
    int count=0;
    for(const auto x: data)
      if(x==val)
        count++;
    return count;
  }

  bool operator==(const Array2D<T> &o){
    if(viewWidth()!=o.viewWidth() || viewHeight()!=o.viewHeight())
      return false;
    if(noData()!=o.noData())
      return false;
    return data==o.data;
  }

  bool isNoData(int x, int y) const {
    return data[y*view_width+x]==no_data;
  }

  void flipVert(){
    for(int y=0;y<view_height/2;y++)
      std::swap_ranges(
        data.begin()+(y+0)*view_width,
        data.begin()+(y+1)*view_width,
        data.begin()+(view_height-1-y)*view_width
      );
  }

  void flipHorz(){
    for(int y=0;y<view_height;y++)
      std::reverse(data.begin()+y*view_width,data.begin()+(y+1)*view_width);
  }

  bool in_grid(int x, int y) const {
    return 0<=x && x<viewWidth() && 0<=y && y<viewHeight();
  }

  void setNoData(const T &ndval){
    no_data = ndval;
  }

  void setAll(const T &val){
    std::fill(data.begin(),data.end(),val);
  }

  void init(T val){
    setAll(val);
  }

  //Destructively resizes the array. All data will die!
  void resize(int width, int height, const T& val = T()){
    data.resize(width*height);
    setAll(val);
    total_height = view_height = height;
    total_width  = view_width  = width;
  }

  void countDataCells(){
    num_data_cells = 0;
    for(const auto x: data)
      if(x!=no_data)
        num_data_cells++;
  }

  int numDataCells(){
    if(num_data_cells==-1)
      countDataCells();
    return num_data_cells;
  }

  T& operator()(int i){
    assert(i>=0);
    assert(i<view_width*view_height);
    return data[i];
  }

  const T& operator()(int i) const {
    assert(i>=0);
    assert(i<view_width*view_height);
    return data[i];
  }

  int getN(int i, int n) const {
    int x = i%view_width+dx[n];
    int y = i/view_width+dy[n];
    if(x<0 || y<0 || x==view_width || y==view_height)
      return -1;
    return y*view_width+x;
  }

  T& operator()(int x, int y){
    assert(x>=0);
    assert(y>=0);
    //std::cerr<<"Width: "<<viewWidth()<<" Height: "<<viewHeight()<<" x: "<<x<<" y: "<<y<<std::endl;
    assert(x<viewWidth());
    assert(y<viewHeight());
    return data[y*view_width+x];
  }

  const T& operator()(int x, int y) const {
    assert(x>=0);
    assert(y>=0);
    assert(x<viewWidth());
    assert(y<viewHeight());
    return data[y*view_width+x];
  }

  // Row&       topRow   ()       { return data.front(); } //TODO
  // Row&       bottomRow()       { return data.back (); } //TODO
  // const Row& topRow   () const { return data.front(); } //TODO
  // const Row& bottomRow() const { return data.back (); } //TODO

  // Row leftColumn() const { //TODO
  //   Row temp(data.size());
  //   for(size_t y=0;y<data.size();y++)
  //     temp[y] = data[y][0];
  //   return temp;
  // }

  // Row rightColumn() const { //TODO
  //   Row temp(data.size());
  //   size_t right = data[0].size()-1;
  //   for(size_t y=0;y<data.size();y++)
  //     temp[y] = data[y][right];
  //   return temp;
  // }

  // Row& rowRef(int rownum){
  //   return data[rownum];
  // }

  // void setRow(int rownum, const T &val){
  //   std::fill(data[rownum].begin(),data[rownum].end(),val);
  // }

  // void setCol(int colnum, const T &val){
  //   for(int y=0;y<viewHeight();y++)
  //     data[y][colnum] = val;
  // }

  // const std::vector<T>& getRowData(int rownum){
  //   return data[rownum].data();
  // }

  void clear(){
    data.clear();
    data.shrink_to_fit();
  }

  template<class U>
  void templateCopy(const Array2D<U> &other){
    geotransform = other.geotransform;
    projection   = other.projection;
    basename     = other.basename;
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

    auto out_geotransform = geotransform;

    //We shift the top-left pixel of hte image eastward to the appropriate
    //coordinate
    out_geotransform[0] += xoffset*geotransform[1];

    //We shift the top-left pixel of the image southward to the appropriate
    //coordinate
    out_geotransform[3] += yoffset*geotransform[5];

    fout->SetGeoTransform(out_geotransform.data());
    fout->SetProjection(projection.c_str());

    #ifdef DEBUG
      std::cerr<<"Filename: "<<std::setw(20)<<filename<<" Xoffset: "<<std::setw(6)<<xoffset<<" Yoffset: "<<std::setw(6)<<yoffset<<" Geotrans0: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotrans[0]<<" Geotrans3: "<<std::setw(10)<<std::setprecision(10)<<std::fixed<<geotrans[3]<< std::endl;
    #endif

    auto temp = oband->RasterIO(GF_Write, 0, 0, view_width, view_height, data.data(), view_width, view_height, myGDALType(), 0, 0);
    if(temp!=CE_None)
      std::cerr<<"Error writing file! Continuing in the hopes that some work can be salvaged."<<std::endl;

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
    bool null_tile  = false;
    bool loaded     = false;
    bool created    = true;
    bool do_set_all = false; //If true, then set all to 'set_all_val' when tile is loaded
    T set_all_val   = 0;
    void lazySetAll(){
      if(do_set_all){
        do_set_all = false;
        setAll(set_all_val);
      }
    }
  };
  std::vector< std::vector< WrappedArray2D > > data;

  LRU< WrappedArray2D* > lru;

  int32_t total_width_in_cells    = 0;
  int32_t total_height_in_cells   = 0;
  int32_t per_tile_width          = -1;
  int32_t per_tile_height         = -1;
  int32_t evictions               = 0;
  int64_t cells_in_not_null_tiles = 0;
  T       no_data_to_set; //Used to disguise null tiles

  bool readonly = true;

  void LoadTile(int tile_x, int tile_y){
    auto& tile = data[tile_y][tile_x];
    if(tile.null_tile)
      return;

    if(tile.loaded){
      lru.insert(&data[tile_y][tile_x]);
      tile.lazySetAll();
      return;
    }

    if(lru.full()){
      auto tile_to_unload = lru.back();

      if(readonly)
        tile_to_unload->clear();
      else
        tile_to_unload->dumpData();

      evictions++;

      tile_to_unload->loaded = false;
      lru.pop_back();
    }

    if(tile.created){
      tile.loadData();
    } else {
      tile.resize(per_tile_width,per_tile_height);
      tile.created = true;
    }
    tile.loaded = true;
    lru.insert(&data[tile_y][tile_x]);

    tile.lazySetAll();
  }

 public:

  A2Array2D(std::string layoutfile, int cachesize){
    lru.setCapacity(cachesize);
    readonly = true;

    int     not_null_tiles = 0;

    LayoutfileReader lf(layoutfile);
    while(lf.next()){
      if(lf.newRow()) //Add a row to the grid of chunks
        data.emplace_back();

      if(lf.isNullTile()){
        data.back().emplace_back();
        data.back().back().null_tile = true;
        continue;
      }

      not_null_tiles++;

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

      auto &this_tile = data.back().back();

      //Get properties of the first file to check against all subsequent ones
      if(per_tile_height==-1){
        per_tile_height = this_tile.viewHeight();
        per_tile_width  = this_tile.viewWidth();
        std::cerr<<"Drawing properties from '"<<this_tile.filename<<"'"<<std::endl;
        std::cerr<<"NoData: "     <<this_tile.noData()    <<std::endl;
        std::cerr<<"Tile Height: "<<this_tile.viewHeight()<<std::endl;
        std::cerr<<"Tile Width: " <<this_tile.viewWidth() <<std::endl;
      }

      cells_in_not_null_tiles += per_tile_width*per_tile_height;

      this_tile.basename = lf.getBasename();

      if(per_tile_width!=this_tile.viewWidth())
        throw std::runtime_error("Tiles were not all the same width!");
      if(per_tile_height!=this_tile.viewHeight())
        throw std::runtime_error("Tiles were not all the same width!");
    }

    total_width_in_cells  = data[0].size()*per_tile_width;
    total_height_in_cells = data.size()*per_tile_height;

    std::cerr<<"Total width: " <<total_width_in_cells<<std::endl;
    std::cerr<<"Total height: "<<total_height_in_cells<<std::endl;

    std::cerr<<"Found "<<not_null_tiles<<" not-null tiles of "<<(data[0].size()*data.size())<<std::endl;
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
  A2Array2D(std::string filename_template, const A2Array2D<U> &other, int cachesize) : A2Array2D(filename_template, other.tileWidth(), other.tileHeight(), other.widthInTiles(), other.heightInTiles(), cachesize) {
    for(int y=0;y<heightInTiles();y++)
    for(int x=0;x<widthInTiles();x++){
      data[y][x].templateCopy(other.data[y][x]);
      data[y][x].filename = filename_template;
      data[y][x].filename.replace(data[y][x].filename.find("%f"), 2, data[y][x].basename);
      data[y][x].null_tile = other.data[y][x].null_tile;
    }
  }

  // T& getn(int tx, int ty, int x, int y, int dx, int dy){
  //   x += dx;
  //   y += dy;
  //   if(x<0){
  //     tx--;
  //     x = per_tile_width-1;
  //   } else if (x==per_tile_width) {
  //     tx++;
  //     x=0;
  //   }
  //   if(y<0){
  //     ty--;
  //     y = per_tile_height-1;
  //   } else if (y==per_tile_height){
  //     ty++;
  //     y=0;
  //   }
  //   assert(x>=0);
  //   assert(y>=0);
  //   assert(x<per_tile_width);
  //   assert(y<per_tile_height);
  //   assert(tx>=0);
  //   assert(ty>=0);
  //   assert(tx<widthInTiles());
  //   assert(ty<heightInTiles());
  //   return data[ty][tx](x,y);
  // }

  // T& operator()(int tx, int ty, int x, int y){
  //   assert(x>=0);
  //   assert(y>=0);
  //   assert(x<per_tile_width);
  //   assert(y<per_tile_height);
  //   assert(tx>=0);
  //   assert(ty>=0);
  //   assert(tx<widthInTiles());
  //   assert(ty<heightInTiles());
  //   LoadTile(tx, ty);
  //   return data[ty][tx](x,y);
  // }

  // const T& operator()(int tx, int ty, int x, int y) const {
  //   assert(x>=0);
  //   assert(y>=0);
  //   assert(x<per_tile_width);
  //   assert(y<per_tile_height);
  //   assert(tx>=0);
  //   assert(ty>=0);
  //   assert(tx<widthInTiles());
  //   assert(ty<heightInTiles());
  //   LoadTile(tx, ty);
  //   return data[ty][tx](x,y);
  // }

  T& operator()(int tx, int ty, int x, int y){
    assert(x>=0);
    assert(y>=0);
    assert(x<per_tile_width);
    assert(y<per_tile_height);

    assert(tx>=0);
    assert(ty>=0);
    assert(tx<data[0].size());
    assert(ty<data.size());

    if(data[ty][tx].null_tile){
      no_data_to_set = data[ty][tx].noData();
      return no_data_to_set;
    }

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

    if(data[tile_y][tile_x].null_tile){
      no_data_to_set = data[tile_y][tile_x].noData();
      return no_data_to_set;
    }

    LoadTile(tile_x, tile_y);

    return data[tile_y][tile_x](x,y);
  }

  // const T& operator()(int x, int y) const {
  //   assert(x>=0);
  //   assert(y>=0);
  //   assert(x<total_width_in_cells);
  //   assert(y<total_height_in_cells);
  //   int tile_x = x/per_tile_width;
  //   int tile_y = y/per_tile_height;
  //   x          = x%per_tile_width;
  //   y          = y%per_tile_height;
  //   if(data[tile_y][tile_x].null_tile)
  //     return no_data;
  //   LoadTile(tile_x, tile_y);
  //   return data[tile_y][tile_x](x,y);
  // }

  int64_t width() const {
    return total_width_in_cells;
  }

  int64_t height() const {
    return total_height_in_cells;
  }

  int64_t widthInTiles() const {
    return data.back().size();
  }

  int64_t heightInTiles() const {
    return data.size();
  }

  int64_t tileWidth() const {
    return per_tile_width;
  }

  int64_t tileHeight() const {
    return per_tile_height;
  }

  void setAll(const T &val){
    for(auto &row: data)
    for(auto &tile: row){
      tile.do_set_all  = true;
      tile.set_all_val = val;
    }
  }

  bool isNoData(int x, int y){
    assert(x>=0);
    assert(y>=0);
    assert(x<total_width_in_cells);
    assert(y<total_height_in_cells);

    int tile_x = x/per_tile_width;
    int tile_y = y/per_tile_height;
    x          = x%per_tile_width;
    y          = y%per_tile_height;

    if(data[tile_y][tile_x].null_tile)
      return true;

    LoadTile(tile_x, tile_y);

    return data[tile_y][tile_x].isNoData(x,y);
  }

  bool isNoData(int tx, int ty, int px, int py){
    assert(px>=0);
    assert(py>=0);
    assert(px<per_tile_width);
    assert(py<per_tile_height);

    assert(tx>=0);
    assert(ty>=0);
    assert(tx<data[0].size());
    assert(ty<data.size());

    if(data[ty][tx].null_tile)
      return true;

    LoadTile(tx, ty);

    return data[ty][tx].isNoData(px,py);
  }

  bool in_grid(int x, int y) const {
    return (x>=0 && y>=0 && x<total_width_in_cells && y<total_height_in_cells);
  }

  bool isReadonly() const {
    return readonly;
  }

  //TODO: Use of checkval here is kinda gross. Is there a good way to get rid of
  //it?
  void saveGDAL(std::string outputname_template) {
    int zero_count      = 0;
    int unvisited_count = 0;
    for(auto &row: data)
    for(auto &tile: row){
      if(tile.null_tile)
        continue;

      //std::cerr<<"Trying to save tile with basename '"<<tile.basename<<"'"<<std::endl;

      if(!tile.loaded)
        tile.loadData();
      //std::cerr<<"\tMin: "<<(int)tile.min()<<" zeros="<<tile.countval(0)<<std::endl;

      if(tile.geotrans[0]<0){
        std::cerr<<"Flip horz"<<std::endl;
        tile.flipHorz();
      }
      if(tile.geotrans[5]<0){
        std::cerr<<"Flip vert"<<std::endl;
        tile.flipVert();
      }

      zero_count      += tile.countval(0);
      unvisited_count += tile.countval(13);

      auto temp = outputname_template;
      temp.replace(temp.find("%f"),2,tile.basename);

      tile.saveGDAL(temp, 0, 0);
      tile.clear();
    }

    std::cerr<<"Found "<<zero_count<<" cells with no flow."<<std::endl;
    std::cerr<<"Found "<<unvisited_count<<" cells that were unvisited."<<std::endl;
  }

  void setNoData(const T &ndval){
    for(auto &row: data)
    for(auto &tile: row)
      if(!tile.null_tile)
        tile.setNoData(ndval); //TODO: Lazy setting: don't set this value until tile is completely loaded
  }

  int32_t getEvictions() const {
    return evictions;
  }

  bool isNullTile(int tx, int ty) const {
    return data[ty][tx].null_tile;
  }
};

#endif