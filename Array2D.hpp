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

//These enable compression in the loadNative() and saveNative() methods
#ifdef WITH_COMPRESSION
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#endif

#ifndef d8flowdirs_dxdy
#define d8flowdirs_dxdy
const int dx[9]={0,-1,-1,0,1,1,1,0,-1};  //TODO: These should be merged with my new dinf_d8 to reflect a uniform and intelligent directional system
///y offsets of D8 neighbours, from a central cell
const int dy[9]={0,0,-1,-1,-1,0,1,1,1};
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


template<class T>
void getGDALHeader(
  const std::string &filename,
  int    &height,
  int    &width,
  T      &no_data,
  double *geotrans
){
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  assert(fin!=NULL);

  GDALRasterBand *band   = fin->GetRasterBand(1);

  height  = band->GetYSize();
  no_data = band->GetNoDataValue();
  width   = band->GetXSize();

  fin->GetGeoTransform(geotrans);

  GDALClose(fin);
}


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





template<class T>
class Array2D {
 public:
  std::string filename;
  std::string basename;
  std::vector<double> geotransform;
  std::string projection;

 private:
  template<typename> friend class Array2D;

  std::vector<T> data;

  GDALDataType data_type;

  size_t header_size = -1;

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

    out.write(reinterpret_cast<const char*>(&total_height),   sizeof(int));
    out.write(reinterpret_cast<const char*>(&total_width),    sizeof(int));
    out.write(reinterpret_cast<const char*>(&view_height),    sizeof(int));
    out.write(reinterpret_cast<const char*>(&view_width),     sizeof(int));
    out.write(reinterpret_cast<const char*>(&view_xoff),      sizeof(int));
    out.write(reinterpret_cast<const char*>(&view_yoff),      sizeof(int));
    out.write(reinterpret_cast<const char*>(&num_data_cells), sizeof(int));
    out.write(reinterpret_cast<const char*>(&no_data),        sizeof(T  ));

    out.write(reinterpret_cast<const char*>(geotransform.data()), 6*sizeof(double));
    std::string::size_type projection_size = projection.size();
    out.write(reinterpret_cast<const char*>(&projection_size), sizeof(std::string::size_type));
    out.write(reinterpret_cast<const char*>(projection.data()), projection.size()*sizeof(const char));

    out.write(reinterpret_cast<const char*>(data.data()), data.size()*sizeof(T));
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

      in.seekg(header_size);

      in.read(reinterpret_cast<char*>(data.data()), view_width*view_height*sizeof(T));
    } else {
      GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
      if(fin==NULL){
        std::cerr<<"Failed to loadData() into tile from '"<<filename<<"'"<<std::endl;
        throw std::runtime_error("Failed to loadData() into tile.");
      }

      GDALRasterBand *band = fin->GetRasterBand(1);

      data.resize(view_width*view_height);
      auto temp = band->RasterIO( GF_Read, view_xoff, view_yoff, view_width, view_height, data.data(), view_width, view_height, myGDALType(), 0, 0 );
      if(temp!=CE_None)
        throw std::runtime_error("Error reading file with GDAL!");

      std::cerr<<"Geotrans: ";
      for(auto x: geotransform)
        std::cerr<<x<<" ";
      std::cerr<<std::endl;

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
    geotransform.resize(6);
    in.read(reinterpret_cast<char*>(geotransform.data()), 6*sizeof(double));

    std::string::size_type projection_size;
    in.read(reinterpret_cast<char*>(&projection_size), sizeof(std::string::size_type));
    projection.resize(projection_size,' ');
    in.read(reinterpret_cast<char*>(&projection[0]), projection.size()*sizeof(char));

    header_size = in.tellg();

    if(load_data)
      loadData();
  }

  //Note: The following functions return signed integers, which make them
  //generally easier to work with. If your DEM has a dimension which exceeds
  //2147483647, some other modifications to this program will probably be
  //necessary.
  int  viewSize   () const { return view_width*view_height; }
  int  totalWidth () const { return total_width;    }
  int  totalHeight() const { return total_height;   }
  size_t viewWidth  () const { return view_width;     }
  size_t viewHeight () const { return view_height;    }
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

  void iToxy(const int i, int &x, int &y) const {
    x = i%view_width;
    y = i/view_width;
  }

  int xyToI(int x, int y) const {
    return y*view_width+x;
  }

  int nToI(int i, int dx, int dy) const {
    int x=i%view_width+dx;
    int y=i/view_width+dy;
    if(x<0 || y<0 || x==view_width || y==view_height)
      return -1;
    return xyToI(x,y);
  }

  template<class U>
  T& operator=(const Array2D<U> &o){
    data = std::vector<T>(o.data.begin(),o.data.end());
    total_height   = o.total_height;
    total_width    = o.total_width;
    view_height    = o.view_height;
    view_width     = o.view_width;
    view_xoff      = o.view_xoff;
    view_yoff      = o.view_yoff;
    num_data_cells = o.num_data_cells;
    geotransform   = o.geotransform;
    no_data        = (T)o.no_data;
    return *this;
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

  void transpose(){
    std::cerr<<"transpose() is an experimental feature."<<std::endl;
    std::vector<T> new_data(view_width*view_height);
    for(int y=0;y<view_height;y++)
    for(int x=0;x<view_width;x++)
      new_data[x*view_height+y] = data[y*view_width+x];
    data = new_data;
    std::swap(view_width,view_height);
    //TODO
  }

  bool in_grid(int x, int y) const {
    return 0<=x && x<view_width && 0<=y && y<view_height;
  }

  bool isEdgeCell(int x, int y) const {
    return x==0 || y==0 || x==view_width-1 || y==view_height-1;
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

  template<class U>
  void resize(const Array2D<U> &other, const T& val = T()){
    resize(other.viewWidth(), other.viewHeight(), val);
    geotransform = other.geotransform;
  }

  void expand(int new_width, int new_height, const T val){
    if(new_width<view_width)
      throw std::runtime_error("expand(): new_width<view_width");
    if(new_height<view_height)
      throw std::runtime_error("expand(): new_height<view_height");
    
    int old_width  = viewWidth();
    int old_height = viewHeight();

    std::vector<T> old_data = std::move(data);

    resize(new_width,new_height,val);

    for(int y=0;y<old_height;y++)
    for(int x=0;x<old_width;x++)
      data[y*new_width+x] = old_data[y*old_width+x];
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

  int numDataCells() const {
    return num_data_cells;
  }

  T& operator()(int i){
    assert(i>=0);
    assert(i<view_width*view_height);
    return data[i];
  }

  T operator()(int i) const {
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

  T& operator()(size_t x, size_t y){
    assert(x>=0);
    assert(y>=0);
    //std::cerr<<"Width: "<<viewWidth()<<" Height: "<<viewHeight()<<" x: "<<x<<" y: "<<y<<std::endl;
    assert(x<viewWidth());
    assert(y<viewHeight());
    return data[y*view_width+x];
  }

  const T& operator()(size_t x, size_t y) const {
    assert(x>=0);
    assert(y>=0);
    assert(x<viewWidth());
    assert(y<viewHeight());
    return data[y*view_width+x];
  }

  std::vector<T> topRow() const {    //TODO: Test
    return std::vector<T>(data.begin(),data.begin()+view_width);
  }

  std::vector<T> bottomRow() const { //TODO: Test
    return std::vector<T>(data.begin()+(view_height-1)*view_width, data.begin()+view_height*view_width);
  }

  std::vector<T> leftColumn() const { //TODO: Test
    return getColData(0); 
  }

  std::vector<T> rightColumn() const { //TODO: Test
    return getColData(view_width-1);
  }

  // Row& rowRef(int rownum){
  //   return data[rownum];
  // }

  void setRow(int y, const T &val){
    std::fill(data.begin()+y*view_width,data.begin()+(y+1)*view_width,val);
  }

  void setCol(int x, const T &val){
    for(int y=0;y<view_height;y++)
      data[y*view_width+x] = val;
  }

  std::vector<T> getRowData(int rownum) const {
    return std::vector<T>(data.begin()+rownum*view_width,data.begin()+(rownum+1)*view_width);
  }

  std::vector<T> getColData(int colnum) const {
    std::vector<T> temp(view_height);
    for(int y=0;y<view_height;y++)
      temp[y]=data[y*view_width+colnum];
    return temp;
  }

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
      std::cerr<<"Could not open GDAL driver!"<<std::endl;
      throw std::runtime_error("Could not open GDAL driver!");
    }
    GDALDataset *fout    = poDriver->Create(filename.c_str(), viewWidth(), viewHeight(), 1, myGDALType(), NULL);
    if(fout==NULL){
      std::cerr<<"Could not open file '"<<filename<<"' for GDAL save!"<<std::endl;
      throw std::runtime_error("Could not open file for GDAL save!");
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

    if(out_geotransform.size()!=6){
      std::cerr<<"Geotransform of output is not the right size. Found "<<out_geotransform.size()<<" expected 6."<<std::endl;
      throw std::runtime_error("saveGDAL(): Invalid output geotransform.");
    }

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

  double getCellArea() const {
    return geotransform[1]*geotransform[5];
  }
};

#endif