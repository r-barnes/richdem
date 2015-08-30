#ifndef _array_2d_hpp_
#define _array_2d_hpp_

#include "gdal_priv.h"
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <functional>

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
void getGDALHeader(const std::string &filename, int &height, int &width, T &no_data, double *geotrans){
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

void getGDALDimensions(const std::string &filename, int &height, int &width){
  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
  assert(fin!=NULL);

  GDALRasterBand *band = fin->GetRasterBand(1);

  height  = band->GetYSize();
  width   = band->GetXSize();

  GDALClose(fin);
}


template<class T>
GDALDataType NativeTypeToGDAL() {
  if(typeid(T)==typeid(int16_t))
    return GDT_Int16;
  else if(typeid(T)==typeid(float))
    return GDT_Float32;
  else if(typeid(T)==typeid(uint8_t))
    return GDT_Byte;
  else
    assert(false);
}


template<class T>
void stackSavedNativeStrips(const std::string &fileA, const std::string &fileB) {
  std::fstream  finA(fileA, std::ios_base::binary | std::ios_base::in | std::ios_base::out);
  std::ifstream finB(fileB);
  assert(finA.good());
  assert(finB.good());

  int heightA, widthA, heightB, widthB;
  finA.read(reinterpret_cast<char*>(&heightA), sizeof(int));
  finA.read(reinterpret_cast<char*>(&widthA ), sizeof(int));
  finB.read(reinterpret_cast<char*>(&heightB), sizeof(int));
  finB.read(reinterpret_cast<char*>(&widthB ), sizeof(int));

  assert(widthA==widthB);

  //Expand A's height to include the other array. Write this to memory.
  heightA += heightB;

  finA.seekp(0);
  finA.write(reinterpret_cast<char*>(&heightA), sizeof(int));

  //Skip to the end of A
  finA.seekp(0, std::ios_base::end);

  //Append the B's saved data to A
  finA << finB.rdbuf();
}

template<class T>
void combineSavedNativeIntoGDAL(const std::string &out_name, const std::string &templ_name, const std::vector<std::string> &filenames, T no_data, std::function<T(T)> trans){
  int width  = -1;
  int height = 0;
  for(auto &in_name: filenames){
    std::ifstream fin(in_name);
    assert(fin.good());
    int this_height, this_width;
    fin.read(reinterpret_cast<char*>(&this_height), sizeof(int));
    fin.read(reinterpret_cast<char*>(&this_width ), sizeof(int));
    if(width==-1)
      width = this_width;
    else
      assert(width==this_width);
    height += this_height;
  }

  GDALDataType data_type = NativeTypeToGDAL<T>();

  GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
  GDALDataset *fout = poDriver->Create(out_name.c_str(), width, height, 1, data_type, NULL);
  assert(fout!=NULL);
  GDALDataset *fintempl = (GDALDataset*)GDALOpen(templ_name.c_str(), GA_ReadOnly);
  assert(fintempl!=NULL);

  double geotrans[6];
  fintempl->GetGeoTransform(geotrans);
  fout->SetGeoTransform(geotrans);

  const char* projection_string=fintempl->GetProjectionRef();
  fout->SetProjection(projection_string);

  GDALRasterBand *oband = fout->GetRasterBand(1);
  oband->SetNoDataValue(no_data);

  std::vector<T> rowdata(width);
  int y_overall = 0;
  for(auto &in_name: filenames){
    std::ifstream fin(in_name);
    assert(fin.good());

    int this_height, this_width;
    fin.read(reinterpret_cast<char*>(&this_height), sizeof(int));
    fin.read(reinterpret_cast<char*>(&this_width ), sizeof(int));

    for(int y=0;y<this_height;y++,y_overall++){
      fin.read(reinterpret_cast<char*>(rowdata.data()), width*sizeof(T));
      //Convert row directions to ArcGIS format for output
      if(trans!=nullptr)
        for(auto &fd: rowdata)
          fd = trans(fd);
      oband->RasterIO(GF_Write, 0, y_overall, width, 1, rowdata.data(), width, 1, data_type, 0, 0);
    }
  }

  GDALClose(fout);
}


template<class T>
class Array2D {
 public:
  typedef std::vector<T>   Row;
 private:
  typedef std::vector<Row> InternalArray;
  InternalArray data;

  std::string ram_name;

  static const int HEADER_SIZE = 2*sizeof(int);

  int total_height;
  int total_width;
  int view_height;
  int view_width;
  int view_xoff;
  int view_yoff;

  T   no_data;

  void loadGDAL(const std::string &filename, int xOffset=0, int yOffset=0, int part_width=0, int part_height=0){
    assert(empty());
    assert(xOffset>=0);
    assert(yOffset>=0);

    GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
    assert(fin!=NULL);

    GDALRasterBand *band = fin->GetRasterBand(1);
    auto data_type       = band->GetRasterDataType();

    assert(
      (data_type==GDT_Int16   && typeid(T)==typeid(int16_t)) ||
      (data_type==GDT_Float32 && typeid(T)==typeid(float))
    );

    total_width  = band->GetXSize();
    total_height = band->GetYSize();
    no_data      = band->GetNoDataValue();

    if(xOffset+part_width>=total_width)
      part_width = total_width-xOffset;
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

    std::cerr<<"Allocating: "<<view_height<<"x"<<view_width<<std::endl;
    data = InternalArray(view_height, Row(view_width));
    for(int y=yOffset;y<yOffset+view_height;y++)
      band->RasterIO( GF_Read, xOffset, y, view_width, 1, data[y-yOffset].data(), view_width, 1, data_type, 0, 0 );

    GDALClose(fin);
  }

  void loadNative(const std::string &filename, int xOffset=0, int yOffset=0, int view_width0=0, int view_height0=0){
    assert(xOffset>=0);
    assert(yOffset>=0);

    std::ifstream fin(filename, std::ios::in | std::ios::binary);
    assert(fin.good());
    fin.read(reinterpret_cast<char*>(&total_height), sizeof(int));
    fin.read(reinterpret_cast<char*>(&total_width ), sizeof(int));
    assert(xOffset+view_width0<total_width || yOffset+view_height0<total_height);

    if(view_width0==0)
      view_width0 = total_width;
    view_width = view_width0;

    if(view_height0==0)
      view_height0 = total_height;
    view_height = view_height0;

    view_xoff = xOffset;
    view_yoff = yOffset;

    data = InternalArray(view_height, Row(view_width));

    for(int y=yOffset;y<yOffset+view_height;y++){
      fin.seekg(HEADER_SIZE+(y*total_width+xOffset)*sizeof(T));
      fin.read(reinterpret_cast<char*>(data[y-yOffset].data()), view_width*sizeof(T));
    }
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
  Array2D(int height, int width, const T& val = T()) : Array2D() {
    data      = InternalArray(height, Row(width, val));
    total_height = view_height = height;
    total_width  = view_width  = width;
  }

  //Create internal array from a file
  Array2D(const std::string &filename, bool use_gdal=true, int xOffset=0, int yOffset=0, int part_width=0, int part_height=0) : Array2D() {
    if(use_gdal)
      loadGDAL(filename, xOffset, yOffset, part_width, part_height);
    else
      loadNative(filename, xOffset, yOffset, part_width, part_height);
  }

  //Create an internal array from a file, dividing it into strips
  Array2D(const std::string &filename, int my_strip_number, int strip_count){
    GDALDataset *fin = (GDALDataset*)GDALOpen(filename.c_str(), GA_ReadOnly);
    assert(fin!=NULL);

    GDALRasterBand *band = fin->GetRasterBand(1);

    int width  = band->GetXSize();
    int height = band->GetYSize();
    GDALClose(fin);

    int segment_first_line = (height/strip_count)*my_strip_number;
    int segment_last_line  = (height/strip_count)*(my_strip_number+1);
    if(my_strip_number==strip_count-1)
      segment_last_line = height;

    int segment_height = segment_last_line - segment_first_line;

    loadGDAL(filename,0,segment_first_line,width,segment_height);
  }

  int  width     () const { return total_width;  }
  int  height    () const { return total_height; }
  int  viewWidth () const { return view_width;   }
  int  viewHeight() const { return view_height;  }
  int  viewXoff  () const { return view_xoff;    }
  int  viewYoff  () const { return view_yoff;    }
  bool empty     () const { return data.empty(); }
  T    noData    () const { return no_data; }

  void setNoData(const T &ndval){
    no_data = ndval;
  }

  void setAll(const T &val){
    for(auto &row: data)
      std::fill(row.begin(),row.end(),val);
  }

  void redimension(){
    total_height=view_height=data.size();
    if(view_height>0)
      total_width=view_width=data.front().size();
    else
      total_width=view_width=0;
  }

  T& operator()(int x, int y){
    assert(x>=0);
    assert(y>=0);
    assert((unsigned int)x<data[0].size());
    assert((unsigned int)y<data.size());
    return data[y][x];
  }

  const T& operator()(int x, int y) const {
    assert(x>=0);
    assert(y>=0);
    assert((unsigned int)x<data[0].size());
    assert((unsigned int)y<data.size());
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

  void setRow(int rownum, const Row &row){
    assert(row.size()==(unsigned int)total_width);
    data[rownum] = row;
  }

  void offloadFromRAM(const std::string &filename){
    ram_name = filename;
    saveNative(filename);
    data.clear();
    data.shrink_to_fit();
  }

  void loadIntoRAM(){
    assert(ram_name!="");

    std::ifstream fin(ram_name, std::ios::in | std::ios::binary);
    assert(fin.good());

    data = InternalArray(view_height, Row(view_width));

    fin.seekg(HEADER_SIZE);

    for(int y=0;y<view_height;y++)
      fin.read(reinterpret_cast<char*>(data[y].data()), view_width*sizeof(T));
  }


  void print() const {
    for(auto &row: data){
      for(auto &col: row)
        std::cout<<col<<" ";
      std::cerr<<std::endl;
    }
  }

  void init(T val){
    for(auto &row: data)
    for(auto &col: row)
      col=val;
  }

  void clear(){
    data.clear();
    data.shrink_to_fit();
  }

  void saveNative(const std::string &filename){
    std::fstream fout;

    fout.open(filename, std::ios_base::binary | std::ios_base::in | std::ios_base::out | std::ios::trunc);
    assert(fout.good());

    fout.write(reinterpret_cast<char*>(&view_height), sizeof(int));
    fout.write(reinterpret_cast<char*>(&view_width ), sizeof(int));
    for(int y=0;y<view_height;y++)
      fout.write(reinterpret_cast<char*>(data[y].data()), view_width*sizeof(T));
  }

  void saveGDAL(const std::string &filename, const std::string &template_name){
    GDALDataset *fintempl = (GDALDataset*)GDALOpen(template_name.c_str(), GA_ReadOnly);
    assert(fintempl!=NULL);

    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    assert(poDriver!=NULL);
    GDALDataset *fout    = poDriver->Create(filename.c_str(), view_width, view_height, 1, myGDALType(), NULL);
    assert(fout!=NULL);

    GDALRasterBand *oband = fout->GetRasterBand(1);
    oband->SetNoDataValue(no_data);

    double geotrans[6];
    fintempl->GetGeoTransform(geotrans);
    fout->SetGeoTransform(geotrans);

    const char* projection_string=fintempl->GetProjectionRef();
    fout->SetProjection(projection_string);

    GDALClose(fintempl);

    for(int y=0;y<view_height;y++)
      oband->RasterIO(GF_Write, 0, y, view_width, 1, data[y].data(), view_width, 1, myGDALType(), 0, 0);

    GDALClose(fout);
  }
};

#endif