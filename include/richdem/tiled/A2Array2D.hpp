/**
  @file
  @brief Experimental tile manager for large datasets (TODO)
  @author Richard Barnes
*/
#ifndef _a2array2d_hpp_
#define _a2array2d_hpp_

#include "richdem/common/logger.hpp"
#include "richdem/common/Layoutfile.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/tiled/lru.hpp"
#include "gdal_priv.h"

namespace richdem {

GDALDataType peekLayoutType(const std::string &layout_filename) {
  LayoutfileReader lf(layout_filename);

  while(lf.next()){
    if(lf.getFilename().size()==0)
      continue;

    GDALAllRegister();
    std::string tile_path = lf.getPath()+lf.getFilename();
    GDALDataset *fin = (GDALDataset*)GDALOpen(tile_path.c_str(), GA_ReadOnly);
    if(fin==NULL)
      throw std::runtime_error("Could not open '"+(lf.getPath()+lf.getFilename())+"' to determine layout type.");

    GDALRasterBand *band   = fin->GetRasterBand(1);
    GDALDataType data_type = band->GetRasterDataType();

    GDALClose(fin);

    return data_type;
  }

  throw std::runtime_error("Empty layout file!");
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
class A2Array2D {
 public:
  bool flipH = false;
  bool flipV = false;
 private:
  template<typename U> friend class A2Array2D;

  std::vector<bool> null_tile_quick;
  int quick_width_in_tiles;
  int quick_height_in_tiles;

  class WrappedArray2D : public Array2D<T> {
   public:
    using Array2D<T>::Array2D;
    bool null_tile         = false;
    bool loaded            = false;
    bool created           = true;
    bool do_set_all        = false; //If true, then set all to 'set_all_val' when tile is loaded
    int create_with_width  = -1;
    int create_with_height = -1;
    T set_all_val          = 0;
    void lazySetAll(){
      if(do_set_all){
        do_set_all = false;
        setAll(set_all_val);
      }
    }
  };
  std::vector< std::vector< WrappedArray2D > > data;

  LRU< WrappedArray2D* > lru;

  int32_t not_null_tiles          = 0;
  int64_t total_width_in_cells    = 0;
  int64_t total_height_in_cells   = 0;
  int32_t per_tile_width          = 0;
  int32_t per_tile_height         = 0;
  int32_t evictions               = 0;
  int64_t cells_in_not_null_tiles = 0;
  T       no_data_to_set; //Used to disguise null tiles

  bool readonly = true;

  void _LoadTile(int tile_x, int tile_y){
    if(isNullTile(tile_x,tile_y))
      return;

    auto& tile = data[tile_y][tile_x];

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
      tile.printStamp(5,"Tile load, before reorientating"); //Print stamp before reorientating since this must match parallel_pf.exe
      if(readonly){
        if((tile.geotransform[1]<0) ^ flipH)
          tile.flipHorz();
        if((tile.geotransform[5]>0) ^ flipV)
          tile.flipVert();
      }
      tile.printStamp(5,"Tile load, after reorientating"); //Print stamp before reorientating since this must match parallel_pf.exe
    } else {
      if(tile.create_with_width!=-1 && tile.create_with_height!=-1)
        tile.resize(tile.create_with_width,tile.create_with_height);
      else
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

      per_tile_height = std::max(per_tile_height,this_tile.height());
      per_tile_width  = std::max(per_tile_width, this_tile.width() );

      cells_in_not_null_tiles += per_tile_width*per_tile_height;

      this_tile.basename = lf.getBasename();
    }

    quick_width_in_tiles  = widthInTiles();
    quick_height_in_tiles = heightInTiles();
    null_tile_quick.resize(quick_width_in_tiles*quick_height_in_tiles, false);

    bool good=true;
    for(int32_t ty=0;ty<heightInTiles();ty++)
    for(int32_t tx=0;tx<widthInTiles();tx++){
      null_tile_quick[ty*quick_width_in_tiles+tx] = data[ty][tx].null_tile;
      if(data[ty][tx].null_tile)
        continue;
      if(data[ty][tx].width()!=per_tile_width){
        RDLOG_WARN<<data[ty][tx].filename<<" has a non-standard width. Found "<<data[ty][tx].width()<<" expected "<<per_tile_width<<".";
        good = false;
      }
      if(data[ty][tx].height()!=per_tile_height){
        RDLOG_WARN<<data[ty][tx].filename<<" has a non-standard height. Found "<<data[ty][tx].height()<<" expected "<<per_tile_height<<".";
        good = false;
      }
    }

    total_width_in_cells  = widthInTiles()*stdTileWidth();
    total_height_in_cells = heightInTiles()*stdTileHeight();

    if(!good)
      throw std::runtime_error("Not all tiles had the same dimensions!");

    RDLOG_MISC<<"Total width = " <<total_width_in_cells;
    RDLOG_MISC<<"Total height = "<<total_height_in_cells;

    RDLOG_MISC<<"Tiles that were not null = "<<not_null_tiles;
    RDLOG_MISC<<"Total tiles = "<<(data[0].size()*data.size());
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
        data.back().back().setCacheFilename(prefix+std::to_string(tile)+".native");
        data.back().back().created  = false;
      }
    }
  }

  template<class U>
  A2Array2D(std::string filename_template, const A2Array2D<U> &other, int cachesize) {
    lru.setCapacity(cachesize);

    readonly = false;

    per_tile_width        = 0;
    per_tile_height       = 0;
    total_width_in_cells  = other.total_width_in_cells;
    total_height_in_cells = other.total_height_in_cells;
    flipV                 = other.flipV;
    flipH                 = other.flipH;

    quick_width_in_tiles  = other.quick_width_in_tiles;
    quick_height_in_tiles = other.quick_height_in_tiles;
    null_tile_quick       = other.null_tile_quick;

    for(int32_t y=0;y<other.heightInTiles();y++){
      data.emplace_back();
      for(int32_t x=0;x<other.widthInTiles();x++){
        data.back().emplace_back();
        auto &this_tile  = data.back().back();
        auto &other_tile = other.data[y][x];

        this_tile.templateCopy(other_tile);
        this_tile.filename = filename_template;
        this_tile.filename.replace(data[y][x].filename.find("%f"), 2, data[y][x].basename);
        per_tile_width               = std::max(per_tile_width, other_tile.width() );
        per_tile_height              = std::max(per_tile_height,other_tile.height());
        this_tile.create_with_width  = other_tile.width();
        this_tile.create_with_height = other_tile.height();
        this_tile.null_tile          = other_tile.null_tile;
        this_tile.created            = false;
      }
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
  //   _LoadTile(tx, ty);
  //   return data[ty][tx](x,y);
  // }

  T& operator()(int32_t tx, int32_t ty, int32_t x, int32_t y){
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<widthInTiles());
    assert(ty<heightInTiles());
    assert(widthInTiles()>0);
    assert(heightInTiles()>0);

    if(isNullTile(tx,ty)){
      no_data_to_set = data[ty][tx].noData();
      return no_data_to_set;
    }

    _LoadTile(tx, ty);

    assert(x>=0);
    assert(y>=0);
    assert(x<data[ty][tx].width() );
    assert(y<data[ty][tx].height());

    return data[ty][tx](x,y);
  }

  T& operator()(int32_t x, int32_t y){
    assert(x>=0);
    assert(y>=0);
    assert(x<total_width_in_cells);
    assert(y<total_height_in_cells);

    int32_t tile_x = x/per_tile_width;
    int32_t tile_y = y/per_tile_height;
    x              = x%per_tile_width;
    y              = y%per_tile_height;

    if(isNullTile(tile_x,tile_y)){
      no_data_to_set = data[tile_y][tile_x].noData();
      return no_data_to_set;
    }

    _LoadTile(tile_x, tile_y);

    return data[tile_y][tile_x](x,y);
  }

  void makeQuadIndex(int32_t x, int32_t y, int32_t &tx, int32_t &ty, int32_t &px, int32_t &py) const {
    assert(x>=0);
    assert(y>=0);
    assert(x<total_width_in_cells);
    assert(y<total_height_in_cells);

    tx = x/per_tile_width;
    px = x%per_tile_width;

    ty = y/per_tile_height;
    py = y%per_tile_height;
  }

  int32_t width() const {
    return total_width_in_cells;
  }

  int32_t height() const {
    return total_height_in_cells;
  }

  int32_t widthInTiles() const {
    return data.back().size();
  }

  int32_t heightInTiles() const {
    return data.size();
  }

  int32_t notNullTiles() const {
    return not_null_tiles;
  }

  int32_t tileWidth(int32_t tx, int32_t ty) const {
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<(int32_t)data[0].size());
    assert(ty<(int32_t)data.size());
    return data[ty][tx].width();
  }

  int32_t tileHeight(int32_t tx, int32_t ty) const {
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<widthInTiles());
    assert(ty<heightInTiles());
    return data[ty][tx].height();
  }

  int32_t stdTileHeight() const {
    return per_tile_height;
  }

  int32_t stdTileWidth() const {
    return per_tile_width;
  }

  void setAll(const T &val){
    for(auto &row: data)
    for(auto &tile: row){
      tile.do_set_all  = true;
      tile.set_all_val = val;
    }
  }

  bool isNoData(int32_t x, int32_t y){
    assert(x>=0);
    assert(y>=0);
    assert(x<total_width_in_cells);
    assert(y<total_height_in_cells);

    int32_t tile_x = x/per_tile_width;
    int32_t tile_y = y/per_tile_height;
    x              = x%per_tile_width;
    y              = y%per_tile_height;

    if(isNullTile(tile_x,tile_y))
      return true;

    _LoadTile(tile_x, tile_y);

    return data[tile_y][tile_x].isNoData(x,y);
  }

  bool isNoData(int32_t tx, int32_t ty, int32_t px, int32_t py){
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<(int32_t)data[0].size());
    assert(ty<(int32_t)data.size());

    assert(px>=0);
    assert(py>=0);
    assert(px<data[ty][tx].width() );
    assert(py<data[ty][tx].height());

    if(isNullTile(tx,ty))
      return true;

    _LoadTile(tx, ty);

    return data[ty][tx].isNoData(px,py);
  }

  bool isReadonly() const {
    return readonly;
  }

  ///Returns the GDAL data type of the A2Array2D template type
  GDALDataType myGDALType() const {
    return NativeTypeToGDAL<T>();
  }

  //TODO: Use of checkval here is kinda gross. Is there a good way to get rid of
  //it?
  void saveGDAL(std::string outputname_template) {
    int zero_count      = 0;

    auto new_layout_name = outputname_template;
    new_layout_name.replace(new_layout_name.find("%f"),2,"layout");

    LayoutfileWriter lfout(new_layout_name);

    for(int32_t ty=0;ty<heightInTiles();ty++){
      lfout.addRow();
      for(int32_t tx=0;tx<widthInTiles();tx++){
        auto& tile = data[ty][tx];

        if(tile.null_tile){
          lfout.addEntry("");
          continue;
        }

        _LoadTile(tx,ty);

        tile.printStamp(5,"Saving, before reorientation");

        if((tile.geotransform[1]<0) ^ flipH)
          tile.flipHorz();
        if((tile.geotransform[5]>0) ^ flipV)
          tile.flipVert();

        tile.printStamp(5,"Saving, after reorientation");

        zero_count += tile.countval(NO_FLOW);

        auto temp = outputname_template;
        temp.replace(temp.find("%f"),2,tile.basename);
        lfout.addEntry(temp);

        tile.saveGDAL(temp, 0, 0);
      }
    }

    RDLOG_MISC<<"Cells with no flow = "<<zero_count;
  }

  void saveUnifiedGDAL(const std::string outputname){
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if(poDriver==NULL)
      throw std::runtime_error("Could not open GDAL driver!");
    GDALDataset *fout    = poDriver->Create(outputname.c_str(), width(), height(), 1, myGDALType(), NULL);
    if(fout==NULL)
      throw std::runtime_error("Could not open file '"+outputname+"' for GDAL save!");

    auto no_data = data[0][0].noData();
    for(int32_t ty=0;ty<heightInTiles();ty++)
    for(int32_t tx=0;tx<widthInTiles();tx++)
      if(data[ty][tx].noData()!=no_data)
        throw std::runtime_error("Files had differing NoData values :-(");

    auto out_geotransform = data[0][0].geotransform;

    if(out_geotransform.size()!=6)
      throw std::runtime_error("Geotransform of output is not the right size. Found "+std::to_string(out_geotransform.size())+" expected 6.");

    fout->SetGeoTransform(out_geotransform.data());
    fout->SetProjection(data[0][0].projection.c_str());

    GDALRasterBand *oband = fout->GetRasterBand(1);
    oband->SetNoDataValue(no_data);
    oband->Fill(no_data);

    for(int32_t ty=0;ty<heightInTiles();ty++)
    for(int32_t tx=0;tx<widthInTiles();tx++){
      auto& tile = data[ty][tx];

      if(tile.null_tile)
        continue;

      _LoadTile(tx,ty);

      auto temp = oband->RasterIO(GF_Write, tx*stdTileWidth(), ty*stdTileHeight(), tileWidth(tx,ty), tileHeight(tx,ty), data[ty][tx].getData(), tileWidth(tx,ty), tileHeight(tx,ty), myGDALType(), 0, 0);
      if(temp!=CE_None)
        throw std::runtime_error("Error writing file!");
    }

    GDALClose(fout);
  }

  void setNoData(const T &ndval){
    for(auto &row: data)
    for(auto &tile: row)
      if(!tile.null_tile)
        tile.setNoData(ndval);
  }

  int32_t getEvictions() const {
    return evictions;
  }

  inline bool isNullTile(int32_t tx, int32_t ty) const {
    return null_tile_quick[ty*quick_width_in_tiles+tx];
  }

  bool isEdgeCell(int32_t x, int32_t y) const {
    return x==0 || y==0 || x==total_width_in_cells-1 || y==total_height_in_cells-1;
  }

  bool in_grid(int32_t x, int32_t y) const {
    return (x>=0 && y>=0 && x<total_width_in_cells && y<total_height_in_cells);
  }

  bool isInteriorCell(int32_t x, int32_t y) const {
    return (1<=x && 1<=y && x<total_width_in_cells-1 && y<total_height_in_cells-1);
  }

  void printStamp(int size){
    for(int32_t ty=0;ty<heightInTiles();ty++)
    for(int32_t tx=0;tx<widthInTiles();tx++){
      if(isNullTile(tx,ty))
        continue;

      _LoadTile(tx, ty);

      data[ty][tx].printStamp(size);
    }
  }

  void loadTile(int tx, int ty){
    _LoadTile(tx,ty);
  }
};

}

#endif
