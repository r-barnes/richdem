#ifndef _a2array2d_hpp_
#define _a2array2d_hpp_

#include "gdal_priv.h"
#include "Layoutfile.hpp"
#include "Array2D.hpp"
#include "lru.hpp"

GDALDataType peekLayoutType(const std::string &layout_filename) {
  LayoutfileReader lf(layout_filename);

  while(lf.next()){
    if(lf.getFilename().size()==0)
      continue;

    GDALAllRegister();
    std::string tile_path = lf.getPath()+lf.getFilename();
    GDALDataset *fin = (GDALDataset*)GDALOpen(tile_path.c_str(), GA_ReadOnly);
    if(fin==NULL){
      std::cerr<<"Could not open '"<<(lf.getPath()+lf.getFilename())<<"' to determine layout type."<<std::endl;
      throw std::runtime_error("Could not open one of the data files!");
    }

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

  size_t total_width_in_cells     = 0;
  size_t total_height_in_cells    = 0;
  size_t  per_tile_width          = 0;
  size_t  per_tile_height         = 0;
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
      if(readonly){
        if((tile.geotransform[0]<0) ^ flipH)
          tile.flipHorz();
        if((tile.geotransform[5]<0) ^ flipV)
          tile.flipVert();
      }
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

      per_tile_height = std::max(per_tile_height,this_tile.viewHeight());
      per_tile_width  = std::max(per_tile_width, this_tile.viewWidth() );

      if(lf.getY()==0)
        total_width_in_cells += this_tile.viewWidth();
      if(lf.getX()==0)
        total_height_in_cells += this_tile.viewHeight();

      cells_in_not_null_tiles += per_tile_width*per_tile_height;

      this_tile.basename = lf.getBasename();
    }

    bool good=true;
    for(size_t y=0;y<heightInTiles()-1;y++)
    for(size_t x=0;x<widthInTiles()-1;x++){
      if(data[y][x].null_tile)
        continue;
      if(data[y][x].viewWidth()!=per_tile_width){
        std::cerr<<data[y][x].filename<<" has a non-standard width. Found "<<data[y][x].viewWidth()<<" expected "<<per_tile_width<<"."<<std::endl;
        good = false;
      }
      if(data[y][x].viewHeight()!=per_tile_height){
        std::cerr<<data[y][x].filename<<" has a non-standard height. Found "<<data[y][x].viewHeight()<<" expected "<<per_tile_height<<"."<<std::endl;
        good = false;
      }
    }

    if(!good){
      throw std::runtime_error("Not all tiles had the same dimensions!");
    }

    std::cerr<<"Total width: " <<total_width_in_cells<<std::endl;
    std::cerr<<"Total height: "<<total_height_in_cells<<std::endl;

    std::cerr<<not_null_tiles<<" of "<<(data[0].size()*data.size())<<" tiles were not null."<<std::endl;
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
  A2Array2D(std::string filename_template, const A2Array2D<U> &other, int cachesize) {
    lru.setCapacity(cachesize);

    readonly = false;

    per_tile_width        = 0;
    per_tile_height       = 0;
    total_width_in_cells  = other.total_width_in_cells;
    total_height_in_cells = other.total_height_in_cells;
    flipV                 = other.flipV;
    flipH                 = other.flipH;

    for(size_t y=0;y<other.heightInTiles();y++){
      data.emplace_back();
      for(size_t x=0;x<other.widthInTiles();x++){
        data.back().emplace_back();
        auto &this_tile  = data.back().back();
        auto &other_tile = other.data[y][x];

        this_tile.templateCopy(other_tile);
        this_tile.filename = filename_template;
        this_tile.filename.replace(data[y][x].filename.find("%f"), 2, data[y][x].basename);
        per_tile_width               = std::max(per_tile_width, other_tile.viewWidth() );
        per_tile_height              = std::max(per_tile_height,other_tile.viewHeight());
        this_tile.create_with_width  = other_tile.viewWidth();
        this_tile.create_with_height = other_tile.viewHeight();
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

  T& operator()(size_t tx, size_t ty, size_t x, size_t y){
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<widthInTiles());
    assert(ty<heightInTiles());
    assert(widthInTiles()>0);
    assert(heightInTiles()>0);

    if(data[ty][tx].null_tile){
      no_data_to_set = data[ty][tx].noData();
      return no_data_to_set;
    }

    LoadTile(tx, ty);

    assert(x>=0);
    assert(y>=0);
    assert(x<data[ty][tx].viewWidth() );
    assert(y<data[ty][tx].viewHeight());

    return data[ty][tx](x,y);
  }

  T& operator()(size_t x, size_t y){
    assert(x>=0);
    assert(y>=0);
    assert(x<total_width_in_cells);
    assert(y<total_height_in_cells);

    //std::cerr<<"From ("<<x<<","<<y<<") derived ";

    size_t tile_x = x/per_tile_width;
    size_t tile_y = y/per_tile_height;
    x             = x%per_tile_width;
    y             = y%per_tile_height;

    //std::cerr<<"tile=("<<tile_x<<","<<tile_y<<") cell=("<<x<<","<<y<<")"<<std::endl;

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

  size_t widthInTiles() const {
    return data.back().size();
  }

  size_t heightInTiles() const {
    return data.size();
  }

  int64_t tileWidth(size_t tx, size_t ty) const {
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<data[0].size());
    assert(ty<data.size());
    return data[ty][tx].viewWidth();
  }

  int64_t tileHeight(size_t tx, size_t ty) const {
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<widthInTiles());
    assert(ty<heightInTiles());
    return data[ty][tx].viewHeight();
  }

  int64_t stdTileHeight() const {
    return per_tile_height;
  }

  int64_t stdTileWidth() const {
    return per_tile_width;
  }

  void setAll(const T &val){
    for(auto &row: data)
    for(auto &tile: row){
      tile.do_set_all  = true;
      tile.set_all_val = val;
    }
  }

  bool isNoData(size_t x, size_t y){
    assert(x>=0);
    assert(y>=0);
    assert(x<total_width_in_cells);
    assert(y<total_height_in_cells);

    size_t tile_x = x/per_tile_width;
    size_t tile_y = y/per_tile_height;
    x             = x%per_tile_width;
    y             = y%per_tile_height;

    if(data[tile_y][tile_x].null_tile)
      return true;

    LoadTile(tile_x, tile_y);

    return data[tile_y][tile_x].isNoData(x,y);
  }

  bool isNoData(size_t tx, size_t ty, size_t px, size_t py){
    assert(tx>=0);
    assert(ty>=0);
    assert(tx<data[0].size());
    assert(ty<data.size());

    assert(px>=0);
    assert(py>=0);
    assert(px<data[ty][tx].viewWidth() );
    assert(py<data[ty][tx].viewHeight());

    if(data[ty][tx].null_tile)
      return true;

    LoadTile(tx, ty);

    return data[ty][tx].isNoData(px,py);
  }

  bool in_grid(size_t x, size_t y) const {
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

      if((tile.geotransform[0]<0) ^ flipH)
        tile.flipHorz();
      if((tile.geotransform[5]<0) ^ flipV)
        tile.flipVert();

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
        tile.setNoData(ndval);
  }

  int32_t getEvictions() const {
    return evictions;
  }

  bool isNullTile(int tx, int ty) const {
    return data[ty][tx].null_tile;
  }

  bool isEdgeCell(size_t x, size_t y){
    return x==0 || y==0 || x==total_width_in_cells-1 || y==total_height_in_cells-1;
  }
};

#endif