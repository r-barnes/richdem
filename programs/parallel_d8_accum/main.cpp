//This algorithm is discussed in the manuscript:
//    Barnes, R., 2017. Parallel non-divergent flow accumulation for trillion 
//    cell digital elevation models on desktops or clusters. Environmental 
//    Modelling & Software 92, 202–212. doi:10.1016/j.envsoft.2017.02.022
#include "gdal_priv.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <stack>
#include <vector>
#include <limits>
#include <fstream> //For reading layout files
#include <sstream> //Used for parsing the <layout_file>
#include "richdem/common/version.hpp"
#include "richdem/common/Layoutfile.hpp"
#include "richdem/common/communication.hpp"
#include "richdem/common/memory.hpp"
#include "richdem/common/timer.hpp"
#include "richdem/common/Array2D.hpp"
#include "richdem/common/grid_cell.hpp"
#include "perimeters.hpp"

using namespace richdem;

const std::string algname  = "Barnes (2017) Parallel Non-divergent Flow Accumulation";
const std::string citation = "Barnes, R., 2017. Parallel non-divergent flow accumulation for trillion cell digital elevation models on desktops or clusters. Environmental Modelling & Software 92, 202-212. doi:10.1016/j.envsoft.2017.02.022";

//We use the cstdint library here to ensure that the program behaves as expected
//across platforms, especially with respect to the expected limits of operation
//for data set sizes and labels. For instance, in C++, a signed integer must be
//at least 16 bits, but not necessarily more. We force a minimum of 32 bits as
//this is, after all, for use with large datasets.
#include <cstdint>

//Define operating system appropriate directory separators
#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
  #define SLASH_CHAR "/"
#elif defined(__WIN32__)
  #define SLASH_CHAR "\\"
#endif

//TODO: Is it possible to run this without mpirun if we specify a single node
//job?

//TODO: What are these for?
const int TAG_WHICH_JOB   = 0;
const int TAG_TILE_DATA   = 1;
const int TAG_DONE_FIRST  = 2;
const int TAG_SECOND_DATA = 3;
const int TAG_DONE_SECOND = 4;

const int SYNC_MSG_KILL = 0;
const int JOB_FIRST     = 2;
const int JOB_SECOND    = 3;

const uint8_t FLIP_VERT   = 1;
const uint8_t FLIP_HORZ   = 2;

//Using int32_t here would allow a single cell to represent flow accumulations
//stemming from a ~46,340^2 cell area. Using uint32_t bumps that to ~65,535^2
//cells. This is easily exceeded in a rather large DEM. Therefore, we would like
//to use int64_t, which allows 3,037,000,499^2 cells. However, GDAL is silly and
//does not allow for 64-bit integers. Therefore, we use double. The IEEE754
//double-precision floating-point type has a 53-bit significand. This is
//sufficient to capture a 94,906,265^2 area.
typedef double accum_t;

//Valid flowdirs are in the range 0-8, inclusive. 0 is the center and 1-8,
//inclusive, are the 8 cells around the central cell. An extra value is needed
//to indicate NoData. Therefore, uint8_t is appropriate.
typedef uint8_t flowdir_t;

//Links need to be able to hold values up to the length of the perimeter of a
//tile of the DEM. uint16_t is therefore acceptable. It allows a perimeter of
//length 65,535 (minus the constants defined below). This allows for sides of
//~16,383.
typedef uint16_t link_t;

//Since, in a grid representation, up to eight cells may flow into a single
//cell, uint8_t is appropriate for appropriate for tracking these dependencies.
typedef uint8_t c_dependency_t;

//In a perimeter representation, any number of cells on the perimeter may
//ultimately flow into a single outlet cell. Therefore, we need this to uint16_t
//for the same reasons discussed for link_t. We will assume that the top value
//(65,535) is not usable in order to raise an assertion error if the dependency
//count goes negative.
typedef uint16_t p_dependency_t;

///Used in the perimeter representation to indicate that a cell has no
///downstream neighbour, for whatever reason
const link_t FLOW_NO_DOWNSTREAM = (link_t)-1;

///Used in the perimeter representation to indicate that flow going into this
///cell will end somewhere internal to the tile. This situation arises if the
///flow path leads to a NoData cell or a cell with NoFlow.
const link_t FLOW_TERMINATES = (link_t)-2;

///Used in the perimeter representation to indicate that this cell passes its
///flow to a neighbouring tile, or off of the DEM entirely.
const link_t FLOW_EXTERNAL   = (link_t)-3;

class TileInfo{
 private:
  friend class cereal::access;
  template<class Archive>
  void serialize(Archive & ar){
    ar(edge,
       flip,
       x,
       y,
       width,
       height,
       gridx,
       gridy,
       nullTile,
       filename,
       outputname,
       retention,
       many,
       analysis);
  }
 public:
  uint8_t     edge;
  uint8_t     flip;
  int32_t     x,y,gridx,gridy,width,height;
  bool        nullTile;
  bool        many;
  std::string filename;
  std::string outputname;
  std::string retention;
  std::string analysis;   //Command line command used to invoke everything
  TileInfo(){
    nullTile = true;
  }
  TileInfo(std::string filename, std::string outputname, std::string retention, int32_t gridx, int32_t gridy, int32_t x, int32_t y, int32_t width, int32_t height, bool many, std::string analysis){
    this->nullTile   = false;
    this->edge       = 0;
    this->x          = x;
    this->y          = y;
    this->width      = width;
    this->height     = height;   
    this->gridx      = gridx;
    this->gridy      = gridy;
    this->filename   = filename;
    this->outputname = outputname;
    this->retention  = retention;
    this->flip       = 0;
    this->many       = many;
    this->analysis   = analysis;
  }
};

typedef std::vector< std::vector< TileInfo > > TileGrid;


class TimeInfo {
 private:
  friend class cereal::access;
  template<class Archive>
  void serialize(Archive & ar){
    ar(calc,overall,io,vmpeak,vmhwm);
  }
 public:
  double calc, overall, io;
  long vmpeak, vmhwm;
  TimeInfo() {
    calc=overall=io=0;
    vmpeak=vmhwm=0;
  }
  TimeInfo(double calc, double overall, double io, long vmpeak, long vmhwm) :
      calc(calc), overall(overall), io(io), vmpeak(vmpeak), vmhwm(vmhwm) {}
  TimeInfo& operator+=(const TimeInfo& o){
    calc    += o.calc;
    overall += o.overall;
    io      += o.io;
    vmpeak   = std::max(vmpeak,o.vmpeak);
    vmhwm    = std::max(vmhwm,o.vmhwm);
    return *this;
  }
};

//TODO: Explain all of the variables
template<class elev_t>
class Job1 {
 private:
  friend class cereal::access;
  template<class Archive>
  void serialize(Archive & ar){
    ar(links,
       accum,
       flowdirs,
       time_info,
       gridy,
       gridx);
  }
 public:
  std::vector<link_t   >      links;
  std::vector<accum_t  >      accum;
  std::vector<flowdir_t>      flowdirs;
  std::vector<accum_t  >      accum_in;     //Used by produce, here for convenience, not communicated, so not serialized.
  std::vector<p_dependency_t> dependencies; //Used by produce, here for convenience, not communicated, so not serialized.
  TimeInfo time_info;
  int gridy, gridx;
  Job1(){}
};



template<class T> using Job1Grid    = std::vector< std::vector< Job1<T> > >;
template<class T> using Job2        = std::vector<accum_t>;
template<class T> using StorageType = std::map<std::pair<int,int>, std::pair< Array2D<flowdir_t>, Array2D<accum_t> > >;



int SuggestTileSize(int selected, int size, int min){
  int best=999999999;
  for(int x=1;x<size;x++)
    if(size%x>min && std::abs(x-selected)<std::abs(x-best))
      best=x;
  return best;
}


















template<class T>
class ConsumerSpecifics {
 public:
  Timer timer_io;
  Timer timer_calc;

 private:
  Array2D<flowdir_t>  flowdirs;
  Array2D<accum_t  >  accum;
  std::vector<link_t> links;


  //TODO: Check this description
  //This function takes a matrix of flow directions and an initial cell (x0,y0).
  //Starting at the initial cell, it follows the path described by the flow
  //direction matrix until it reaches an edge of the grid, a no_data cell, or a
  //cell which does not flow into any neighbour.

  //The initial cell is assumed to be on a top or bottom edge. If the flow path
  //terminates at a top or bottom edge, the initial cell is marked as flowing into
  //the terminal cell; otherwise, the initial cell is marked as terminating at an
  //unimportant location.

  //After running this function on all the top and bottom edge cells, we can
  //quickly determine the ultimate destination of any initial cell.
  template<class flowdir_t>
  void FollowPath(
    const int x0,                       //x-coordinate of initial cell
    const int y0,                       //y-coordinate of initial cell
    const TileInfo          &tile,      //Used to determine which tile we are in
    const Array2D<flowdir_t> &flowdirs, //Flow directions matrix
    std::vector<link_t>      &links
  ){
    int x = x0;
    int y = y0;

    int path_len = 0;

    int x0y0serial = xyToSerial(x0,y0,flowdirs.width(),flowdirs.height());

    const int max_path_length = flowdirs.width()*flowdirs.height();

    //Follow the flow path until it terminates
    while(path_len++<max_path_length){ //Follow the flow path until we reach its end
      const int n = flowdirs(x,y);     //Neighbour the current cell flows towards

      //Show the final part of the loop path
      if(path_len>max_path_length-10)
        std::cerr<<"E Path: "<<x<<","<<y<<" with flowdir "<<n<<std::endl;


      //If the neighbour this cell flows into is a no_data cell or this cell does
      //not flow into any neighbour, then mark the initial cell from which we
      //began this flow path as terminating somewhere unimportant: its flow cannot
      //pass to neighbouring segments/nodes for further processing.
      if(flowdirs.isNoData(x,y) || n==NO_FLOW){
        links.at(x0y0serial) = FLOW_TERMINATES;
        return;
      }

      //Flow direction was valid. Where does it lead?
      const int nx = x+dx[n]; //Get neighbour's x-coordinate.
      const int ny = y+dy[n]; //Get neighbour's y-coordinate.

      //The neighbour cell is off one of the sides of the tile. Therefore, its
      //flow may be passed on to a neighbouring tile. Thus, we need to link this
      //flow path's initial cell to this terminal cell.
      if(!flowdirs.inGrid(nx,ny)) {
        if(x==x0 && y==y0) {
          //If (x,y)==(x0,y0), then this is the first cell and it points off of
          //the grid. Mark it as being FLOW_EXTERNAL.
          links.at(xyToSerial(x0,y0,flowdirs.width(),flowdirs.height())) = FLOW_EXTERNAL;
        } else {
          //The flow path entered the grid one place and left another. We mark
          //only the start of the flow path since the end of the flow path will be
          //handled by another call to this function
          links.at(xyToSerial(x0,y0,flowdirs.width(),flowdirs.height())) = xyToSerial(x,y,flowdirs.width(),flowdirs.height());
        }
        return;
      }

      //The flow path has not yet terminated. Continue following it.
      x = nx;
      y = ny;
    }

    //The loop breaks with a return. This is only reached if more cells are
    //visited than are in the tile, which implies that a loop must exist.
    std::cerr<<"E File '"<<tile.filename<<"' contains a loop!"<<std::endl;
    throw std::logic_error("FollowPath() found a loop in the flow path!");
  }

  //As in the function above, we will start at an initial cell and follow its flow
  //direction to its neighbour. Then we will follow the neighbour's flow direction
  //to the neighbour's neighbour, and so on, until we can go no farther (as
  //indicated by running of one of the edges of the segment or hitting a no_data
  //cell).

  //However, at this point we know how much flow accumulation to add to each cell.
  //Therefore, we add this additional accumulation to each cell as we pass it.
  void FollowPathAdd(
    int x,                              //Initial x-coordinate
    int y,                              //Initial y-coordinate
    const Array2D<flowdir_t> &flowdirs, //Flow directions matrix
    Array2D<accum_t>         &accum,    //Output: Accumulation matrix
    const accum_t additional_accum      //Add this to every cell in the flow path
  ){
    //Follow the flow path until it terminates
    while(true){
      if(!flowdirs.inGrid(x,y))
        return;

      //Break when we reach a no_data cell
      if(flowdirs.isNoData(x,y))
        return;

      //Add additional flow accumulation to this cell
      accum(x,y) += additional_accum;

      int n = flowdirs(x,y); //Get neighbour
      if(n==NO_FLOW)          //This cell doesn't flow to a neighbour
        return;              

      x += dx[n];
      y += dy[n];
    }
  }


  void FlowAccumulation(
    const Array2D<flowdir_t> &flowdirs,
    Array2D<accum_t>         &accum
  ){
    typedef Array2D<flowdir_t>::i_t i_t;
    auto NO_I = Array2D<flowdir_t>::NO_I;

    //Each cell that flows points to a neighbouring cell. But not every cell
    //is pointed at. Cells which are not pointed at are the peaks from which
    //flow originates. In order to calculate the flow accumulation we begin at
    //peaks and pass flow downwards. Once flow has been passed downwards, the
    //cell receiving the flow no longer needs to wait for the cell which
    //passed the flow. When the receiving cell has no cells on which it is
    //waiting, it then becomes a peak itself. The number of cells pointing at
    //a cell is its "dependency count". In this section of the code we find
    //each cell's dependency count.

    accum.resize(flowdirs,0); //TODO: Is the initialization needed?
    accum.setNoData(ACCUM_NO_DATA);

    std::vector<c_dependency_t> dependencies(flowdirs.size(),0);

    for(i_t i=0;i<flowdirs.size();i++){
      if(flowdirs.isNoData(i)){    //This cell is a no_data cell
        accum(i) = ACCUM_NO_DATA;
        continue;                
      }         

      int n = flowdirs(i);         //The neighbour this cell flows into
      if(n==NO_FLOW)               //This cell does not flow into a neighbour
        continue;

      auto ni = flowdirs.getN(i,n);
        
      //Neighbour is not on the grid
      if(ni==NO_I)
        continue;

      //Neighbour is valid and is part of the grid. The neighbour depends on this
      //cell, so increment its dependency count.
      dependencies[ni]++;
    }

    //Now that we know how many dependencies each cell has, we can determine which
    //cells are the peaks: the sources of flow. We make a note of where the peaks
    //are for later use.
    std::stack<i_t> sources;
    for(i_t i=0;i<dependencies.size();i++)
      //Valid cell with no dependencies: a peak!
      if(dependencies[i]==0 && !flowdirs.isNoData(i))
        sources.emplace(i);

    //Now that we know where the sources of flow are, we can start at this. Each
    //cell will have at least an accumulation of 1: itself. It then passes this
    //and any other accumulation it has gathered along its flow path to its
    //neighbour and decrements the neighbours dependency count. When a neighbour
    //has no more dependencies, it becomes a source.
    while(!sources.empty()){         //There are sources remaining
      auto i = sources.top();        //Grab a source. Order is not important here.
      sources.pop();                 //We've visited this source. Discard it.

      if(flowdirs.isNoData(i))       //Oh snap! This isn't a real cell!
        continue;

      accum(i)++;                    //This is a real cell, and it accumulates
                                     //one cell's worth of flow automatically.

      int n = flowdirs(i);           //Who is this source's neighbour?

      if(n==NO_FLOW)                 //This cell doesn't flow anywhere.
        continue;                    //Move on to the next source.

      auto ni = flowdirs.getN(i,n);  //Okay, this cell is going somewhere.
                                     //Make a note of where

      //This cell flows of the edge of the grid. Move on to next source.
      if(ni==NO_I)
        continue;
      //This cell flows into a no_data cell. Move on to next source.
      if(flowdirs.isNoData(ni))
        continue;

      //This cell has a neighbour it flows into. Add to its accumulation.
      accum(ni) += accum(i);
      //Decrement the neighbour's dependencies.
      dependencies[ni]--;

      //The neighbour has no more dependencies, so it has become a source
      if(dependencies[ni]==0)
        sources.emplace(ni);
    }
  }

 public:
  void LoadFromEvict(const TileInfo &tile){
    #ifdef DEBUG
      std::cerr<<"d Grid tile: "<<tile.gridx<<","<<tile.gridy<<std::endl;
      std::cerr<<"d Opening "<<tile.filename<<" as flowdirs."<<std::endl;
    #endif

    timer_io.start();
    flowdirs = Array2D<flowdir_t>(tile.filename, false, tile.x, tile.y, tile.width, tile.height);

    //TODO: Figure out a clever way to allow tiles of different widths/heights
    if(flowdirs.width()!=tile.width){
      std::cerr<<"E Tile '"<<tile.filename<<"' had unexpected width. Found "<<flowdirs.width()<<" expected "<<tile.width<<std::endl;
      throw std::runtime_error("Unexpected width.");
    }

    if(flowdirs.height()!=tile.height){
      std::cerr<<"E Tile '"<<tile.filename<<"' had unexpected height. Found "<<flowdirs.height()<<" expected "<<tile.height<<std::endl;
      throw std::runtime_error("Unexpected height.");
    }

    flowdirs.printStamp(5,"LoadFromEvict() before reorientation");

    if(tile.flip & FLIP_VERT)
      flowdirs.flipVert();
    if(tile.flip & FLIP_HORZ)
      flowdirs.flipHorz();
    timer_io.stop();

    flowdirs.printStamp(5,"LoadFromEvict() after reorientation");

    timer_calc.start();
    FlowAccumulation(flowdirs,accum);
    timer_calc.stop();
  }

  void VerifyInputSanity(){
    //Let's double-check that the flowdirs are valid
    timer_calc.start();
    for(int32_t y=0;y<flowdirs.height();y++)
    for(int32_t x=0;x<flowdirs.width();x++)
      if(!flowdirs.isNoData(x,y) && !(1<=flowdirs(x,y) && flowdirs(x,y)<=8) && !(flowdirs(x,y)==NO_FLOW))
        throw std::domain_error("Invalid flow direction found: "+std::to_string(flowdirs(x,y)));
    timer_calc.stop();
  }

  void FirstRound(const TileInfo &tile, Job1<T> &job1){
    //-2 removes duplicate cells on vertical edges which would otherwise
    //overlap horizontal edges
    links.resize(2*flowdirs.width()+2*(flowdirs.height()-2), FLOW_TERMINATES);

    //The following considers corner cells more than once. That's okay, though:
    //the repeated effort merely produces the same results in the same places
    //twice and the loss in efficiency is small since at most 4 cells are
    //repeated.

    //If we are the top segment, nothing can flow into us, so we do not need
    //to know where flow paths originating at the top go to. On the otherhand,
    //if we are not the top segment, then consider each cell of the top row
    //and find out where its flow goes to.
    timer_calc.start();
    if(!(tile.edge & GRID_TOP))
      for(int32_t x=0;x<flowdirs.width();x++)
        FollowPath(x,0,tile,flowdirs,links);

    //If we are the bottom segment, nothing can flow into us, so we do not
    //need to know where flow paths originating at the bottom go to. On the
    //otherhand, if we are not the bottom segment, then consider each cell of
    //the bottom row and find out where its flow goes to.
    if(!(tile.edge & GRID_BOTTOM))
      for(int32_t x=0;x<flowdirs.width();x++)
        FollowPath(x,flowdirs.height()-1,tile,flowdirs,links);

    if(!(tile.edge & GRID_LEFT))
      for(int32_t y=0;y<flowdirs.height();y++)
        FollowPath(0,y,tile,flowdirs,links);

    if(!(tile.edge & GRID_RIGHT))
      for(int32_t y=0;y<flowdirs.height();y++)
        FollowPath(flowdirs.width()-1,y,tile,flowdirs,links);

    job1.links = std::move(links);

    //Construct output arrays
    GridPerimToArray(flowdirs, job1.flowdirs);
    GridPerimToArray(accum,    job1.accum   );
    timer_calc.stop();
  }

  void SecondRound(const TileInfo &tile, Job2<T> &job2){
    #ifdef DEBUG
      std::cerr<<"d SECOND ROUND"<<std::endl;
      std::cerr<<"d Grid tile: "<<tile.gridx<<","<<tile.gridy<<std::endl;
    #endif

    auto &accum_offset = job2;

    timer_calc.start();
    for(int s=0;s<(int)accum_offset.size();s++){
      if(accum_offset.at(s)==0)
        continue;
      int x,y;
      serialToXY(s, x, y, accum.width(), accum.height());
      FollowPathAdd(x,y,flowdirs,accum,accum_offset.at(s));
    }
    timer_calc.stop();

    //At this point we're done with the calculation! Boo-yeah!

    accum.printStamp(5,"Saving output before reorientation");

    timer_io.start();
    if(tile.flip & FLIP_HORZ)
      accum.flipHorz();
    if(tile.flip & FLIP_VERT)
      accum.flipVert();
    timer_io.stop();

    accum.printStamp(5,"Saving output after reorientation");

    timer_io.start();
    accum.saveGDAL(tile.outputname, tile.analysis, tile.x, tile.y);
    timer_io.stop();
  }

  void SaveToCache(const TileInfo &tile){
    timer_io.start();
    flowdirs.setCacheFilename(tile.retention+"-flowdirs.dat");
    accum.setCacheFilename(tile.retention+"-accum.dat");
    flowdirs.dumpData();
    accum.dumpData();
    timer_io.stop();
  }

  void LoadFromCache(const TileInfo &tile){
    timer_io.start();
    flowdirs = Array2D<flowdir_t>(tile.retention+"-flowdirs.dat", true);
    accum    = Array2D<accum_t  >(tile.retention+"-accum.dat",    true);
    timer_io.stop();
  }

  void SaveToRetain(TileInfo &tile, StorageType<T> &storage){
    timer_io.start();
    auto &temp  = storage[std::make_pair(tile.gridy,tile.gridx)];
    temp.first  = std::move(flowdirs);
    temp.second = std::move(accum);
    timer_io.stop();
  }

  void LoadFromRetain(TileInfo &tile, StorageType<T> &storage){
    timer_io.start();
    auto &temp = storage.at(std::make_pair(tile.gridy,tile.gridx));
    flowdirs   = std::move(temp.first);
    accum      = std::move(temp.second);
    timer_io.stop();
  }
};















template<class T>
class ProducerSpecifics {
 public:
  Timer timer_calc;
  Timer timer_io;

 private:
  Job1Grid<T> job2s_to_dist;

  void DownstreamCell(
    const Job1Grid<T> &jobs,
    const TileGrid   &tiles,
    const int gx,     //Grid tile x of cell we wish to find downstream cell of
    const int gy,     //Grid tile y of cell we wish to find downstream cell of
    const link_t s,   //Array index of the cell in the tile in question
    link_t &ns,       //Array index of downstream cell in the tile
    int &gnx,         //Grid tile x of downstream cell
    int &gny          //Grid tile y of downstream cell
  ){
    const auto &j = jobs.at(gy).at(gx);

    assert(j.links.size()==j.flowdirs.size());

    ns = FLOW_NO_DOWNSTREAM;

    gnx = gx;
    gny = gy;

    if(j.links.at(s)==FLOW_TERMINATES || j.flowdirs.at(s)==NO_FLOW){
      //Flow ends somewhere internal to the tile or this particular cell has no
      //flow direction. The TERMINATES should also take care of NoData flowdirs
      return;

    } else if(j.links.at(s)==FLOW_EXTERNAL){
      //Flow goes into a valid neighbouring tile
      const auto &c = tiles.at(gy).at(gx);

      int x,y;
      serialToXY(s,x,y,c.width,c.height);

      const int nfd = j.flowdirs.at(s);
      int nx        = x+dx[nfd];
      int ny        = y+dy[nfd];

      //Since this is a FLOW_EXTERNAL, this cell must point off of the grid
      assert(nx==-1 || ny==-1 || nx==c.width || ny==c.height);

      //Identify which neighbouring tile the flow will be going to
      if(nx==c.width)
        gnx += 1;
      else if(nx==-1)
        gnx -= 1;
      
      if(ny==c.height)
        gny += 1;
      else if(ny==-1)
        gny -= 1;

      //Ensure that we are within the bounds of the grid
      //NOTE: gridwidth=jobs.front().size() and gridheight=jobs.size()
      if(gnx<0 || gny<0 || gnx==(int)jobs.front().size() || gny==(int)jobs.size()){
        gnx = -1;
        gny = -1;
        ns  = FLOW_NO_DOWNSTREAM;
        return;
      }

      const auto &nc = tiles.at(gny).at(gnx);

      //Ensure that the neighbouring tile is not null
      if(nc.nullTile){
        gnx = -1;
        gny = -1;
        ns  = FLOW_NO_DOWNSTREAM;
        return;
      }

      //Now that we know the neighbouring tile, set the next-x and next-y
      //coordinates with reference to the bounds of that tile
      if(nx==c.width)
        nx = 0;
      else if(nx==-1)
        nx = nc.width-1;
      
      if(ny==c.height)
        ny = 0;
      else if(ny==-1)
        ny = nc.height-1;

      //Serialize the coordinates
      ns = xyToSerial(nx,ny,nc.width,nc.height);

      assert(ns<(int)jobs.at(gny).at(gnx).flowdirs.size());
    } else {
      //Flow goes to somewhere else on the perimeter of the same tile
      ns = j.links.at(s);
    }
  }
 public:
  void Calculations(
    TileGrid   &tiles,
    Job1Grid<T> &jobs1
  ){
    timer_calc.start();
    const int gridheight = tiles.size();
    const int gridwidth  = tiles[0].size();

    //Set initial values for all dependencies to zero
    for(auto &row: jobs1)
    for(auto &this_job: row){
      this_job.dependencies.resize(this_job.links.size(),0);      
      this_job.accum_in.resize(this_job.links.size(),0);      
    }

    //TODO: Used for manuscript
    // for(size_t y=0;y<jobs1.size();y++)
    // for(size_t x=0;x<jobs1[0].size();x++){
    //   std::cerr<<"AccumOrig: "<<x<<" "<<y<<" ";
    //   for(const auto x: jobs1[y][x].accum)
    //     std::cerr<<x<<" ";
    //   std::cerr<<"\n";
    // }

    std::cerr<<"p Calculating dependencies..."<<std::endl;
    for(int y=0;y<gridheight;y++)
    for(int x=0;x<gridwidth;x++){
      if(tiles.at(y).at(x).nullTile)
        continue;

      const int flowdir_size = jobs1.at(y).at(x).flowdirs.size();
      for(link_t s=0;s<flowdir_size;s++){
        //Using -1 will hopefully cause problems if these variables are misused
        link_t ns = FLOW_NO_DOWNSTREAM; //Next cell in a tile, using a serialized coordinate
        int gnx   = -1; //Next tile, x-coordinate
        int gny   = -1; //Next tile, y-coordinate
        DownstreamCell(jobs1, tiles, x, y, s, ns, gnx, gny);
        if(ns==FLOW_NO_DOWNSTREAM || tiles.at(gny).at(gnx).nullTile) 
          continue;

        jobs1.at(gny).at(gnx).dependencies.at(ns)++;
      }
    }

    //Used for storing the coordinates of cells in this perimeter representation
    //of the grid
    class atype {
     public:
      int gx;   //X-coordinate of tile
      int gy;   //Y-coordinate of tile
      link_t s; //Serialized coordinate of cell's position on the tile's perimeter
      atype(int gx, int gy, link_t s) : gx(gx), gy(gy), s(s) {};
    };

    //Stack of cells which are known to have no dependencies
    std::stack<atype> q;

    //Search for cells without dependencies, these will constitute the set from
    //which we will begin push accumulation to downstream cells
    for(int y=0;y<gridheight;y++)
    for(int x=0;x<gridwidth;x++){
      if(tiles[y][x].nullTile)
        continue;

      auto &this_job = jobs1.at(y).at(x);

      for(link_t s=0;s<this_job.dependencies.size();s++){
        if(this_job.dependencies.at(s)==0)
          q.emplace(x,y,s);

        //We are just receiving information about the tiles' perimeters now. Any
        //flow at an input cell has already been transferred to output cells.
        //So, to prevent double-counting, we zero the inputs.
        if(this_job.links.at(s)!=FLOW_EXTERNAL)  //TODO: NO_FLOW
          this_job.accum.at(s) = 0;
      }
    }

    //TODO: Used for manuscript
    // for(size_t y=0;y<jobs1.size();y++)
    // for(size_t x=0;x<jobs1[0].size();x++){
    //   std::cerr<<"Accum: "<<x<<" "<<y<<" ";
    //   for(const auto x: jobs1[y][x].accum)
    //     std::cerr<<x<<" ";
    //   std::cerr<<"\nLinks: "<<x<<" "<<y<<" ";
    //   for(const auto x: jobs1[y][x].links)
    //     std::cerr<<x<<" ";
    //   std::cerr<<"\n";
    // }

    std::cerr<<"m Peaks found in aggregated problem = "<<q.size()<<std::endl;

    //TODO: Detect loops!
    int processed_cells = 0;
    while(!q.empty()){
      atype    c  = q.top();
      Job1<T> &j  = jobs1.at(c.gy).at(c.gx);
      q.pop();
      processed_cells++;

      assert(!tiles.at(c.gy).at(c.gx).nullTile);

      link_t ns = FLOW_NO_DOWNSTREAM; //Initial value designed to cause devastation if misused
      int gnx   = -1; //Initial value designed to cause devastation if misused
      int gny   = -1; //Initial value designed to cause devastation if misused
      DownstreamCell(jobs1, tiles, c.gx, c.gy, c.s, ns, gnx, gny);
      if(ns==FLOW_NO_DOWNSTREAM || tiles.at(gny).at(gnx).nullTile)
        continue;

      jobs1.at(gny).at(gnx).accum.at(ns) += j.accum.at(c.s);
      if(gny!=c.gy || gnx!=c.gx)
        jobs1.at(gny).at(gnx).accum_in.at(ns) += j.accum.at(c.s);

      if( (--jobs1.at(gny).at(gnx).dependencies.at(ns))==0 )
        q.emplace(gnx,gny,ns);

      //This ensures that the dependency count is >=0 for both sign and unsigned
      //types.
      assert(jobs1.at(gny).at(gnx).dependencies.at(ns)!=(p_dependency_t)-1); 
    }

    std::cerr<<"m Perimeter cells processed = "<<processed_cells<<std::endl;

    job2s_to_dist = std::move(jobs1);
    timer_calc.stop();
  }

  Job2<T> DistributeJob2(const TileGrid &tiles, int tx, int ty){
    auto &this_job = job2s_to_dist.at(ty).at(tx);
    // for(size_t s=0;s<this_job.accum.size();s++)
    //   this_job.accum[s] -= this_job.accum_orig[s];

    //TODO: Used for manuscript
    // std::cerr<<"AccumEnd: "<<tx<<" "<<ty<<" "; //TODO
    // for(const auto x: this_job.accum_in)
    //   std::cerr<<x<<" ";
    // std::cerr<<"\n";

    return this_job.accum_in;
  }
};















































template<class T>
void Consumer(){
  TileInfo      tile;
  StorageType<T> storage;

  //Have the consumer process messages as long as they are coming using a
  //blocking receive to wait.
  while(true){
    // When probe returns, the status object has the size and other attributes
    // of the incoming message. Get the message size. TODO
    int the_job = CommGetTag(0);

    //This message indicates that everything is done and the Consumer should shut
    //down.
    if(the_job==SYNC_MSG_KILL){
      return;

    //This message indicates that the consumer should prepare to perform the
    //first part of the distributed Priority-Flood algorithm on an incoming job
    } else if (the_job==JOB_FIRST){
      Timer timer_overall;
      timer_overall.start();
      
      CommRecv(&tile, nullptr, 0);

      ConsumerSpecifics<T> consumer;
      Job1<T>              job1;

      job1.gridy = tile.gridy;
      job1.gridx = tile.gridx;

      consumer.LoadFromEvict(tile);
      consumer.VerifyInputSanity();

      consumer.FirstRound(tile, job1);

      if(tile.retention=="@evict"){
        //Nothing to do: it will all get overwritten
      } else if(tile.retention=="@retain"){
        consumer.SaveToRetain(tile,storage);
      } else {
        consumer.SaveToCache(tile);
      }

      timer_overall.stop();

      long vmpeak, vmhwm;
      ProcessMemUsage(vmpeak,vmhwm);

      job1.time_info = TimeInfo(consumer.timer_calc.accumulated(),timer_overall.accumulated(),consumer.timer_io.accumulated(),vmpeak,vmhwm);

      CommSend(&job1,nullptr,0,TAG_DONE_FIRST);
    } else if (the_job==JOB_SECOND){
      Timer timer_overall;
      timer_overall.start();

      ConsumerSpecifics<T> consumer;
      Job2<T>              job2;

      CommRecv(&tile, &job2, 0);

      //These use the same logic as the analogous lines above
      if(tile.retention=="@evict")
        consumer.LoadFromEvict(tile);
      else if(tile.retention=="@retain")
        consumer.LoadFromRetain(tile,storage);
      else
        consumer.LoadFromCache(tile);

      consumer.SecondRound(tile, job2);

      timer_overall.stop();

      long vmpeak, vmhwm;
      ProcessMemUsage(vmpeak,vmhwm);

      TimeInfo temp(consumer.timer_calc.accumulated(), timer_overall.accumulated(), consumer.timer_io.accumulated(),vmpeak,vmhwm);
      CommSend(&temp, nullptr, 0, TAG_DONE_SECOND);
    }
  }
}







//Producer takes a collection of Jobs and delegates them to Consumers. Once all
//of the jobs have received their initial processing, it uses that information
//to compute the global properties necessary to the solution. Each Job, suitably
//modified, is then redelegated to a Consumer which ultimately finishes the
//processing.
template<class T>
void Producer(TileGrid &tiles){
  Timer timer_overall;
  timer_overall.start();

  ProducerSpecifics<T> producer;

  const int gridheight = tiles.size();
  const int gridwidth  = tiles.front().size();

  //How many processes to send to
  const int active_consumer_limit = CommSize()-1;
  //Used to hold message buffers while non-blocking sends are used
  std::vector<msg_type> msgs;
  //Number of jobs for which we are waiting for a return
  int jobs_out=0;

  ////////////////////////////////////////////////////////////
  //SEND JOBS

  //Distribute jobs to the consumers. Since this is non-blocking, all of the
  //jobs will be sent and then we will wait to hear back below.
  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++){
    if(tiles[y][x].nullTile)
      continue;

    msgs.push_back(CommPrepare(&tiles.at(y).at(x),nullptr));
    CommISend(msgs.back(), (jobs_out%active_consumer_limit)+1, JOB_FIRST);
    jobs_out++;
  }

  std::cerr<<"m Jobs created = "<<jobs_out<<std::endl;

  //Grid to hold returned jobs
  Job1Grid<T> jobs1(tiles.size(), std::vector< Job1<T> >(tiles[0].size()));
  while(jobs_out--){
    std::cerr<<"p Jobs remaining = "<<jobs_out<<std::endl;
    Job1<T> temp;
    CommRecv(&temp, nullptr, -1);
    jobs1.at(temp.gridy).at(temp.gridx) = temp;
  }

  std::cerr<<"n First stage Tx = "<<CommBytesSent()<<" B"<<std::endl;
  std::cerr<<"n First stage Rx = "<<CommBytesRecv()<<" B"<<std::endl;
  CommBytesReset();

  //Get timing info
  TimeInfo time_first_total;
  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++)
    time_first_total += jobs1[y][x].time_info;


  ////////////////////////////////////////////////////////////
  //PRODUCER NODE PERFORMS PROCESSING ON ALL THE RETURNED DATA

  producer.Calculations(tiles,jobs1);

  ////////////////////////////////////////////////////////////
  //SEND OUT JOBS TO FINALIZE GLOBAL SOLUTION

  //Reset these two variables
  jobs_out = 0; 
  msgs     = std::vector<msg_type>();

  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++){
    if(tiles[y][x].nullTile)
      continue;

    auto job2 = producer.DistributeJob2(tiles, x, y);

    msgs.push_back(CommPrepare(&tiles.at(y).at(x),&job2));
    CommISend(msgs.back(), (jobs_out%active_consumer_limit)+1, JOB_SECOND);
    jobs_out++;
  }

  //There's no further processing to be done at this point, but we'll gather
  //timing and memory statistics from the consumers.
  TimeInfo time_second_total;

  while(jobs_out--){
    std::cerr<<"p Jobs left to receive = "<<jobs_out<<std::endl;
    TimeInfo temp;
    CommRecv(&temp, nullptr, -1);
    time_second_total += temp;
  }

  //Send out a message to tell the consumers to politely quit. Their job is
  //done.
  for(int i=1;i<CommSize();i++){
    int temp;
    CommSend(&temp,nullptr,i,SYNC_MSG_KILL);
  }

  timer_overall.stop();

  std::cerr<<"t First stage total overall time = "<<time_first_total.overall<<" s"<<std::endl;
  std::cerr<<"t First stage total IO time = "     <<time_first_total.io     <<" s"<<std::endl;
  std::cerr<<"t First stage total calc time = "   <<time_first_total.calc   <<" s"<<std::endl;
  std::cerr<<"r First stage peak child VmPeak = " <<time_first_total.vmpeak <<std::endl;
  std::cerr<<"r First stage peak child VmHWM = "  <<time_first_total.vmhwm  <<std::endl;

  std::cerr<<"n Second stage Tx = "<<CommBytesSent()<<" B"<<std::endl;
  std::cerr<<"n Second stage Rx = "<<CommBytesRecv()<<" B"<<std::endl;

  std::cerr<<"t Second stage total overall time = "<<time_second_total.overall<<" s"<<std::endl;
  std::cerr<<"t Second stage total IO time = "     <<time_second_total.io     <<" s"<<std::endl;
  std::cerr<<"t Second stage total calc time = "   <<time_second_total.calc   <<" s"<<std::endl;
  std::cerr<<"r Second stage peak child VmPeak = " <<time_second_total.vmpeak <<std::endl;
  std::cerr<<"r Second stage peak child VmHWM = "  <<time_second_total.vmhwm  <<std::endl;

  std::cerr<<"t Producer overall time = "<<timer_overall.accumulated()       <<" s"<<std::endl;
  std::cerr<<"t Producer calc time = "   <<producer.timer_calc.accumulated() <<" s"<<std::endl;

  long vmpeak, vmhwm;
  ProcessMemUsage(vmpeak,vmhwm);
  std::cerr<<"r Producer's VmPeak = "   <<vmpeak <<std::endl;
  std::cerr<<"r Producer's VmHWM = "    <<vmhwm  <<std::endl;
}



//Preparer divides up the input raster file into tiles which can be processed
//independently by the Consumers. Since the tileing may be done on-the-fly or
//rely on preparation the user has done, the Preparer routine knows how to deal
//with both. Once assemebled, the collection of jobs is passed off to Producer,
//which is agnostic as to the original form of the jobs and handles
//communication and solution assembly.
void Preparer(
  std::string many_or_one,
  const std::string retention,
  const std::string input_file,
  const std::string output_name,
  int bwidth,
  int bheight,
  int flipH,
  int flipV,
  std::string analysis
){
  Timer timer_overall;
  timer_overall.start();

  TileGrid tiles;
  std::string  filename;
  GDALDataType file_type;        //All tiles must have a common file_type
  TileInfo *reptile = nullptr; //Pointer to a representative tile

  std::string output_layout_name = output_name;
  if(output_name.find("%f")!=std::string::npos){
    output_layout_name.replace(output_layout_name.find("%f"), 2, "layout");
  } else if(output_name.find("%n")!=std::string::npos){
    output_layout_name.replace(output_layout_name.find("%n"), 2, "layout");
  } else { //Should never happen
    std::cerr<<"E Outputname for mode-many must contain '%f' or '%n'!"<<std::endl;
    throw std::runtime_error("Outputname for mode-many must contain '%f' or '%n'!");
  }
  LayoutfileWriter lfout(output_layout_name);

  if(many_or_one=="many"){
    int32_t tile_width     = -1; //Width of 1st tile. All tiles must equal this
    int32_t tile_height    = -1; //Height of 1st tile, all tiles must equal this
    long    cell_count     = 0;
    int     not_null_tiles = 0;
    std::vector<double> tile_geotransform(6);

    LayoutfileReader lf(input_file);

    while(lf.next()){
      if(lf.newRow()){ //Add a row to the grid of tiles
        tiles.emplace_back();
        lfout.addRow();
      }

      if(lf.isNullTile()){
        tiles.back().emplace_back();
        lfout.addEntry(""); //Add a null tile to the output
        continue;
      }

      not_null_tiles++;

      if(tile_height==-1){
        //Retrieve information about this tile. All tiles must have the same
        //dimensions, which we could check here, but opening and closing
        //thousands of files is expensive. Therefore, we rely on the user to
        //check this beforehand if they want to. We will, however, verify that
        //things are correct in Consumer() as we open the files for reading.
        try{
          getGDALDimensions(lf.getFullPath(),tile_height,tile_width,file_type,tile_geotransform.data());
        } catch (...) {
          std::cerr<<"E Error getting file information from '"<<lf.getFullPath()<<"'!"<<std::endl;
          CommAbort(-1); //TODO
        }
      }

      cell_count += tile_width*tile_height;

      std::string this_retention = retention;
      if(retention.find("%f")!=std::string::npos){
        this_retention.replace(this_retention.find("%f"), 2, lf.getBasename());          
      } else if(retention.find("%n")!=std::string::npos){
        this_retention.replace(this_retention.find("%n"), 2, lf.getGridLocName());
      } else if(retention[0]=='@') {
        this_retention = retention;
      } else { //Should never happen
        std::cerr<<"E Outputname for mode-many must contain '%f' or '%n'!"<<std::endl;
        throw std::runtime_error("Outputname for mode-many must contain '%f' or '%n'!");
      }

      std::string this_output_name = output_name;
      if(output_name.find("%f")!=std::string::npos){
        this_output_name.replace(this_output_name.find("%f"), 2, lf.getBasename());          
      } else if(output_name.find("%n")!=std::string::npos){
        this_output_name.replace(this_output_name.find("%n"), 2, lf.getGridLocName());
      } else { //Should never happen
        std::cerr<<"E Outputname for mode-many must contain '%f' or '%n'!"<<std::endl;
        throw std::runtime_error("Outputname for mode-many must contain '%f' or '%n'!");
      }

      tiles.back().emplace_back(
        lf.getFullPath(),
        this_output_name,
        this_retention,
        lf.getX(),
        lf.getY(),
        0,
        0,
        tile_width,
        tile_height,
        true,
        analysis
      );

      //Get a representative tile, if we don't already have one
      if(reptile==nullptr)
        reptile = &tiles.back().back();

      lfout.addEntry(this_output_name);

      //Flip tiles if the geotransform demands it
      if(tile_geotransform[1]<0)
        tiles.back().back().flip ^= FLIP_HORZ;
      if(tile_geotransform[5]>0)
        tiles.back().back().flip ^= FLIP_VERT;

      //Flip (or reverse the above flip!) if the user demands it
      if(flipH)
        tiles.back().back().flip ^= FLIP_HORZ;
      if(flipV)
        tiles.back().back().flip ^= FLIP_VERT;
    }

    std::cerr<<"c Loaded "<<tiles.size()<<" rows each of which had "<<tiles[0].size()<<" columns."<<std::endl;
    std::cerr<<"m Total cells to be processed = "<<cell_count<<std::endl;
    std::cerr<<"m Number of tiles which were not null = "<<not_null_tiles<<std::endl;

    //nullTiles imply that the tiles around them have edges, as though they
    //are on the edge of the raster.
    for(int y=0;y<(int)tiles.size();y++)
    for(int x=0;x<(int)tiles[0].size();x++){
      if(tiles[y][x].nullTile)
        continue;
      if(y-1>0 && x>0 && tiles[y-1][x].nullTile)
        tiles[y][x].edge |= GRID_TOP;
      if(y+1<(int)tiles.size() && x>0 && tiles[y+1][x].nullTile)
        tiles[y][x].edge |= GRID_BOTTOM;
      if(y>0 && x-1>0 && tiles[y][x-1].nullTile)
        tiles[y][x].edge |= GRID_LEFT;
      if(y>0 && x+1<(int)tiles[0].size() && tiles[y][x+1].nullTile)
        tiles[y][x].edge |= GRID_RIGHT;
    }

  } else if(many_or_one=="one") {
    int32_t total_height;
    int32_t total_width;

    //Get the total dimensions of the input file
    try {
      getGDALDimensions(input_file, total_height, total_width, file_type, NULL);
    } catch (...) {
      std::cerr<<"E Error getting file information from '"<<input_file<<"'!"<<std::endl;
      CommAbort(-1); //TODO
    }

    //If the user has specified -1, that implies that they want the entire
    //dimension of the raster along the indicated axis to be processed within a
    //single job.
    if(bwidth==-1)
      bwidth  = total_width;
    if(bheight==-1)
      bheight = total_height;

    std::cerr<<"m Total width =  "<<total_width <<"\n";
    std::cerr<<"m Total height = "<<total_height<<"\n";
    std::cerr<<"m Block width =  "<<bwidth      <<"\n";
    std::cerr<<"m Block height = "<<bheight     <<std::endl;
    std::cerr<<"m Total cells to be processed = "<<(total_width*total_height)<<std::endl;

    //Create a grid of jobs
    for(int32_t y=0,gridy=0;y<total_height; y+=bheight, gridy++){
      tiles.emplace_back(std::vector<TileInfo>());
      lfout.addRow();
      for(int32_t x=0,gridx=0;x<total_width;x+=bwidth,  gridx++){
        //if(total_height-y<100){
        //  std::cerr<<"E At least one tile is <100 cells in height. Please change rectangle size to avoid this!"<<std::endl;
        //  std::cerr<<"E I suggest you use bheight="<<SuggestTileSize(bheight, total_height, 100)<<std::endl;
        //  throw std::logic_error("Tile height too small!");
        //}
        //if(total_width -x<100){
        //  std::cerr<<"E At least one tile is <100 cells in width. Please change rectangle size to avoid this!"<<std::endl;
        //  std::cerr<<"E I suggest you use bwidth="<<SuggestTileSize(bwidth, total_width, 100)<<std::endl;
        //  throw std::logic_error("Tile width too small!");
        //}

        if(retention[0]!='@' && retention.find("%n")==std::string::npos){
          std::cerr<<"E In <one> mode '%n' must be present in the retention path."<<std::endl;
          throw std::invalid_argument("'%n' not found in retention path!");
        }

        if(output_name.find("%n")==std::string::npos){
          std::cerr<<"E In <one> mode '%n' must be present in the output path."<<std::endl;
          throw std::invalid_argument("'%n' not found in output path!");
        }

        //Used for '%n' formatting
        std::string coord_string = std::to_string(gridx)+"_"+std::to_string(gridy);

        std::string this_retention = retention;
        if(this_retention[0]!='@')
          this_retention.replace(this_retention.find("%n"), 2, coord_string);
        std::string this_output_name = output_name;
        if(this_output_name.find("%n")==std::string::npos){
          std::cerr<<"E Outputname must include '%n' for <one> mode."<<std::endl;
          throw std::runtime_error("Outputname must include '%n' for <one> mode.");
        }
        this_output_name.replace(this_output_name.find("%n"),2,coord_string);

        lfout.addEntry(this_output_name);

        tiles.back().emplace_back(
          input_file,
          this_output_name,
          this_retention,
          gridx,
          gridy,
          x,
          y,
          (total_width-x >=bwidth )?bwidth :total_width -x, //TODO: Check
          (total_height-y>=bheight)?bheight:total_height-y,
          false,
          analysis
        );
      }
    }

  } else {
    std::cout<<"Unrecognised option! Must be 'many' or 'one'!"<<std::endl;
    CommAbort(-1);
  }

  //If a job is on the edge of the raster, mark it as having this property so
  //that it can be handled with elegance later.
  for(auto &e: tiles.front())
    e.edge |= GRID_TOP;
  for(auto &e: tiles.back())
    e.edge |= GRID_BOTTOM;
  for(size_t y=0;y<tiles.size();y++){
    tiles[y].front().edge |= GRID_LEFT;
    tiles[y].back().edge  |= GRID_RIGHT;
  }

  CommBroadcast(&file_type,0);
  timer_overall.stop();
  std::cerr<<"t Preparer time = "<<timer_overall.accumulated()<<" s"<<std::endl;

  if(reptile!=nullptr){
    std::cerr<<"c Flip horizontal = "<<((reptile->flip & FLIP_HORZ)?"YES":"NO")<<std::endl;
    std::cerr<<"c Flip vertical =   "<<((reptile->flip & FLIP_VERT)?"YES":"NO")<<std::endl;
  }
  std::cerr<<"c Input data type = "<<GDALGetDataTypeName(file_type)<<std::endl;

  switch(file_type){
    case GDT_Byte:
      return Producer<uint8_t >(tiles);
    case GDT_UInt16:
      return Producer<uint16_t>(tiles);
    case GDT_Int16:
      return Producer<int16_t >(tiles);
    case GDT_UInt32:
      return Producer<uint32_t>(tiles);
    case GDT_Int32:
      return Producer<int32_t >(tiles);
    case GDT_Float32:
      return Producer<float   >(tiles);
    case GDT_Float64:
      return Producer<double  >(tiles);
    case GDT_CInt16:
    case GDT_CInt32:
    case GDT_CFloat32:
    case GDT_CFloat64:
      std::cerr<<"E Complex types are not supported. Sorry!"<<std::endl;
      CommAbort(-1); //TODO
    case GDT_Unknown:
    default:
      std::cerr<<"E Unrecognised data type: "<<GDALGetDataTypeName(file_type)<<std::endl;
      CommAbort(-1); //TODO
  }
}



int main(int argc, char **argv){
  CommInit(&argc,&argv);

  if(CommRank()==0){
    std::string many_or_one;
    std::string retention;
    std::string input_file;
    std::string output_name;
    int         bwidth    = -1;
    int         bheight   = -1;
    int         flipH     = false;
    int         flipV     = false;

    Timer timer_master;
    timer_master.start();

    std::string analysis = PrintRichdemHeader(argc,argv);

    std::cerr<<"A "<<algname <<std::endl;
    std::cerr<<"C "<<citation<<std::endl;

    std::string help=
    #include "help.txt"
    ;

    try{
      for(int i=1;i<argc;i++){
        if(strcmp(argv[i],"--bwidth")==0 || strcmp(argv[i],"-w")==0){
          if(i+1==argc)
            throw std::invalid_argument("-w followed by no argument.");
          bwidth = std::stoi(argv[i+1]);
          // if(bwidth<300 && bwidth!=-1) //TODO
          //   throw std::invalid_argument("Width must be at least 500.");
          i++;
          continue;
        } else if(strcmp(argv[i],"--bheight")==0 || strcmp(argv[i],"-h")==0){
          if(i+1==argc)
            throw std::invalid_argument("-h followed by no argument.");
          bheight = std::stoi(argv[i+1]);
          // if(bheight<300 && bheight!=-1) //TODO
          //   throw std::invalid_argument("Height must be at least 500.");
          i++;
          continue;
        } else if(strcmp(argv[i],"--help")==0){
          std::cerr<<help<<std::endl;
          int good_to_go=0;
          CommBroadcast(&good_to_go,0);
          CommFinalize();
          return -1;
        } else if(strcmp(argv[i],"--flipH")==0 || strcmp(argv[i],"-H")==0){
          flipH = true;
        } else if(strcmp(argv[i],"--flipV")==0 || strcmp(argv[i],"-V")==0){
          flipV = true;
        } else if(argv[i][0]=='-'){
          throw std::invalid_argument("Unrecognised flag: "+std::string(argv[i]));
        } else if(many_or_one==""){
          many_or_one = argv[i];
        } else if(retention==""){
          retention = argv[i];
        } else if(input_file==""){
          input_file = argv[i];
        } else if(output_name==""){
          output_name = argv[i];
        } else {
          throw std::invalid_argument("Too many arguments.");
        }
      }
      if(many_or_one=="" || retention=="" || input_file=="" || output_name=="")
        throw std::invalid_argument("Too few arguments.");
      if(retention[0]=='@' && !(retention=="@evict" || retention=="@retain"))
        throw std::invalid_argument("Retention must be @evict or @retain or a path.");
      if(many_or_one!="many" && many_or_one!="one")
        throw std::invalid_argument("Must specify many or one.");
      if(CommSize()==1) //TODO: Invoke a special "one-process mode?"
        throw std::invalid_argument("Must run program with at least two processes!");
      if( !((output_name.find("%f")==std::string::npos) ^ (output_name.find("%n")==std::string::npos)) )
        throw std::invalid_argument("Output filename must indicate either file number (%n) or name (%f).");
      if(retention[0]!='@' && retention.find("%n")==std::string::npos && retention.find("%f")==std::string::npos)
        throw std::invalid_argument("Retention filename must indicate file number with '%n' or '%f'.");
      if(retention==output_name)
        throw std::invalid_argument("Retention and output filenames must differ.");
    } catch (const std::invalid_argument &ia){
      std::string output_err;
      if(ia.what()==std::string("stoi"))
        output_err = "Invalid width or height.";
      else
        output_err = ia.what();

      std::cerr<<"parallel_d8_accum.exe [--flipV] [--flipH] [--bwidth #] [--bheight #] <many/one> <retention> <input> <output>"<<std::endl;
      std::cerr<<"\tUse '--help' to show help."<<std::endl;

      std::cerr<<"E "<<output_err<<std::endl;

      int good_to_go=0;
      CommBroadcast(&good_to_go,0);
      CommFinalize();
      return -1;
    }

    int good_to_go = 1;
    std::cerr<<"c Processes = "              <<CommSize()<<std::endl;
    std::cerr<<"c Many or one = "            <<many_or_one<<std::endl;
    std::cerr<<"c Input file = "             <<input_file<<std::endl;
    std::cerr<<"c Retention strategy = "     <<retention <<std::endl;
    std::cerr<<"c Block width = "            <<bwidth    <<std::endl;
    std::cerr<<"c Block height = "           <<bheight   <<std::endl;
    std::cerr<<"c Flip horizontal = "        <<flipH     <<std::endl;
    std::cerr<<"c Flip vertical = "          <<flipV     <<std::endl;

    #ifdef WITH_COMPRESSION
      std::cerr<<"c Cache compression = TRUE"<<std::endl;
    #else
      std::cerr<<"c Cache compression = FALSE"<<std::endl;
    #endif

    CommBroadcast(&good_to_go,0);
    Preparer(many_or_one, retention, input_file, output_name, bwidth, bheight, flipH, flipV, analysis);

    timer_master.stop();
    std::cerr<<"t Total wall-time = "<<timer_master.accumulated()<<" s"<<std::endl;

  } else {
    int good_to_go;
    CommBroadcast(&good_to_go,0);
    if(good_to_go){
      GDALDataType file_type;
      CommBroadcast(&file_type,0);
      switch(file_type){
        case GDT_Byte:
          Consumer<uint8_t >();break;
        case GDT_UInt16:
          Consumer<uint16_t>();break;
        case GDT_Int16:
          Consumer<int16_t >();break;
        case GDT_UInt32:
          Consumer<uint32_t>();break;
        case GDT_Int32:
          Consumer<int32_t >();break;
        case GDT_Float32:
          Consumer<float   >();break;
        case GDT_Float64:
          Consumer<double  >();break;
        default:
          return -1;
      }
    }
  }

  CommFinalize();

  return 0;
}