#include "gdal_priv.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <queue>
#include <vector>
#include <limits>
#include <fstream> //For reading layout files
#include <sstream> //Used for parsing the <layout_file>
#include "Layoutfile.hpp"
#include "communication.hpp"
#include "memory.hpp"
#include "timer.hpp"
#include "Array2D.hpp"
#include "grid_cell.hpp"

const char* program_version = "1";

//We use the cstdint library here to ensure that the program behaves as expected
//across platforms, especially with respect to the expected limits of operation
//for data set sizes and labels. For instance, in C++, a signed integer must be
//at least 16 bits, but not necessarily more. We force a minimum of 32 bits as
//this is, after all, for use with large datasets.
#include <cstdint>
//#define DEBUG 1

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
const int TAG_CHUNK_DATA  = 1;
const int TAG_DONE_FIRST  = 2;
const int TAG_SECOND_DATA = 3;
const int TAG_DONE_SECOND = 4;

const int SYNC_MSG_KILL = 0;
const int JOB_FIRST     = 2;
const int JOB_SECOND    = 3;

const uint8_t FLIP_VERT   = 1;
const uint8_t FLIP_HORZ   = 2;

#define FLOW_TERMINATES -3 //TODO: Explain and ensure it fits in link_t
#define FLOW_EXTERNAL   -4 //TODO: Explain and ensure it fits in link_t

typedef int32_t accum_t;
typedef uint8_t dependency_t;
typedef uint8_t flowdir_t;
typedef int32_t link_t;

class ChunkInfo{
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
       nullChunk,
       filename,
       outputname,
       retention,
       many);
  }
 public:
  uint8_t     edge;
  uint8_t     flip;
  int32_t     x,y,width,height,gridx,gridy;
  bool        nullChunk;
  bool        many;
  std::string filename;
  std::string outputname;
  std::string retention;
  ChunkInfo(){
    nullChunk = true;
  }
  ChunkInfo(std::string filename, std::string outputname, std::string retention, int32_t gridx, int32_t gridy, int32_t x, int32_t y, int32_t width, int32_t height, bool many){
    this->nullChunk  = false;
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
  }
};

typedef std::vector< std::vector< ChunkInfo > > ChunkGrid;


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
  std::vector<link_t   >    links;
  std::vector<accum_t  >    accum;
  std::vector<flowdir_t>    flowdirs;
  std::vector<accum_t  >    accum_in;     //Not communicated, so not serialized. Here for convenience.
  std::vector<dependency_t> dependencies; //Not communicated, so not serialized. Here for convenience.
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








//TODO: Explain
int xyToSerial(const int x, const int y, const int width, const int height){
  //Ensure cell is on the perimeter
  assert( (x==0 || x==width-1 || y==0 || y==height-1) && x>=0 && y>=0 && x<width && y<height);

  if(y==0)                         //Top row
    return x;

  if(x==width-1)                   //Right hand side
    return (width-1)+y;

  if(y==height-1)                  //Bottom-row
    return (width-1)+(height)+x;   

  return 2*(width-1)+(height-1)+y; //Left-hand side
}

//TODO: Explain
void serialToXY(const int serial, int &x, int &y, const int width, const int height){
  if(serial<width){                        //Top row
    x = serial;
    y = 0;
  } else if(serial<(width-1)+height){     //Right-hand side
    x = width-1;
    y = serial-(width-1);
  } else if(serial<2*(width-1)+(height)){ //Bottom row
    x = serial-(width-1)-(height-1)-1;
    y = height-1;
  } else {                                //Left-hand side
    x = 0;
    y = serial-2*(width-1)-(height-1);
  }

  //Ensure cell is on the perimeter
  assert( (x==0 || x==width-1 || y==0 || y==height-1) && x>=0 && y>=0 && x<width && y<height);
}


//TODO: Remove
template<class T>
void print2d(const Array2D<T> &arr){
  for(size_t y=0;y<arr.height();y++){
    for(size_t x=0;x<arr.width();x++)
      std::cerr<<std::setw(5)<<(int)arr(x,y);
    std::cerr<<std::endl;
  }
}

//TODO: Remove
// template<class T>
// void print2dradius(const Array2D<T> &arr, int xcen, int ycen, int radius){
//   int minx  = std::max(xcen-radius,0);
//   int maxx  = std::min(xcen+radius,arr.width()-1);
//   int miny  = std::max(ycen-radius,0);
//   int maxy  = std::min(ycen+radius,arr.height()-1);
//   for(int y=ymin;y<=ymax;y++){
//     for(int x=xmin;x<=xmax;x++)
//       std::cerr<<std::setw(5)<<(int)arr(x,y);
//     std::cerr<<std::endl;
//   }
// }

//TODO: Remove
template<class T>
void print1d(const std::vector<T> &v){
  for(int x=0;x<(int)v.size();x++)
    std::cerr<<std::setw(5)<<(int)v[x];
  std::cerr<<std::endl;
}

//TODO: Remove
template<class T>
void print1das2d(const std::vector<T> &v, const int width, const int height){
  for(int y=0;y<height;y++){
    for(int x=0;x<width;x++){
      if(! (x==0 || y==0 || x==width-1 || y==height-1))
        std::cerr<<std::setw(5)<<" ";
      else
        std::cerr<<std::setw(5)<<(int)v.at(xyToSerial(x,y,width,height));
    }
    std::cerr<<std::endl;
  }
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
    const ChunkInfo          &chunk,    //Used to determine which tile we are in
    const Array2D<flowdir_t> &flowdirs, //Flow directions matrix
    std::vector<link_t>      &links
  ){
    int x = x0;
    int y = y0;

    int path_len = 0;

    int x0y0serial = xyToSerial(x0,y0,flowdirs.width(),flowdirs.height());

    #ifdef DEBUG
      std::cerr<<std::endl<<"FP: "<<chunk.gridx<<","<<chunk.gridy<<": ("<<x0<<","<<y0<<")";
    #endif

    const int max_path_length = flowdirs.width()*flowdirs.height(); //TODO: Should this have +1?

    //Follow the flow path until it terminates
    while(path_len++<max_path_length){ //Follow the flow path until we reach its end
      const int n = flowdirs(x,y);     //Neighbour the current cell flows towards

      //Show the final part of the loop path (TODO)
      if(path_len>max_path_length-10)
        std::cerr<<"Path: "<<x<<","<<y<<" with flowdir "<<n<<std::endl;

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

      #ifdef DEBUG
        std::cerr<<" -- ("<<nx<<","<<ny<<")";
      #endif

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
    std::cerr<<"File "<<chunk.filename<<"("<<chunk.gridx<<","<<chunk.gridy<<") dimensions=("<<flowdirs.width()<<","<<flowdirs.height()<<") contains a loop!"<<std::endl;
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
    //Each cell that flows points to a neighbouring cell. But not every cell
    //is pointed at. Cells which are not pointed at are the peaks from which
    //flow originates. In order to calculate the flow accumulation we begin at
    //peaks and pass flow downwards. Once flow has been passed downwards, the
    //cell receiving the flow no longer needs to wait for the cell which
    //passed the flow. When the receiving cell has no cells on which it is
    //waiting, it then becomes a peak itself. The number of cells pointing at
    //a cell is its "dependency count". In this section of the code we find
    //each cell's dependency count.

    accum.resize(flowdirs,0);
    accum.setNoData(ACCUM_NO_DATA);

    Array2D<dependency_t> dependencies(flowdirs,0);
    for(size_t y=0;y<flowdirs.height();y++)
    for(size_t x=0;x<flowdirs.width();x++){
      if(flowdirs.isNoData(x,y)){  //This cell is a no_data cell
        accum(x,y) = ACCUM_NO_DATA;
        continue;                
      }         

      int n = flowdirs(x,y);       //The neighbour this cell flows into
      if(n==NO_FLOW)               //This cell does not flow into a neighbour
        continue;

      int nx = x+dx[n];            //x-coordinate of the neighbour
      int ny = y+dy[n];            //y-coordinate of the neighbour
        
      //Neighbour is not on the grid
      if(!flowdirs.inGrid(nx,ny))
        continue;

      //Neighbour is valid and is part of the grid. The neighbour depends on this
      //cell, so increment its dependency count.
      dependencies(nx,ny)++;
    }

    //Now that we know how many dependencies each cell has, we can determine which
    //cells are the peaks: the sources of flow. We make a note of where the peaks
    //are for later use.
    std::queue<GridCell> sources;
    for(size_t y=0;y<dependencies.height();y++)
    for(size_t x=0;x<dependencies.width();x++)
      //Valid cell with no dependencies: a peak!
      if(dependencies(x,y)==0 && !flowdirs.isNoData(x,y))
        sources.emplace(x,y);

    //Now that we know where the sources of flow are, we can start at this. Each
    //cell will have at least an accumulation of 1: itself. It then passes this
    //and any other accumulation it has gathered along its flow path to its
    //neighbour and decrements the neighbours dependency count. When a neighbour
    //has no more dependencies, it becomes a source.
    while(!sources.empty()){         //There are sources remaining
      GridCell c = sources.front();  //Grab a source. Order is not important here.
      sources.pop();                 //We've visited this source. Discard it.

      if(flowdirs.isNoData(c.x,c.y)) //Oh snap! This isn't a real cell!
        continue;

      accum(c.x,c.y)++;              //This is a real cell, and it accumulates
                                     //one cell's worth of flow automatically.

      int n = flowdirs(c.x,c.y);     //Who is this source's neighbour?

      if(n==NO_FLOW)                 //This cell doesn't flow anywhere.
        continue;                    //Move on to the next source.

      int nx = c.x+dx[n];            //Okay, this cell is going somewhere.
      int ny = c.y+dy[n];            //Make a note of where

      //This cell flows of the edge of the grid. Move on to next source.
      if(!flowdirs.inGrid(nx,ny))
        continue;
      //This cell flows into a no_data cell. Move on to next source.
      if(flowdirs.isNoData(nx,ny))
        continue;

      //This cell has a neighbour it flows into. Add to its accumulation.
      accum(nx,ny) += accum(c.x,c.y);
      //Decrement the neighbour's dependencies.
      dependencies(nx,ny)--;

      //The neighbour has no more dependencies, so it has become a source
      if(dependencies(nx,ny)==0)
        sources.emplace(nx,ny);
    }
  }

  template<class U>
  void GridPerimToArray(const Array2D<U> &grid, std::vector<U> &vec){
    assert(vec.size()==0); //Ensure receiving array is empty

    std::vector<U> vec2copy;

    vec2copy = grid.getRowData(0);                         //Top
    vec.insert(vec.end(),vec2copy.begin(),vec2copy.end());

    vec2copy = grid.getColData(grid.width()-1);        //Right
    vec.insert(vec.end(),vec2copy.begin()+1,vec2copy.end());
    
    vec2copy = grid.getRowData(grid.height()-1);       //Bottom
    vec.insert(vec.end(),vec2copy.begin(),vec2copy.end()-1);
    
    vec2copy = grid.getColData(0);                         //Left
    vec.insert(vec.end(),vec2copy.begin()+1,vec2copy.end()-1);
  }

 public:
  void LoadFromEvict(const ChunkInfo &chunk){
    #ifdef DEBUG
      std::cerr<<"Grid tile: "<<chunk.gridx<<","<<chunk.gridy<<std::endl;
      std::cerr<<"Opening "<<chunk.filename<<" as flowdirs."<<std::endl;
    #endif

    timer_io.start();
    flowdirs = Array2D<flowdir_t>(chunk.filename, false, chunk.x, chunk.y, chunk.width, chunk.height);

    //TODO: Figure out a clever way to allow tiles of different widths/heights
    if(flowdirs.width()!=chunk.width){
      std::cerr<<"Tile '"<<chunk.filename<<"' had unexpected width. Found "<<flowdirs.width()<<" expected "<<chunk.width<<std::endl;
      throw std::runtime_error("Unexpected width.");
    }

    if(flowdirs.height()!=chunk.height){
      std::cerr<<"Tile '"<<chunk.filename<<"' had unexpected height. Found "<<flowdirs.height()<<" expected "<<chunk.height<<std::endl;
      throw std::runtime_error("Unexpected height.");
    }

    std::cerr<<"Loading with width="<<chunk.width<<", height="<<chunk.height<<std::endl;

    #ifdef DEBUG
     std::cerr<<"Flowdirs raw: "<<std::endl;
     print2d(flowdirs);
    #endif

    flowdirs.printStamp(5,"LoadFromEvict() before reorientation");

    if(chunk.flip & FLIP_VERT)
      flowdirs.flipVert();
    if(chunk.flip & FLIP_HORZ)
      flowdirs.flipHorz();
    timer_io.stop();

    flowdirs.printStamp(5,"LoadFromEvict() after reorientation");

    timer_calc.start();
    FlowAccumulation(flowdirs,accum);
    timer_calc.stop();

    #ifdef DEBUG
      std::cerr<<"Accum first: "<<std::endl;
      print2d(accum);
    #endif
  }

  void VerifyInputSanity(){
    //Let's double-check that the flowdirs are valid
    for(size_t y=0;y<flowdirs.height();y++)
    for(size_t x=0;x<flowdirs.width();x++)
      if(!flowdirs.isNoData(x,y) && !(1<=flowdirs(x,y) && flowdirs(x,y)<=8) && !(flowdirs(x,y)==NO_FLOW))
        throw std::domain_error("Invalid flow direction found: "+std::to_string(flowdirs(x,y)));
  }

  void FirstRound(const ChunkInfo &chunk, Job1<T> &job1){
    std::cerr<<"Grid tile: "<<chunk.gridx<<","<<chunk.gridy<<std::endl;

    //-2 removes duplicate cells on vertical edges which would otherwise
    //overlap horizontal edges
    links.resize(2*flowdirs.width()+2*(flowdirs.height()-2), FLOW_TERMINATES);

    //TODO: Although the following may consider a cell more than once, the
    //repeated effort merely produces the same results in the same places
    //twice, so it's okay.

    //If we are the top segment, nothing can flow into us, so we do not need
    //to know where flow paths originating at the top go to. On the otherhand,
    //if we are not the top segment, then consider each cell of the top row
    //and find out where its flow goes to.
    timer_calc.start();
    if(!(chunk.edge & GRID_TOP))
      for(size_t x=0;x<flowdirs.width();x++)
        FollowPath(x,0,chunk,flowdirs,links);

    //If we are the bottom segment, nothing can flow into us, so we do not
    //need to know where flow paths originating at the bottom go to. On the
    //otherhand, if we are not the bottom segment, then consider each cell of
    //the bottom row and find out where its flow goes to.
    if(!(chunk.edge & GRID_BOTTOM))
      for(size_t x=0;x<flowdirs.width();x++)
        FollowPath(x,flowdirs.height()-1,chunk,flowdirs,links);

    if(!(chunk.edge & GRID_LEFT))
      for(size_t y=0;y<flowdirs.height();y++)
        FollowPath(0,y,chunk,flowdirs,links);

    if(!(chunk.edge & GRID_RIGHT))
      for(size_t y=0;y<flowdirs.height();y++)
        FollowPath(flowdirs.width()-1,y,chunk,flowdirs,links);

    job1.links = std::move(links);

    //Construct output arrays
    GridPerimToArray(flowdirs, job1.flowdirs);
    GridPerimToArray(accum,    job1.accum   );
    timer_calc.stop();

    #ifdef DEBUG
      std::cerr<<job1.gridx<<","<<job1.gridy<<" Links: ";
      for(auto const x: job1.links)
        std::cerr<<x<<" ";
      std::cerr<<std::endl;

      std::cerr<<job1.gridx<<","<<job1.gridy<<" Flowdirs: ";
      for(auto const x: job1.flowdirs)
        std::cerr<<x<<" ";
      std::cerr<<std::endl;

      std::cerr<<job1.gridx<<","<<job1.gridy<<" Accum: ";
      for(auto const x: job1.accum)
        std::cerr<<x<<" ";
      std::cerr<<std::endl;
    #endif

    //TODO: Remove
    #ifdef DEBUG
      std::cerr<<"Flowdirs: "<<std::endl;
      print2d(flowdirs);
      std::cerr<<"Accum: "<<std::endl;
      print2d(accum);
      std::cerr<<"Flowdirs Perim: "<<std::endl;
      print1d(job1.flowdirs);
      std::cerr<<"Accum perim: "<<std::endl;
      print1das2d(job1.accum,flowdirs.width(),flowdirs.height());
      std::cerr<<"Links: "<<std::endl;
      print1das2d(job1.links,flowdirs.width(),flowdirs.height());
    #endif
  }

  void SecondRound(const ChunkInfo &chunk, Job2<T> &job2){
    #ifdef DEBUG
      std::cerr<<"SECOND ROUND"<<std::endl;
      std::cerr<<"Grid tile: "<<chunk.gridx<<","<<chunk.gridy<<std::endl;
    #endif

    auto &accum_offset = job2;

    #ifdef DEBUG
      std::cerr<<"Received increments: "<<std::endl;
      print1das2d(accum_offset,flowdirs.width(),flowdirs.height());

      std::cerr<<"Accumulation"<<std::endl;
      print2d(accum);

      std::cerr<<"Adjusted accumulation"<<std::endl;
      print2d(accum);
    #endif

    for(int s=0;s<(int)accum_offset.size();s++){
      if(accum_offset.at(s)==0)
        continue;
      int x,y;
      serialToXY(s, x, y, accum.width(), accum.height());
      FollowPathAdd(x,y,flowdirs,accum,accum_offset.at(s));
    }

    //TODO: Remove
    #ifdef DEBUG
      std::cerr<<"Final accumulation"<<std::endl;
      print2d(accum);
    #endif

    //At this point we're done with the calculation! Boo-yeah!

    accum.printStamp(5,"Saving output before reorientation");

    timer_io.start();
    if(chunk.flip & FLIP_HORZ)
      accum.flipHorz();
    if(chunk.flip & FLIP_VERT)
      accum.flipVert();
    timer_io.stop();

    accum.printStamp(5,"Saving output after reorientation");

    timer_io.start();
    accum.saveGDAL(chunk.outputname, chunk.x, chunk.y);
    timer_io.stop();
  }

  void SaveToCache(const ChunkInfo &chunk){
    timer_io.start();
    flowdirs.setCacheFilename(chunk.retention+"-flowdirs.dat");
    accum.setCacheFilename(chunk.retention+"-accum.dat");
    flowdirs.dumpData();
    accum.dumpData();
    timer_io.stop();
  }

  void LoadFromCache(const ChunkInfo &chunk){
    timer_io.start();
    flowdirs = Array2D<flowdir_t>(chunk.retention+"-flowdirs.dat", true);
    accum    = Array2D<accum_t  >(chunk.retention+"-accum.dat",    true);
    timer_io.stop();
  }

  void SaveToRetain(ChunkInfo &chunk, StorageType<T> &storage){
    auto &temp  = storage[std::make_pair(chunk.gridy,chunk.gridx)];
    temp.first  = std::move(flowdirs);
    temp.second = std::move(accum);
  }

  void LoadFromRetain(ChunkInfo &chunk, StorageType<T> &storage){
    auto &temp = storage.at(std::make_pair(chunk.gridy,chunk.gridx));
    flowdirs   = std::move(temp.first);
    accum      = std::move(temp.second);
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
    const ChunkGrid   &chunks,
    const int gx,     //Grid tile x of cell we wish to find downstream cell of
    const int gy,     //Grid tile y of cell we wish to find downstream cell of
    const int s,      //Array index of the cell in the tile in question
    int &ns,          //Array index of downstream cell in the tile
    int &gnx,         //Grid tile x of downstream cell
    int &gny          //Grid tile y of downstream cell
  ){
    const auto &j = jobs.at(gy).at(gx);

    assert(j.links.size()==j.flowdirs.size());

    ns = -1;

    gnx = gx;
    gny = gy;

    if(j.links.at(s)==FLOW_TERMINATES || j.flowdirs.at(s)==NO_FLOW){
      //Flow ends somewhere internal to the tile or this particular cell has no
      //flow direction. The TERMINATES should also take care of NoData flowdirs
      return;

    } else if(j.links.at(s)==FLOW_EXTERNAL){
      //Flow goes into a valid neighbouring tile
      const auto &c = chunks.at(gy).at(gx);

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

      //Ensure that the neighbouring tile is valid (TODO: nullChunk here?)
      //NOTE: gridwidth=jobs.front().size() and gridheight=jobs.size()
      if(gnx<0 || gny<0 || gnx==(int)jobs.front().size() || gny==(int)jobs.size()){
        gnx=gny=ns=-1;
        return;
      }

      const auto &nc = chunks.at(gny).at(gnx);

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

      #ifdef DEBUG
        //Ensure that the serialization is valid
        if(ns>=(int)jobs.at(gny).at(gnx).flowdirs.size()){
          std::cerr<<nx<<" "<<ny<<" "<<nc.width<<" "<<nc.height<<std::endl;
          std::cerr<<ns<<std::endl;
          std::cerr<<(int)jobs.at(gny).at(gnx).flowdirs.size()<<std::endl;
        }
      #endif
      assert(ns<(int)jobs.at(gny).at(gnx).flowdirs.size());
    } else {
      //Flow goes to somewhere else on the perimeter of the same tile
      ns = j.links.at(s);
    }
  }
 public:
  void Calculations(
    ChunkGrid   &chunks,
    Job1Grid<T> &jobs1
  ){
    const int gridheight = chunks.size();
    const int gridwidth  = chunks[0].size();

    #ifdef DEBUG
      std::cerr<<std::endl;
      std::cerr<<std::endl;
      std::cerr<<"========================="<<std::endl;
      std::cerr<<"========================="<<std::endl;
      std::cerr<<"======PRODUCER==========="<<std::endl;
      std::cerr<<"========================="<<std::endl;
      std::cerr<<"========================="<<std::endl;
      std::cerr<<std::endl;
      std::cerr<<std::endl;
    #endif

    //Set initial values for all dependencies to zero
    for(auto &row: jobs1)
    for(auto &this_job: row){
      this_job.dependencies.resize(this_job.links.size(),0);      
      this_job.accum_in.resize(this_job.links.size(),0);      
    }

    std::cerr<<"Calculating dependencies..."<<std::endl;
    for(int y=0;y<gridheight;y++)
    for(int x=0;x<gridwidth;x++){
      auto &this_chunk = chunks.at(y).at(x);
      if(this_chunk.nullChunk)
        continue;

      auto &this_job = jobs1.at(y).at(x);

      for(int s=0;s<(int)this_job.flowdirs.size();s++){
        int ns  = -1; //Initial values designed to cause devastation if misused
        int gnx = -1;
        int gny = -1;
        DownstreamCell(jobs1, chunks, x, y, s, ns, gnx, gny);
        if(ns==-1 || chunks.at(gny).at(gnx).nullChunk) 
          continue;

        jobs1.at(gny).at(gnx).dependencies.at(ns)++;
      }
    }

    //TODO: Remove
    #ifdef DEBUG
      for(int y=0;y<gridheight;y++)
      for(int x=0;x<gridwidth;x++){
        std::cerr<<"Dependencies of chunk ("<<x<<","<<y<<")"<<std::endl;
        print1das2d(jobs1.at(y).at(x).dependencies,chunks.at(y).at(x).width,chunks.at(y).at(x).height);
      }
    #endif

    class atype {
     public:
      int gx,gy,s;
      atype(int gx, int gy, int s) : gx(gx), gy(gy), s(s) {};
    };

    std::queue<atype> q;

    //Search for cells without dependencies
    for(int y=0;y<gridheight;y++)
    for(int x=0;x<gridwidth;x++){
      if(chunks[y][x].nullChunk)
        continue;

      auto &this_job = jobs1.at(y).at(x);

      for(size_t s=0;s<this_job.dependencies.size();s++){
        if(this_job.dependencies.at(s)==0) // && jobs1.at(y).at(x).links.at(s)==FLOW_EXTERNAL) //TODO
          q.emplace(x,y,s);

        //Accumulated flow at an input will be transfered to an output resulting
        //in double-counting (TODO: More logical place to put this?)
        if(this_job.links.at(s)!=FLOW_EXTERNAL)  //TODO: NO_FLOW
          this_job.accum.at(s) = 0;
      }

      //this_job.accum_orig = this_job.accum;
    }

    std::cerr<<q.size()<<" peaks founds in aggregated problem."<<std::endl;

    int processed_peaks = 0;
    while(!q.empty()){
      atype      c  = q.front();
      Job1<T>   &j  = jobs1.at(c.gy).at(c.gx);
      q.pop();
      processed_peaks++;

      //TODO: Cut
      ChunkInfo &ci = chunks.at(c.gy).at(c.gx);
      assert(!ci.nullChunk);

      int ns  = -1; //Initial value designed to cause devastation if misused
      int gnx = -1; //Initial value designed to cause devastation if misused
      int gny = -1; //Initial value designed to cause devastation if misused
      DownstreamCell(jobs1, chunks, c.gx, c.gy, c.s, ns, gnx, gny);
      if(ns==-1 || chunks.at(gny).at(gnx).nullChunk)
        continue;

      jobs1.at(gny).at(gnx).accum.at(ns) += j.accum.at(c.s);
      if(gny!=c.gy || gnx!=c.gx)
        jobs1.at(gny).at(gnx).accum_in.at(ns) += j.accum.at(c.s);

      if( (--jobs1.at(gny).at(gnx).dependencies.at(ns))==0 )
        q.emplace(gnx,gny,ns);
    }

    std::cerr<<processed_peaks<<" peaks processed."<<std::endl;

    //TODO: Remove this block as this should always be zero once I have things
    //right.
    #ifdef DEBUG
      int unprocessed_dependencies=0;
      for(int y=0;y<gridheight;y++)
      for(int x=0;x<gridwidth;x++)
      for(size_t s=0;s<jobs1[y][x].dependencies.size();s++)
        unprocessed_dependencies+=jobs1[y][x].dependencies[s];
      std::cerr<<unprocessed_dependencies<<" unprocessed dependencies."<<std::endl;
    #endif

    job2s_to_dist = std::move(jobs1);

    #ifdef DEBUG
      std::cerr<<std::endl;
      std::cerr<<std::endl;
      std::cerr<<"========================="<<std::endl;
      std::cerr<<"========================="<<std::endl;
      std::cerr<<"==========DONE==========="<<std::endl;
      std::cerr<<"========================="<<std::endl;
      std::cerr<<"========================="<<std::endl;
      std::cerr<<std::endl;
      std::cerr<<std::endl;
    #endif
  }

  Job2<T> DistributeJob2(const ChunkGrid &chunks, int tx, int ty){
    auto &this_job = job2s_to_dist.at(ty).at(tx);
    // for(size_t s=0;s<this_job.accum.size();s++)
    //   this_job.accum[s] -= this_job.accum_orig[s];

    return this_job.accum_in;
  }
};















































template<class T>
void Consumer(){
  ChunkInfo      chunk;
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
      
      CommRecv(&chunk, nullptr, 0);

      ConsumerSpecifics<T> consumer;
      Job1<T>              job1;

      job1.gridy = chunk.gridy;
      job1.gridx = chunk.gridx;

      consumer.LoadFromEvict(chunk);
      consumer.VerifyInputSanity();

      consumer.FirstRound(chunk, job1);

      if(chunk.retention=="@evict"){
        //Nothing to do: it will all get overwritten
      } else if(chunk.retention=="@retain"){
        consumer.SaveToRetain(chunk,storage);
      } else {
        consumer.SaveToCache(chunk);
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

      CommRecv(&chunk, &job2, 0);

      //These use the same logic as the analogous lines above
      if(chunk.retention=="@evict")
        consumer.LoadFromEvict(chunk);
      else if(chunk.retention=="@retain")
        consumer.LoadFromRetain(chunk,storage);
      else
        consumer.LoadFromCache(chunk);

      consumer.SecondRound(chunk, job2);

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
void Producer(ChunkGrid &chunks){
  Timer timer_overall;
  timer_overall.start();

  ProducerSpecifics<T> producer;

  const int gridheight = chunks.size();
  const int gridwidth  = chunks.front().size();

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
    if(chunks[y][x].nullChunk)
      continue;

    msgs.push_back(CommPrepare(&chunks.at(y).at(x),nullptr));
    CommISend(msgs.back(), (jobs_out%active_consumer_limit)+1, JOB_FIRST);
    jobs_out++;
  }

  std::cerr<<jobs_out<<" jobs out."<<std::endl;

  //Grid to hold returned jobs
  Job1Grid<T> jobs1(chunks.size(), std::vector< Job1<T> >(chunks[0].size()));
  while(jobs_out--){
    std::cerr<<jobs_out<<" jobs left to receive."<<std::endl;
    Job1<T> temp;
    CommRecv(&temp, nullptr, -1);
    jobs1.at(temp.gridy).at(temp.gridx) = temp;
  }

  std::cerr<<"!First stage: "<<CommBytesSent()<<"B sent."<<std::endl;
  std::cerr<<"!First stage: "<<CommBytesRecv()<<"B received."<<std::endl;
  CommBytesReset();

  //Get timing info
  TimeInfo time_first_total;
  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++)
    time_first_total += jobs1[y][x].time_info;


  ////////////////////////////////////////////////////////////
  //PRODUCER NODE PERFORMS PROCESSING ON ALL THE RETURNED DATA

  producer.Calculations(chunks,jobs1);

  ////////////////////////////////////////////////////////////
  //SEND OUT JOBS TO FINALIZE GLOBAL SOLUTION

  //Reset these two variables
  jobs_out = 0; 
  msgs     = std::vector<msg_type>();

  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++){
    if(chunks[y][x].nullChunk)
      continue;

    auto job2 = producer.DistributeJob2(chunks, x, y);

    msgs.push_back(CommPrepare(&chunks.at(y).at(x),&job2));
    CommISend(msgs.back(), (jobs_out%active_consumer_limit)+1, JOB_SECOND);
    jobs_out++;
  }

  //There's no further processing to be done at this point, but we'll gather
  //timing and memory statistics from the consumers.
  TimeInfo time_second_total;

  while(jobs_out--){
    std::cerr<<jobs_out<<" jobs left to receive."<<std::endl;
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

  std::cout<<"!TimeInfo: First stage total overall time="<<time_first_total.overall<<std::endl;
  std::cout<<"!TimeInfo: First stage total io time="     <<time_first_total.io     <<std::endl;
  std::cout<<"!TimeInfo: First stage total calc time="   <<time_first_total.calc   <<std::endl;
  std::cout<<"!TimeInfo: First stage peak child VmPeak=" <<time_first_total.vmpeak <<std::endl;
  std::cout<<"!TimeInfo: First stage peak child VmHWM="  <<time_first_total.vmhwm  <<std::endl;

  std::cout<<"!Second stage: "<<CommBytesSent()<<"B sent."<<std::endl;
  std::cout<<"!Second stage: "<<CommBytesRecv()<<"B received."<<std::endl;

  std::cout<<"!TimeInfo: Second stage total overall time="<<time_second_total.overall<<std::endl;
  std::cout<<"!TimeInfo: Second stage total IO time="     <<time_second_total.io     <<std::endl;
  std::cout<<"!TimeInfo: Second stage total calc time="   <<time_second_total.calc   <<std::endl;
  std::cout<<"!TimeInfo: Second stage peak child VmPeak=" <<time_second_total.vmpeak <<std::endl;
  std::cout<<"!TimeInfo: Second stage peak child VmHWM="  <<time_second_total.vmhwm  <<std::endl;

  std::cout<<"!TimeInfo: Producer overall="<<timer_overall.accumulated()       <<std::endl;
  std::cout<<"!TimeInfo: Producer calc="   <<producer.timer_calc.accumulated() <<std::endl;

  long vmpeak, vmhwm;
  ProcessMemUsage(vmpeak,vmhwm);
  std::cout<<"!TimeInfo: Producer's VmPeak="   <<vmpeak <<std::endl;
  std::cout<<"!TimeInfo: Producer's VmHWM="    <<vmhwm  <<std::endl;
}



//Preparer divides up the input raster file into chunks which can be processed
//independently by the Consumers. Since the chunking may be done on-the-fly or
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
  int flipV
){
  Timer overall;
  overall.start();

  ChunkGrid chunks;
  std::string  filename;
  GDALDataType file_type; //All chunks must have a common file_type

  std::string output_layout_name = output_name;
  if(output_name.find("%f")!=std::string::npos){
    output_layout_name.replace(output_layout_name.find("%f"), 2, "layout");
  } else if(output_name.find("%n")!=std::string::npos){
    output_layout_name.replace(output_layout_name.find("%n"), 2, "layout");
  } else { //Should never happen
    std::cerr<<"Outputname for mode-many must contain '%f' or '%n'!"<<std::endl;
    throw std::runtime_error("Outputname for mode-many must contain '%f' or '%n'!");
  }
  LayoutfileWriter lfout(output_layout_name);

  if(many_or_one=="many"){
    int32_t chunk_width    = -1; //Width of 1st chunk. All chunks must equal this
    int32_t chunk_height   = -1; //Height of 1st chunk, all chunks must equal this
    long    cell_count     = 0;
    int     not_null_tiles = 0;
    std::vector<double> chunk_geotransform(6);

    LayoutfileReader lf(input_file);

    while(lf.next()){
      if(lf.newRow()){ //Add a row to the grid of chunks
        chunks.emplace_back();
        lfout.addRow();
      }

      if(lf.isNullTile()){
        chunks.back().emplace_back();
        continue;
      }

      not_null_tiles++;

      if(chunk_height==-1){
        //Retrieve information about this chunk. All chunks must have the same
        //dimensions, which we could check here, but opening and closing
        //thousands of files is expensive. Therefore, we rely on the user to
        //check this beforehand if they want to. We will, however, verify that
        //things are correct in Consumer() as we open the files for reading.
        if(getGDALDimensions(
            lf.getFullPath(),
            chunk_height,
            chunk_width,
            file_type,
            chunk_geotransform.data()
        )!=0){
          std::cerr<<"Error getting file information from '"<<lf.getFullPath()<<"'!"<<std::endl;
          CommAbort(-1); //TODO
        }
      }

      cell_count += chunk_width*chunk_height;

      std::string this_retention = retention;
      if(retention.find("%f")!=std::string::npos){
        this_retention.replace(this_retention.find("%f"), 2, lf.getBasename());          
      } else if(retention.find("%n")!=std::string::npos){
        this_retention.replace(this_retention.find("%n"), 2, lf.getGridLocName());
      } else if(retention[0]=='@') {
        this_retention = retention;
      } else { //Should never happen
        std::cerr<<"Outputname for mode-many must contain '%f' or '%n'!"<<std::endl;
        throw std::runtime_error("Outputname for mode-many must contain '%f' or '%n'!");
      }

      std::string this_output_name = output_name;
      if(output_name.find("%f")!=std::string::npos){
        this_output_name.replace(this_output_name.find("%f"), 2, lf.getBasename());          
      } else if(output_name.find("%n")!=std::string::npos){
        this_output_name.replace(this_output_name.find("%n"), 2, lf.getGridLocName());
      } else { //Should never happen
        std::cerr<<"Outputname for mode-many must contain '%f' or '%n'!"<<std::endl;
        throw std::runtime_error("Outputname for mode-many must contain '%f' or '%n'!");
      }

      chunks.back().emplace_back(
        lf.getFullPath(),
        this_output_name,
        this_retention,
        lf.getX(),
        lf.getY(),
        0,
        0,
        chunk_width,
        chunk_height,
        true
      );

      lfout.addEntry(this_output_name);

      //Flip tiles if the geotransform demands it
      if(chunk_geotransform[0]<0)
        chunks.back().back().flip ^= FLIP_HORZ;
      if(chunk_geotransform[5]<0)
        chunks.back().back().flip ^= FLIP_VERT;

      //Flip (or reverse the above flip!) if the user demands it
      if(flipH)
        chunks.back().back().flip ^= FLIP_HORZ;
      if(flipV)
        chunks.back().back().flip ^= FLIP_VERT;
    }

    std::cerr<<"!Loaded "<<chunks.size()<<" rows each of which had "<<chunks[0].size()<<" columns."<<std::endl;
    std::cerr<<"!Total cells to be processed: "<<cell_count<<std::endl;
    std::cerr<<"!Number of tiles which were not null: "<<not_null_tiles<<std::endl;

    //nullChunks imply that the chunks around them have edges, as though they
    //are on the edge of the raster.
    for(int y=0;y<(int)chunks.size();y++)
    for(int x=0;x<(int)chunks[0].size();x++){
      if(chunks[y][x].nullChunk)
        continue;
      if(y-1>0 && x>0 && chunks[y-1][x].nullChunk)
        chunks[y][x].edge |= GRID_TOP;
      if(y+1<(int)chunks.size() && x>0 && chunks[y+1][x].nullChunk)
        chunks[y][x].edge |= GRID_BOTTOM;
      if(y>0 && x-1>0 && chunks[y][x-1].nullChunk)
        chunks[y][x].edge |= GRID_LEFT;
      if(y>0 && x+1<(int)chunks[0].size() && chunks[y][x+1].nullChunk)
        chunks[y][x].edge |= GRID_RIGHT;
    }

  } else if(many_or_one=="one") {
    int32_t total_height;
    int32_t total_width;

    //Get the total dimensions of the input file
    if(getGDALDimensions(input_file, total_height, total_width, file_type, NULL)!=0){
      std::cerr<<"Error getting file information from '"<<input_file<<"'!"<<std::endl;
      CommAbort(-1); //TODO
    }

    //If the user has specified -1, that implies that they want the entire
    //dimension of the raster along the indicated axis to be processed within a
    //single job.
    if(bwidth==-1)
      bwidth  = total_width;
    if(bheight==-1)
      bheight = total_height;

    std::cerr<<"!Total width:  "<<total_width <<"\n";
    std::cerr<<"!Total height: "<<total_height<<"\n";
    std::cerr<<"!Block width:  "<<bwidth      <<"\n";
    std::cerr<<"!Block height: "<<bheight     <<std::endl;
    std::cerr<<"!Total cells to be processed: "<<(total_width*total_height)<<std::endl;

    //Create a grid of jobs
    //TODO: Avoid creating extremely narrow or small strips
    for(int32_t y=0,gridy=0;y<total_height; y+=bheight, gridy++){
      chunks.emplace_back(std::vector<ChunkInfo>());
      for(int32_t x=0,gridx=0;x<total_width;x+=bwidth,  gridx++){
        //if(total_height-y<100){
        //  std::cerr<<"At least one tile is <100 cells in height. Please change rectangle size to avoid this!"<<std::endl;
        //  std::cerr<<"I suggest you use bheight="<<SuggestTileSize(bheight, total_height, 100)<<std::endl;
        //  throw std::logic_error("Tile height too small!");
        //}
        //if(total_width -x<100){
        //  std::cerr<<"At least one tile is <100 cells in width. Please change rectangle size to avoid this!"<<std::endl;
        //  std::cerr<<"I suggest you use bwidth="<<SuggestTileSize(bwidth, total_width, 100)<<std::endl;
        //  throw std::logic_error("Tile width too small!");
        //}

        if(retention[0]!='@' && retention.find("%n")==std::string::npos){
          std::cerr<<"In <one> mode '%n' must be present in the retention path."<<std::endl;
          throw std::invalid_argument("'%n' not found in retention path!");
        }

        if(output_name.find("%n")==std::string::npos){
          std::cerr<<"In <one> mode '%n' must be present in the output path."<<std::endl;
          throw std::invalid_argument("'%n' not found in output path!");
        }

        //Used for '%n' formatting
        std::string coord_string = std::to_string(gridx)+"_"+std::to_string(gridy);

        std::string this_retention = retention;
        if(this_retention[0]!='@')
          this_retention.replace(this_retention.find("%n"), 2, coord_string);
        std::string this_output_name = output_name;
        if(this_output_name.find("%n")==std::string::npos){
          std::cerr<<"Outputname must include '%n' for <one> mode."<<std::endl;
          throw std::runtime_error("Outputname must include '%n' for <one> mode.");
        }
        this_output_name.replace(this_output_name.find("%n"),2,coord_string);

        chunks.back().emplace_back(
          input_file,
          this_output_name,
          this_retention,
          gridx,
          gridy,
          x,
          y,
          (total_width-x >=bwidth )?bwidth :total_width -x, //TODO: Check
          (total_height-y>=bheight)?bheight:total_height-y,
          false
        );
      }
    }

  } else {
    std::cout<<"Unrecognised option! Must be 'many' or 'one'!"<<std::endl;
    CommAbort(-1);
  }

  //If a job is on the edge of the raster, mark it as having this property so
  //that it can be handled with elegance later.
  for(auto &e: chunks.front())
    e.edge |= GRID_TOP;
  for(auto &e: chunks.back())
    e.edge |= GRID_BOTTOM;
  for(size_t y=0;y<chunks.size();y++){
    chunks[y].front().edge |= GRID_LEFT;
    chunks[y].back().edge  |= GRID_RIGHT;
  }

  CommBroadcast(&file_type,0);
  overall.stop();
  std::cerr<<"!Preparer time: "<<overall.accumulated()<<"s."<<std::endl;

  std::cerr<<"!Input data type: "<<GDALGetDataTypeName(file_type)<<std::endl;

  switch(file_type){
    case GDT_Unknown:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(file_type)<<std::endl;
      CommAbort(-1); //TODO
    case GDT_Byte:
      return Producer<uint8_t >(chunks);
    case GDT_UInt16:
      return Producer<uint16_t>(chunks);
    case GDT_Int16:
      return Producer<int16_t >(chunks);
    case GDT_UInt32:
      return Producer<uint32_t>(chunks);
    case GDT_Int32:
      return Producer<int32_t >(chunks);
    case GDT_Float32:
      return Producer<float   >(chunks);
    case GDT_Float64:
      return Producer<double  >(chunks);
    case GDT_CInt16:
    case GDT_CInt32:
    case GDT_CFloat32:
    case GDT_CFloat64:
      std::cerr<<"Complex types are not supported. Sorry!"<<std::endl;
      CommAbort(-1); //TODO
    default:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(file_type)<<std::endl;
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

    Timer master_time;
    master_time.start();

    std::cerr<<"!Running program version: "<<program_version<<std::endl;

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
      if(retention[0]!='@' && retention.find("%n")==std::string::npos)
        throw std::invalid_argument("Retention filename must indicate file number with '%n'.");
      if(retention==output_name)
        throw std::invalid_argument("Retention and output filenames must differ.");
    } catch (const std::invalid_argument &ia){
      std::string output_err;
      if(ia.what()==std::string("stoi"))
        output_err = "Invalid width or height.";
      else
        output_err = ia.what();

      std::cerr<<"parallel_pflood.exe [--flipV] [--flipH] [--bwidth #] [--bheight #] <many/one> <retention> <input> <output>"<<std::endl;
      std::cerr<<"\tUse '--help' to show help."<<std::endl;

      std::cerr<<"###Error: "<<output_err<<std::endl;

      int good_to_go=0;
      CommBroadcast(&good_to_go,0);
      CommFinalize();
      return -1;
    }

    int good_to_go = 1;
    std::cerr<<"!Running with "            <<CommSize()<<" processes."<<std::endl;
    std::cerr<<"!Many or one: "            <<many_or_one<<std::endl;
    std::cerr<<"!Input file: "             <<input_file<<std::endl;
    std::cerr<<"!Retention strategy: "     <<retention <<std::endl;
    std::cerr<<"!Block width: "            <<bwidth    <<std::endl;
    std::cerr<<"!Block height: "           <<bheight   <<std::endl;
    std::cerr<<"!Flip horizontal: "        <<flipH     <<std::endl;
    std::cerr<<"!Flip vertical: "          <<flipV     <<std::endl;
    std::cerr<<"!World Size: "             <<CommSize()<<std::endl;
    CommBroadcast(&good_to_go,0);
    Preparer(many_or_one, retention, input_file, output_name, bwidth, bheight, flipH, flipV);

    master_time.stop();
    std::cerr<<"!TimeInfo: Total wall-time was "<<master_time.accumulated()<<"s."<<std::endl;

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