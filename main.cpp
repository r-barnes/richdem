//Compile with
// mpic++ -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
// mpirun -n 3 ./a.out ~/projects/watershed/data/beauford03.flt
// TODO: MPI abort
// TODO: See optimization notes at "http://www.boost.org/doc/libs/1_56_0/doc/html/mpi/tutorial.html"
// For memory usage see: http://info.prelert.com/blog/stl-container-memory-usage
// abs("single_proc@1"-"merged@1")>0
// SRTM data: https://dds.cr.usgs.gov/srtm/version2_1/SRTM1/Region_03/
#include "gdal_priv.h"
#include <iostream>
#include <boost/mpi.hpp>
#include <string>
#include <queue>
#include <vector>
#include <fstream> //For reading layout files
#include <sstream> //Used for parsing the <layout_file>
#include <boost/filesystem.hpp>

//We use the cstdint library here to ensure that the program behaves as expected
//across platforms, especially with respect to the expected limits of operation
//for data set sizes and labels. For instance, in C++, a signed integer must be
//at least 16 bits, but not necessarily more. We force a minimum of 32 bits as
//this is, after all, for use with large datasets.
#include <cstdint>
#include "Array2D.hpp"
#include "common.hpp"
//#define DEBUG 1

//TODO: Is it possible to run this without mpirun if we specify a single node
//job?

#define NO_FLOW         37 //TODO: Expalin and ensure it fits in flowdir_t
#define FLOW_TERMINATES -3 //TODO: Expalin and ensure it fits in flowdir_t
#define FLOW_EXTERNAL   -4 //TODO: Expalin and ensure it fits in flowdir_t

//TODO: What are these for?
const int TAG_WHICH_JOB   = 0;
const int TAG_CHUNK_DATA  = 1;
const int TAG_DONE_FIRST  = 2;
const int TAG_SECOND_DATA = 3;
const int TAG_DONE_SECOND = 4;

std::vector<std::string> tag_trans = {{"TAG_WHICH_JOB","TAG_CHUNK_DATA","TAG_DONE_FIRST","TAG_SECOND_DATA","TAG_DONE_SECOND"}};

const int SYNC_MSG_KILL = 0;
const int JOB_CHUNK     = 1;
const int JOB_FIRST     = 2;
const int JOB_SECOND    = 3;

const uint8_t FLIP_VERT   = 1;
const uint8_t FLIP_HORZ   = 2;

const int ACCUM_NO_DATA = -1;

//TODO
template<typename T> boost::mpi::status WorldRecv(int source, int tag, T & value) {
  boost::mpi::communicator world;
  assert(world.rank()!=source);
  //std::cerr<<"I am "<<world.rank()<<" receiving from "<<source<<" with tag "<<tag_trans[tag]<<std::endl;
  boost::mpi::status status = world.recv(source,tag,value);
  //std::cerr<<"Received."<<std::endl;
  return status;
}

//TODO
template<typename T> void WorldSend(int dest, int tag, const T & value) {
  boost::mpi::communicator world;
  assert(world.rank()!=dest);
  //std::cerr<<"I am "<<world.rank()<<" sending to "<<dest<<" with tag "<<tag_trans[tag]<<std::endl;
  world.send(dest,tag,value);
  //std::cerr<<"Sent"<<std::endl;
}

typedef uint8_t dependency_t;
typedef int32_t link_t;
typedef int32_t accum_t;

class ChunkInfo{
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & edge;
    ar & flip;
    ar & x;
    ar & y;
    ar & width;
    ar & height;
    ar & gridx;
    ar & gridy;
    ar & id;
    ar & nullChunk;
    ar & filename;
    ar & outputname;
    ar & retention;
  }
 public:
  uint8_t     edge;
  uint8_t     flip;
  int32_t     x,y,width,height,gridx,gridy;
  int32_t     id;
  bool        nullChunk;
  std::string filename;
  std::string outputname;
  std::string retention;
  ChunkInfo(){
    nullChunk = true;
  }
  ChunkInfo(int32_t id, std::string filename, std::string outputname, std::string retention, int32_t gridx, int32_t gridy, int32_t x, int32_t y, int32_t width, int32_t height){
    this->nullChunk    = false;
    this->edge         = 0;
    this->x            = x;
    this->y            = y;
    this->width        = width;
    this->height       = height;   
    this->gridx        = gridx;
    this->gridy        = gridy;
    this->id           = id;
    this->filename     = filename;
    this->outputname   = outputname;
    this->retention    = retention;
    this->flip         = 0;
  }
};

typedef std::vector< std::vector< ChunkInfo > > ChunkGrid;


class TimeInfo {
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & calc;
    ar & overall;
    ar & io;
  }
 public:
  double calc, overall, io;
  TimeInfo() {
    calc=overall=io=0;
  }
  TimeInfo(double calc, double overall, double io) :
      calc(calc), overall(overall), io(io) {}
  TimeInfo& operator+=(const TimeInfo& o){
    calc    += o.calc;
    overall += o.overall;
    io      += o.io;
    return *this;
  }
};

template<class flowdir_t>
class Job1 {
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & links;
    ar & accum;
    ar & flowdirs;
    ar & dependencies;
    ar & time_info;
  }
 public:
  std::vector<link_t   >    links;
  std::vector<accum_t  >    accum;
  std::vector<flowdir_t>    flowdirs;
  std::vector<dependency_t> dependencies;
  TimeInfo time_info;
  Job1(){}
};

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
  const Array2D<flowdir_t> &flowdirs, //Flow directions matrix
  std::vector<link_t>      &links
){
  int x = x0;
  int y = y0;

  int path_len = 0;

  //Follow the flow path until it terminates
  while(path_len++<10000000){      //Follow the flow path until we reach its end
    const int n = flowdirs(x,y);   //Neighbour the current cell flows towards

    //Show the final part of the loop path
    if(path_len>10000000-200)
      std::cerr<<"Path: "<<x<<","<<y<<" with flowdir "<<n<<std::endl;

    //If the neighbour this cell flows into is a no_data cell or this cell does
    //not flow into any neighbour, then mark the initial cell from which we
    //began this flow path as terminating somewhere unimportant: its flow cannot
    //pass to neighbouring segments/nodes for further processing.
    if(flowdirs.isNoData(x,y) || n==NO_FLOW){
      links.at(xyToSerial(x0,y0,flowdirs.viewWidth(),flowdirs.viewHeight())) = FLOW_TERMINATES;
      return;
    }

    //Flow direction was valid. Where does it lead?
    const int nx = x+dx[n]; //Get neighbour's x-coordinate.
    const int ny = y+dy[n]; //Get neighbour's y-coordinate.

    //The neighbour cell is off one of the sides of the tile. Therefore, its
    //flow may be passed on to a neighbouring tile. Thus, we need to link this
    //flow path's initial cell to this terminal cell.
    if(!flowdirs.in_grid(nx,ny)) {
      if(x==x0 && y==y0) {
        //If (x,y)==(x0,y0), then this is the first cell and it points off of
        //the grid. Mark it as being FLOW_EXTERNAL.
        links.at(xyToSerial(x0,y0,flowdirs.viewWidth(),flowdirs.viewHeight())) = FLOW_EXTERNAL;
      } else {
        //The flow path entered the grid one place and left another. We mark
        //only the start of the flow path since the end of the flow path will be
        //handled by another call to this function
        links.at(xyToSerial(x0,y0,flowdirs.viewWidth(),flowdirs.viewHeight())) = xyToSerial(x,y,flowdirs.viewWidth(),flowdirs.viewHeight());
      }
      return;
    }

    //The flow path has not yet terminated. Continue following it.
    x = nx;
    y = ny;
  }

  //The loop breaks with a return, so this is only reached for very long paths
  throw std::logic_error("A flow path was more than 10M cells long - there's probably a loop!");
}

//As in the function above, we will start at an initial cell and follow its flow
//direction to its neighbour. Then we will follow the neighbour's flow direction
//to the neighbour's neighbour, and so on, until we can go no farther (as
//indicated by running of one of the edges of the segment or hitting a no_data
//cell).

//However, at this point we know how much flow accumulation to add to each cell.
//Therefore, we add this additional accumulation to each cell as we pass it.
template<class flowdir_t>
void FollowPathAdd(
  int x,                              //Initial x-coordinate
  int y,                              //Initial y-coordinate
  const Array2D<flowdir_t> &flowdirs, //Flow directions matrix
  Array2D<accum_t>         &accum,    //Output: Accumulation matrix
  const accum_t additional_accum      //Add this to every cell in the flow path
){
  int count = 0;

  //Follow the flow path until it terminates
  while(count++<10000000){
    //Break when we reach a no_data cell
    if(flowdirs.isNoData(x,y))
      return;
    
    //Add additional flow accumulation to this cell
    accum(x,y) += additional_accum;

    int n = flowdirs(x,y); //Get neighbour
    if(n==NO_FLOW)         //This cell doesn't flow to a neighbour
      return;              

    //Move to neighbour
    x += dx[n];
    y += dy[n];

    //The neighbour is off the edge of the grid. Flow path has terminated
    if(!flowdirs.in_grid(x,y))
      return;
  }

  //The loop breaks with a return, so this is only reached for very long paths
  throw std::logic_error("A flow path was more than 10M cells long - there's probably a loop!");
}

template<class flowdir_t>
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

  std::cerr<<"Calculating dependencies..."<<std::endl;
  Array2D<dependency_t> dependencies;
  dependencies.resize(flowdirs,0);
  for(int y=0;y<flowdirs.viewHeight();y++)
  for(int x=0;x<flowdirs.viewWidth();x++){
    int n = flowdirs(x,y);       //The neighbour this cell flows into

    
    if(flowdirs.isNoData(x,y)){  //This cell is a no_data cell
      accum(x,y) = ACCUM_NO_DATA;
      continue;                
    }         
    if(n==NO_FLOW)               //This cell does not flow into a neighbour
      continue;

    int nx = x+dx[n];            //x-coordinate of the neighbour
    int ny = y+dy[n];            //y-coordinate of the neighbour
      
    //Neighbour is not on the grid
    if(!flowdirs.in_grid(nx,ny))
      continue;

    //Neighbour is valid and is part of the grid. The neighbour depends on this
    //cell, so increment its dependency count.
    dependencies(nx,ny)++;
  }

  //Now that we know how many dependencies each cell has, we can determine which
  //cells are the peaks: the sources of flow. We make a note of where the peaks
  //are for later use.
  std::queue<GridCell> sources;
  for(int y=0;y<dependencies.viewHeight();y++)
  for(int x=0;x<dependencies.viewWidth();x++)
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
    if(!flowdirs.in_grid(nx,ny))
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


template<class T>
void GridPerimToArray(const Array2D<T> &grid, std::vector<T> &vec){
  assert(vec.size()==0); //Ensure receiving array is empty

  std::vector<T> vec2copy;

  vec2copy = grid.getRowData(0);                         //Top
  vec.insert(vec.end(),vec2copy.begin(),vec2copy.end());

  vec2copy = grid.getColData(grid.viewWidth()-1);        //Right
  vec.insert(vec.end(),vec2copy.begin()+1,vec2copy.end());
  
  vec2copy = grid.getRowData(grid.viewHeight()-1);       //Bottom
  vec.insert(vec.end(),vec2copy.begin(),vec2copy.end()-1);
  
  vec2copy = grid.getColData(0);                         //Left
  vec.insert(vec.end(),vec2copy.begin()+1,vec2copy.end()-1);
}

template<class flowdir_t>
void DownstreamCell(
  const std::vector< std::vector< Job1<flowdir_t> > > &jobs,
  const ChunkGrid &chunks,
  const int gx,
  const int gy,
  const int s,
  int &ns,
  int &gnx, //Input is the tile of this cell
  int &gny  //Input is the tile of this cell
){
  const auto &j = jobs.at(gy).at(gx);

  assert(j.links.size()==j.flowdirs.size());

  ns=-1;

  gnx = gx;
  gny = gy;

  //Flow ends somewhere internal to the tile or this particular cell has no
  //flow direction. The TERMINATES should also take care of NoData flowdirs
  if(j.links.at(s)==FLOW_TERMINATES || j.flowdirs.at(s)==NO_FLOW){
    return;
  } else if(j.links.at(s)==FLOW_EXTERNAL){ //Flow goes into a valid neighbouring tile
    const auto &c = chunks.at(gy).at(gx);

    int x,y;
    serialToXY(s,x,y,c.width,c.height);

    const int nfd = j.flowdirs.at(s);
    int nx        = x+dx[nfd];
    int ny        = y+dy[nfd];

    //Since this is a FLOW_EXTERNAL, this cell must point off of the grid
    assert(nx<0 || ny<0 || nx==c.width || ny==c.height);

    if(nx==c.width){
      gnx += 1;
      nx   = 0;
    } else if(nx==-1){
      gnx -= 1;
      nx   = c.width-1;
    } 
    
    if(ny==c.height){
      gny += 1;
      ny   = 0;
    } else if(ny==-1) {
      gny -= 1;
      ny   = c.height-1;
    }

    //NOTE: gridwidth=jobs.front().size() and gridheight=jobs.size()
    if(gnx<0 || gny<0 || gnx==(int)jobs.front().size() || gny==(int)jobs.size()){
      gnx=gny=ns=-1;
      return;
    }

    const auto &nc = chunks.at(gny).at(gnx);

    ns = xyToSerial(nx,ny,nc.width,nc.height);
  } else { //Flow goes to somewhere else on the perimeter of the same tile
    ns = j.links.at(s);
  }
}


template<class flowdir_t>
void Consumer(){
  boost::mpi::environment env;
  boost::mpi::communicator world;

  int the_job;   //Which job should consumer perform?

  ChunkInfo chunk;

  Array2D<flowdir_t> flowdirs;
  Array2D<accum_t  > accum;

  //Have the consumer process messages as long as they are coming using a
  //blocking receive to wait.
  while(true){
    WorldRecv(0, TAG_WHICH_JOB, the_job); //Blocking wait

    //This message indicates that everything is done and the Consumer should shut
    //down.
    if(the_job==SYNC_MSG_KILL){
      return;

    } else if (the_job==JOB_CHUNK){
      WorldRecv(0, TAG_CHUNK_DATA, chunk);

    //This message indicates that the consumer should prepare to perform the
    //first part of the distributed Priority-Flood algorithm on an incoming job
    } else if (the_job==JOB_FIRST){
      Timer timer_calc,timer_io,timer_overall;
      timer_overall.start();

      Job1<flowdir_t> job1;

      //Read in the data associated with the job
      timer_io.start();
      flowdirs = Array2D<flowdir_t>(chunk.filename, false, chunk.x, chunk.y, chunk.width, chunk.height);
      if(chunk.flip & FLIP_VERT)
        flowdirs.flipVert();
      if(chunk.flip & FLIP_HORZ)
        flowdirs.flipHorz();
      timer_io.stop();

      //Let's double-check that the flowdirs are valid
      for(int y=0;y<flowdirs.viewHeight();y++)
      for(int x=0;x<flowdirs.viewWidth();x++)
        if(!flowdirs.isNoData(x,y) && !(1<=flowdirs(x,y) && flowdirs(x,y)<=8) && !(flowdirs(x,y)==NO_FLOW))
          throw std::domain_error("Invalid flow direction found: "+std::to_string(flowdirs(x,y)));


      timer_calc.start();
      FlowAccumulation(flowdirs,accum);

      //-2 removes duplicate cells on vertical edges which would otherwise
      //overlap horizontal edges
      job1.links.resize(2*flowdirs.viewWidth()+2*(flowdirs.viewHeight()-2), FLOW_TERMINATES);

      //TODO: Although the following may consider a cell more than once, the
      //repeated effort merely produces the same results in the same places
      //twice, so it's okay.

      //If we are the top segment, nothing can flow into us, so we do not need
      //to know where flow paths originating at the top go to. On the otherhand,
      //if we are not the top segment, then consider each cell of the top row
      //and find out where its flow goes to.
      if(!(chunk.edge & GRID_TOP))
        for(int x=0;x<flowdirs.viewWidth();x++)
          FollowPath(x,0,flowdirs,job1.links);

      //If we are the bottom segment, nothing can flow into us, so we do not
      //need to know where flow paths originating at the bottom go to. On the
      //otherhand, if we are not the bottom segment, then consider each cell of
      //the bottom row and find out where its flow goes to.
      if(!(chunk.edge & GRID_BOTTOM))
        for(int x=0;x<flowdirs.viewWidth();x++)
          FollowPath(x,flowdirs.viewHeight()-1,flowdirs,job1.links);

      if(!(chunk.edge & GRID_LEFT))
        for(int y=0;y<flowdirs.viewHeight();y++)
          FollowPath(0,y,flowdirs,job1.links);

      if(!(chunk.edge & GRID_RIGHT))
        for(int y=0;y<flowdirs.viewHeight();y++)
          FollowPath(flowdirs.viewWidth()-1,y,flowdirs,job1.links);

      //Construct output arrays
      GridPerimToArray(flowdirs, job1.flowdirs);
      GridPerimToArray(accum,    job1.accum   );

      timer_calc.stop();

      if(chunk.retention=="@offloadall"){
        //Nothing to do: it will all get overwritten
      } else if(chunk.retention=="@retainall"){
        //Nothing to do: it will all be kept because this process won't get
        //called again
      } else {
        timer_io.start();
        flowdirs.saveNative(chunk.retention+"dem.dat"   );
        accum.saveNative   (chunk.retention+"labels.dat");
        timer_io.stop();
      }

      timer_overall.stop();
      std::cerr<<"Node "<<world.rank()<<" finished with Calc="<<timer_calc.accumulated()<<"s. timer_Overall="<<timer_overall.accumulated()<<"s. timer_IO="<<timer_io.accumulated()<<"s."<<std::endl;

      job1.time_info = TimeInfo(timer_calc.accumulated(), timer_overall.accumulated(), timer_io.accumulated());

      WorldSend(0, TAG_DONE_FIRST, job1);
    } else if (the_job==JOB_SECOND){
      Timer timer_calc,timer_io,timer_overall;
      timer_overall.start();
      std::vector<accum_t> accum_offset; //TODO

      WorldRecv(0, TAG_SECOND_DATA, accum_offset);

      //These use the same logic as the analogous lines above
      if(chunk.retention=="@offloadall"){
        //Read in the data associated with the job
        timer_io.start();
        flowdirs = Array2D<flowdir_t>(chunk.filename, false, chunk.x, chunk.y, chunk.width, chunk.height);
        if(chunk.flip & FLIP_VERT)
          flowdirs.flipVert();
        if(chunk.flip & FLIP_HORZ)
          flowdirs.flipHorz();
        timer_io.stop();

        timer_calc.start();
        FlowAccumulation(flowdirs,accum);
        timer_calc.stop();
      } else if(chunk.retention=="@retainall"){
        //Nothing to do: we have it all in memory
      } else {
        timer_io.start();
        flowdirs = Array2D<flowdir_t>(chunk.retention+"dem.dat"    ,true); //TODO: There should be an exception if this fails
        accum    = Array2D<accum_t  >(chunk.retention+"labels.dat",true);
        timer_io.stop();
      }

      for(int s=0;s<(int)accum_offset.size();s++){
        if(accum_offset[s]==0)
          continue;
        int x,y;
        serialToXY(s, x, y, accum.viewWidth(), accum.viewHeight());
        FollowPathAdd(x,y,flowdirs,accum,accum_offset[s]);
      }

      //At this point we're done with the calculation! Boo-yeah!

      timer_io.start();
      if(chunk.flip & FLIP_HORZ)
        accum.flipHorz();
      if(chunk.flip & FLIP_VERT)
        accum.flipVert();
      timer_io.stop();

      timer_io.start();
      Timer timer_save_gdal;
      timer_save_gdal.start(); //TODO: Remove
      accum.saveGDAL(chunk.outputname, chunk.filename, chunk.x, chunk.y);
      timer_save_gdal.stop();
      timer_io.stop();
      timer_overall.stop();
      std::cerr<<"GDAL save took "<<timer_save_gdal.accumulated()<<"s."<<std::endl;

      std::cerr<<"Node "<<world.rank()<<" finished ("<<chunk.gridx<<","<<chunk.gridy<<") with Calc="<<timer_calc.accumulated()<<"s. timer_Overall="<<timer_overall.accumulated()<<"s. timer_IO="<<timer_io.accumulated()<<"s."<<std::endl;

      TimeInfo temp(timer_calc.accumulated(), timer_overall.accumulated(), timer_io.accumulated());
      WorldSend(0, TAG_DONE_SECOND, temp);
    }
  }
}








//Producer takes a collection of Jobs and delegates them to Consumers. Once all
//of the jobs have received their initial processing, it uses that information
//to compute the global properties necessary to the solution. Each Job, suitably
//modified, is then redelegated to a Consumer which ultimately finishes the
//processing.
template<class flowdir_t>
void Producer(ChunkGrid &chunks){
  boost::mpi::environment env;
  boost::mpi::communicator world;
  Timer timer_overall,timer_calc;
  timer_overall.start();
  int active_nodes = 0;

  TimeInfo time_first_total;
  TimeInfo time_second_total;
  int      time_first_count  = 0;
  int      time_second_count = 0;

  std::map<int,ChunkInfo> rank_to_chunk;

  std::vector< std::vector< Job1<flowdir_t> > > jobs1(chunks.size(), std::vector< Job1<flowdir_t> >(chunks[0].size()));

  //Loop through all of the jobs, delegating them to Consumers
  active_nodes=0;
  for(size_t y=0;y<chunks.size();y++)
  for(size_t x=0;x<chunks[0].size();x++){
    std::cerr<<"Sending job "<<(y*chunks[0].size()+x+1)<<"/"<<(chunks.size()*chunks[0].size())<<" ("<<(x+1)<<"/"<<chunks[0].size()<<","<<(y+1)<<"/"<<chunks.size()<<")"<<std::endl;
    if(chunks.at(y).at(x).nullChunk){
      std::cerr<<"\tNull chunk: skipping."<<std::endl;
      continue;
    }

    //If fewer jobs have been delegated than there are Consumers available,
    //delegate the job to a new Consumer.
    if(active_nodes<world.size()-1){
      active_nodes++;
      // std::cerr<<"Sending init to "<<(active_nodes+1)<<std::endl;
      WorldSend(active_nodes,TAG_WHICH_JOB,JOB_CHUNK);
      WorldSend(active_nodes,TAG_CHUNK_DATA,chunks.at(y).at(x));

      rank_to_chunk[active_nodes] = chunks.at(y).at(x);
      WorldSend(active_nodes,TAG_WHICH_JOB,JOB_FIRST);

    //Once all of the consumers are active, wait for them to return results. As
    //each Consumer returns a result, pass it the next unfinished Job until
    //there are no jobs left.
    } else {
      Job1<flowdir_t> finished_job;

      //TODO: Note here about how the code below could be incorporated

      //Execute a blocking receive until some consumer finishes its work.
      //Receive that work.
      boost::mpi::status status = WorldRecv(boost::mpi::any_source,TAG_DONE_FIRST,finished_job);

      //NOTE: This could be used to implement more robust handling of lost nodes
      ChunkInfo received_chunk = rank_to_chunk[status.source()];
      jobs1.at(received_chunk.gridy).at(received_chunk.gridx) = finished_job;

      //Delegate new work to that consumer
      WorldSend(status.source(),TAG_WHICH_JOB,JOB_CHUNK);
      WorldSend(status.source(),TAG_CHUNK_DATA,chunks.at(y).at(x));

      rank_to_chunk[status.source()] = chunks.at(y).at(x);
      WorldSend(status.source(),TAG_WHICH_JOB,JOB_FIRST);
    }
  }

  while(active_nodes>0){
    Job1<flowdir_t> finished_job;

    //Execute a blocking receive until some consumer finishes its work.
    //Receive that work
    boost::mpi::status status = WorldRecv(boost::mpi::any_source,TAG_DONE_FIRST,finished_job);
    ChunkInfo received_chunk = rank_to_chunk[status.source()];
    jobs1.at(received_chunk.gridy).at(received_chunk.gridx) = finished_job;

    //Decrement the number of consumers we are waiting on. When this hits 0 all
    //of the jobs have been completed and we can move on
    active_nodes--;
  }


  //TODO: Note here about how the code above code be incorporated

  const int gridwidth  = jobs1.front().size();
  const int gridheight = jobs1.size();

  //Set initial values for all dependencies to zero
  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++){
    auto &this_job = jobs1.at(y).at(x);
    this_job.dependencies.resize(this_job.links.size(),0);
  }

  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++){
    auto &this_chunk = chunks.at(y).at(x);
    if(this_chunk.nullChunk)
      continue;

    auto &this_job = jobs1.at(y).at(x);

    time_first_total += this_job.time_info;
    time_first_count++;

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

  class atype {
   public:
    int gx,gy,s;
    atype(int gx, int gy, int s) : gx(gx), gy(gy), s(s) {};
  };

  std::queue<atype> q;

  int accum_to_dist=0; //TODO
  //Search for cells without dependencies
  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++){
    if(chunks.at(y).at(x).nullChunk)
      continue;

    auto &this_job = jobs1.at(y).at(x);

    for(size_t s=0;s<this_job.dependencies.size();s++){
      if(this_job.dependencies.at(s)==0) // && jobs1.at(y).at(x).links.at(s)==FLOW_EXTERNAL) //TODO
        q.emplace(x,y,s);

      //Accumulated flow at an input will be transfered to an output resulting
      //in double-counting
      if(this_job.links.at(s)!=FLOW_EXTERNAL) //TODO: NO_FLOW
        this_job.accum.at(s) = 0;
      else
        accum_to_dist += this_job.accum.at(s); //TODO
    }
  }

  std::cerr<<accum_to_dist<<" accumulation to distribute."<<std::endl;

  std::cerr<<q.size()<<" peaks founds in aggregated problem."<<std::endl;

  int processed_peaks = 0;
  while(!q.empty()){
    atype            c  = q.front();
    Job1<flowdir_t> &j  = jobs1.at(c.gy).at(c.gx);
    ChunkInfo       &ci = chunks.at(c.gy).at(c.gx);
    q.pop();
    processed_peaks++;

    assert(!ci.nullChunk);

    int ns  = -1; //Initial value designed to cause devastation if misused
    int gnx = -1; //Initial value designed to cause devastation if misused
    int gny = -1; //Initial value designed to cause devastation if misused
    DownstreamCell(jobs1, chunks, c.gx, c.gy, c.s, ns, gnx, gny);
    if(ns==-1 || chunks.at(gny).at(gnx).nullChunk)
      continue;

    jobs1.at(gny).at(gnx).accum.at(ns) += j.accum.at(c.s);

    if( (--jobs1.at(gny).at(gnx).dependencies.at(ns))==0 )
      q.emplace(gnx,gny,ns);
  }

  std::cerr<<processed_peaks<<" peaks processed."<<std::endl;

  //TODO: Remove this block as this should always be zero once I have things
  //right.
  int unprocessed_dependencies=0;
  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++)
  for(size_t s=0;s<jobs1[y][x].dependencies.size();s++)
    unprocessed_dependencies+=jobs1[y][x].dependencies[s];
  std::cerr<<unprocessed_dependencies<<" unprocessed dependencies."<<std::endl;

  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++){
    if(chunks.at(y).at(x).nullChunk)
      continue;

    for(size_t s=0;s<jobs1.at(y).at(x).accum.size();s++) //TODO: NO_FLOW
      if(jobs1.at(y).at(x).links.at(s)==FLOW_EXTERNAL)
        jobs1.at(y).at(x).accum.at(s) = 0;

    //TODO: Could start clearing memory here
  }

  std::cerr<<"Sending out final jobs..."<<std::endl;
  //Loop through all of the jobs, delegating them to Consumers
  active_nodes = 0;
  for(size_t y=0;y<chunks.size();y++)
  for(size_t x=0;x<chunks[0].size();x++){
    std::cerr<<"Sending job "<<(y*chunks[0].size()+x+1)<<"/"<<(chunks.size()*chunks[0].size())<<" ("<<(x+1)<<"/"<<chunks[0].size()<<","<<(y+1)<<"/"<<chunks.size()<<")"<<std::endl;
    if(chunks.at(y).at(x).nullChunk){
      std::cerr<<"\tNull chunk: skipping."<<std::endl;
      continue;
    }

    int destination_node = -1;
    if(active_nodes<world.size()-1){
      //If fewer jobs have been delegated than there are Consumers available,
      //delegate the job to a new Consumer
      destination_node = ++active_nodes;
    } else {
      //Once all of the consumers are active, wait for them to return results.
      //As each Consumer returns a result, pass it the next unfinished Job until
      //there are no jobs left.
      TimeInfo temp;
      boost::mpi::status status = WorldRecv(boost::mpi::any_source,TAG_DONE_SECOND,temp);

      time_second_total += temp;
      time_second_count++;

      destination_node = status.source();
    }

    WorldSend(destination_node,TAG_WHICH_JOB,JOB_CHUNK);
    WorldSend(destination_node,TAG_CHUNK_DATA,chunks.at(y).at(x));

    rank_to_chunk[destination_node] = chunks.at(y).at(x);
    WorldSend(destination_node,TAG_WHICH_JOB,JOB_SECOND);
    WorldSend(destination_node,TAG_SECOND_DATA,jobs1.at(y).at(x).accum);
  }

  while(active_nodes>0){
    TimeInfo temp;
    //Execute a blocking receive until some consumer finishes its work.
    //Receive that work
    WorldRecv(boost::mpi::any_source,TAG_DONE_SECOND,temp);

    time_second_total += temp;
    time_second_count++;

    //Decrement the number of consumers we are waiting on. When this hits 0 all
    //of the jobs have been completed and we can move on
    active_nodes--;
  }

  for(int i=1;i<world.size();i++)
    WorldSend(i,TAG_WHICH_JOB,SYNC_MSG_KILL);

  timer_overall.stop();

  std::cout<<"TimeInfo: First stage total overall time="<<time_first_total.overall<<std::endl;
  std::cout<<"TimeInfo: First stage total io time="     <<time_first_total.io     <<std::endl;
  std::cout<<"TimeInfo: First stage total calc time="   <<time_first_total.calc   <<std::endl;
  std::cout<<"TimeInfo: First stage block count="       <<time_first_count        <<std::endl;

  std::cout<<"TimeInfo: Second stage total overall time="<<time_second_total.overall<<std::endl;
  std::cout<<"TimeInfo: Second stage total IO time="     <<time_second_total.io     <<std::endl;
  std::cout<<"TimeInfo: Second stage total calc time="   <<time_second_total.calc   <<std::endl;
  std::cout<<"TimeInfo: Second stage block count="       <<time_second_count        <<std::endl;

  std::cout<<"TimeInfo: Producer overall="<<timer_overall.accumulated()<<std::endl;
  std::cout<<"TimeInfo: Producer calc="   <<timer_calc.accumulated()   <<std::endl;
}



std::string trimStr(std::string const& str){
  if(str.empty())
      return str;

  std::size_t firstScan = str.find_first_not_of(' ');
  std::size_t first     = firstScan == std::string::npos ? str.length() : firstScan;
  std::size_t last      = str.find_last_not_of(' ');
  return str.substr(first, last-first+1);
}


//Preparer divides up the input raster file into chunks which can be processed
//independently by the Consumers. Since the chunking may be done on-the-fly or
//rely on preparation the user has done, the Preparer routine knows how to deal
//with both. Once assemebled, the collection of jobs is passed off to Producer,
//which is agnostic as to the original form of the jobs and handles
//communication and solution assembly.
void Preparer(
  std::string many_or_one,
  std::string retention_base,
  std::string input_file,
  std::string output_prefix,
  int bwidth,
  int bheight,
  int flipH,
  int flipV
){
  boost::mpi::environment  env;
  boost::mpi::communicator world;
  int chunkid             = 0;
  Timer overall;
  overall.start();

  ChunkGrid chunks;
  std::string  filename;
  GDALDataType file_type; //All chunks must have a common file_type

  if(many_or_one=="many"){
    std::cerr<<"Multi file mode"<<std::endl;

    int32_t gridx        =  0; //Current x coordinate in the chunk grid
    int32_t gridy        = -1; //Current y coordinate in the chunk grid
    int32_t row_width    = -1; //Width of 1st row. All rows must equal this
    int32_t chunk_width  = -1; //Width of 1st chunk. All chunks must equal this
    int32_t chunk_height = -1; //Height of 1st chunk, all chunks must equal this
    
    boost::filesystem::path layout_path_and_name = input_file;
    auto path = layout_path_and_name.parent_path();

    //Read each line of the layout file
    std::ifstream fin_layout(input_file);
    while(fin_layout.good()){
      gridy++;
      std::string line;
      std::getline(fin_layout,line); //Read layout file line

      //If this line is a blank row, we're done
      if(line.find_first_not_of("\t\n ")==std::string::npos)
        break;

      //Add a row to the grid of chunks
      chunks.emplace_back(std::vector<ChunkInfo>());

      std::stringstream cells(line);
      std::string       filename;
      gridx = -1;

      //Split the row of file names at the commas
      while(std::getline(cells,filename,',')){
        gridx++;
        filename = trimStr(filename);

        //If the comma delimits only whitespace, then this is a null chunk
        if(filename==""){
          chunks.back().emplace_back();
          continue;
        }

        //Okay, the file exists. Let's check it out.
        auto path_and_filename = path / filename;
        auto path_and_filestem = path / path_and_filename.stem();
        auto outputname        = output_prefix+path_and_filename.stem().string()+"-fill.tif";

        std::string retention = retention_base;
        if(retention[0]!='@')
          retention = retention_base+path_and_filestem.stem().string()+"-int-";

        //Place to store information about this chunk
        int          this_chunk_width;
        int          this_chunk_height;
        GDALDataType this_chunk_type;
        double       this_geotransform[6];

        //Retrieve information about this chunk
        if(getGDALDimensions(
            path_and_filename.string(),
            this_chunk_height,
            this_chunk_width,
            this_chunk_type,
            this_geotransform
        )!=0){
          std::cerr<<"Error getting file information from '"<<path_and_filename.string()<<"'!"<<std::endl;
          env.abort(-1); //TODO
        }

        //If this is the first chunk we've seen, it determines what every other
        //chunk should be. Note it.
        if(chunk_height==-1){
          chunk_height = this_chunk_height;
          chunk_width  = this_chunk_width;
          file_type    = this_chunk_type;
        }

        //Check that this chunk conforms to what we've seen previously so that
        //all chunks are the same.
        if( this_chunk_height!=chunk_height || 
            this_chunk_width !=chunk_width  ||
            this_chunk_type  !=file_type   ){
          std::cerr<<"All of the files specified by <layout_file> must be the same dimensions and type!"<<std::endl;
          env.abort(-1); //TODO: Set error code
        }

        //Add the chunk to the grid
        chunks.back().emplace_back(chunkid++, path_and_filename.string(), outputname, retention, gridx, gridy, 0, 0, chunk_width, chunk_height);

        //Flip tiles if the geotransform demands it
        if(this_geotransform[0]<0)
          chunks.back().back().flip ^= FLIP_HORZ;
        if(this_geotransform[5]<0)
          chunks.back().back().flip ^= FLIP_VERT;

        //Flip (or reverse the above flip!) if the user demands it
        if(flipH)
          chunks.back().back().flip ^= FLIP_HORZ;
        if(flipV)
          chunks.back().back().flip ^= FLIP_VERT;
      }

      if(row_width==-1){ //This is the first row
        row_width = gridx;
      } else if(row_width!=gridx){
        std::cerr<<"All rows of <layout_file> must specify the same number of files! First row="<<(row_width+1)<<", this row="<<(gridx+1)<<"."<<std::endl;
        env.abort(-1); //TODO: Set error code
      }
    }

    std::cerr<<"Loaded "<<chunks.size()<<" rows of "<<chunks[0].size()<<" columns."<<std::endl;

    //nullChunks imply that the chunks around them have edges, as though they
    //are on the edge of the raster.
    for(int y=0;y<(int)chunks.size();y++)
    for(int x=0;x<(int)chunks[0].size();x++){
      if(chunks.at(y).at(x).nullChunk)
        continue;
      if(y-1>0 && x>0 && chunks[y-1][x].nullChunk)
        chunks.at(y).at(x).edge |= GRID_TOP;
      if(y+1<(int)chunks.size() && x>0 && chunks[y+1][x].nullChunk)
        chunks.at(y).at(x).edge |= GRID_BOTTOM;
      if(y>0 && x-1>0 && chunks[y][x-1].nullChunk)
        chunks.at(y).at(x).edge |= GRID_LEFT;
      if(y>0 && x+1<(int)chunks[0].size() && chunks[y][x+1].nullChunk)
        chunks.at(y).at(x).edge |= GRID_RIGHT;
    }

  } else if(many_or_one=="one") {
    std::cerr<<"Single file mode"<<std::endl;
    int32_t total_height;
    int32_t total_width;

    filename = input_file;

    auto filepath   = boost::filesystem::path(filename);
    filepath        = filepath.parent_path() / filepath.stem();

    //Get the total dimensions of the input file
    if(getGDALDimensions(filename, total_height, total_width, file_type, NULL)!=0){
      std::cerr<<"Error getting file information from '"<<filename<<"'!"<<std::endl;
      env.abort(-1); //TODO
    }

    //If the user has specified -1, that implies that they want the entire
    //dimension of the raster along the indicated axis to be processed within a
    //single job.
    if(bwidth==-1)
      bwidth  = total_width;
    if(bheight==-1)
      bheight = total_height;

    std::cerr<<"Total width:  "<<total_width <<"\n";
    std::cerr<<"Total height: "<<total_height<<"\n";
    std::cerr<<"Block width:  "<<bwidth      <<"\n";
    std::cerr<<"Block height: "<<bheight     <<std::endl;

    //Create a grid of jobs
    //TODO: Avoid creating extremely narrow or small strips
    for(int32_t y=0,gridy=0;y<total_height; y+=bheight, gridy++){
      chunks.emplace_back(std::vector<ChunkInfo>());
      for(int32_t x=0,gridx=0;x<total_width;x+=bwidth,  gridx++){
        if(total_height-y<100 || total_width-x<100){
          throw std::logic_error("At least one tile is <100 cells in at least one dimensions. Please change rectangle size to avoid this!");
        }
        auto outputname = output_prefix+filepath.stem().string()+"-"+std::to_string(chunkid)+"-fill.tif";
        std::string retention = retention_base;
        if(retention[0]!='@')
          retention = retention_base+filepath.stem().string()+"-int-"+std::to_string(chunkid)+"-";
        chunks.back().emplace_back(chunkid++, filename, outputname, retention, gridx, gridy, x, y, bwidth, bheight);
        //Adjust the label_offset by the total number of perimeter cells of this
        //chunk plus one (to avoid another chunk's overlapping the last label of
        //this chunk). Obviously, the right and bottom edges of the global grid
        //may not be a perfect multiple of bwidth and bheight; therefore, labels
        //could be dolled out more conservatively. But this would require a
        //correctness proof and introduces the potential for truly serious bugs
        //into the code. Since we are unlikely to run out of labels (see
        //associated manuscript), it is better to waste a few and make this
        //section obviously correct. Although we do not expect to run out of
        //labels, it is possible to rigorously check for this condition here,
        //before we have used much time or I/O.
      }
    }

    if(retention_base=="@retainall" && (int)(chunks.size()*chunks[0].size())>=(int)(world.size()-1)){
      std::cerr<<"This job requires "<<(chunks.size()*chunks[0].size()+1)<<" processes. Only "<<world.size()<<" are available."<<std::endl;
      env.abort(-1);
    }

  } else {
    std::cout<<"Unrecognised option! Must be 'many' or 'one'!"<<std::endl;
    env.abort(-1);
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

  boost::mpi::broadcast(world,file_type,0);
  overall.stop();
  std::cerr<<"Preparer took "<<overall.accumulated()<<"s."<<std::endl;

  switch(file_type){
    case GDT_Unknown:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(file_type)<<std::endl;
      env.abort(-1); //TODO
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
      env.abort(-1); //TODO
    default:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(file_type)<<std::endl;
      env.abort(-1); //TODO
  }
}



int main(int argc, char **argv){
  boost::mpi::environment env;
  boost::mpi::communicator world;

  if(world.rank()==0){
    std::string many_or_one;
    std::string retention;
    std::string input_file;
    std::string output_prefix;
    int         bwidth  = -1;
    int         bheight = -1;
    int         flipH   = false;
    int         flipV   = false;

    std::string help=
    #include "help.txt"
    ;

    try{
      for(int i=1;i<argc;i++){
        if(strcmp(argv[i],"--bwidth")==0 || strcmp(argv[i],"-w")==0){
          if(i+1==argc)
            throw std::invalid_argument("-w followed by no argument.");
          bwidth = std::stoi(argv[i+1]);
          if(bwidth<300 && bwidth!=-1)
            throw std::invalid_argument("Width must be at least 500.");
          i++;
          continue;
        } else if(strcmp(argv[i],"--bheight")==0 || strcmp(argv[i],"-h")==0){
          if(i+1==argc)
            throw std::invalid_argument("-h followed by no argument.");
          bheight = std::stoi(argv[i+1]);
          if(bheight<300 && bheight!=-1)
            throw std::invalid_argument("Height must be at least 500.");
          i++;
          continue;
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
        } else if(output_prefix==""){
          output_prefix = argv[i];
        } else {
          throw std::invalid_argument("Too many arguments.");
        }
      }
      if(many_or_one=="" || retention=="" || input_file=="" || output_prefix=="")
        throw std::invalid_argument("Too few arguments.");
      if(retention[0]=='@' && !(retention=="@offloadall" || retention=="@retainall"))
        throw std::invalid_argument("Retention must be @offloadall or @retainall or a path.");
      if(many_or_one!="many" && many_or_one!="one")
        throw std::invalid_argument("Must specify many or one.");
    } catch (const std::invalid_argument &ia){
      std::string output_err;
      if(ia.what()==std::string("stoi"))
        output_err = "Invalid width or height.";
      else
        output_err = ia.what();
      std::cerr<<"###Error: "<<output_err<<std::endl;
      std::cerr<<help<<std::endl;
      std::cerr<<"###Error: "<<output_err<<std::endl;

      int good_to_go=0;
      boost::mpi::broadcast(world,good_to_go,0);

      return -1;
    }

    int good_to_go = 1;
    boost::mpi::broadcast(world,good_to_go,0);
    Preparer(many_or_one, retention, input_file, output_prefix, bwidth, bheight, flipH, flipV);

  } else {
    int good_to_go;
    boost::mpi::broadcast(world, good_to_go, 0);
    if(!good_to_go)
      return -1;

    GDALDataType file_type;
    boost::mpi::broadcast(world, file_type, 0);
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

  return 0;
}