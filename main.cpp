//Compile with
// mpic++ -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization
// mpirun -n 3 ./a.out ~/projects/watershed/data/beauford03.flt
// TODO: MPI abort
// TODO: See optimization notes at "http://www.boost.org/doc/libs/1_56_0/doc/html/mpi/tutorial.html"
// For memory usage see: http://info.prelert.com/blog/stl-container-memory-usage
#include "gdal_priv.h"
#include <iostream>
#include <boost/mpi.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <queue>
#include "Array2D.hpp"
#include "common.hpp"
//#define DEBUG 1



/*
  For reference, this is the definition of the RasterIO() function
  CPLErr GDALRasterBand::RasterIO( GDALRWFlag eRWFlag,
                                   int nXOff, int nYOff, int nXSize, int nYSize,
                                   void * pData, int nBufXSize, int nBufYSize,
                                   GDALDataType eBufType,
                                   int nPixelSpace,
                                   int nLineSpace )
*/

class Job {
 private:
  friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version){
      ar & edge;
      ar & gridx;
      ar & gridy;
      ar & filename;
      ar & x;
      ar & y;
      ar & width;
      ar & height;
    }
 public:
  char edge;
  int x,y,width,height,gridx,gridy;
  std::string filename;
  Job(){}
  Job(std::string filename, int gridx, int gridy, int x, int y, int width, int height){
    this->edge     = 0;
    this->gridx    = gridx;
    this->gridy    = gridy;
    this->filename = filename;
    this->x        = x;
    this->y        = y;
    this->width    = width;
    this->height   = height;
  }

};

typedef std::vector<Job> JobRow;
typedef std::vector<JobRow> JobGrid;
typedef int label_t;



const int TAG_JOB      = 0;
const int TAG_SYNC_MSG = 1;

const int SYNC_MSG_KILL = 0;
const int SYNC_MSG_JOB  = 1;

const char GRID_LEFT   = 1;
const char GRID_TOP    = 2;
const char GRID_RIGHT  = 4;
const char GRID_BOTTOM = 8;



template<class elev_t, GDALDataType gdt_t>
void PriorityFlood(Array2D<elev_t> &dem, Array2D<label_t> &labels, char edge){
  std::priority_queue<GridCellZ<elev_t>, std::vector<GridCellZ<elev_t> >, std::greater<GridCellZ<elev_t> > > pq;
  std::queue< GridCellZ<elev_t> > pit;
  typedef std::map<label_t, std::map<label_t, elev_t> > Graph;

  labels.init(0);

  Graph my_graph;

  int current_label = 2;

  for(int x=0;x<dem.width();x++){
    pq.emplace_back(x,0,dem(x,0));
    pq.emplace_back(x,dem.height()-1,dem(x,dem.height()-1));
    if(edge & GRID_TOP)
      labels(x,0) = -1;
    if(edge & GRID_BOTTOM)
      labels(x,dem.height()-1) = -1;
  }
  for(int y=0;y<dem.height();y++){
    pq.emplace_back(0,y,dem(0,y));
    pq.emplace_back(dem.width()-1,y,dem(dem.width()-1,y));
    if(edge & GRID_LEFT)
      labels(0,y) = -1;
    if(edge & GRID_RIGHT)
      labels(dem.width()-1,y) = -1;
  }

  while(!pq.empty() || !pit.empty()){
    GridCellZ<elev_t> c;
    if(pit.size()>0){
      c=pit.front();
      c.pop();
    } else {
      c=pq.top();
      pq.pop();
    }

    //Since labels are inherited from parent cells we need to be able to process
    //previously labeled cells. But all the edge cells are already in the
    //my_open queue and may also be added to that queue by a parent cell. So
    //they could be processed twice! The solution is to assign the negative of a
    //label and then make the label positive when we actually process the cell.
    //This way, if we ever see we are processing a cell with a positive label,
    //we will know it was already processed.
    if(labels(c.x,c.y)>0)             //Positive label. Cell already processed.
      continue;
    else if(labels(c.x,c.y)==0)           //Label=0. Cell has never been seen.
      labels(c.x,c.y) = current_label++;  //current_label is the next new label.
    else if(labels(c.x,c.y)<0)            //Cell's parent added it.
      labels(c.x,c.y) = -labels(c.x,c.y); //Mark cell as visited

    //At this point the cell's label is guaranteed to be positive and in the
    //range [1,MAX_INTEGER] (unless we overflow).
    auto my_label = labels(c.x,c.y);

    for(int n=1;n<=8;n++){
      auto nx = c.x+dx[n]; //Neighbour's x-coordinate
      auto ny = c.y+dy[n]; //Neighbour's y-coordinate

      //Check to see if the neighbour coordinates are valid
      if(nx<0 || ny<0 || nx==dem.width() || ny==dem.viewHeight()) continue;
      auto n_label = std::abs(labels(nx,ny)); //Neighbour's label
      //Does the neighbour have a label? If so, it is part of the edge, has
      //already been assigned a label by a parent cell which must be of lower or
      //equal elevation to the current cell, or has already been processed, in
      //which case its elevation is lower or equal to this cell.
      if(n_label!=0){
        //If the neighbour's label were the same as the current cell's, then the
        //current cell's flow and the neighbour's flow eventually comingle. If
        //the neighbour's label is different it has been added by a cell whose
        //flow drains the opposite side of a watershed from this cell. Here, we
        //make a note of the height of that watershed.
        if(n_label!=my_label){
          auto elev_over = std::max(dem(nx,ny),dem(c.x,c.y)); //TODO: I think this should always be the neighbour.
          //Haven't seen this watershed before
          if(my_graph[my_label].count(n_label)==0){
            my_graph[my_label][n_label] = elev_over;
            my_graph[n_label][my_label] = elev_over;
          //We've seen this watershed before, so only make a note of the
          //spill-over elevation if it is lower than what we've seen before.
          } else if(elev_over<my_graph[my_label][n_label]){
            my_graph[my_label][n_label] = elev_over;
            my_graph[n_label][my_label] = elev_over;
          }
        }
        continue;
      }

      //The neighbour is not one we've seen before, so mark it as being part of
      //our watershed and add it as an unprocessed item to the queue.
      my_labels(nx,ny) = -my_labels(c.x,c.y);


      if(dem(nx,ny)!=dem.noData())        //The neighbour is part of the DEM's data
        my_flowdirs(nx,ny) = d8_inverse[n]; //and flows into this cell

      //If the neighbour is lower than this cell, elevate it to the level of
      //this cell so that a depression is not formed. The flow directions will
      //be fixed later, after all the depressions have been filled.
      if(dem(nx,ny)<=c.z){
        dem(nx,ny) = c.z;
        pit.emplace(nx,ny,c.z);
      //Otherwise, if the neighbour is higher, do not adjust its elevation.
      } else
        pq.emplace(nx,ny,dem(nx,ny));
    }
  }

}


template<class elev_t, GDALDataType gdt_t>
void Serf(){
  boost::mpi::environment env;
  boost::mpi::communicator world;
  int sync_message;

  while(true){
    world.recv(0, TAG_SYNC_MSG, sync_message);
    if(sync_message==SYNC_MSG_KILL)
      return;
    else if (sync_message==SYNC_MSG_JOB){
      Job job;
      world.recv(0, TAG_JOB, job);
      std::cout<<"Node "<<world.rank()<<" has job "<<job.x<<"-"<<(job.x+job.width)<<" "<< job.y<<"-"<<(job.y+job.height)<<std::endl;
      Array2D<elev_t> grid(job.filename, true, job.x, job.y, job.x+job.width, job.y+job.height);


      world.send(0, TAG_JOB, job);
    }
  }
}



template<class elev_t, GDALDataType gdt_t>
void Master(JobGrid &jobs){
  boost::mpi::environment env;
  boost::mpi::communicator world;
  int active_nodes  = 0;

  for(size_t y=0;y<jobs.size();y++)
  for(size_t x=0;x<jobs[0].size();x++){
    std::cerr<<"Sending job "<<y<<"/"<<jobs.size()<<", "<<x<<"/"<<jobs[0].size()<<std::endl;
    if(active_nodes<world.size()-1){
      // std::cerr<<"Sending init to "<<(active_nodes+1)<<std::endl;
      world.send(active_nodes+1,TAG_SYNC_MSG,SYNC_MSG_JOB);
      world.isend(active_nodes+1,TAG_JOB,jobs[y][x]);
      active_nodes++;
    } else {
      Job temp;
      boost::mpi::status status    = world.recv(boost::mpi::any_source,TAG_JOB,temp);
      jobs[temp.gridy][temp.gridx] = temp;
      world.send (status.source(),TAG_SYNC_MSG,SYNC_MSG_JOB);
      world.isend(status.source(),TAG_JOB     ,jobs[y][x]  );
    }
  }

  while(active_nodes>0){
    Job temp;
    boost::mpi::status status    = world.recv(boost::mpi::any_source,TAG_JOB,temp);
    jobs[temp.gridy][temp.gridx] = temp;
    active_nodes--;
    (void)status; //TODO: Suppress non-usage warning
  }

  for(int i=1;i<world.size();i++)
    world.isend(i,TAG_SYNC_MSG,SYNC_MSG_KILL);

}





int main(int argc, char **argv){
  boost::mpi::environment env;
  boost::mpi::communicator world;
  std::string filename;

  if(world.rank()==0){
    if(argc<5){
      std::cerr<<"Syntax: "<<argv[0]<<" <many> <retain/offload> <dir> <extension>\n";
      std::cerr<<"Syntax: "<<argv[0]<<" <one>  <retain/offload> <file> <bwidth> <bheight>\n";
      std::cerr<<"    many      - Implies that the data has already been tiled and all tiles are\n";
      std::cerr<<"                in <dir> and end with <extension>. Files are assumed to be small\n";
      std::cerr<<"                enough that each individual file can fit into RAM. File must be\n";
      std::cerr<<"                non-overlapping square blocks.\n";
      std::cerr<<"    one       - Implies that the data is in a single file and must be split up\n";
      std::cerr<<"                into blocks of size <bwidth> and <bheight>\n";
      std::cerr<<"    retain    - Threads should retain all relevant information in memory\n";
      std::cerr<<"                throughout the life of the program. Good for smaller datasets or\n";
      std::cerr<<"                high RAM environments.\n";
      std::cerr<<"    offload   - Threads will write intermediate products to memory. Good for\n";
      std::cerr<<"                large datasets, low RAM environments, or unstable environments.\n";
      std::cerr<<"    file      - Single data file to be processed by <one>\n";
      std::cerr<<"    extension - Final characters of file name, such as: tif, flt\n";
      std::cerr<<"    dir       - Directory in which files to be processed by <many> are located\n";
      std::cerr<<"    bwidth    - Block width in cells. Should be >1000\n";
      std::cerr<<"    bheight   - Block height in cells. Should be >1000\n";
      return -1;
    }

    JobGrid     jobs;
    std::string filename;

    if(argv[1]==std::string("many")){
      std::cout<<"This feature is still experimental!\n";
      return -1;

      std::cout<<"Locating files...\n";
      //Convert argv[4], the extension argument, into a string for easy
      //manipulation
      std::string extension(argv[4]);
      //Recursively loop through each file in the indicated directory and its
      //subdirectories. Check to see if it ends with the appropriate extension
      //and, if it does, add it to the list of files to be processed.
      for (boost::filesystem::recursive_directory_iterator i(argv[3]), end; i != end; ++i)
        if (!is_directory(i->path())){
          std::string this_filename=i->path().filename().string();
          if (0==this_filename.compare(this_filename.length() - extension.length(), extension.length(), extension))
            //files_to_process.push_back(this_filename);
            std::cout << "   " << this_filename << "\n";
        }
      std::cout<<std::flush; //Ensure all file names are written before we continue
    } else if(argv[1]==std::string("one")) {
      std::cerr<<"Single file mode"<<std::endl;
      int bwidth  = std::stoi(argv[4]);
      int bheight = std::stoi(argv[5]);
      int total_height;
      int total_width;

      filename = argv[3];

      getGDALDimensions(filename, total_height, total_width);

      if(bwidth==-1)
        bwidth  = total_width;
      if(bheight==-1)
        bheight = total_height;

      std::cerr<<"Total width:  "<<total_width <<"\n";
      std::cerr<<"Total height: "<<total_height<<"\n";
      std::cerr<<"Block width:  "<<bwidth      <<"\n";
      std::cerr<<"Block height: "<<bheight     <<std::endl;

      for(int y=0,gridy=0;y<total_height; y+=bheight, gridy++){
        jobs.push_back(JobRow());
        for(int x=0,gridx=0;x<total_width;x+=bwidth,  gridx++){
          std::cerr<<"Constructing job ("<<x<<","<<y<<")"<<std::endl;
          jobs.back().emplace_back(filename, gridx, gridy, x, y, bwidth, bheight);
        }
      }

      //Indicate edges
      for(int x=0;x<(int)(total_width/bwidth);x++){
        jobs[0][x].edge             |= GRID_TOP;
        jobs[jobs.size()-1][x].edge |= GRID_BOTTOM;
      }
      for(int y=0;y<(int)(total_height/bheight);y++){
        jobs[y][0].edge                |= GRID_LEFT;
        jobs[y][jobs[y].size()-1].edge |= GRID_RIGHT;
      }

    } else {
      std::cout<<"Unrecognised option! Must be 'many' or 'one'!"<<std::endl;
      return -1;
    }

    switch(peekGDALType(filename)){
      case GDT_Int16:
        Master<int16_t, GDT_Int16>(jobs);
        break;
      case GDT_Float32:
        Master<float, GDT_Float32>(jobs);
        break;
      default:
        if(world.rank()==0){
          std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[3]))<<std::endl;
          return -1;
        }
        break;
    }

  } else {
    if(argc<5)
      return -1;

    std::string filename(argv[3]);

    switch(peekGDALType(filename)){
      case GDT_Int16:
        Serf<int16_t, GDT_Int16>();
        break;
      case GDT_Float32:
        Serf<float, GDT_Float32>();
        break;
      default:
        break;
    }
  }

  return 0;
}