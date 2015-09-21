//Compile with
// mpic++ -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization
// mpirun -n 3 ./a.out ~/projects/watershed/data/beauford03.flt
// TODO: MPI abort
// TODO: See optimization notes at "http://www.boost.org/doc/libs/1_56_0/doc/html/mpi/tutorial.html"
// For memory usage see: http://info.prelert.com/blog/stl-container-memory-usage
#include "gdal_priv.h"
#include <iostream>
#include <iomanip>
#include <boost/mpi.hpp>
#include <boost/filesystem.hpp>
#include <boost/serialization/map.hpp>
#include <string>
#include <queue>
#include <vector>
#include <limits>

//We use the cstdint library here to ensure that the program behaves as expected
//across platforms, especially with respect to the expected limits of operation
//for data set sizes and labels. For instance, in C++, a signed integer must be
//at least 16 bits, but not necessarily more. We force a minimum of 32 bits as
//this is, after all, for use with large datasets.
#include <cstdint>
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

class ChunkInfo{
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & edge;
    ar & x;
    ar & y;
    ar & width;
    ar & height;
    ar & gridx;
    ar & gridy;
    ar & label_offset;
    ar & max_label;
    ar & filename;
    ar & id;
  }
 public:
  uint8_t edge;
  int32_t x,y,width,height,gridx,gridy;
  int32_t label_offset,max_label;
  int32_t id;
  std::string filename;
  ChunkInfo(){}
  ChunkInfo(std::string filename, int32_t id, int32_t label_offset, int32_t max_label, int32_t gridx, int32_t gridy, int32_t x, int32_t y, int32_t width, int32_t height){
    this->edge         = 0;
    this->x            = x;
    this->y            = y;
    this->width        = width;
    this->height       = height;   
    this->gridx        = gridx;
    this->gridy        = gridy;
    this->label_offset = label_offset;
    this->max_label    = max_label;
    this->id           = id;
    this->filename     = filename;
  }
};

template<class elev_t>
class Job1 {
 private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & top_elev;
    ar & bot_elev;
    ar & left_elev;
    ar & right_elev;
    ar & top_label;
    ar & bot_label;
    ar & left_label;
    ar & right_label;
    ar & graph;
  }
 public:
  std::vector<elev_t > top_elev,  bot_elev,  left_elev,  right_elev;  //TODO: Consider using std::array instead
  std::vector<label_t> top_label, bot_label, left_label, right_label; //TODO: Consider using std::array instead
  std::map<label_t, std::map<label_t, elev_t> > graph;
  std::map<label_t,elev_t> graph_elev;
  Job1(){}
};

typedef int32_t label_t;


//TODO: What are these for?
const int TAG_WHICH_JOB   = 0;
const int TAG_CHUNK_DATA  = 1;
const int TAG_DONE_FIRST  = 2;
const int TAG_SECOND_DATA = 3;
const int TAG_DONE_SECOND = 4;

const int SYNC_MSG_KILL = 0;
const int JOB_CHUNK     = 1;
const int JOB_FIRST     = 2;
const int JOB_SECOND    = 3;

const uint8_t GRID_LEFT   = 1;
const uint8_t GRID_TOP    = 2;
const uint8_t GRID_RIGHT  = 4;
const uint8_t GRID_BOTTOM = 8;



template<class elev_t>
void PriorityFlood(Array2D<elev_t> &dem, Array2D<label_t> &labels, std::map<label_t, std::map<label_t, elev_t> > &my_graph, uint8_t edge){
  typedef std::priority_queue<GridCellZ<elev_t>, std::vector<GridCellZ<elev_t> >, std::greater<GridCellZ<elev_t> > > GridCellZ_pq;
  GridCellZ_pq pq;
  std::queue< GridCellZ<elev_t> > pit;

  labels.init(0);

  int32_t current_label = 2;

  for(int x=0;x<dem.viewWidth();x++){
    pq.emplace(x,0,dem(x,0));
    pq.emplace(x,dem.viewHeight()-1,dem(x,dem.viewHeight()-1));
  }

  for(int y=1;y<dem.viewHeight()-1;y++){
    pq.emplace(0,y,dem(0,y));
    pq.emplace(dem.viewWidth()-1,y,dem(dem.viewWidth()-1,y));
  }

  if(edge & GRID_TOP)
    for(int x=0;x<dem.viewWidth();x++)
      labels(x,0) = -1;

  if(edge & GRID_BOTTOM)
    for(int x=0;x<dem.viewWidth();x++)
      labels(x,dem.viewHeight()-1) = -1;

  if(edge & GRID_LEFT)
    for(int y=0;y<dem.viewHeight();y++)
      labels(0,y) = -1;

  if(edge & GRID_RIGHT)
    for(int y=0;y<dem.viewHeight();y++)
      labels(dem.viewWidth()-1,y) = -1;

  while(!pq.empty() || !pit.empty()){
    GridCellZ<elev_t> c;
    if(pit.size()>0){
      c=pit.front();
      pit.pop();
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
      if(nx<0 || ny<0 || nx==dem.viewWidth() || ny==dem.viewHeight()) continue;
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
      labels(nx,ny) = -labels(c.x,c.y);


//      if(dem(nx,ny)!=dem.noData())        //The neighbour is part of the DEM's data
//        my_flowdirs(nx,ny) = d8_inverse[n]; //and flows into this cell

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
void Consumer(){
  boost::mpi::environment env;
  boost::mpi::communicator world;

  int the_job;   //Which job should consumer perform?

  ChunkInfo chunk;

  //Have the consumer process messages as long as they are coming using a
  //blocking receive to wait.
  while(true){
    world.recv(0, TAG_WHICH_JOB, the_job); //Blocking wait


    //This message indicates that everything is done and the Consumer should shut
    //down.
    if(the_job==SYNC_MSG_KILL){
      std::cerr<<"Killing "<<world.rank()<<std::endl;
      MPI_Finalize();
      return;

    } else if (the_job==JOB_CHUNK){
      world.recv(0, TAG_CHUNK_DATA, chunk);

    //This message indicates that the consumer should prepare to perform the
    //first part of the distributed Priority-Flood algorithm on an incoming job
    } else if (the_job==JOB_FIRST){
      Job1<elev_t> job1;
      std::cout<<"Node "<<world.rank()<<" has job "<<chunk.x<<"-"<<(chunk.x+chunk.width)<<" "<< chunk.y<<"-"<<(chunk.y+chunk.height)<<std::endl;

      //Read in the data associated with the job
      Array2D<elev_t>   dem(chunk.filename, true, chunk.x, chunk.y, chunk.width, chunk.height);

      //These variables are needed by Priority-Flood. The internal
      //interconnections of labeled regions (named "graph") are also needed to
      //solve the problem, but that can be passed directly from the job object.
      Array2D<label_t>  labels(dem.viewWidth(),dem.viewHeight(),0);

      //Perform the usual Priority-Flood algorithm on the chunk. TODO: Use the
      //faster algorithm by Zhou, Sun, Fu
      PriorityFlood(dem,labels,job1.graph,chunk.edge);

      //The chunk's edge info is needed to solve the global problem. Collect it.
      job1.top_elev    = dem.topRow     ();
      job1.bot_elev    = dem.bottomRow  ();
      job1.left_elev   = dem.leftColumn ();
      job1.right_elev  = dem.rightColumn();


      job1.top_label   = labels.topRow     ();
      job1.bot_label   = labels.bottomRow  ();
      job1.left_label  = labels.leftColumn ();
      job1.right_label = labels.rightColumn();

      std::cerr<<"TE: "<<job1.top_elev.size()<<" TL: "<<job1.top_label.size()<<std::endl;
      //job.graph = graph;

      world.send(0, TAG_DONE_FIRST, job1);
    } else if (the_job==JOB_SECOND){
      std::map<label_t, elev_t> graph_elev;
      world.recv(0, TAG_SECOND_DATA, graph_elev);

      std::map<label_t, std::map<label_t, elev_t> > graph;

      //These use the same logic as the analogous lines above
      Array2D<elev_t>   dem(chunk.filename, true, chunk.x, chunk.y, chunk.width, chunk.height);
      Array2D<label_t>  labels(dem.viewWidth(),dem.viewHeight(),0);
      PriorityFlood(dem,labels,graph,chunk.edge);
      for(int y=0;y<dem.viewHeight();y++)
      for(int x=0;x<dem.viewWidth();x++)
        if(graph_elev.count(labels(x,y))){
          if(dem(x,y)<graph_elev[labels(x,y)])
            dem(x,y)=graph_elev[labels(x,y)];
        }

      std::cerr<<"Finished: "<<chunk.gridx<<" "<<chunk.gridy<<std::endl;

      //At this point we're done with the calculation! Boo-yeah!

      std::string output_name = std::string("/z/output")+std::to_string(chunk.id)+std::string(".tif");
      dem.saveGDAL(output_name, chunk.filename, chunk.x, chunk.y);

      world.send(0, TAG_DONE_SECOND, 0);
    }
  }
}




template<class elev_t>
void HandleEdge(
  const std::vector<elev_t>  &elev_a,
  const std::vector<elev_t>  &elev_b,
  const std::vector<label_t> &label_a,
  const std::vector<label_t> &label_b,
  std::map<label_t, std::map<label_t, elev_t> > &mastergraph//,
  //const elev_t no_data
){
  //Guarantee that all vectors are of the same length
  assert(elev_a.size ()==elev_b.size ());
  assert(label_a.size()==label_b.size());
  assert(elev_a.size ()==label_b.size());

  int len = elev_a.size();

  for(size_t i=0;i<len;i++){
    //if(elev_a[i]==no_data) //TODO: Does no_data really matter here?
    //  continue;

    auto my_label = label_a[i];

    for(int di=-1;i<=1;i++){
      int ni = i+di;
      if(ni<0 || ni==len) // || elev_b[ni]==no_data) //TODO: Does no_data really matter here?
        continue;
      auto other_label = label_b[ni];
      //TODO: Does this really matter? We could just ignore these entries
      if(my_label==other_label) //Only happens when labels are both 1
        continue;

      auto elev_over = std::max(elev_a[i],elev_b[ni]);
      if(mastergraph[my_label].count(other_label)==0){
        mastergraph[my_label][other_label] = elev_over;
        mastergraph[other_label][my_label] = elev_over;
      } else if(elev_over<mastergraph[my_label][other_label]){
        mastergraph[my_label][other_label] = elev_over;
        mastergraph[other_label][my_label] = elev_over;
      }
    }
  }
}

template<class elev_t>
void HandleCorner(
  const elev_t  elev_a,
  const elev_t  elev_b,
  const label_t label_a,
  const label_t label_b,
  std::map<label_t, std::map<label_t, elev_t> > &mastergraph//,
  //const elev_t no_data
){
//  if(elev_a==no_data || elev_b==no_data)
//    return;

  auto elev_over = std::max(elev_a,elev_b);
  if(mastergraph[label_a].count(label_b)==0){
    mastergraph[label_a][label_b] = elev_over;
    mastergraph[label_b][label_a] = elev_over;
  } else if(elev_over<mastergraph[label_a][label_b]){
    mastergraph[label_a][label_b] = elev_over;
    mastergraph[label_b][label_a] = elev_over;
  }
}







//Producer takes a collection of Jobs and delegates them to Consumers. Once all
//of the jobs have received their initial processing, it uses that information
//to compute the global properties necessary to the solution. Each Job, suitably
//modified, is then redelegated to a Consumer which ultimately finishes the
//processing.
template<class elev_t, GDALDataType gdt_t>
void Producer(std::vector< std::vector< ChunkInfo > > &chunks){
  boost::mpi::environment env;
  boost::mpi::communicator world;
  int active_nodes = 0;

  std::map<int,ChunkInfo> rank_to_chunk;

  std::vector< std::vector< Job1<elev_t> > > jobs1(chunks.size(), std::vector< Job1<elev_t> >(chunks[0].size()));


  //Loop through all of the jobs, delegating them to Consumers
  active_nodes=0;
  for(size_t y=0;y<chunks.size();y++)
  for(size_t x=0;x<chunks[0].size();x++){
    std::cerr<<"Sending job "<<y<<"/"<<chunks.size()<<", "<<x<<"/"<<chunks[0].size()<<std::endl;

    //If fewer jobs have been delegated than there are Consumers available,
    //delegate the job to a new Consumer.
    if(active_nodes<world.size()-1){
      // std::cerr<<"Sending init to "<<(active_nodes+1)<<std::endl;
      world.send(active_nodes+1,TAG_WHICH_JOB,JOB_CHUNK);
      world.send(active_nodes+1,TAG_CHUNK_DATA,chunks[y][x]);

      rank_to_chunk[active_nodes+1] = chunks[y][x];
      world.isend(active_nodes+1,TAG_WHICH_JOB,JOB_FIRST);
      active_nodes++;

    //Once all of the consumers are active, wait for them to return results. As
    //each Consumer returns a result, pass it the next unfinished Job until
    //there are no jobs left.
    } else {
      Job1<elev_t> finished_job;

      //Execute a blocking receive until some consumer finishes its work.
      //Receive that work.
      boost::mpi::status status = world.recv(boost::mpi::any_source,TAG_DONE_FIRST,finished_job);

      //NOTE: This could be used to implement more robust handling of lost nodes
      ChunkInfo received_chunk = rank_to_chunk[status.source()];
      jobs1[received_chunk.gridy][received_chunk.gridx] = finished_job;

      //Delegate new work to that consumer
      world.send(status.source(),TAG_WHICH_JOB,JOB_CHUNK);
      world.send(status.source(),TAG_CHUNK_DATA,chunks[y][x]);

      rank_to_chunk[status.source()] = chunks[y][x];
      world.isend (status.source(),TAG_WHICH_JOB,JOB_FIRST);
    }
  }

  while(active_nodes>0){
    Job1<elev_t> finished_job;

    //Execute a blocking receive until some consumer finishes its work.
    //Receive that work
    boost::mpi::status status = world.recv(boost::mpi::any_source,TAG_DONE_FIRST,finished_job);
    ChunkInfo received_chunk = rank_to_chunk[status.source()];
    jobs1[received_chunk.gridy][received_chunk.gridx] = finished_job;

    //Decrement the number of consumers we are waiting on. When this hits 0 all
    //of the jobs have been completed and we can move on
    active_nodes--;
  }


  //Merge all of the graphs together into one very big graph. Clear information
  //as we go in order to save space, though I am not sure if the map::clear()
  //method is not guaranteed to release space.
  std::map<label_t, std::map<label_t, elev_t> > mastergraph;
  for(auto &row: jobs1)
  for(auto &cell: row){
    for(auto const &fkey: cell.graph)
    for(auto const &skey: fkey.second)
      mastergraph[fkey.first][skey.first] = skey.second;
    cell.graph.clear();
  }

  //Now look at the chunks around a central chunk c.
  const int gridwidth  = jobs1.front().size();
  const int gridheight = jobs1.size();

  for(size_t y=0;y<jobs1.size();y++)
  for(size_t x=0;x<jobs1.front().size();x++){
    auto &c = jobs1[y][x];

    if(y>0)
      HandleEdge(c.top_elev,   jobs1[y-1][x].bot_elev,   c.top_label,   jobs1[y-1][x].bot_label,   mastergraph);
    if(y<gridheight-1)
      HandleEdge(c.bot_elev,   jobs1[y+1][x].top_elev,   c.bot_label,   jobs1[y+1][x].top_label,   mastergraph);
    if(x>0)
      HandleEdge(c.left_elev,  jobs1[y][x-1].right_elev, c.left_label,  jobs1[y][x-1].right_label, mastergraph);
    if(x<gridwidth-1)
      HandleEdge(c.right_elev, jobs1[y][x+1].left_elev,  c.right_label, jobs1[y][x+1].left_label,  mastergraph);


    //I wish I had wrote it all in LISP.
    if(y>0 && x>0)                      //Top left
      HandleCorner(c.top_elev.front(), jobs1[y-1][x-1].bot_elev.back(),  c.top_label.front(), jobs1[y-1][x-1].bot_label.back(),  mastergraph);
    if(y<gridheight-1 && x<gridwidth-1) //Bottom right
      HandleCorner(c.bot_elev.back(),  jobs1[y+1][x+1].top_elev.front(), c.bot_label.back(),  jobs1[y+1][x+1].top_label.front(), mastergraph);
    if(y>0 && x<gridwidth-1)            //Top right
      HandleCorner(c.top_elev.back(),  jobs1[y-1][x+1].bot_elev.front(), c.top_label.back(),  jobs1[y-1][x+1].bot_label.front(), mastergraph);
    if(x>0 && y<gridheight-1)           //Bottom left
      HandleCorner(c.bot_elev.front(), jobs1[y+1][x-1].top_elev.back(),  c.bot_label.front(), jobs1[y+1][x-1].top_label.back(),  mastergraph);
  }

  jobs1.clear();
  jobs1.shrink_to_fit();


  std::cerr<<"Mastergraph constructed!"<<std::endl;


  std::cerr<<"Performing aggregated priority flood"<<std::endl;
  typedef std::pair<elev_t, label_t>  graph_node;
  std::priority_queue<graph_node, std::vector<graph_node>, std::greater<graph_node> > open;
  std::map<label_t,bool>              visited;
  std::map<label_t,elev_t>            graph_elev;

  open.emplace(std::numeric_limits<elev_t>::min(),1);

  while(open.size()>0){
    graph_node c=open.top();
    open.pop();
    auto my_elev       = c.first;
    auto my_vertex_num = c.second;
    #ifdef DEBUG
      std::cerr<<"Popped "<<my_vertex_num<<std::endl;
    #endif
    if(visited[my_vertex_num])
      continue;

    graph_elev[my_vertex_num] = my_elev;
    visited[my_vertex_num]    = true;

    for(auto &n: mastergraph[my_vertex_num]){
      auto n_vertex_num = n.first;
      auto n_elev       = n.second;
      if(visited.count(n_vertex_num))
        continue;
      #ifdef DEBUG
        std::cerr<<"Proposing going to "<<n_vertex_num<<" with "<<std::max(n_elev,my_elev)<<std::endl;
      #endif
      open.emplace(std::max(n_elev,my_elev),n_vertex_num);
    }
  }

  #ifdef DEBUG
    std::cerr<<"Graph elevations"<<std::endl;
    for(const auto &ge: graph_elev)
      std::cerr<<std::setw(3)<<ge.first<<" = "<<std::setw(3)<<ge.second<<std::endl;
  #endif

  std::cerr<<"Sending out final jobs..."<<std::endl;
  //Loop through all of the jobs, delegating them to Consumers
  active_nodes=0;
  for(size_t y=0;y<chunks.size();y++)
  for(size_t x=0;x<chunks[0].size();x++){
    std::cerr<<"Sending job "<<y<<"/"<<chunks.size()<<", "<<x<<"/"<<chunks[0].size()<<std::endl;

    std::map<elev_t, label_t> job2;
    for(const auto &ge: graph_elev)
      if(chunks[y][x].label_offset<=ge.first && ge.first<=chunks[y][x].max_label)
        job2[ge.first] = ge.second;

    //If fewer jobs have been delegated than there are Consumers available,
    //delegate the job to a new Consumer.
    if(active_nodes<world.size()-1){
      // std::cerr<<"Sending init to "<<(active_nodes+1)<<std::endl;
      world.send(active_nodes+1,TAG_WHICH_JOB,JOB_CHUNK);
      world.send(active_nodes+1,TAG_CHUNK_DATA,chunks[y][x]);

      rank_to_chunk[active_nodes+1] = chunks[y][x];
      world.send (active_nodes+1,TAG_WHICH_JOB,JOB_SECOND);
      world.isend(active_nodes+1,TAG_SECOND_DATA,job2);
      active_nodes++;

    //Once all of the consumers are active, wait for them to return results. As
    //each Consumer returns a result, pass it the next unfinished Job until
    //there are no jobs left.
    } else {
      int msg;

      //Execute a blocking receive until some consumer finishes its work.
      //Receive that work.
      boost::mpi::status status = world.recv(boost::mpi::any_source,TAG_DONE_SECOND,msg);

      //Delegate new work to that consumer
      world.send(status.source(),TAG_WHICH_JOB,JOB_CHUNK);
      world.send(status.source(),TAG_CHUNK_DATA,chunks[y][x]);

      rank_to_chunk[status.source()] = chunks[y][x];
      world.send (status.source(),TAG_WHICH_JOB,JOB_SECOND);
      world.isend(status.source(),TAG_SECOND_DATA,job2);
    }
  }

  while(active_nodes>0){
    int msg;
    //Execute a blocking receive until some consumer finishes its work.
    //Receive that work
    world.recv(boost::mpi::any_source,TAG_DONE_SECOND,msg);


    //Decrement the number of consumers we are waiting on. When this hits 0 all
    //of the jobs have been completed and we can move on
    active_nodes--;

    std::cerr<<"Nodes left: "<<active_nodes<<std::endl;
  }

  for(int i=1;i<world.size();i++)
    world.isend(i,TAG_WHICH_JOB,SYNC_MSG_KILL);
  MPI_Finalize();
}




//Preparer divides up the input raster file into chunks which can be processed
//independently by the Consumers. Since the chunking may be done on-the-fly or
//rely on preparation the user has done, the Preparer routine knows how to deal
//with both. Once assemebled, the collection of jobs is passed off to Producer,
//which is agnostic as to the original form of the jobs and handles
//communication and solution assembly.
template<class elev_t, GDALDataType gdt_t>
int Preparer(int argc, char **argv){
  boost::mpi::environment env;
  int chunkid = 0;

  std::vector< std::vector< ChunkInfo > >     chunks;
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
    int32_t bwidth  = std::stoi(argv[4]);
    int32_t bheight = std::stoi(argv[5]);
    int32_t total_height;
    int32_t total_width;

    filename = argv[3];

    //Get the total dimensions of the input file
    getGDALDimensions(filename, total_height, total_width);

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
    int32_t label_offset = 2;
    const int32_t label_increment = 2*bheight+2*bwidth+1;
    for(int32_t y=0,gridy=0;y<total_height; y+=bheight, gridy++){
      chunks.push_back( std::vector< ChunkInfo >() );
      for(int32_t x=0,gridx=0;x<total_width;x+=bwidth,  gridx++){
        std::cerr<<"Constructing job ("<<x<<","<<y<<")"<<std::endl;
        chunks.back().emplace_back(filename, chunkid++, label_offset, label_offset+label_increment-1, gridx, gridy, x, y, bwidth, bheight);
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
        if(std::numeric_limits<int>::max()-label_offset<label_increment){
          std::cerr<<"Ran out of labels. Cannot proceed."<<std::endl;
          env.abort(-1); //TODO: Set error code
        }
        label_offset += label_increment;
      }
    }

    //#ifdef DEBUG
      for(int y=0;y<chunks.size();y++)
      for(int x=0;x<chunks[0].size();x++){
        std::cerr<<"Chunk offset: "<<y<<" "<<x<<" "<<chunks[y][x].label_offset<<" "<<chunks[y][x].max_label<< std::endl;
      }
    //#endif
  
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

  } else {
    std::cout<<"Unrecognised option! Must be 'many' or 'one'!"<<std::endl;
    return -1;
  }

  //Pass off the list of jobs to Producer.
  Producer<elev_t, gdt_t>(chunks);

  return 0;
}



int main(int argc, char **argv){
  boost::mpi::environment env;
  boost::mpi::communicator world;

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

    switch(peekGDALType(argv[3])){
      case GDT_Int16:
        return Preparer<int16_t, GDT_Int16>(argc,argv);
      case GDT_Float32:
        return Preparer<float, GDT_Float32>(argc,argv);
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

    switch(peekGDALType(argv[3])){
      case GDT_Int16:
        Consumer<int16_t, GDT_Int16>();
        break;
      case GDT_Float32:
        Consumer<float, GDT_Float32>();
        break;
      default:
        break;
    }
  }

  return 0;
}