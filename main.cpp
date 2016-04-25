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
#include <iomanip>
#include "Zhou2015pf.hpp"
#include "Barnes2014pf.hpp"
#include <string>
#include <queue>
#include <vector>
#include <limits>
#include <fstream> //For reading layout files
#include <sstream> //Used for parsing the <layout_file>
#include <boost/filesystem.hpp>
#include "communication.hpp"
#include "memory.hpp"

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

typedef uint32_t label_t;

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
       id,
       nullChunk,
       filename,
       outputname,
       retention);
  }
 public:
  uint8_t     edge;
  uint8_t     flip;
  int32_t     x,y,width,height,gridx,gridy;
  int32_t     id;
  bool        nullChunk;
  label_t     label_offset,label_increment; //Used for convenience in Producer()
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

template<class elev_t>
class Job1 {
 private:
  friend class cereal::access;
  template<class Archive>
  void serialize(Archive & ar){
    ar(top_elev,
       bot_elev,
       left_elev,
       right_elev,
       top_label,
       bot_label,
       left_label,
       right_label,
       graph,
       time_info);
  }
 public:
  std::vector<elev_t > top_elev,  bot_elev,  left_elev,  right_elev;  //TODO: Consider using std::array instead
  std::vector<label_t> top_label, bot_label, left_label, right_label; //TODO: Consider using std::array instead
  std::vector< std::map<label_t, elev_t> > graph;
  TimeInfo time_info;
  Job1(){}
};

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

const uint8_t FLIP_VERT   = 1;
const uint8_t FLIP_HORZ   = 2;


template<class elev_t>
void Consumer(){
  ChunkInfo chunk;

  Array2D<elev_t>  dem;
  Array2D<label_t> labels;

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

    } else if (the_job==JOB_CHUNK){
      assert(false);
    //This message indicates that the consumer should prepare to perform the
    //first part of the distributed Priority-Flood algorithm on an incoming job
    } else if (the_job==JOB_FIRST){
      CommRecv(&chunk, nullptr, 0);

      std::cerr<<"Rec name: "<<chunk.filename<<std::endl;

      Timer timer_calc,timer_io,timer_overall;
      timer_overall.start();

      Job1<elev_t> job1;

      //The upper limit on unique watersheds is the number of edge cells. Resize
      //the graph to this number. The Priority-Flood routines will shrink it to
      //the actual number needed.
      job1.graph.resize(2*chunk.width+2*chunk.height);

      //Read in the data associated with the job
      timer_io.start();
      dem = Array2D<elev_t>(chunk.filename, false, chunk.x, chunk.y, chunk.width, chunk.height);
      if(chunk.flip & FLIP_VERT)
        dem.flipVert();
      if(chunk.flip & FLIP_HORZ)
        dem.flipHorz();
      timer_io.stop();

      //These variables are needed by Priority-Flood. The internal
      //interconnections of labeled regions (named "graph") are also needed to
      //solve the problem, but that can be passed directly from the job object.
      labels = Array2D<label_t>(dem.viewWidth(),dem.viewHeight(),0);

      //Perform the watershed Priority-Flood algorithm on the chunk. The variant
      //by Zhou, Sun, and Fu (2015) is used for this; however, I have modified
      //their algorithm to label watersheds similarly to what is described in
      //Barnes, Lehman, and Mulla (2014).

      timer_calc.start();
      Zhou2015Labels(dem,labels,job1.graph,chunk.edge);
      timer_calc.stop();

      //The chunk's edge info is needed to solve the global problem. Collect it.
      job1.top_elev    = dem.topRow     ();
      job1.bot_elev    = dem.bottomRow  ();
      job1.left_elev   = dem.leftColumn ();
      job1.right_elev  = dem.rightColumn();

      job1.top_label   = labels.topRow     ();
      job1.bot_label   = labels.bottomRow  ();
      job1.left_label  = labels.leftColumn ();
      job1.right_label = labels.rightColumn();

      if(chunk.retention=="@offloadall"){
        //Nothing to do: it will all get overwritten
      } else if(chunk.retention=="@retainall"){
        //Nothing to do: it will all be kept because this process won't get
        //called again
      } else {
        timer_io.start();
        dem.saveNative   (chunk.retention+"dem.dat"   );
        labels.saveNative(chunk.retention+"labels.dat");
        timer_io.stop();
      }

      timer_overall.stop();
      std::cerr<<"Node "<<CommRank()<<" finished with Calc="<<timer_calc.accumulated()<<"s. timer_Overall="<<timer_overall.accumulated()<<"s. timer_IO="<<timer_io.accumulated()<<"s. Labels used="<<job1.graph.size()<<std::endl;

      long vmpeak, vmhwm;
      ProcessMemUsage(vmpeak,vmhwm);
      job1.time_info = TimeInfo(timer_calc.accumulated(),timer_overall.accumulated(),timer_io.accumulated(),vmpeak,vmhwm);

      CommSend(&job1,nullptr,0,TAG_DONE_FIRST);
    } else if (the_job==JOB_SECOND){
      Timer timer_calc,timer_io,timer_overall;
      timer_overall.start();
      std::vector<elev_t> graph_elev; //TODO

      CommRecv(&chunk, &graph_elev, 0);

      std::vector<std::map<label_t, elev_t> > graph(2*chunk.width+2*chunk.height); //TODO

      //These use the same logic as the analogous lines above
      if(chunk.retention=="@offloadall"){
        timer_io.start();
        dem = Array2D<elev_t>(chunk.filename, false, chunk.x, chunk.y, chunk.width, chunk.height);
        if(chunk.flip & FLIP_VERT)
          dem.flipVert();
        if(chunk.flip & FLIP_HORZ)
          dem.flipHorz();
        timer_io.stop();
  
        labels = Array2D<label_t>(dem.viewWidth(),dem.viewHeight(),0);
        timer_calc.start();
        Zhou2015Labels(dem,labels,graph,chunk.edge);
        timer_calc.stop();
      } else if(chunk.retention=="@retainall"){
        //Nothing to do: we have it all in memory
      } else {
        timer_io.start();
        dem    = Array2D<elev_t>(chunk.retention+"dem.dat"    ,true); //TODO: There should be an exception if this fails
        labels = Array2D<label_t>(chunk.retention+"labels.dat",true);
        timer_io.stop();
      }

      timer_calc.start();
      for(int y=0;y<dem.viewHeight();y++)
      for(int x=0;x<dem.viewWidth();x++)
        if(labels(x,y)>1 && dem(x,y)<graph_elev.at(labels(x,y)))
          dem(x,y) = graph_elev.at(labels(x,y));
      timer_calc.stop();

      //At this point we're done with the calculation! Boo-yeah!

      timer_io.start();
      if(chunk.flip & FLIP_HORZ)
        dem.flipHorz();
      if(chunk.flip & FLIP_VERT)
        dem.flipVert();
      timer_io.stop();

      timer_io.start();
      dem.saveGDAL(chunk.outputname, chunk.filename, chunk.x, chunk.y);
      timer_io.stop();
      timer_overall.stop();

      std::cerr<<"Node "<<CommRank()<<" finished ("<<chunk.gridx<<","<<chunk.gridy<<") with Calc="<<timer_calc.accumulated()<<"s. timer_Overall="<<timer_overall.accumulated()<<"s. timer_IO="<<timer_io.accumulated()<<"s."<<std::endl;

      long vmpeak, vmhwm;
      ProcessMemUsage(vmpeak,vmhwm);
      TimeInfo temp(timer_calc.accumulated(), timer_overall.accumulated(), timer_io.accumulated(),vmpeak,vmhwm);
      CommSend(&temp, nullptr, 0, TAG_DONE_SECOND);
    }
  }
}




template<class elev_t>
void HandleEdge(
  const std::vector<elev_t>  &elev_a,
  const std::vector<elev_t>  &elev_b,
  const std::vector<label_t> &label_a,
  const std::vector<label_t> &label_b,
  std::vector< std::map<label_t, elev_t> > &mastergraph,
  const label_t label_a_offset,
  const label_t label_b_offset
){
  //Guarantee that all vectors are of the same length
  assert(elev_a.size ()==elev_b.size ());
  assert(label_a.size()==label_b.size());
  assert(elev_a.size ()==label_b.size());

  int len = elev_a.size();

  for(int i=0;i<len;i++){
    auto c_l = label_a[i];
    if(c_l>1) c_l+=label_a_offset;

    for(int ni=i-1;ni<=i+1;ni++){
      if(ni<0 || ni==len)
        continue;
      auto n_l = label_b[ni];
      if(n_l>1) n_l+=label_b_offset;
      //TODO: Does this really matter? We could just ignore these entries
      if(c_l==n_l) //Only happens when labels are both 1
        continue;

      auto elev_over = std::max(elev_a[i],elev_b[ni]);
      if(mastergraph.at(c_l).count(n_l)==0 || elev_over<mastergraph.at(c_l)[n_l]){
        mastergraph[c_l][n_l] = elev_over;
        mastergraph[n_l][c_l] = elev_over;
      }
    }
  }
}

template<class elev_t>
void HandleCorner(
  const elev_t  elev_a,
  const elev_t  elev_b,
  label_t       l_a,
  label_t       l_b,
  std::vector< std::map<label_t, elev_t> > &mastergraph,
  const label_t l_a_offset,
  const label_t l_b_offset
){
  if(l_a>1) l_a += l_a_offset;
  if(l_b>1) l_b += l_b_offset;
  auto elev_over = std::max(elev_a,elev_b);
  if(mastergraph.at(l_a).count(l_b)==0 || elev_over<mastergraph.at(l_a)[l_b]){
    mastergraph[l_a][l_b] = elev_over;
    mastergraph[l_b][l_a] = elev_over;
  }
}







//Producer takes a collection of Jobs and delegates them to Consumers. Once all
//of the jobs have received their initial processing, it uses that information
//to compute the global properties necessary to the solution. Each Job, suitably
//modified, is then redelegated to a Consumer which ultimately finishes the
//processing.
template<class elev_t>
void Producer(std::vector< std::vector< ChunkInfo > > &chunks){
  Timer timer_overall,timer_calc;
  timer_overall.start();
  int active_nodes = 0;

  std::string omsg;

  TimeInfo time_first_total;
  TimeInfo time_second_total;
  int      time_first_count  = 0;
  int      time_second_count = 0;

  std::map<int,ChunkInfo> rank_to_chunk;

  std::vector< std::vector< Job1<elev_t> > > jobs1(chunks.size(), std::vector< Job1<elev_t> >(chunks[0].size()));


  //Loop through all of the jobs, delegating them to Consumers
  active_nodes=0;
  for(size_t y=0;y<chunks.size();y++)
  for(size_t x=0;x<chunks[0].size();x++){
    std::cerr<<"Sending job "<<(y*chunks[0].size()+x+1)<<"/"<<(chunks.size()*chunks[0].size())<<" ("<<(x+1)<<"/"<<chunks[0].size()<<","<<(y+1)<<"/"<<chunks.size()<<")"<<std::endl;
    if(chunks[y][x].nullChunk){
      std::cerr<<"\tNull chunk: skipping."<<std::endl;
      continue;
    }

    //If fewer jobs have been delegated than there are Consumers available,
    //delegate the job to a new Consumer.
    if(active_nodes<CommSize()-1){
      active_nodes++;

      std::cerr<<chunks.at(y).at(x).filename<<std::endl;

      rank_to_chunk[active_nodes] = chunks.at(y).at(x);
      CommSend(&chunks.at(y).at(x),nullptr,active_nodes,JOB_FIRST);

    //Once all of the consumers are active, wait for them to return results. As
    //each Consumer returns a result, pass it the next unfinished Job until
    //there are no jobs left.
    } else {
      //Execute a blocking receive until some consumer finishes its work.
      //Receive that work.
      int source = CommGetSource();
      ChunkInfo received_chunk = rank_to_chunk[source];
      CommRecv(&jobs1.at(received_chunk.gridy).at(received_chunk.gridx), nullptr, source);

      //Delegate new work to that consumer
      rank_to_chunk[source] = chunks.at(y).at(x);
      CommSend(&chunks.at(y).at(x),nullptr,source,JOB_FIRST);
    }
  }

  while(active_nodes>0){
    //Execute a blocking receive until some consumer finishes its work.
    //Receive that work
    int source = CommGetSource();
    ChunkInfo received_chunk = rank_to_chunk[source];
    CommRecv(&jobs1.at(received_chunk.gridy).at(received_chunk.gridx), nullptr, source);

    //Decrement the number of consumers we are waiting on. When this hits 0 all
    //of the jobs have been completed and we can move on
    active_nodes--;
  }


  //TODO: Note here about how the code above code be incorporated

  const int gridwidth  = jobs1.front().size();
  const int gridheight = jobs1.size();

  std::cerr<<"First stage: "<<CommBytesSent()<<"B sent."<<std::endl;
  std::cerr<<"First stage: "<<CommBytesRecv()<<"B received."<<std::endl;
  CommBytesReset();

  //Merge all of the graphs together into one very big graph. Clear information
  //as we go in order to save space, though I am not sure if the map::clear()
  //method is not guaranteed to release space.
  std::cerr<<"Constructing mastergraph..."<<std::endl;
  std::cerr<<"Merging graphs..."<<std::endl;
  timer_calc.start();
  Timer timer_mg_construct;
  timer_mg_construct.start();
  
  //Get a chunk size so we can calculate the max label
  label_t maxlabel = 0;
  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++)
    maxlabel+=jobs1[y][x].graph.size();
  std::cerr<<"Total labels required: "<<maxlabel<<std::endl;

  std::vector< std::map<label_t, elev_t> > mastergraph(maxlabel);

  label_t label_offset = 0;
  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++){
    if( (y*gridwidth+x)%10==0 )
      std::cerr<<"\tcmg: "<<(y*gridwidth+x)<<"/"<<(gridheight*gridwidth)<<std::endl;

    if(chunks[y][x].nullChunk)
      continue;

    Job1<elev_t> &this_job = jobs1[y][x];

    chunks[y][x].label_offset = label_offset;

    for(int l=0;l<(int)this_job.graph.size();l++)
    for(auto const &skey: this_job.graph[l]){
      label_t first_label  = l;
      label_t second_label = skey.first;
      if(first_label >1) first_label +=label_offset;
      if(second_label>1) second_label+=label_offset;
      //We insert both ends of the bidirectional edge because in the watershed
      //labeling process, we only inserted one. We need both here because we
      //don't know which end of the edge we will approach from as we traverse
      //the spillover graph.
      mastergraph.at(first_label)[second_label] = skey.second;
      mastergraph.at(second_label)[first_label] = skey.second;
    }
    chunks[y][x].label_increment = this_job.graph.size();
    label_offset                += this_job.graph.size();
    this_job.graph.clear();
  }

  std::cerr<<"Handling adjacent edges and corners..."<<std::endl;
  for(int y=0;y<gridheight;y++)
  for(int x=0;x<gridwidth;x++){
    if( (y*gridwidth+x)%10==0 )
      std::cerr<<"\tha: "<<(y*gridwidth+x)<<"/"<<(gridheight*gridwidth)<<"\n";

    if(chunks[y][x].nullChunk)
      continue;

    auto &c = jobs1[y][x];

    time_first_total += c.time_info;
    time_first_count++;

    if(y>0            && !chunks[y-1][x].nullChunk)
      HandleEdge(c.top_elev,   jobs1[y-1][x].bot_elev,   c.top_label,   jobs1[y-1][x].bot_label,   mastergraph, chunks[y][x].label_offset, chunks[y-1][x].label_offset);

    if(y<gridheight-1 && !chunks[y+1][x].nullChunk)
      HandleEdge(c.bot_elev,   jobs1[y+1][x].top_elev,   c.bot_label,   jobs1[y+1][x].top_label,   mastergraph, chunks[y][x].label_offset, chunks[y+1][x].label_offset);
    
    if(x>0            && !chunks[y][x-1].nullChunk)
      HandleEdge(c.left_elev,  jobs1[y][x-1].right_elev, c.left_label,  jobs1[y][x-1].right_label, mastergraph, chunks[y][x].label_offset, chunks[y][x-1].label_offset);
    
    if(x<gridwidth-1  && !chunks[y][x+1].nullChunk)
      HandleEdge(c.right_elev, jobs1[y][x+1].left_elev,  c.right_label, jobs1[y][x+1].left_label,  mastergraph, chunks[y][x].label_offset, chunks[y][x+1].label_offset);


    //I wish I had wrote it all in LISP.
    //Top left
    if(y>0 && x>0                      && !chunks[y-1][x-1].nullChunk)   
      HandleCorner(c.top_elev.front(), jobs1[y-1][x-1].bot_elev.back(),  c.top_label.front(), jobs1[y-1][x-1].bot_label.back(),  mastergraph, chunks[y][x].label_offset, chunks[y-1][x-1].label_offset);
    
    //Bottom right
    if(y<gridheight-1 && x<gridwidth-1 && !chunks[y+1][x+1].nullChunk) 
      HandleCorner(c.bot_elev.back(),  jobs1[y+1][x+1].top_elev.front(), c.bot_label.back(),  jobs1[y+1][x+1].top_label.front(), mastergraph, chunks[y][x].label_offset, chunks[y+1][x+1].label_offset);
    
    //Top right
    if(y>0 && x<gridwidth-1            && !chunks[y-1][x+1].nullChunk)            
      HandleCorner(c.top_elev.back(),  jobs1[y-1][x+1].bot_elev.front(), c.top_label.back(),  jobs1[y-1][x+1].bot_label.front(), mastergraph, chunks[y][x].label_offset, chunks[y-1][x+1].label_offset);
    
    //Bottom left
    if(x>0 && y<gridheight-1           && !chunks[y+1][x-1].nullChunk) 
      HandleCorner(c.bot_elev.front(), jobs1[y+1][x-1].top_elev.back(),  c.bot_label.front(), jobs1[y+1][x-1].top_label.back(),  mastergraph, chunks[y][x].label_offset, chunks[y+1][x-1].label_offset);
  }
  timer_mg_construct.stop();

  std::cerr<<"Mastergraph constructed in "<<timer_mg_construct.accumulated()<<"s. "<<std::endl;


  jobs1.clear();
  jobs1.shrink_to_fit();


  std::cerr<<"Performing aggregated priority flood"<<std::endl;
  Timer agg_pflood_timer;
  agg_pflood_timer.start();
  typedef std::pair<elev_t, label_t>  graph_node;
  std::priority_queue<graph_node, std::vector<graph_node>, std::greater<graph_node> > open;
  std::queue<graph_node> pit;
  std::vector<bool>   visited   (maxlabel,false); //TODO
  std::vector<elev_t> graph_elev(maxlabel);       //TODO

  open.emplace(std::numeric_limits<elev_t>::lowest(),1);

  while(open.size()>0 || pit.size()>0){
    graph_node c;
    if(pit.size()>0){
      c = pit.front();
      pit.pop();
    } else {
      c = open.top();
      open.pop();
    }

    auto my_elev       = c.first;
    auto my_vertex_num = c.second;
    if(visited[my_vertex_num])
      continue;

    graph_elev[my_vertex_num] = my_elev;
    visited   [my_vertex_num] = true;

    for(auto &n: mastergraph[my_vertex_num]){
      auto n_vertex_num = n.first;
      auto n_elev       = n.second;
      if(visited[n_vertex_num])
        continue;
      open.emplace(std::max(my_elev,n_elev),n_vertex_num);
      // if(n_elev<=my_elev){
      //   pit.emplace(my_elev,n_vertex_num);
      // } else {
      //   open.emplace(n_elev,n_vertex_num);
      // }
    }
  }
  agg_pflood_timer.stop();
  std::cerr<<"Aggregated priority flood took "<<agg_pflood_timer.accumulated()<<"s."<<std::endl;
  timer_calc.stop();

  std::cerr<<"Sending out final jobs..."<<std::endl;
  //Loop through all of the jobs, delegating them to Consumers
  active_nodes = 0;
  for(size_t y=0;y<chunks.size();y++)
  for(size_t x=0;x<chunks[0].size();x++){
    std::cerr<<"Sending job "<<(y*chunks[0].size()+x+1)<<"/"<<(chunks.size()*chunks[0].size())<<" ("<<(x+1)<<"/"<<chunks[0].size()<<","<<(y+1)<<"/"<<chunks.size()<<")"<<std::endl;
    if(chunks[y][x].nullChunk){
      std::cerr<<"\tNull chunk: skipping."<<std::endl;
      continue;
    }

    timer_calc.start();
    std::vector<elev_t> job2(graph_elev.begin()+chunks[y][x].label_offset,graph_elev.begin()+chunks[y][x].label_offset+chunks[y][x].label_increment);
    timer_calc.stop();

    //If fewer jobs have been delegated than there are Consumers available,
    //delegate the job to a new Consumer.
    if(active_nodes<CommSize()-1){
      // std::cerr<<"Sending init to "<<(active_nodes+1)<<std::endl;
      active_nodes++;
      CommSend(&chunks.at(y).at(x), &job2, active_nodes, JOB_SECOND);

    //Once all of the consumers are active, wait for them to return results. As
    //each Consumer returns a result, pass it the next unfinished Job until
    //there are no jobs left.
    } else {
      //Execute a blocking receive until some consumer finishes its work.
      //Receive that work.
      int source = CommGetSource();

      TimeInfo temp;
      CommRecv(&temp, nullptr, source);
      
      time_second_total += temp;
      time_second_count++;

      //Delegate new work to that consumer
      CommSend(&chunks.at(y).at(x), &job2, source, JOB_SECOND);
    }
  }

  while(active_nodes>0){
    //Execute a blocking receive until some consumer finishes its work.
    //Receive that work
    int source = CommGetSource();

    TimeInfo temp;
    CommRecv(&temp, nullptr, source);

    time_second_total += temp;
    time_second_count++;

    //Decrement the number of consumers we are waiting on. When this hits 0 all
    //of the jobs have been completed and we can move on
    active_nodes--;
  }

  for(int i=1;i<CommSize();i++){
    int temp;
    CommSend(&temp,nullptr,i,SYNC_MSG_KILL);
  }

  timer_overall.stop();

  std::cout<<"TimeInfo: First stage total overall time="<<time_first_total.overall<<std::endl;
  std::cout<<"TimeInfo: First stage total io time="     <<time_first_total.io     <<std::endl;
  std::cout<<"TimeInfo: First stage total calc time="   <<time_first_total.calc   <<std::endl;
  std::cout<<"TimeInfo: First stage block count="       <<time_first_count        <<std::endl;
  std::cout<<"TimeInfo: First stage peak child RSS="    <<time_first_total.vmpeak <<std::endl;
  std::cout<<"TimeInfo: First stage peak child HWM="    <<time_first_total.vmhwm  <<std::endl;

  std::cout<<"Second stage: "<<CommBytesSent()<<"B sent."<<std::endl;
  std::cout<<"Second stage: "<<CommBytesRecv()<<"B received."<<std::endl;

  std::cout<<"TimeInfo: Second stage total overall time="<<time_second_total.overall<<std::endl;
  std::cout<<"TimeInfo: Second stage total IO time="     <<time_second_total.io     <<std::endl;
  std::cout<<"TimeInfo: Second stage total calc time="   <<time_second_total.calc   <<std::endl;
  std::cout<<"TimeInfo: Second stage block count="       <<time_second_count        <<std::endl;
  std::cout<<"TimeInfo: Second stage peak child RSS="    <<time_second_total.vmpeak <<std::endl;
  std::cout<<"TimeInfo: Second stage peak child HWM="    <<time_second_total.vmhwm  <<std::endl;

  std::cout<<"TimeInfo: Producer overall="<<timer_overall.accumulated()<<std::endl;
  std::cout<<"TimeInfo: Producer calc="   <<timer_calc.accumulated()   <<std::endl;

  long vmpeak, vmhwm;
  ProcessMemUsage(vmpeak,vmhwm);
  std::cout<<"TimeInfo: Producer child RSS="    <<vmpeak <<std::endl;
  std::cout<<"TimeInfo: Producer child HWM="    <<vmhwm  <<std::endl;
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
  int chunkid             = 0;
  Timer overall;
  overall.start();

  std::vector< std::vector< ChunkInfo > > chunks;
  std::string  filename;
  GDALDataType file_type; //All chunks must have a common file_type

  if(many_or_one=="many"){
    std::cerr<<"Multi file mode"<<std::endl;

    int32_t gridx        =  0; //Current x coordinate in the chunk grid
    int32_t gridy        = -1; //Current y coordinate in the chunk grid
    int32_t row_width    = -1; //Width of 1st row. All rows must equal this
    int32_t chunk_width  = -1; //Width of 1st chunk. All chunks must equal this
    int32_t chunk_height = -1; //Height of 1st chunk, all chunks must equal this
    long    cell_count   = 0;

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
          CommAbort(-1); //TODO
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
          CommAbort(-1); //TODO: Set error code
        }

        cell_count += this_chunk_width*this_chunk_height;

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
        CommAbort(-1); //TODO: Set error code
      }
    }

    std::cerr<<"Loaded "<<chunks.size()<<" rows of "<<chunks[0].size()<<" columns."<<std::endl;
    std::cerr<<"Total cells to be processed: "<<cell_count<<std::endl;

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
    std::cerr<<"Single file mode"<<std::endl;
    int32_t total_height;
    int32_t total_width;

    filename = input_file;

    auto filepath   = boost::filesystem::path(filename);
    filepath        = filepath.parent_path() / filepath.stem();

    //Get the total dimensions of the input file
    if(getGDALDimensions(filename, total_height, total_width, file_type, NULL)!=0){
      std::cerr<<"Error getting file information from '"<<filename<<"'!"<<std::endl;
      CommAbort(-1); //TODO
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
    std::cerr<<"Total cells to be processed: "<<(total_width*total_height)<<std::endl;

    //Create a grid of jobs
    //TODO: Avoid creating extremely narrow or small strips
    for(int32_t y=0,gridy=0;y<total_height; y+=bheight, gridy++){
      chunks.emplace_back(std::vector<ChunkInfo>());
      for(int32_t x=0,gridx=0;x<total_width;x+=bwidth,  gridx++){
        if(total_height-y<100)
          throw std::logic_error("At least one tile is <100 cells in height. Please change rectangle size to avoid this!");
        if(total_width -x<100)
          throw std::logic_error("At least one tile is <100 cells in width. Please change rectangle size to avoid this!");
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

    if(retention_base=="@retainall" && ((int)(chunks.size()*chunks[0].size()))>=world.size()-1){
      std::cerr<<"This job requires "<<(chunks.size()*chunks[0].size()+1)<<" processes. Only "<<world.size()<<" are available."<<std::endl;
      env.abort(-1);
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
  std::cerr<<"Preparer took "<<overall.accumulated()<<"s."<<std::endl;

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
      if(CommSize()==1) //TODO: Invoke a special "one-process mode?"
        throw std::invalid_argument("Must run program with at least two processes!");
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
      CommBroadcast(&good_to_go,0);
      CommFinalize();
    }

    int good_to_go = 1;
    std::cerr<<"gtg broadcast"<<std::endl;
    CommBroadcast(&good_to_go,0);
    Preparer(many_or_one, retention, input_file, output_prefix, bwidth, bheight, flipH, flipV);

  } else {
    int good_to_go;
    std::cerr<<CommRank()<<" gtg broadcast"<<std::endl;
    CommBroadcast(&good_to_go,0);
    if(!good_to_go){
      CommFinalize();
      return -1;
    }

    GDALDataType file_type;
    std::cerr<<CommRank()<<" filetype broadcast"<<std::endl;
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

  CommFinalize();

  return 0;
}