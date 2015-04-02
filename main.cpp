//Compile with
// mpic++ -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization
// mpirun -n 3 ./a.out ~/projects/watershed/data/beauford03.flt
// TODO: MPI abort
// TODO: See optimization notes at "http://www.boost.org/doc/libs/1_56_0/doc/html/mpi/tutorial.html"
#include "gdal_priv.h"
#include <iostream>
#include <map>
#include <queue>
#include <boost/mpi.hpp>
#include <string>
#include <iomanip>
#include <queue>
#include <limits>
#include <tuple>
#include <ctime> //For generating a basename
#include <fstream>
#include "DisjointDenseIntSet.hpp"
#include "Array2D.hpp"
#include "FillExtent.hpp"
#include "common.hpp"
//#include <unordered_map>
//#define DEBUG 1

#ifdef DEBUG
  #include <fstream>
  #include <string>
#endif

using namespace std;


#define TOP_ELEVATIONS_TAG   1
#define BOT_ELEVATIONS_TAG   2
#define TOP_LABELS_TAG       3
#define BOT_LABELS_TAG       4
#define GRAPH_TAG            5
#define LABEL_OFFSETS        6
#define SYNC_SIG             7
#define FLATS_SIG            8

/*
  For reference, this is the definition of the RasterIO() function
  CPLErr GDALRasterBand::RasterIO( GDALRWFlag eRWFlag,
                                   int nXOff, int nYOff, int nXSize, int nYSize,
                                   void * pData, int nBufXSize, int nBufYSize,
                                   GDALDataType eBufType,
                                   int nPixelSpace,
                                   int nLineSpace )
*/



std::string ElevFilename(const std::string &basename, int node_number){
  return basename+"_elev_"+std::to_string(node_number)+".tmp";
}

std::string FlowdirsFilename(const std::string &basename, int node_number){
  return basename+"_flow_"+std::to_string(node_number)+".tmp";
}





template<class elev_t, GDALDataType gdt_t>
void DoNode(int ram, char *dem_filename){
  typedef std::priority_queue<GridCellZ<elev_t>, std::vector<GridCellZ<elev_t> >, std::greater<GridCellZ<elev_t> > > GridCellZ_pq;
  typedef std::map<label_t, std::map<label_t, elev_t> > Graph;
  boost::mpi::communicator world;

  const int my_node_number        = world.rank()-1;
  const int total_number_of_nodes = world.size()-1; //Don't count master node

  //Synchronize Threads TODO
  std::cerr.flush();
  if(my_node_number>0){
    int msg;
    world.recv(my_node_number,SYNC_SIG,msg);
  }

  std::string basename;
  boost::mpi::broadcast(world,basename,0);

  Array2D<elev_t> elev(dem_filename, my_node_number, total_number_of_nodes);

  //Setup flowdir matrix
  Array2D<flowdir_t> flowdirs; //TODO: Keep a separate line so we can disable this operation if the user wants
  flowdirs = Array2D<flowdir_t>(elev.viewHeight(),elev.width(),FLOWDIR_NO_DATA);

  //////////////////////////////////////////
  //Use the improved priority flood from Barnes 2014 TODO
  GridCellZ_pq open;
  std::queue< GridCellZ<elev_t> > pit;
  Array2D<label_t> labels(elev.viewHeight(), elev.width(), 0);

  //A node value of 1 will be used to indicate areas into which we want flow to
  //ultimately drain. By assumption, this includes any cell on the edge of the
  //DEM. We therefore initialize the edge cells to -1, since a negative label
  //value indicates a label that has been assigned, but whose cell has not yet
  //been visitied.

  std::cerr<<"Adding edges..."<<std::endl;
  //Add the top and bottom row
  for(int x=1;x<elev.width()-1;x++){
    open.emplace(x,0,elev(x,0));
    open.emplace(x,elev.viewHeight()-1,elev(x,elev.viewHeight()-1));
    //My node has the top row of the DEM
    if(my_node_number==0)
      labels(x,0) = -1;
    //My node has the bottom row of the DEM
    else if(my_node_number==total_number_of_nodes-1)
      labels(x,elev.viewHeight()-1) = -1;
  }
  //Add the left and right edge
  for(int y=0;y<elev.viewHeight();y++){
    open.emplace(0,y,elev(0,y));
    open.emplace(elev.width()-1,y,elev(elev.width()-1,y));
    labels(0,y)       = -1;
    labels(elev.width()-1,y) = -1;
  }

  Graph graph;

  //We use an initial current label of 2. This is because a value of 1 is any
  //area of the DEM, universally across all nodes, to which we want flow
  //ultimately to go. Therefore, 2 is the first unoccupied label to use.
  int current_label=2;
  std::cerr<<"Performing Priority-Flood..."<<std::endl;
  while(open.size()>0 || pit.size()>0){
    GridCellZ<elev_t> c;
    if(pit.size()>0){
      c=pit.front();
      pit.pop();
    } else {
      c=open.top();
      open.pop();
    }

    //Since labels are inherited from parent cells we need to be able to process
    //previously labeled cells. But all the edge cells are already in the open
    //queue and may also be added to that queue by a parent cell. So they could
    //be processed twice! The solution is to assign the negative of a label and
    //then make the label positive when we actually process the cell. This way,
    //if we ever see we are processing a cell with a positive label, we will
    //know it was already processed.
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
      if(nx<0 || ny<0 || nx==elev.width() || ny==elev.viewHeight()) continue;
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
          auto elev_over = std::max(elev(nx,ny),elev(c.x,c.y)); //TODO: I think this should always be the neighbour.
          //Haven't seen this watershed before
          if(graph[my_label].count(n_label)==0){
            graph[my_label][n_label] = elev_over;
            graph[n_label][my_label] = elev_over;
          //We've seen this watershed before, so only make a note of the
          //spill-over elevation if it is lower than what we've seen before.
          } else if(elev_over<graph[my_label][n_label]){
            graph[my_label][n_label] = elev_over;
            graph[n_label][my_label] = elev_over;
          }
        }
        continue;
      }

      //The neighbour is not one we've seen before, so mark it as being part of
      //our watershed and add it as an unprocessed item to the queue.
      labels(nx,ny) = -labels(c.x,c.y);


      if(elev(nx,ny)!=elev.noData())        //The neighbour is part of the DEM's data
        flowdirs(nx,ny) = d8_inverse[n]; //and flows into this cell

      //If the neighbour is lower than this cell, elevate it to the level of
      //this cell so that a depression is not formed. The flow directions will
      //be fixed later, after all the depressions have been filled.
      if(elev(nx,ny)<=c.z){
        elev(nx,ny) = c.z;
        pit.emplace(nx,ny,c.z);
      //Otherwise, if the neighbour is higher, do not adjust its elevation.
      } else
        open.emplace(nx,ny,elev(nx,ny));
    }
  }

  //At this point all labels are positive. All cells have passed through the
  //queue and that process ensures that each cell's label is positive.

  #ifdef DEBUG
    std::cerr<<"Accumulation grid"<<std::endl;
    for(int y=0;y<segment_height;y++){
      for(int x=0;x<width;x++)
        std::cerr<<setw(3)<<elev(x,y)<<" ";
      std::cerr<<std::endl;
    }
    std::cerr<<"Labels grid"<<std::endl;
    for(int y=0;y<segment_height;y++){
      for(int x=0;x<width;x++)
        std::cerr<<setw(3)<<labels(x,y)<<" ";
      std::cerr<<std::endl;
    }
  #endif

  if(total_number_of_nodes>1){
    //Synchronize Threads TODO
    std::cerr.flush();
    if(my_node_number<total_number_of_nodes-1)
      world.send(my_node_number+2,SYNC_SIG, 0);

    #ifdef DEBUG
      for(int x=0;x<width;x++){
        assert(labels.front()[x]>0);
        assert(labels.back()[x]>0);
      }
    #endif

    //TODO: Could have master load elevations while nodes are processing. This
    //would speed things up slightly and reduce communication overhead.
    world.send(0, TOP_ELEVATIONS_TAG, elev.topRow     ());
    world.send(0, BOT_ELEVATIONS_TAG, elev.bottomRow  ());
    world.send(0, TOP_LABELS_TAG,     labels.topRow   ());
    world.send(0, BOT_LABELS_TAG,     labels.bottomRow());
    world.send(0, GRAPH_TAG,          graph);

    //Synchronize Threads TODO
    std::cerr.flush();
    if(my_node_number>0){
      int msg;
      world.recv(my_node_number,SYNC_SIG,msg);
    }

    std::cerr<<"=========="<<my_node_number<<std::endl;
    std::cerr<<"Receiving label offsets..."<<std::endl;
    std::map<label_t,elev_t> label_offsets;
    world.recv(0,LABEL_OFFSETS,label_offsets);

    std::map<label_t,FillExtent<elev_t> > flats;
    for(auto const &l: label_offsets){
      #ifdef DEBUG
        assert(l.first>0);
      #endif
      flats.emplace(l.first,l.second);
    }

    #ifdef DEBUG
      for(auto const &lo: label_offsets)
        std::cerr<<lo.first<<"->"<<lo.second<<std::endl;
    #endif

    std::cerr<<"Applying label offsets..."<<std::endl;
    for(int y=0;y<elev.viewHeight();y++)
    for(int x=0;x<elev.width();x++){
      auto my_label = labels(x,y);
      if(elev(x,y)==elev.noData())
        continue;
      auto my_offset = label_offsets[my_label];
      if(elev(x,y)<=my_offset){
        elev(x,y) = my_offset;
        flats.at(my_label).expand(x,y);
      }
    }

    elev.saveNative    (ElevFilename    (basename,my_node_number));
    flowdirs.saveNative(FlowdirsFilename(basename,my_node_number));

    // std::cerr<<"Flats: "<<std::endl;
    // for(auto const &f: flats){
    //   std::cerr<<"("<<f.second.xmin<<","<<f.second.ymin<<")-("<<f.second.xmax<<","<<f.second.ymax<<") "<<(f.second.init?"":"UNUSED")<<std::endl;
    // }

    std::cerr<<my_node_number<<" sending flats."<<std::endl;
    world.send(0,FLATS_SIG,flats);
  } else {
    elev.saveGDAL    (basename+"_elev.tif",  dem_filename);
    flowdirs.saveGDAL(basename+"_flows.tif", dem_filename);
  }

  //Synchronize Threads TODO
  std::cerr.flush();
  if(my_node_number<total_number_of_nodes-1)
    world.send(my_node_number+2,SYNC_SIG, 0);
}











template<class elev_t, GDALDataType gdt_t>
void DoMaster(int ram, char *dem_filename){
  typedef std::map<label_t, std::map<label_t, elev_t> > Graph;
  boost::mpi::communicator world;
  GDALAllRegister();

  const int total_number_of_nodes = world.size()-1; //Don't count master node

  std::string basename = std::to_string(std::time(nullptr));
  boost::mpi::broadcast(world,basename,0);

  if(total_number_of_nodes==1)
    return;

  GDALDataset *fin = (GDALDataset*)GDALOpen(dem_filename, GA_ReadOnly);
  if(fin==NULL){
    cerr<<"Master could not open file: "<<dem_filename<<endl;
    return;
  }

  GDALRasterBand *demband = fin->GetRasterBand(1);
  elev_t no_data          = demband->GetNoDataValue();
  int width               = demband->GetXSize();
  int height              = demband->GetYSize();

  Array2D<elev_t>  elev  (total_number_of_nodes*2,width);
  Array2D<label_t> labels(total_number_of_nodes*2,width);
  std::vector<Graph> graphs(total_number_of_nodes);


  for(int i=1;i<=total_number_of_nodes;i++){
    int n=i-1;
    world.recv(i, TOP_ELEVATIONS_TAG, elev.rowRef(2*n));
    world.recv(i, BOT_ELEVATIONS_TAG, elev.rowRef(2*n+1));
    world.recv(i, TOP_LABELS_TAG, labels.rowRef(2*n));
    world.recv(i, BOT_LABELS_TAG, labels.rowRef(2*n+1));
    world.recv(i, GRAPH_TAG, graphs[n]);
  }

  std::cerr<<"\n\n====================\nMASTER\n===================="<<std::endl;

  #ifdef DEBUG
    std::cerr<<"Elevations"<<std::endl;
    for(size_t y=0;y<elev.size();y++){
      if(y%2==0)
        std::cerr<<"------------"<<std::endl;
      for(auto x=0;x<width;x++)
        std::cerr<<setw(3)<<elev[y][x]<<" ";
      std::cerr<<std::endl;
    }

    std::cerr<<"Labels"<<std::endl;
    for(size_t y=0;y<labels.size();y++){
      if(y%2==0)
        std::cerr<<"------------"<<std::endl;
      for(auto x=0;x<width;x++)
        std::cerr<<setw(3)<<labels[y][x]<<" ";
      std::cerr<<std::endl;
    }
  #endif

  //Merge graphs
  std::cerr<<"Merging graphs"<<std::endl;
  Graph mastergraph;
  std::map<int, label_t> strip_to_max_label;
  std::map<label_t, int> label_to_strip;
  label_t maxlabel = 0;
  for(size_t i=0;i<graphs.size();i++){
    label_t newmaxlabel   = 0;
    strip_to_max_label[i] = maxlabel;
    //Recall that graphs have the form "std::map<label_t, std::map<label_t, elev_t> >"
    for(auto const &fkey: graphs[i])
    for(auto const &skey: fkey.second){
      auto flabel = fkey.first;
      auto slabel = skey.first;
      assert(flabel>0); //TODO
      assert(slabel>0); //TODO
      if(flabel!=1)
        flabel += maxlabel;
      if(slabel!=1)
        slabel += maxlabel;
      label_to_strip[flabel]      = i;
      label_to_strip[slabel]      = i;
      mastergraph[flabel][slabel] = skey.second;
      newmaxlabel                 = std::max(flabel,newmaxlabel);
      newmaxlabel                 = std::max(slabel,newmaxlabel);
    }
    for(auto &s: labels.rowRef(2*i))
      if(s!=1)
        s += maxlabel;
    for(auto &s: labels.rowRef(2*i+1))
      if(s!=1)
        s += maxlabel;
    maxlabel = newmaxlabel; //TODO: Should I be adding 1 to newmaxlabel?
  }

  graphs.clear();
  graphs.shrink_to_fit();

  #ifdef DEBUG
    std::cerr<<"Merged labels"<<std::endl;
    for(size_t y=0;y<labels.size();y++){
      if(y%2==0)
        std::cerr<<"------------"<<std::endl;
      for(int x=0;x<width;x++)
        std::cerr<<setw(3)<<labels[y][x]<<" ";
      std::cerr<<std::endl;
    }
  #endif

  std::cerr<<"Joining graphs"<<std::endl;
  //Consider the interior adjoining edge of each pair of strips
  for(int y=1;y<2*total_number_of_nodes-1;y+=2){
    for(int x=0;x<width;x++){
      if(elev(x,y)==no_data)
        continue;
      auto my_label = labels(x,y);
      //Consider the situation from the upper strip looking towards the lower
      //strip. We will implicitly handle the case in which the lower strip looks
      //towards the upper strip.
      for(int n=6;n<=8;n++){
        auto nx = x+dx[n];
        auto ny = y+dy[n];
        #ifdef DEBUG
          std::cerr<<"Considering ("<<x<<","<<y<<") with n="<<n<<" pointing to ("<<nx<<","<<ny<<")"<<std::endl;
        #endif
        if(nx<0 || nx==width || elev(nx,ny)==no_data)
          continue;
        auto other_label = labels(nx,ny);
        if(my_label==other_label) //Only when both labels=1
          continue;

          auto elev_over = std::max(elev(nx,ny),elev(x,y));
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

  elev.clear();

  DisjointDenseIntSet unified_labels(maxlabel);

  #ifdef DEBUG
    std::cerr<<"Master Graph"<<std::endl;
    for(auto const &vertex: mastergraph)
    for(auto const &edges: vertex.second)
      std::cerr<<setw(3)<<vertex.first<<" -> "<<setw(3)<<edges.first<<" = "<<setw(3)<<edges.second<<std::endl;
  #endif


  std::cerr<<"Performing aggregated priority flood"<<std::endl;
  //Keep track of elevation of the current node, the current node's label, and
  //the parent label from which the node was explored
  typedef std::tuple<elev_t, int, int> GraphNode;
  std::priority_queue<GraphNode, std::vector<GraphNode>, std::greater<GraphNode> > open;
  std::map<label_t,bool>              visited;
  std::map<label_t,elev_t>            graph_elev;

  open.emplace(std::numeric_limits<elev_t>::min(),1,-1);

  //Consider the labels as graph unto themselves. What elevation should each
  //label have in order to form a meta-graph without depressions?
  while(open.size()>0){
    GraphNode c=open.top();
    open.pop();

    auto my_elev   = std::get<0>(c);
    auto my_label  = std::get<1>(c);
    auto my_parent = std::get<2>(c);

    #ifdef DEBUG
      std::cerr<<"Popped "<<my_label<<std::endl;
    #endif
    if(visited[my_label])
      continue;


    graph_elev[my_label] = my_elev;
    visited[my_label]    = true;

    //If I have a parent, and I am now the same elevation as my parent, then my
    //parent and I are part of a common depression and should be filled together
    if(my_parent!=-1 && my_elev==graph_elev[my_parent])
      unified_labels.unionSet(my_label,my_parent);

    for(auto const &n: mastergraph[my_label]){
      auto n_label = n.first;
      auto n_elev  = n.second;
      if(visited.count(n_label))
        continue;
      #ifdef DEBUG
        std::cerr<<"Proposing going to "<<n_label<<" with "<<std::max(n_elev,my_elev)<<std::endl;
      #endif
      open.emplace(std::max(n_elev,my_elev),n_label,my_label);
    }
  }

  #ifdef DEBUG
    std::cerr<<"Graph elevations"<<std::endl;
    for(auto const &ge: graph_elev)
      std::cerr<<setw(3)<<ge.first<<" = "<<setw(3)<<ge.second<<std::endl;
  #endif

  std::vector< std::map<label_t,elev_t> > strip_label_elevations(total_number_of_nodes);
  for(auto const &ge: graph_elev){
    auto vertex_num   = ge.first;
    if(vertex_num==1) //1 indicates the edge of the DEM, which is filled with
      continue;       //no_data values or other things we don't want to touch
    auto elevation    = ge.second;
    auto strip_num    = label_to_strip[vertex_num];
    auto label_offset = strip_to_max_label[strip_num];
    vertex_num       -= label_offset;
    assert(vertex_num>0); //TODO
    strip_label_elevations[strip_num][vertex_num] = elevation;
  }

  std::cerr.flush();
  boost::mpi::request reqs[total_number_of_nodes];
  for(int i=1;i<=total_number_of_nodes;i++){
    int n   = i-1;
    reqs[n] = world.isend(i, LABEL_OFFSETS, strip_label_elevations[n]); //TODO
    //world.send(i, LABEL_OFFSETS, strip_label_elevations[n]); //TODO
  }
  boost::mpi::wait_all(reqs, reqs + total_number_of_nodes);




  std::cerr<<"Gathering flats..."<<std::endl;
  std::cerr<<"Maximum element: "<<unified_labels.maxElement()<<std::endl;
  std::map<label_t, FillExtent<elev_t> > flats;
  //For each flat, store the node and that node's label for its part of the flat
  for(int i=1;i<=total_number_of_nodes;i++){
    std::map<label_t,FillExtent<elev_t> > tempflats;
    std::cerr<<"Gathering from "<<i<<std::endl;
    world.recv(i,FLATS_SIG,tempflats);
    int segment_first_line = (height/total_number_of_nodes)*(i-1);
    for(auto &f: tempflats){
      auto my_set    = unified_labels.findSet(f.first+strip_to_max_label[i-1]);
      f.second.ymin += segment_first_line;
      f.second.ymax += segment_first_line;
      flats[my_set].expand(f.second);
    }
  }
  std::cerr<<"Gathered flats."<<std::endl;

  double geotrans[6];
  // In case of north up images, the GT(2) and GT(4) coefficients are zero, and
  // the GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3))
  // position is the top left corner of the top left pixel of the raster.
  fin->GetGeoTransform(geotrans);

  std::cerr<<"Saving flats."<<std::endl;
  std::ofstream flatpolys(basename+"_flat_polys.tmp");
  for(auto &f: flats){
    int xmin = geotrans[0]+geotrans[1]*f.second.xmin;
    int ymax = geotrans[3]+geotrans[5]*f.second.ymin;
    int xmax = geotrans[0]+geotrans[1]*f.second.xmax;
    int ymin = geotrans[3]+geotrans[5]*f.second.ymax;
    flatpolys<<"POLYGON(("<<xmin<<" "<<ymin<<","
                         <<xmin<<" "<<ymax<<","
                         <<xmax<<" "<<ymax<<","
                         <<xmax<<" "<<ymin<<","
                         <<xmin<<" "<<ymin<<"))"<<std::endl;
  }
  flatpolys.close();

  std::ofstream flatsout(basename+"_flats.tmp");
  flatsout<<flats.size()<<std::endl;
  for(auto &f: flats){
    assert(f.second.ymin<=f.second.ymax);
    flatsout<<f.second.elevation<<" "
        <<f.second.xmin<<" "
        <<f.second.ymin<<" "
        <<f.second.xmax<<" "
        <<f.second.ymax<<std::endl;
  }
  flatsout.close();

  std::cerr<<"Writing out elevation data."<<std::endl;
  GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
  if(poDriver==NULL){
    std::cerr<<"Could not open GDAL driver."<<std::endl;
    return;
  }

  std::cerr<<"Writing out elevations."<<std::endl;
  {
    std::vector<std::string> saved_elevs;
    for(int i=1;i<=total_number_of_nodes;i++)
      saved_elevs.push_back(ElevFilename(basename,i-1));

    combineSavedNativeIntoGDAL<elev_t>(basename+"_elev.tif",dem_filename,saved_elevs, no_data, nullptr);
  }

  std::cerr<<"Writing out flow directions."<<std::endl;
  {
    std::vector<std::string> saved_flows;
    for(int i=1;i<=total_number_of_nodes;i++)
      saved_flows.push_back(FlowdirsFilename(basename,i-1));

    combineSavedNativeIntoGDAL<flowdir_t>(basename+"_flows.tif",dem_filename,saved_flows,FLOWDIR_NO_DATA,[](flowdir_t fd){return d8_arcgis[fd];});
  }

  std::cerr<<"Combining flow outputs in temporary file for flat resolution."<<std::endl;
  // Array2D<elev_t> elev_combined (ElevFilename    (basename,0));
  for(int i=2;i<=total_number_of_nodes;i++){
    int n = i-1;
    stackSavedNativeStrips<flowdir_t>(FlowdirsFilename(basename,0),FlowdirsFilename(basename,n));
  }
}








template<class elev_t, GDALDataType gdt_t>
void DoRouter(int ram, char *dem_filename){
  boost::mpi::communicator world;
  if(world.rank()==0)
    DoMaster<elev_t,gdt_t>(ram, dem_filename);
  else
    DoNode<elev_t,gdt_t>(ram, dem_filename);
}

int main(int argc, char **argv){
  boost::mpi::environment env;
  boost::mpi::communicator world;

  if(argc!=3){
    std::cerr<<"Syntax: "<<argv[0]<<" <RAM> <DEM>"<<std::endl;
    std::cerr<<"RAM is an integer value estimating the number of MEGABYTES of available RAM."<<std::endl;
    std::cerr<<"0 implies UNLIMITED RAM. All nodes will retain their data in RAM at all times."<<std::endl;
    std::cerr<<"Anything greater than 0 implies LIMITED RAM. Sequential processing with temporary files will be used limit RAM usage."<<std::endl;
    return -1;
  }

  int ram = std::stoi(argv[1]);

  switch(peekGDALType(argv[2])){
    case GDT_Int16:
      DoRouter<int16_t, GDT_Int16>(ram,argv[2]);
      break;
    case GDT_Float32:
      DoRouter<float, GDT_Float32>(ram,argv[2]);
      break;
    default:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[2]))<<std::endl;
      break;
  }

  return 0;
}