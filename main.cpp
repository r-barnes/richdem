//Compile with
// mpic++ -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization
// mpirun -n 3 ./a.out ~/projects/watershed/data/beauford03.flt
// TODO: MPI abort
// TODO: See optimization notes at "http://www.boost.org/doc/libs/1_56_0/doc/html/mpi/tutorial.html"
// For memory usage see: http://info.prelert.com/blog/stl-container-memory-usage
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

template<class T>
void build2DFromInterleavedVectors(Array2D<T> &combined, std::vector< std::vector<T> > &a, std::vector< std::vector<T> > &b){
  assert(a.size()==b.size());
  for(unsigned int i=0;i<a.size();i++){
    combined.emplaceRow(a[i]);
    combined.emplaceRow(b[i]);
  }
  combined.redimension();
}



template<class elev_t, GDALDataType gdt_t>
void FillDepressions(int ram, char *dem_filename){
  typedef std::priority_queue<GridCellZ<elev_t>, std::vector<GridCellZ<elev_t> >, std::greater<GridCellZ<elev_t> > > GridCellZ_pq;
  typedef std::map<label_t, std::map<label_t, elev_t> > Graph;
  boost::mpi::communicator world;

  int dummy_sync_sig_msg = 0;

  std::string basename;
  if(world.rank()==0)
     basename = std::to_string(std::time(nullptr));
  boost::mpi::broadcast(world,basename,0);


  int    dem_height, dem_width;
  elev_t dem_no_data;
  double geotrans[6];
  if(world.rank()==0)
    getGDALHeader(dem_filename, dem_height, dem_width, dem_no_data, geotrans);


  //Synchronize Threads TODO
  std::cerr.flush();
  if(world.rank()>0)
    world.recv(world.rank()-1,SYNC_SIG,dummy_sync_sig_msg);


  Array2D<elev_t> my_elev(dem_filename, world.rank(), world.size());

  std::map<label_t,elev_t> my_label_offsets;

  DisjointDenseIntSet unified_labels;
  std::map<int, label_t> strip_to_max_label;

  //Setup flowdir matrix
  Array2D<flowdir_t> my_flowdirs; //TODO: Keep a separate line so we can disable this operation if the user wants
  my_flowdirs = Array2D<flowdir_t>(my_elev.viewHeight(),my_elev.width(),FLOWDIR_NO_DATA);
  my_flowdirs.setNoData(FLOWDIR_NO_DATA);

  //////////////////////////////////////////
  //Use the improved priority flood from Barnes 2014 TODO
  GridCellZ_pq my_open;
  std::queue< GridCellZ<elev_t> > my_pit;
  Array2D<label_t> my_labels(my_elev.viewHeight(), my_elev.width(), 0);

  //A node value of 1 will be used to indicate areas into which we want flow to
  //ultimately drain. By assumption, this includes any cell on the edge of the
  //DEM. We therefore initialize the edge cells to -1, since a negative label
  //value indicates a label that has been assigned, but whose cell has not yet
  //been visitied.

  std::cerr<<"Adding edges..."<<std::endl;
  //Add the top and bottom row
  for(int x=1;x<my_elev.width()-1;x++){
    my_open.emplace(x,0,my_elev(x,0));
    my_open.emplace(x,my_elev.viewHeight()-1,my_elev(x,my_elev.viewHeight()-1));
    //My node has the top row of the DEM
    if(world.rank()==0)
      my_labels(x,0) = -1;
    //My node has the bottom row of the DEM
    else if(world.rank()==world.size()-1)
      my_labels(x,my_elev.viewHeight()-1) = -1;
  }
  //Add the left and right edge
  for(int y=0;y<my_elev.viewHeight();y++){
    my_open.emplace(0,y,my_elev(0,y));
    my_open.emplace(my_elev.width()-1,y,my_elev(my_elev.width()-1,y));
    my_labels(0,y)       = -1;
    my_labels(my_elev.width()-1,y) = -1;
  }

  Graph my_graph;

  //We use an initial current label of 2. This is because a value of 1 is any
  //area of the DEM, universally across all nodes, to which we want flow
  //ultimately to go. Therefore, 2 is the first unoccupied label to use.
  int current_label=2;
  std::cerr<<"Performing Priority-Flood..."<<std::endl;
  while(my_open.size()>0 || my_pit.size()>0){
    GridCellZ<elev_t> c;
    if(my_pit.size()>0){
      c=my_pit.front();
      my_pit.pop();
    } else {
      c=my_open.top();
      my_open.pop();
    }

    //Since labels are inherited from parent cells we need to be able to process
    //previously labeled cells. But all the edge cells are already in the
    //my_open queue and may also be added to that queue by a parent cell. So
    //they could be processed twice! The solution is to assign the negative of a
    //label and then make the label positive when we actually process the cell.
    //This way, if we ever see we are processing a cell with a positive label,
    //we will know it was already processed.
    if(my_labels(c.x,c.y)>0)             //Positive label. Cell already processed.
      continue;
    else if(my_labels(c.x,c.y)==0)           //Label=0. Cell has never been seen.
      my_labels(c.x,c.y) = current_label++;  //current_label is the next new label.
    else if(my_labels(c.x,c.y)<0)            //Cell's parent added it.
      my_labels(c.x,c.y) = -my_labels(c.x,c.y); //Mark cell as visited

    //At this point the cell's label is guaranteed to be positive and in the
    //range [1,MAX_INTEGER] (unless we overflow).
    auto my_label = my_labels(c.x,c.y);

    for(int n=1;n<=8;n++){
      auto nx = c.x+dx[n]; //Neighbour's x-coordinate
      auto ny = c.y+dy[n]; //Neighbour's y-coordinate

      //Check to see if the neighbour coordinates are valid
      if(nx<0 || ny<0 || nx==my_elev.width() || ny==my_elev.viewHeight()) continue;
      auto n_label = std::abs(my_labels(nx,ny)); //Neighbour's label
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
          auto elev_over = std::max(my_elev(nx,ny),my_elev(c.x,c.y)); //TODO: I think this should always be the neighbour.
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


      if(my_elev(nx,ny)!=my_elev.noData())        //The neighbour is part of the DEM's data
        my_flowdirs(nx,ny) = d8_inverse[n]; //and flows into this cell

      //If the neighbour is lower than this cell, elevate it to the level of
      //this cell so that a depression is not formed. The flow directions will
      //be fixed later, after all the depressions have been filled.
      if(my_elev(nx,ny)<=c.z){
        my_elev(nx,ny) = c.z;
        my_pit.emplace(nx,ny,c.z);
      //Otherwise, if the neighbour is higher, do not adjust its elevation.
      } else
        my_open.emplace(nx,ny,my_elev(nx,ny));
    }
  }

  //At this point all my_labels are positive. All cells have passed through the
  //queue and that process ensures that each cell's label is positive.

  if(world.size()==1){
    my_elev.saveGDAL    (basename+"_elev.tif",  dem_filename);
    my_flowdirs.saveGDAL(basename+"_flows.tif", dem_filename);
    return;
  }

  //Synchronize Threads TODO
  std::cerr.flush();
  if(world.rank()<world.size()-1)
    world.send(world.rank()+1,SYNC_SIG,0);



  //TODO: Could have master load elevations while nodes are processing. This
  //would speed things up slightly and reduce communication overhead.

  if(world.rank()>0){
    boost::mpi::gather(world, my_elev.topRow     (), 0);
    boost::mpi::gather(world, my_elev.bottomRow  (), 0);
    boost::mpi::gather(world, my_labels.topRow   (), 0);
    boost::mpi::gather(world, my_labels.bottomRow(), 0);
    boost::mpi::gather(world, my_graph,              0);
  }

  if(world.rank()==0){
    std::vector< std::vector<elev_t > > agg_elev_tops,  agg_elev_bottoms;
    std::vector< std::vector<label_t> > agg_label_tops, agg_label_bottoms;
    Array2D<elev_t>    agg_elev;
    agg_elev.setNoData(dem_no_data);
    Array2D<label_t>   agg_label;
    std::vector<Graph> graphs;
    boost::mpi::gather(world, my_elev.topRow     (), agg_elev_tops,     0);
    boost::mpi::gather(world, my_elev.bottomRow  (), agg_elev_bottoms,  0);
    boost::mpi::gather(world, my_labels.topRow   (), agg_label_tops,    0);
    boost::mpi::gather(world, my_labels.bottomRow(), agg_label_bottoms, 0);
    boost::mpi::gather(world, my_graph,              graphs,            0);
    build2DFromInterleavedVectors(agg_elev,  agg_elev_tops,  agg_elev_bottoms);
    build2DFromInterleavedVectors(agg_label, agg_label_tops, agg_label_bottoms);


    std::cerr<<"\n\n====================\nMASTER\n===================="<<std::endl;

    std::cerr<<"Agg_elev height: "<<agg_elev.height()<<", "<<agg_elev.width()<<std::endl;
    //Merge graphs
    std::cerr<<"Merging graphs"<<std::endl;
    Graph mastergraph;
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
      for(auto &s: agg_label.rowRef(2*i))
        if(s!=1)
          s += maxlabel;
      for(auto &s: agg_label.rowRef(2*i+1))
        if(s!=1)
          s += maxlabel;
      maxlabel = newmaxlabel; //TODO: Should I be adding 1 to newmaxlabel?
    }

    graphs.clear();
    graphs.shrink_to_fit();

    std::cerr<<"Joining graphs"<<std::endl;
    //Consider the interior adjoining edge of each pair of strips
    for(int y=1;y<agg_elev.height()-1;y+=2){
      std::cerr<<"y: "<<y<<std::endl;
      for(int x=0;x<agg_elev.width();x++){
        if(agg_elev(x,y)==agg_elev.noData())
          continue;
        auto my_label = agg_label(x,y);
        //Consider the situation from the upper strip looking towards the lower
        //strip. We will implicitly handle the case in which the lower strip looks
        //towards the upper strip.
        for(int n=6;n<=8;n++){
          auto nx = x+dx[n];
          auto ny = y+dy[n];
          if(nx<0 || nx==agg_elev.width() || agg_elev(nx,ny)==agg_elev.noData())
            continue;
          auto other_label = agg_label(nx,ny);
          if(my_label==other_label) //Only when both agg_label=1
            continue;

            auto elev_over = std::max(agg_elev(nx,ny),agg_elev(x,y));
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

    agg_elev.clear();

    std::cerr<<"Making labels: "<<maxlabel<<std::endl;
    for(int i=0;i<=maxlabel;i++)
      unified_labels.makeSet(maxlabel);

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

      auto c_elev   = std::get<0>(c);
      auto c_label  = std::get<1>(c);
      auto c_parent = std::get<2>(c);

      if(visited[c_label])
        continue;


      graph_elev[c_label] = c_elev;
      visited[c_label]    = true;

      //If I have a parent, and I am now the same elevation as my parent, then my
      //parent and I are part of a common depression and should be filled together
      if(c_parent!=-1 && c_elev==graph_elev[c_parent])
        unified_labels.unionSet(c_label,c_parent);

      for(auto const &n: mastergraph[c_label]){
        auto n_label = n.first;
        auto n_elev  = n.second;
        if(visited.count(n_label))
          continue;
        open.emplace(std::max(n_elev,c_elev),n_label,c_label);
      }
    }

    std::vector< std::map<label_t,elev_t> > strip_label_elevations(world.size());
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

    boost::mpi::scatter(world, strip_label_elevations, my_label_offsets, 0);
  }

  if(world.rank()>0)
    boost::mpi::scatter(world, my_label_offsets,0);

  //Synchronize Threads TODO
  std::cerr.flush();
  if(world.rank()>0)
    world.recv(world.rank()-1,SYNC_SIG,dummy_sync_sig_msg);

  std::cerr<<"=========="<<world.rank()<<std::endl;
  std::cerr<<"Receiving label offsets..."<<std::endl;

  std::map<label_t,FillExtent<elev_t> > my_flats;
  for(auto const &l: my_label_offsets)
    my_flats.emplace(l.first,l.second);

  std::cerr<<"Applying label offsets..."<<std::endl;
  for(int y=0;y<my_elev.viewHeight();y++)
  for(int x=0;x<my_elev.width();x++){
    auto my_label = my_labels(x,y);
    if(my_elev(x,y)==my_elev.noData())
      continue;
    auto my_offset = my_label_offsets[my_label];
    if(my_elev(x,y)<=my_offset){
      my_elev(x,y) = my_offset;
      my_flats.at(my_label).expand(x,y);
    }
  }

  my_elev.saveNative    (ElevFilename    (basename,world.rank()));
  my_flowdirs.saveNative(FlowdirsFilename(basename,world.rank()));

  //Synchronize Threads TODO
  std::cerr.flush();
  if(world.rank()<world.size()-1)
    world.send(world.rank()+1,SYNC_SIG, dummy_sync_sig_msg);



  if(world.rank()>0)
    world.send(0,FLATS_SIG,my_flats);

  if(world.rank()==0){
    std::cerr<<"Gathering flats..."<<std::endl;
    std::cerr<<"Maximum element: "<<unified_labels.maxElement()<<std::endl;
    std::map<label_t, FillExtent<elev_t> > agg_flats;
    //For each flat, store the node and that node's label for its part of the flat
    for(int i=0;i<world.size();i++){
      std::cerr<<"Gathering from "<<i<<std::endl;
      if(i>0)
        world.recv(i,FLATS_SIG,my_flats);
      int segment_first_line = (dem_height/world.size())*i;
      for(auto &f: my_flats){
        std::cerr<<"###"<<i<<" "<<f.first<<" "<<strip_to_max_label[i]<<std::endl;
        auto my_set    = unified_labels.findSet(f.first+strip_to_max_label[i]);
        f.second.ymin += segment_first_line;
        f.second.ymax += segment_first_line;
        agg_flats[my_set].expand(f.second);
      }
    }
    std::cerr<<"Gathered flats."<<std::endl;

    std::cerr<<"Saving flats."<<std::endl;
    std::ofstream flatpolys(basename+"_flat_polys.tmp");
    for(auto &f: agg_flats){
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
    flatsout<<agg_flats.size()<<std::endl;
    for(auto &f: agg_flats){
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
      for(int i=0;i<world.size();i++)
        saved_elevs.push_back(ElevFilename(basename,i));

      combineSavedNativeIntoGDAL<elev_t>(basename+"_elev.tif",dem_filename,saved_elevs,dem_no_data, nullptr);
    }

    std::cerr<<"Writing out flow directions."<<std::endl;
    {
      std::vector<std::string> saved_flows;
      for(int i=0;i<world.size();i++)
        saved_flows.push_back(FlowdirsFilename(basename,i));

      //combineSavedNativeIntoGDAL<flowdir_t>(basename+"_flows.tif",dem_filename,saved_flows,FLOWDIR_NO_DATA,[](flowdir_t fd){if(fd==FLOWDIR_NO_DATA) return FLOWDIR_NO_DATA; else return d8_arcgis[fd];});
      combineSavedNativeIntoGDAL<flowdir_t>(basename+"_flows.tif",dem_filename,saved_flows,FLOWDIR_NO_DATA, nullptr);
    }

    std::cerr<<"Combining flow outputs in temporary file for flat resolution."<<std::endl;
    // Array2D<elev_t> elev_combined (ElevFilename    (basename,0));
    for(int i=1;i<world.size();i++)
      stackSavedNativeStrips<flowdir_t>(FlowdirsFilename(basename,0),FlowdirsFilename(basename,i));
  }
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
      FillDepressions<int16_t, GDT_Int16>(ram,argv[2]);
      break;
    case GDT_Float32:
      FillDepressions<float, GDT_Float32>(ram,argv[2]);
      break;
    default:
      if(world.rank()==0)
        std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(peekGDALType(argv[2]))<<std::endl;
      break;
  }

  return 0;
}