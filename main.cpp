//Compile with
// mpic++ -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization
// mpirun -n 3 ./a.out ~/projects/watershed/data/beauford03.flt
#include "gdal_priv.h"
#include <iostream>
#include <map>
#include <queue>
#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp>
#include <string>
#include <iomanip>
#include <queue>
#include <limits>
//#include <unordered_map>
//#define DEBUG 1

namespace mpi = boost::mpi;

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

/*
  For reference, this is the definition of the RasterIO() function
  CPLErr GDALRasterBand::RasterIO( GDALRWFlag eRWFlag,
                                   int nXOff, int nYOff, int nXSize, int nYSize,
                                   void * pData, int nBufXSize, int nBufYSize,
                                   GDALDataType eBufType,
                                   int nPixelSpace,
                                   int nLineSpace )
*/

//D8 Directions
///x offsets of D8 neighbours, from a central cell
const int dx[9]={0,-1,-1,0,1,1,1,0,-1};  //TODO: These should be merged with my new dinf_d8 to reflect a uniform and intelligent directional system
///y offsets of D8 neighbours, from a central cell
const int dy[9]={0,0,-1,-1,-1,0,1,1,1};
//234
//105
//876

/// Stores the (x,y) coordinates of a grid cell
class GridCell {
 public:
  int x; ///< Grid cell's x-coordinate
  int y; ///< Grid cell's y-coordinate
  /// Initiate the grid cell without coordinates; should generally be avoided
  GridCell(){}
  /// Initiate the grid cell to the coordinates (x0,y0)
  GridCell(int x, int y):x(x),y(y){}
};

template<class elev_t>
class GridCellZ : public GridCell {
 public:
  elev_t z;         ///< Grid cell's z-coordinate
  GridCellZ(int x, int y, elev_t z): GridCell(x,y), z(z) {}
  GridCellZ(){}
  bool operator> (const GridCellZ& a) const { return z>a.z; }
};

typedef int label_t;
typedef std::vector<label_t> LabelRow;
typedef std::vector<LabelRow> Labels;




template<class elev_t, GDALDataType gdt_t>
void DoNode(char *dem_filename){
  typedef std::priority_queue<GridCellZ<elev_t>, std::vector<GridCellZ<elev_t> >, std::greater<GridCellZ<elev_t> > > GridCellZ_pq;
  typedef std::vector<elev_t>       ElevationRow;
  typedef std::vector<ElevationRow> Elevations;
  typedef std::map<label_t, std::map<label_t, elev_t> > Graph;
  mpi::communicator world;

  const int my_node_number        = world.rank()-1;
  const int total_number_of_nodes = world.size()-1; //Don't count master node

  //Synchronize Threads TODO
  std::cerr.flush();
  if(my_node_number>0){
    int msg;
    world.recv(my_node_number,SYNC_SIG,msg);
  }

  GDALAllRegister();

  GDALDataset *fin = (GDALDataset*)GDALOpen(dem_filename, GA_ReadOnly);
  if(fin==NULL){
    cerr<<"Node "<<my_node_number<<" could not open file: "<<dem_filename<<endl;
    return;
  }

  GDALRasterBand *demband = fin->GetRasterBand(1);
  int width               = demband->GetXSize();
  int height              = demband->GetYSize();
  elev_t no_data          = demband->GetNoDataValue();

  int segment_first_line = (height/total_number_of_nodes)*my_node_number;
  int segment_last_line  = (height/total_number_of_nodes)*(my_node_number+1);
  if(my_node_number==total_number_of_nodes-1)
    segment_last_line = height;
  int segment_height = segment_last_line - segment_first_line;

  //Read in DEM data
  Elevations elev(segment_height, ElevationRow(width));
  for(int y=segment_first_line;y<segment_last_line;y++){
    ElevationRow temp(width);
    demband -> RasterIO( GF_Read, 0, y, width, 1, temp.data(), width, 1, gdt_t, 0, 0 );
    elev[y-segment_first_line]=temp;
  }
  //TODO: Why? "y-segment_first_line" in the above?

  //////////////////////////////////////////
  //Use the improved priority flood from Barnes 2014 TODO
  GridCellZ_pq open;
  std::queue< GridCellZ<elev_t> > pit;
  Labels labels(segment_height, LabelRow(width,0));

  std::cerr<<"Adding edges..."<<std::endl;
  for(int x=1;x<width-1;x++){
    open.emplace(x,0,elev[0][x]);
    open.emplace(x,segment_height-1,elev[segment_height-1][x]);
    if(my_node_number==0)
      labels[0][x] = -1;
    else if(my_node_number==total_number_of_nodes-1)
      labels[segment_height-1][x] = -1;
  }
  for(int y=0;y<segment_height;y++){
    open.emplace(0,y,elev[y][0]);
    open.emplace(width-1,y,elev[y][width-1]);
    labels[y][0]       = -1;
    labels[y][width-1] = -1;
  }

  Graph graph;

  int current_label=2;  //TODO: Thought this was more clear than zero in the results.
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
    if(labels[c.y][c.x]>0)
      continue;
    else if(labels[c.y][c.x]==0)
      labels[c.y][c.x] = current_label++;   //Use current_label and increment it
    else if(labels[c.y][c.x]<0)
      labels[c.y][c.x] = -labels[c.y][c.x]; //Mark cell as visited
    auto my_label = labels[c.y][c.x];

    for(int n=1;n<=8;n++){
      auto nx = c.x+dx[n];
      auto ny = c.y+dy[n];
      if(nx<0 || ny<0 || nx==width || ny==segment_height) continue;
      auto other_label = std::abs(labels[ny][nx]);
      if(other_label!=0){
        if(other_label!=my_label){
          auto elev_over = std::max(elev[ny][nx],elev[c.y][c.x]);
          if(graph[my_label].count(other_label)==0){
            graph[my_label][other_label] = elev_over;
            graph[other_label][my_label] = elev_over;
          } else if(elev_over<graph[my_label][other_label]){
            graph[my_label][other_label] = elev_over;
            graph[other_label][my_label] = elev_over;
          }
        }
        continue;
      }

      labels[ny][nx] = -labels[c.y][c.x];

      if(elev[ny][nx]<=c.z){
        elev[ny][nx] = c.z;
        pit.emplace(nx,ny,c.z);
      } else
        open.emplace(nx,ny,elev[ny][nx]);
    }
  }

  #ifdef DEBUG
    std::cerr<<"Accumulation grid"<<std::endl;
    for(int y=0;y<segment_height;y++){
      for(int x=0;x<width;x++)
        std::cerr<<setw(3)<<elev[y][x]<<" ";
      std::cerr<<std::endl;
    }
    std::cerr<<"Labels grid"<<std::endl;
    for(int y=0;y<segment_height;y++){
      for(int x=0;x<width;x++)
        std::cerr<<setw(3)<<labels[y][x]<<" ";
      std::cerr<<std::endl;
    }
  #endif

  if(total_number_of_nodes>1){
    //Synchronize Threads TODO
    std::cerr.flush();
    if(my_node_number<total_number_of_nodes-1)
      world.send(my_node_number+2,SYNC_SIG, 0);

    //TODO: Could have master load elevations while nodes are processing. This
    //would speed things up slightly and reduce communication overhead.
    world.send(0, TOP_ELEVATIONS_TAG, elev.front  ());
    world.send(0, BOT_ELEVATIONS_TAG, elev.back   ());
    world.send(0, TOP_LABELS_TAG,     labels.front());
    world.send(0, BOT_LABELS_TAG,     labels.back ());
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

    #ifdef DEBUG
      for(auto &lo: label_offsets)
        std::cerr<<lo.first<<"->"<<lo.second<<std::endl;
    #endif

    std::cerr<<"Applying label offsets..."<<std::endl;
    for(int y=0;y<segment_height;y++)
    for(int x=0;x<width;x++){
      auto my_label = labels[y][x];
      if(elev[y][x]==no_data)
        continue;
      auto my_offset = label_offsets[my_label];
      if(elev[y][x]<=my_offset)
        elev[y][x] = my_offset;
    }
  }

  std::cerr<<"Writing out from "<<(my_node_number)<<std::endl;
  GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
  if(poDriver==NULL){
    std::cerr<<"Could not open GDAL driver."<<std::endl;
    return;
  }

  std::string output_name = std::string("output")+std::to_string(my_node_number)+std::string(".tif");
  GDALDataset *fout       = poDriver->Create(output_name.c_str(), width, segment_height, 1, gdt_t, NULL);

  if(fout==NULL){
    std::cerr<<"could not create output file."<<std::endl;
    return;
  }

  double geotrans[6];
  fin->GetGeoTransform(geotrans);

  // Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
  // Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
  // In case of north up images, the GT(2) and GT(4) coefficients are zero, and
  // the GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3))
  // position is the top left corner of the top left pixel of the raster.
  geotrans[3] += my_node_number*segment_height*geotrans[5];
  fout->SetGeoTransform(geotrans);

  const char* projection_string=fin->GetProjectionRef();
  fout->SetProjection(projection_string);

  GDALRasterBand *oband = fout->GetRasterBand(1);
  oband->SetNoDataValue(no_data);
  //poBand->RasterIO( GF_Write, 0, 0, no2output.shape()[0], no2output.shape()[1], no2output.origin(), no2output.shape()[0], no2output.shape()[1], GDT_Float32, 0, 0 );

  std::cerr<<"Writing out."<<std::endl;
  #ifdef DEBUG
  std::ofstream foutasc( std::string("output") + std::to_string(my_node_number) + std::string(".asc") );
  #endif
  for(int y=0;y<segment_height;y++){
    oband->RasterIO(GF_Write, 0, y, width, 1, elev[y].data(), width, 1, gdt_t, 0, 0);
    #ifdef DEBUG
      for(int x=0;x<width;x++){
        foutasc<<setw(3)<<elev[y][x]<<" ";
        cerr<<setw(3)<<elev[y][x]<<" ";
      }
      foutasc<<std::endl;
      std::cerr<<std::endl;
    #endif
  }

  GDALClose(fin);
  GDALClose(fout); //TODO

  //Synchronize Threads TODO
  std::cerr.flush();
  if(my_node_number<total_number_of_nodes-1)
    world.send(my_node_number+2,SYNC_SIG, 0);
}











template<class elev_t, GDALDataType gdt_t>
void DoMaster(char *dem_filename){
  typedef std::vector<elev_t>       ElevationRow;
  typedef std::vector<ElevationRow> Elevations;
  typedef std::map<label_t, std::map<label_t, elev_t> > Graph;
  mpi::communicator world;
  GDALAllRegister();

  const int total_number_of_nodes = world.size()-1; //Don't count master node

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
  GDALClose(fin);

  Elevations elev  (total_number_of_nodes*2,ElevationRow(width));
  Labels     labels(total_number_of_nodes*2,LabelRow    (width));
  std::vector<Graph> graphs(total_number_of_nodes);


  for(int i=1;i<=total_number_of_nodes;i++){
    int n=i-1;
    world.recv(i, TOP_ELEVATIONS_TAG, elev[2*n]);
    world.recv(i, BOT_ELEVATIONS_TAG, elev[2*n+1]);
    world.recv(i, TOP_LABELS_TAG, labels[2*n]);
    world.recv(i, BOT_LABELS_TAG, labels[2*n+1]);
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
    label_t newmaxlabel = 0;
    for(auto &fkey: graphs[i])
    for(auto &skey: fkey.second){
      auto flabel = fkey.first;
      auto slabel = skey.first;
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
    strip_to_max_label[i] = maxlabel;
    for(auto &s: labels[2*i])
      if(s!=1)
        s += maxlabel;
    for(auto &s: labels[2*i+1])
      if(s!=1)
        s += maxlabel;
    maxlabel = newmaxlabel;
  }

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
  for(int y=1;y<2*total_number_of_nodes-1;y+=2){
    for(int x=0;x<width;x++){
      if(elev[y][x]==no_data)
        continue;
      auto my_label = labels[y][x];
      for(int n=6;n<=8;n++){
        auto nx = x+dx[n];
        auto ny = y+dy[n];
        #ifdef DEBUG
          std::cerr<<"Considering ("<<x<<","<<y<<") with n="<<n<<" pointing to ("<<nx<<","<<ny<<")"<<std::endl;
        #endif
        if(nx<0 || nx==width || elev[ny][nx]==no_data)
          continue;
        auto other_label = labels[ny][nx];
        if(my_label==other_label) //Only when both are one
          continue;

          auto elev_over = std::max(elev[ny][nx],elev[y][x]);
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

  #ifdef DEBUG
    std::cerr<<"Master Graph"<<std::endl;
    for(auto &vertex: mastergraph)
    for(auto &edges: vertex.second)
      std::cerr<<setw(3)<<vertex.first<<" -> "<<setw(3)<<edges.first<<" = "<<setw(3)<<edges.second<<std::endl;
  #endif


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
    for(auto &ge: graph_elev)
      std::cerr<<setw(3)<<ge.first<<" = "<<setw(3)<<ge.second<<std::endl;
  #endif

  std::vector< std::map<label_t,elev_t> > strip_label_elevations(total_number_of_nodes);
  for(auto &ge: graph_elev){
    auto vertex_num   = ge.first;
    auto elevation    = ge.second;
    auto strip_num    = label_to_strip[vertex_num];
    auto label_offset = strip_to_max_label[strip_num];
    vertex_num       -= label_offset;
    strip_label_elevations[strip_num][vertex_num] = elevation;
  }

  std::cerr.flush();
  mpi::request reqs[total_number_of_nodes];
  for(int i=1;i<=total_number_of_nodes;i++){
    int n   = i-1;
    reqs[n] = world.isend(i, LABEL_OFFSETS, strip_label_elevations[n]); //TODO
    //world.send(i, LABEL_OFFSETS, strip_label_elevations[n]); //TODO
  }
  mpi::wait_all(reqs, reqs + total_number_of_nodes);
}

template<class elev_t, GDALDataType gdt_t>
void DoRouter(char *dem_filename){
  mpi::communicator world;
  if(world.rank()==0)
    DoMaster<elev_t,gdt_t>(dem_filename);
  else
    DoNode<elev_t,gdt_t>(dem_filename);
}

int main(int argc, char **argv){
  mpi::environment env;
  mpi::communicator world;

  if(argc!=2){
    std::cerr<<"Syntax: "<<argv[0]<<" <DEM>"<<std::endl;
    return -1;
  }

  GDALAllRegister();
  GDALDataset *fin = (GDALDataset*)GDALOpen(argv[1], GA_ReadOnly);
  if(fin==NULL){
    if(world.rank()==0)
      cerr<<"Could not peek at file: "<<argv[1]<<endl;
    return -1;
  }

  GDALRasterBand *band = fin->GetRasterBand(1);

  if(band->GetRasterDataType()==GDT_Int16)
    DoRouter<int16_t, GDT_Int16>(argv[1]);
  else if(band->GetRasterDataType()==GDT_Float32)
    DoRouter<float, GDT_Float32>(argv[1]);
  else if(world.rank()==0)
    std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(band->GetRasterDataType())<<std::endl;

  GDALClose(fin);

  return 0;
}