//Compile with
// mpic++ -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall
// mpirun -n 3 ./a.out ~/projects/watershed/data/beauford03.flt
// gdalbuildvrt -input_file_list ifiles gdalbuildvrt -input_file_list ifiles output.vrt
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <iostream>
#include <cmath>
#include <stdint.h>
#include <map>
#include <queue>
#include <mpi.h>
#include <string>
#include <iomanip>
#include <cassert>

#define INTERNODE_DELAY 2000

#include <chrono>
#include <thread>

using namespace std;

#define TOP_LINKS_TAG        1
#define BOT_LINKS_TAG        2
#define TOP_ACCUMULATION_TAG 3
#define BOT_ACCUMULATION_TAG 4
#define TOP_FLOWDIRS_TAG     5
#define BOT_FLOWDIRS_TAG     6

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

class GridCell{
 public:
  int x,y;
  GridCell(int x, int y) : x(x), y(y) {}
};

void FollowPath(const int x0, const int y0, const int width, const int height, const std::vector< std::vector<char > > &flowdirs, const char no_data, std::vector<int> &top_row_links, std::vector<int> &bottom_row_links){
  int x = x0;
  int y = y0;
  //std::cerr<<"----------------------------"<<std::endl;
  while(true){
    int n  = flowdirs[y][x];

    int nx,ny;
    if(n<=0 || n==no_data){
      n=0;
    } else {
      nx = x+dx[n];
      ny = y+dy[n];
    }

    //std::cerr<<x<<" "<<y<<std::endl;

    if(n==0 || nx<0 || nx==width){
      if(y0==0)
        top_row_links[x0]    = 99999999;
      else
        bottom_row_links[x0] = 99999999;
      return;
    } else if(ny<0 || ny==height){
      if(ny<0)
        x = -x;

      if(y0==0)
        top_row_links[x0]    = x;
      else
        bottom_row_links[x0] = x;
      return;
    } else {
      x = nx;
      y = ny;
    }
  }
}

void FollowPathAdd(const int x0, const int y0, const int width, const int height, const std::vector< std::vector<char > > &flowdirs, const char no_data, std::vector< std::vector<int> > &accum, const int additional_accum){
  int x = x0;
  int y = y0;
  //std::cerr<<"----------------------------"<<std::endl;
  while(true){
    accum[y][x] += additional_accum;

    int n  = flowdirs[y][x];
    if(n<=0 || n==no_data)
      return;

    int nx = x+dx[n];
    int ny = y+dy[n];


    if(nx<0 || nx==width || ny<0 || ny==height)
      return;

    x = nx;
    y = ny;
  }
}





void doNode(int my_node_number, int total_number_of_nodes, char *flowdir_fname){
  std::this_thread::sleep_for(std::chrono::milliseconds(INTERNODE_DELAY*my_node_number));
  std::cerr<<std::endl<<std::endl;
  std::cerr<<"Node #"<<my_node_number<<std::endl;
  std::cerr<<"================="<<std::endl;

  GDALAllRegister();

  GDALDataset *fin = (GDALDataset*)GDALOpen(flowdir_fname, GA_ReadOnly);
  if(fin==NULL){
    cerr<<"Could not open file: "<<flowdir_fname<<endl;
    return;
  }

  GDALRasterBand *flowband = fin->GetRasterBand(1);

  char no_data = flowband->GetNoDataValue();
  std::cerr<<"No data value: "<<(int)no_data<<std::endl;

  int width  = flowband->GetXSize();
  int height = flowband->GetYSize();

  int segment_first_line = (height/total_number_of_nodes)*my_node_number;
  int segment_last_line  = (height/total_number_of_nodes)*(my_node_number+1);
  if(my_node_number==total_number_of_nodes-1)
    segment_last_line = height;
  int segment_height = segment_last_line - segment_first_line;

  //Read in DEM data
  std::vector< std::vector<char> > flowdirs(segment_height, std::vector<char>(width));
  for(int y=segment_first_line;y<segment_last_line;y++){
    std::vector<int> temp(width);
//    flowband -> RasterIO( GF_Read, 0, y, width, 1, temp.data(), width, 1, GDT_Float64, 0, 0 );
    flowband -> RasterIO( GF_Read, 0, y, width, 1, temp.data(), width, 1, GDT_Int32, 0, 0 );
    flowdirs[y-segment_first_line]=std::vector<char>(temp.begin(),temp.end());
  }
  //TODO: Why? "y-segment_first_line" in the above?

  // for(int y=0;y<segment_height;y++){
  //   for(int x=0;x<width;x++)
  //     std::cout<<setw(3)<<(int)flowdirs[y][x]<<" ";
  //   std::cout<<std::endl;
  // }

  ////////////////////////
  //Find dependencies
  ////////////////////////
  std::cerr<<"Calculating dependencies..."<<std::endl;
  std::vector< std::vector<char> > dependencies(segment_height, std::vector<char>(width,0));
  for(int y=0;y<segment_height;y++)
  for(int x=0;x<width;x++){
    int n = flowdirs[y][x];

    if(n<=0 || n==no_data) continue; //TODO: Make this compatible with ArcGIS

    int nx = x+dx[n];
    int ny = y+dy[n];
    if(nx<0 || ny<0 || nx==width || ny==segment_height) //TODO: Check ArcGIS compatibility
      continue;

    dependencies[ny][nx]++;
  }

  ////////////////////////
  //Find peaks
  ////////////////////////
  std::cerr<<"Finding peaks"<<std::endl;
  std::queue<GridCell> sources;
  for(int y=0;y<segment_height;y++)
  for(int x=0;x<width;x++)
    if(dependencies[y][x]==0 && flowdirs[y][x]!=no_data)
      sources.emplace(x,y);

  ////////////////////////
  //Calculate accumulation
  ////////////////////////
  std::cerr<<"Calculating accumulation"<<std::endl;
  std::vector< std::vector<int> > accum(segment_height,std::vector<int>(width,0));

  while(!sources.empty()){
    GridCell c = sources.front();
    sources.pop();

    int n = flowdirs[c.y][c.x];
    if(n==no_data)
      continue;

    accum[c.y][c.x]++;

    if(n<=0)
      continue;

    int nx = c.x+dx[n];
    int ny = c.y+dy[n];
    if(nx<0 || nx==width || ny<0 || ny==segment_height)
      continue;
    if(flowdirs[ny][nx]==no_data)
      continue;

    accum[ny][nx] += accum[c.y][c.x];
    dependencies[ny][nx]--;

    if(dependencies[ny][nx]==0)
      sources.emplace(nx,ny);
  }

  std::vector<int> top_row_links   (width,-1);
  std::vector<int> bottom_row_links(width,-1);

  ////////////////////////
  //Load path beginnings for incoming paths
  ////////////////////////
  std::cerr<<"Connecting top edges"<<std::endl;
  if(my_node_number!=0)
    for(int x=0;x<width;x++){
      // std::cerr<<x<<" ";
      switch(flowdirs[0][x]){
        case 1:
        case 5:
        case 8:
        case 7:
        case 6:
          FollowPath(x,0,width,segment_height,flowdirs,no_data,top_row_links,bottom_row_links);
      }
    }
  // std::cerr<<std::endl;

  std::cerr<<"Connecting bottom edges"<<std::endl;
  if(my_node_number!=total_number_of_nodes-1)
    for(int x=0;x<width;x++)
      switch(flowdirs[segment_height-1][x]){
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
          FollowPath(x,segment_height-1,width,segment_height,flowdirs,no_data,top_row_links,bottom_row_links);
      }

  //int MPI_Send(void* buf,int count,MPI_Datatype datatype,int dest,int tag,MPI_Comm comm);
  std::cerr<<"Sending "<<(my_node_number+1)<<std::endl;
  MPI_Send(top_row_links   .data(), top_row_links   .size(), MPI_INT, 0,TOP_LINKS_TAG,        MPI_COMM_WORLD);
  MPI_Send(bottom_row_links.data(), bottom_row_links.size(), MPI_INT, 0,BOT_LINKS_TAG,        MPI_COMM_WORLD);
  MPI_Send(accum.front()   .data(), accum.front()   .size(), MPI_INT, 0,TOP_ACCUMULATION_TAG, MPI_COMM_WORLD);
  MPI_Send(accum.back ()   .data(), accum.back ()   .size(), MPI_INT, 0,BOT_ACCUMULATION_TAG, MPI_COMM_WORLD);
  MPI_Send(flowdirs.front().data(), flowdirs.front().size(), MPI_BYTE,0,TOP_FLOWDIRS_TAG,     MPI_COMM_WORLD);
  MPI_Send(flowdirs.back ().data(), flowdirs.back ().size(), MPI_BYTE,0,BOT_FLOWDIRS_TAG,     MPI_COMM_WORLD);
  std::cerr<<"Sent "<<(my_node_number+1)<<std::endl;

  //std::cerr<<"Result receiving "<<(my_node_number+1)<<std::endl;
  MPI_Status status;
  std::vector<int> accum_top(width), accum_bot(width);
  MPI_Recv(accum_top.data(), accum_top.size(), MPI_INT, 0, TOP_ACCUMULATION_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(accum_bot.data(), accum_bot.size(), MPI_INT, 0, BOT_ACCUMULATION_TAG, MPI_COMM_WORLD, &status);
  //std::cerr<<"Result received "<<(my_node_number+1)<<std::endl;

  std::this_thread::sleep_for(std::chrono::milliseconds(INTERNODE_DELAY*my_node_number));

  // std::cerr<<std::endl<<std::endl;
  // std::cerr<<"Node #"<<my_node_number<<std::endl;
  // std::cerr<<"================="<<std::endl;
/*
  std::cerr<<"Received top: ";
  for(auto &x: accum_top)
    std::cerr<<setw(4)<<x<<" ";
  std::cerr<<endl;

  std::cerr<<"Received bot: ";
  for(auto &x: accum_bot)
    std::cerr<<setw(4)<<x<<" ";
  std::cerr<<endl;
*/
  if(my_node_number!=0){
  //  std::cerr<<"Now adding to paths on top."<<std::endl;
    for(int x=0;x<width;x++)
      if(accum_top[x]){
        // std::cerr<<"Path at ("<<x<<",0) had "<<accum_top[x]<<" accumulation pointing towards "<<(int)flowdirs[0][x]<<std::endl;
        switch(flowdirs[0][x]){
          case 1:
          case 5:
          case 8:
          case 7:
          case 6:
            // std::cerr<<"Adding to path at "<<x<<std::endl;
            FollowPathAdd(x, 0, width, segment_height, flowdirs, no_data, accum, accum_top[x]);//-accum[0][x]);
        }
      }
  }

  //std::cerr<<"Connecting bottom edges"<<std::endl;
  if(my_node_number!=total_number_of_nodes-1){
    //std::cerr<<"Now adding to paths on bot."<<std::endl;
    for(int x=0;x<width;x++)
      if(accum_bot[x]){
        // std::cerr<<"Path at ("<<x<<",Ym) had "<<accum_bot[x]<<" accumulation pointing towards "<<(int)flowdirs[segment_height-1][x]<<std::endl;
        switch(flowdirs[segment_height-1][x]){
          case 1:
          case 2:
          case 3:
          case 4:
          case 5:
            // std::cerr<<"Adding to path at "<<x<<std::endl;
            FollowPathAdd(x, segment_height-1, width, segment_height, flowdirs, no_data, accum, accum_bot[x]);//-accum[segment_height-1][x]);
        }
      }
  }


  std::cerr<<"Writing out from "<<(my_node_number)<<std::endl;
  GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
  if(poDriver==NULL){
    std::cerr<<"Could not open GDAL driver."<<std::endl;
    return;
  }

  std::string output_name = std::string("output")+std::to_string(my_node_number)+std::string(".tif");
  GDALDataset *fout       = poDriver->Create(output_name.c_str(), width, segment_height, 1, GDT_Int32, NULL);

  if(fout==NULL){
    std::cerr<<"could not create output file."<<std::endl;
    return;
  }

  std::cerr<<"Files opened successfully."<<std::endl;
  double geotrans[6];
  fin->GetGeoTransform(geotrans);

  // Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
  // Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
  // In case of north up images, the GT(2) and GT(4) coefficients are zero, and
  // the GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3))
  // position is the top left corner of the top left pixel of the raster.
  geotrans[3] += my_node_number*segment_height*geotrans[5];
  fout->SetGeoTransform(geotrans);

  std::cerr<<"Got geotransform!"<<std::endl;

  const char* projection_string=fin->GetProjectionRef();
  fout->SetProjection(projection_string);

  GDALRasterBand *oband = fout->GetRasterBand(1);
  oband->SetNoDataValue(0);
  //poBand->RasterIO( GF_Write, 0, 0, no2output.shape()[0], no2output.shape()[1], no2output.origin(), no2output.shape()[0], no2output.shape()[1], GDT_Float32, 0, 0 );

  std::cerr<<"Writing out."<<std::endl;
  for(int y=0;y<segment_height;y++){
    for(int x=0;x<width;x++){
      //int temp = accum[y][x];
      //oband->RasterIO( GF_Write, x, y, 1, 1, &temp, 1, 1, GDT_Int32, 0, 0 );
      std::cout<<setw(3)<<accum[y][x]<<" ";
    }
    std::cout<<std::endl;
  }

  GDALClose(fin);
  //GDALClose(fout); //TODO
}












void DoMaster(int my_node_number, int total_number_of_nodes, char *flowdir_fname){
  GDALAllRegister();

  GDALDataset *fin = (GDALDataset*)GDALOpen(flowdir_fname, GA_ReadOnly);
  if(fin==NULL){
    cerr<<"Could not open file: "<<flowdir_fname<<endl;
    return;
  }

  GDALRasterBand *flowband = fin->GetRasterBand(1);

  int no_data = flowband->GetNoDataValue();

  int width  = flowband->GetXSize();
  GDALClose(fin);

  std::vector< std::vector<int> >  links        (total_number_of_nodes*2,std::vector<int> (width  ));
  std::vector< std::vector<int> >  accum        (total_number_of_nodes*2,std::vector<int> (width  ));
  std::vector< std::vector<char> > flowdirs     (total_number_of_nodes*2,std::vector<char>(width  ));
  std::vector< std::vector<char> > dependencies (total_number_of_nodes*2,std::vector<char>(width,0));

  //int MPI_Recv(void* buf,int count,MPI_Datatype datatype,int source,int tag, MPI_Comm comm, MPI_Status *status);

  //Gather top links
  std::cerr<<"Total nodes: "<<(total_number_of_nodes)<<std::endl;
  for(int i=1;i<=total_number_of_nodes;i++){
    int n=i-1;
    MPI_Status status;
    std::cerr<<"Receiving "<<i<<std::endl;
    MPI_Recv(links   [2*n]  .data(), links   [2*n]  .size(),MPI_INT, i,TOP_LINKS_TAG,        MPI_COMM_WORLD,&status);
    MPI_Recv(links   [2*n+1].data(), links   [2*n+1].size(),MPI_INT, i,BOT_LINKS_TAG,        MPI_COMM_WORLD,&status);
    MPI_Recv(accum   [2*n]  .data(), accum   [2*n]  .size(),MPI_INT, i,TOP_ACCUMULATION_TAG, MPI_COMM_WORLD,&status);
    MPI_Recv(accum   [2*n+1].data(), accum   [2*n+1].size(),MPI_INT, i,BOT_ACCUMULATION_TAG, MPI_COMM_WORLD,&status);
    MPI_Recv(flowdirs[2*n]  .data(), flowdirs[2*n]  .size(),MPI_BYTE,i,TOP_FLOWDIRS_TAG,     MPI_COMM_WORLD,&status);
    MPI_Recv(flowdirs[2*n+1].data(), flowdirs[2*n+1].size(),MPI_BYTE,i,BOT_FLOWDIRS_TAG,     MPI_COMM_WORLD,&status);
    std::cerr<<"Received "<<i<<std::endl;
  }

  std::cerr<<"Received accumulation grid"<<std::endl;
  for(size_t y=0;y<accum.size();y++){
    if(y%2==0)
      std::cerr<<"--------"<<std::endl;
    for(size_t x=0;x<accum[0].size();x++)
      std::cerr<<setw(3)<<accum[y][x]<<" ";
    std::cerr<<std::endl;
  }

  std::cerr<<"Flowdirs"<<std::endl;
  for(size_t y=0;y<flowdirs.size();y++){
    if(y%2==0)
      std::cerr<<"--------"<<std::endl;
    for(size_t x=0;x<flowdirs[0].size();x++)
      std::cerr<<setw(3)<<(int)flowdirs[y][x]<<" ";
    std::cerr<<std::endl;
  }

  std::cerr<<"Links"<<std::endl;
  for(size_t y=0;y<links.size();y++){
    if(y%2==0)
      std::cerr<<"--------"<<std::endl;
    for(size_t x=0;x<links[0].size();x++)
      std::cerr<<setw(3)<<(int)links[y][x]<<" ";
    std::cerr<<std::endl;
  }

  std::cerr<<"Finding dependencies.."<<std::endl;
  for(size_t y=0;y<links.size();y++)
  for(int x=0;x<width;x++){
    int n = flowdirs[y][x];

    if(n<=0 || n==no_data)
      continue;

    int nx = x+dx[n];
    int ny = y+dy[n];
    if(nx<0 || ny<0 || nx==width || ny==flowdirs.size()) //TODO: Consider using some kind of height-esque variable here
      continue;

    //Part of the same strip. Use the links
    if(y%2==ny%2){
      nx = links[y][x];
      if(nx==99999999)
        continue;
      if( (nx<0 && ny>y) || (ny>0 && ny<y) )
        ny = y;
      if(nx<0)
        nx = -nx;
    }

    dependencies[ny][nx]++;
  }

  std::cerr<<"Convert from accumulation to delta-accumulation..."<<std::endl;
  for(size_t y=0;y<links.size();y++)
  for(int x=0;x<width;x++){
    int fd = flowdirs[y][x];
    //On the bottom going up into the strip
    if( y%2==1 && (fd==1 || fd==2 || fd==3 || fd==4 || fd==5) )
      accum[y][x] = 0;
    //On the top going down into the strip
    if( y%2==0 && (fd==1 || fd==5 || fd==8 || fd==7 || fd==6) )
      accum[y][x] = 0;
  }

  std::cerr<<"Accumulation grid following delta accumulation"<<std::endl;
  for(size_t y=0;y<accum.size();y++){
    if(y%2==0)
      std::cerr<<"--------"<<std::endl;
    for(size_t x=0;x<accum[0].size();x++)
      std::cerr<<setw(3)<<accum[y][x]<<" ";
    std::cerr<<std::endl;
  }

  std::cerr<<"Finding sources..."<<std::endl;
  std::queue<GridCell> sources;
  for(size_t y=1;y<links.size()-1;y++) //Don't need to worry about top and bottom strips
  for(int x=0;x<width;x++)
    if(dependencies[y][x]==0){
      // std::cerr<<"Source at ("<<x<<","<<y<<")"<<std::endl;
      sources.emplace(x,y);
    }

  std::cerr<<"Accumulating across aggregated grid"<<std::endl;
  while(!sources.empty()){
    GridCell c = sources.front();
    sources.pop();

    int n = flowdirs[c.y][c.x];
    if(n<=0 || n==no_data)
      continue;

    int nx = c.x+dx[n];
    int ny = c.y+dy[n];
    if(nx<0 || ny<0 || nx==width || ny==flowdirs.size())
      continue;

    //Part of the same strip. Use the links
    if(c.y%2==ny%2){
      nx = links[c.y][c.x];
      if(nx==99999999)
        continue;
      if( (nx<0 && ny>c.y) || (ny>0 && ny<c.y) )
        ny = c.y;
      if(nx<0)
        nx = -nx;
    }

    std::cout<<"Pouring ("<<c.x<<","<<c.y<<") into ("<<nx<<","<<ny<<") with "<<accum[c.y][c.x]<<std::endl;

    accum[ny][nx] += accum[c.y][c.x];
    dependencies[ny][nx]--;

    if(dependencies[ny][nx]==0)
      sources.emplace(nx,ny);
  }


  // for(auto &r: accum){
  //   for(auto &x: r)
  //     std::cerr<<setw(4)<<x<<" ";
  //   std::cerr<<std::endl;
  // }

  std::cerr<<"Net accumulation grid"<<std::endl;
  for(size_t y=0;y<accum.size();y++){
    if(y%2==0)
      std::cerr<<"--------"<<std::endl;
    for(size_t x=0;x<accum[0].size();x++)
      std::cerr<<setw(3)<<accum[y][x]<<" ";
    std::cerr<<std::endl;
  }




  std::cerr<<"Dispersing accumulation..."<<std::endl;
  for(int i=1;i<=total_number_of_nodes;i++){
    int n=i-1;
    // accum[2*n]=std::vector<int>(width,400);
    // accum[2*n+1]=std::vector<int>(width,400);
    MPI_Send(accum[2*n]  .data(), accum[2*n]  .size(),MPI_INT, i, TOP_ACCUMULATION_TAG, MPI_COMM_WORLD);
    MPI_Send(accum[2*n+1].data(), accum[2*n+1].size(),MPI_INT, i, BOT_ACCUMULATION_TAG, MPI_COMM_WORLD);
  }
}

int main(int argc, char **argv){
  MPI_Init(&argc,&argv);

  if(argc!=2){
    std::cerr<<"Syntax: "<<argv[0]<<" <DEM>"<<std::endl;
    MPI_Finalize();
    return -1;
  }

  int tid, nthreads;
  MPI_Comm_rank(MPI_COMM_WORLD, &tid);
  MPI_Comm_size(MPI_COMM_WORLD, &nthreads);

  if(tid>0)
    doNode(tid-1,nthreads-1,argv[1]);
  else
    DoMaster(tid-1,nthreads-1,argv[1]);

  MPI_Finalize();

  return 0;
}