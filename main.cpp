//Compile with
// g++ -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal
// ./a.out ~/projects/watershed/data/beauford03.flt
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <iostream>
#include <cmath>
#include <stdint.h>
#include <map>
#include <queue>
using namespace std;

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

void doNode(int my_node_number, int total_number_of_nodes, char *filename){
  GDALAllRegister();

  GDALDataset *demdata = (GDALDataset*)GDALOpen(filename, GA_ReadOnly);
  if(demdata==NULL){
    cerr<<"Could not open file: "<<filename<<endl;
    return;
  }

  GDALRasterBand *demband = demdata->GetRasterBand(1);

  int width  = demband->GetXSize();
  int height = demband->GetYSize();

  std::cerr<<"width line: "<<width<<std::endl;
  std::cerr<<"height line: "<<height<<std::endl;

  int first_line = (height/total_number_of_nodes)*my_node_number-1;
  int last_line  = (height/total_number_of_nodes)*(my_node_number+1)+1;

  if(first_line<0)     first_line = 0;
  if(last_line>height) last_line  = height;
  int lines_to_read = last_line-first_line;

  std::cerr<<"First line: "<<first_line<<std::endl;
  std::cerr<<"Last line: "<<last_line<<std::endl;

  if(my_node_number==total_number_of_nodes)
    last_line = height;

  std::vector< std::vector<float> > dem         (lines_to_read, std::vector<float>(width));
  std::vector< std::vector<char > > flowdirs    (lines_to_read, std::vector<char >(width));
  std::vector< std::vector<char > > dependencies(lines_to_read, std::vector<char >(width,0));
  std::vector< std::vector<int> >   accum       (lines_to_read, std::vector<int  >(width,0));

  for(int y=first_line;y<last_line;y++){
    //float *vals = (float*) CPLMalloc(sizeof(float)*width);
    //if(vals==NULL){
    //  cerr<<"Could not allocate vals."<<endl;
    //  return -1;
    //}
    //if(wprox_band->GetXSize()!=nlcd_band->GetXSize() ||

    demband -> RasterIO( GF_Read, 0, y, width, 1, &dem[y-first_line][0], width, 1, GDT_Float32, 0, 0 );

    //CPLFree(vals);
  }

  ////////////////////////
  //Find flow directions
  ////////////////////////

  //Set outside to point out
  for(int y=0;y<lines_to_read;y++){
    flowdirs[y][0]       = 1;
    flowdirs[y][width-1] = 5;
  }
  for(int x=0;x<width;x++){
    flowdirs[0][x]               = 3;
    flowdirs[lines_to_read-1][x] = 7;
  }

  std::cerr<<dem.size()<<" "<<dem[0].size()<<std::endl;

  //TODO: Add no data

  for(int y=0;y<lines_to_read;y++)
  for(int x=0;x<width;x++){
    int lowest      = 0;
    int elev_lowest = 99999999; //TODO

    for(int n=1;n<=8;n++){
      int nx=x+dx[n];
      int ny=y+dy[n];
      if(nx<0 || nx==width || ny<0 || ny==lines_to_read) continue;
      if(dem[ny][nx]>=dem[y][x]) continue;
      if(lowest==0 || elev_lowest>dem[ny][nx]){
        lowest      = n;
        elev_lowest = dem[ny][nx];
      }
    }
    flowdirs[y][x] = lowest;
    dependencies[y+dy[lowest]][x+dx[lowest]]++;
  }

  ////////////////////////
  //Find peaks
  ////////////////////////
  std::queue<GridCell> sources;
  for(int y=0;y<lines_to_read;y++)
  for(int x=0;x<width;x++)
    if(dependencies[y][x]==0)
      sources.push(GridCell(x,y));

  while(!sources.empty()){
    GridCell c = sources.front();
    sources.pop();

    accum[c.y][c.x]++;

    if(flowdirs[c.y][c.x]==0) continue;

    int nx = c.x+dx[flowdirs[c.y][c.x]];
    int ny = c.y+dy[flowdirs[c.y][c.x]];
    if(nx<0 || nx==width || ny<0 || ny==lines_to_read)
      continue;

    dependencies[ny][nx]--;

    if(dependencies[ny][nx]==0)
      sources.push(GridCell(nx,ny));
  }

  GDALClose(demdata);
}

int main(int argc, char **argv){
  int nodes = 2;
  doNode(0,nodes,argv[1]);
}