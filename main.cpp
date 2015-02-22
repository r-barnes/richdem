//Parallel D8 Flow Accumulation
//Author: Richard Barnes (rbarnes@umn.edu)
//Date:   2015-02-21
#include "gdal_priv.h"
#include <iostream>
#include <queue>
#include <boost/mpi.hpp>
#include <string>
#include <limits>
#include <cstdint>
//#define DEBUG 1

#ifdef DEBUG
  #include <fstream>
#endif

//These tags are used to identify the data being sent around by MPI. Each piece
//of information being sent should have a unique code.
#define TOP_LINKS_TAG        1
#define BOT_LINKS_TAG        2
#define TOP_ACCUMULATION_TAG 3
#define BOT_ACCUMULATION_TAG 4
#define TOP_FLOWDIRS_TAG     5
#define BOT_FLOWDIRS_TAG     6
#define SYNC_SIG             7


//This program assumes its input is in the form of a D8 flow direction matrix.
//Each cell's accumulated flow will be directed to one neighbouring cell. That
//neighbouring cell is identified with a neighbour code. If no flow is sent to
//any neighbouring cell, then the code `0` is used. The following codes are used
//to identify the eight neighbours:
//    234
//    105
//    876

//Using the above definition of neighbours, we can generate an x- and y-offset
//for each neighbour allowing us to succinctly address neighbours using:
//    x0+dx[n]
//    y0+dy[n]
//Note that the neighbours will be in the range [1,8] whereas index 0 refers to
//the central cell.

///x offsets of D8 neighbours, from a central cell
const int dx[9]={0,-1,-1,0,1,1,1,0,-1};
///y offsets of D8 neighbours, from a central cell
const int dy[9]={0,0,-1,-1,-1,0,1,1,1};

//Later in the program we will generate queues of grid cells to visit. The
//GridCell class is used to store the x and y values of these cells.
class GridCell{
 public:
  int x,y;
  GridCell(int x, int y) : x(x), y(y) {}
};

//The flow directions matrix stores the flow directions, as defined above. Since
//valid flow directions are in the range [0,8] with an additional value for
//no_data, we can store everything in the space of a single byte.
typedef int8_t flowdir_t;
typedef std::vector<flowdir_t> FlowdirsRow;
typedef std::vector< FlowdirsRow > Flowdirs;

//Since each cell drains entirely into a neighbouring cell, the flow
//accumulation is an integer. In the worst case, the maximal accumulation may be
//width*height of the entire input. This calls for at least a 32-bit integer. At
//the time of writing this usually permits values up to 2,147,483,647.
typedef int32_t accum_t;
typedef std::vector<accum_t> AccumRow;
typedef std::vector< AccumRow > Accum;

//Most cells flow into another cell. The number of cells flow into a cell
//represent its dependency count, which can be a value in the range [0,8]. Thus,
//a single byte is appropriate for storing dependency information.
typedef int8_t dependency_t;
typedef std::vector<dependency_t> DependenciesRow;
typedef std::vector< DependenciesRow > Dependencies;

//The links array maps the originating cells of flow paths that begin on a top
//or bottom edge of a segment to their terminal cells. It does this by storing
//the x-location of the terminal cell and using the integer value's negative
//sign to indicate whether the terminal cell is in the top or bottom row. Hence,
//a signed integer capable of storing values up to the width of the input is
//necessary.
typedef int32_t links_t;
typedef std::vector<links_t> LinksRow;
typedef std::vector< LinksRow > Links;

//When we are following flow paths sometimes a flow path will terminate in a
//place which guarnatees its flow cannot reach another segment and, therefore,
//that no other segment can be dependent on that flow path for its processing.
//We mark all such flow paths with the FLOW_TERMINATES value.
#define FLOW_TERMINATES std::numeric_limits<int>::max()




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
void FollowPath(
  const int x0,              //x-coordinate of initial cell
  const int y0,              //y-coordinate of initial cell
  const Flowdirs &flowdirs,  //Flow directions matrix
  const flowdir_t no_data,   //no_data value for flow directions matrix
  LinksRow &top_row_links,   //Output! Where the top initial cells go
  LinksRow &bottom_row_links //Output! Where the bottom initial cells go
){
  const int width  = flowdirs[0].size();
  const int height = flowdirs.size();
  int x            = x0;
  int y            = y0;

  while(true){                     //Follow the flow path until we reach its end
    int n  = flowdirs[y][x];       //Neighbour the current cell flows towards

    int nx,ny;                     //Neighbour's coordinates
    if(n<=0 || n==no_data){        //n<=0 indicates internal drainage or an
      n = 0;                       //invalid flow direction. Make a note of this
    } else {                       //with n=0.
      nx = x+dx[n];       //Valid flow direction. Mark neighbour's x-coordinate.
      ny = y+dy[n];       //Valid flow direction. Mark neighbour's y-coordinate.
    }

    //If the neighbour this cell flows into is off one of the sides of the grid,
    //or a no_data cell, or this cell does not flow into any neighbour, then
    //mark the initial cell from which we began this flow path as terminating
    //somewhere unimportant: its flow cannot pass to neighbouring segments/nodes
    //for further processing.
    if(n==0 || nx<0 || nx==width){
      if(y0==0) //Initial cell was on the top row of the segment
        top_row_links[x0]    = FLOW_TERMINATES;
      else      //Initial cell was on the bottom row of the segment
        bottom_row_links[x0] = FLOW_TERMINATES;

      return;

    //The neighbour cell is off the top or bottom of the segment. Therefore, its
    //flow may be passed on to a neighbouring segment. Thus, we need to link
    //this flow path's initial cell to this terminal cell. We want to minimize
    //data storage and transfer requirements, so we wish to indicate where the
    //terminal cell is using a single integer. Therefore, we store the terminal
    //cell's x-coordinate and use the negative sign to indicate that the
    //terminal cell was in the top row and the positive sign to indicate the
    //terminal cell was in the bottom row.
    } else if(ny<0 || ny==height){
      if(ny<0)  //If ny<0, the terminal cell is in the top row.
        x = -x; //Indicate this by making the x-coordinate negative.
      else      //If ny==height, terminal cell is in bottom row. Take x+1.
        x++;    //Thus (-Inf,0] implies ny<0 while [1,Inf) implies ny==height.
                //Ignoring the value INT_MAX that we use to indicate a terminal
                //flow path.

      if(y0==0) //Initial cell was in the top row
        top_row_links[x0]    = x;
      else      //Initial cell was in the bottom row
        bottom_row_links[x0] = x;

      return;

    //The flow path has not yet terminated. Continue following it.
    } else {
      x = nx;
      y = ny;
    }
  }
}


//As in the function above, we will start at an initial cell and follow its flow
//direction to its neighbour. Then we will follow the neighbour's flow direction
//to the neighbour's neighbour, and so on, until we can go no farther (as
//indicated by running of one of the edges of the segment or hitting a no_data
//cell).

//However, at this point we know how much flow accumulation to add to each cell.
//Therefore, we add this additional accumulation to each cell as we pass it.
void FollowPathAdd(
  int x,                     //Initial x-coordinate
  int y,                     //Initial y-coordinate
  const Flowdirs &flowdirs,  //Flow directions matrix
  const flowdir_t no_data,   //no_data value for  the flow directions matrix
  Accum &accum,              //Flow accumulation matrix
  const int additional_accum //Accumulation to add to each cell on the flow path
){
  int height = flowdirs.size();
  int width  = flowdirs[0].size();

  //Break when we reach a no_data cell
  if(flowdirs[y][x]==no_data)
    return;

  //Follow the flow path until it terminates
  while(true){
    //Add additional flow accumulation to this cell
    accum[y][x] += additional_accum;

    int n = flowdirs[y][x];  //Get neighbour
    if(n<=0 || n==no_data)   //This cell doesn't flow to a neighbour or is a
      return;                //no_data cell

    //Move to neighbour
    x += dx[n];
    y += dy[n];

    //The neighbour is off the edge of the grid. Flow path has terminated
    if(x<0 || x==width || y<0 || y==height)
      return;
  }
}




//This function handles all the functions of non-master nodes. The non-master
//nodes are labeled [0,total_number_of_nodes). Input is the name of the flow
//directions file to be read and processed.

//Each non-master node reads in a portion of the total input data and performs
//as many calculations with it as it can. It then passes a small amount of data
//to the master node which returns a small amount of data that enables each
//non-master node to complete the calculation as though everything had been done
//on a single node.
void doNode(int my_node_number, int total_number_of_nodes, char *flowdir_fname){
  GDALAllRegister();              //Load GDAL drivers so we can read and write
  boost::mpi::communicator world; //Open MPI communication line

  //Open dataset
  GDALDataset *fin = (GDALDataset*)GDALOpen(flowdir_fname, GA_ReadOnly);
  if(fin==NULL){
    std::cerr<<"Could not open file: "<<flowdir_fname<<std::endl;
    return;
  }

  //The first band is assumed to contain the flow directions
  GDALRasterBand *flowband = fin->GetRasterBand(1);
  flowdir_t no_data        = flowband->GetNoDataValue();
  int width                = flowband->GetXSize();
  int height               = flowband->GetYSize();

  //Each non-master node reads in a segment of the total input data. Here we
  //figure out what segment is assigned to this node.
  int segment_first_line = (height/total_number_of_nodes)*my_node_number;
  int segment_last_line  = (height/total_number_of_nodes)*(my_node_number+1);

  //If this is the last node, then read in all of the reamining data
  if(my_node_number==total_number_of_nodes-1)
    segment_last_line = height;
  int segment_height = segment_last_line - segment_first_line;

  //Read in the D8 flow directions data, one row at a time
  Flowdirs flowdirs(segment_height, FlowdirsRow(width));
  for(int y=segment_first_line;y<segment_last_line;y++){
    //Read data using its native format
    std::vector<int> temp(width);
    flowband -> RasterIO(GF_Read, 0, y, width, 1, temp.data(), width, 1, GDT_Int32, 0, 0);

    //Cast the data to our internal flowdir_t format.

    //Using `y-segment_first_line` converts from where we are reading in the
    //data to where that data is stored internally.
    flowdirs[y-segment_first_line] = FlowdirsRow(temp.begin(),temp.end());
  }


  //Each cell that flows points to a neighbouring cell. But not every cell is
  //pointed at. Cells which are not pointed at are the peaks from which flow
  //originates. In order to calculate the flow accumulation we begin at peaks
  //and pass flow downwards. Once flow has been passed downwards, the cell
  //receiving the flow no longer needs to wait for the cell which passed the
  //flow. When the receiving cell has no cells on which it is waiting, it then
  //becomes a peak itself. The number of cells pointing at a cell is its
  //"dependency count". In this section of the code we find each cell's
  //dependency count.
  std::cerr<<"Calculating dependencies..."<<std::endl;
  Dependencies dependencies(segment_height, DependenciesRow(width,0));
  for(int y=0;y<segment_height;y++)
  for(int x=0;x<width;x++){
    int n = flowdirs[y][x]; //The neighbour this cell flows into

    if(n<=0 || n==no_data)  //This cell does not flow into a neighbour
      continue;             //Or this cell is a no_data cell

    int nx = x+dx[n];       //x-coordinate of the neighbour
    int ny = y+dy[n];       //y-coordinate of the neighbour

    //Neighbour is not on the grid
    if(nx<0 || ny<0 || nx==width || ny==segment_height)
      continue;

    //Neighbour is valid and is part of the grid. The neighbour depends on this
    //cell, so increment its dependency count.
    dependencies[ny][nx]++;
  }

  //Now that we know how many dependencies each cell has, we can determine which
  //cells are the peaks: the sources of flow. We make a note of where the peaks
  //are for later use.
  std::queue<GridCell> sources;
  for(int y=0;y<segment_height;y++)
  for(int x=0;x<width;x++)
    //Valid cell with no dependencies: a peak!
    if(dependencies[y][x]==0 && flowdirs[y][x]!=no_data)
      sources.emplace(x,y);

  //Now that we know where the sources of flow are, we can start at this. Each
  //cell will have at least an accumulation of 1: itself. It then passes this
  //and any other accumulation it has gathered along its flow path to its
  //neighbour and decrements the neighbours dependency count. When a neighbour
  //has no more dependencies, it becomes a source.
  Accum accum(segment_height,AccumRow(width,0));
  while(!sources.empty()){         //There are sources remaining
    GridCell c = sources.front();  //Grab a source. Order is not important here.
    sources.pop();                 //We've visited this source. Discard it.

    int n = flowdirs[c.y][c.x];    //Who is this source's neighbour?
    if(n==no_data)                 //Oh snap! This isn't a real cell!
      continue;

    accum[c.y][c.x]++;             //This is a real cell, and it accumulates one
                                   //cell's worth of flow automatically.

    if(n<=0)                       //This cell doesn't flow anywhere.
      continue;                    //Move on to the next source.

    int nx = c.x+dx[n];            //Okay, this cell is going somewhere.
    int ny = c.y+dy[n];            //Make a note of where

    //This cell flows of the edge of the grid. Move on to next source.
    if(nx<0 || nx==width || ny<0 || ny==segment_height)
      continue;
    //This cell flows into a no_data cell. Move on to next source.
    if(flowdirs[ny][nx]==no_data)
      continue;

    //This cell has a neighbour it flows into. Add to its accumulation.
    accum[ny][nx] += accum[c.y][c.x];
    //Decrement the neighbour's dependencies.
    dependencies[ny][nx]--;

    //The neighbour has no more dependencies, so it has become a source
    if(dependencies[ny][nx]==0)
      sources.emplace(nx,ny);
  }

  //Some flow paths originate at the top and bottom of the segment this node is
  //processing. These arrays indicate where flow paths beginning at the top or
  //bottom of the segment terminate.
  LinksRow top_row_links   (width,FLOW_TERMINATES);
  LinksRow bottom_row_links(width,FLOW_TERMINATES);

  //If we are the top segment, nothing can flow into us, so we do not need to
  //know where flow paths originating at the top go to. On the otherhand, if we
  //are not the top segment, then consider each cell of the top row and find out
  //where its flow goes to.
  if(my_node_number!=0)
    for(int x=0;x<width;x++)
      //Only consider cells whose flow is directed down into the segment
      if( (5<=flowdirs[0][x] && flowdirs[0][x]<=8) || flowdirs[0][x]==1)
        FollowPath(x,0,flowdirs,no_data,top_row_links,bottom_row_links);

  //If we are the bottom segment, nothing can flow into us, so we do not need to
  //know where flow paths originating at the bottom go to. On the otherhand, if
  //we are not the bottom segment, then consider each cell of the bottom row and
  //find out where its flow goes to.
  if(my_node_number!=total_number_of_nodes-1)
    for(int x=0;x<width;x++)
      //Only consider cells whose flow is directed up into the segment
      if(1<=flowdirs[segment_height-1][x] && flowdirs[segment_height-1][x]<=5)
        FollowPath(x,segment_height-1,flowdirs,no_data,top_row_links,bottom_row_links);

  //Send our partial computation to the master node
  world.send(0, TOP_LINKS_TAG,        top_row_links);
  world.send(0, BOT_LINKS_TAG,        bottom_row_links);
  world.send(0, TOP_ACCUMULATION_TAG, accum.front());
  world.send(0, BOT_ACCUMULATION_TAG, accum.back());
  world.send(0, TOP_FLOWDIRS_TAG,     flowdirs.front());
  world.send(0, BOT_FLOWDIRS_TAG,     flowdirs.back());

  //Receive results of master node's computation. The master node will return
  //two arrays indicating how much additional flow must be added to each flow
  //cell's flow path.
  AccumRow accum_top(width), accum_bot(width);
  world.recv(0, TOP_ACCUMULATION_TAG, accum_top);
  world.recv(0, BOT_ACCUMULATION_TAG, accum_bot);

  //If we are not the top segment, considering each cell in the top row. If that
  //cell and its flow path should have additional accumulation, add it.
  if(my_node_number!=0)
    for(int x=0;x<width;x++)
      if(accum_top[x])
        FollowPathAdd(x, 0, flowdirs, no_data, accum, accum_top[x]);

  //If we are not the bottom segment, considering each cell in the bottom row.
  //If that cell and its flow path should have additional accumulation, add it.
  if(my_node_number!=total_number_of_nodes-1)
    for(int x=0;x<width;x++)
      if(accum_bot[x])
        FollowPathAdd(x, segment_height-1, flowdirs, no_data, accum, accum_bot[x]);


  //At this point we're done with the calculation! Boo-yeah!

  //Load an output driver
  GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
  if(poDriver==NULL){
    std::cerr<<"Could not open GDAL driver."<<std::endl;
    return;
  }

  //Open an output file
  std::string output_name = std::string("output")+std::to_string(my_node_number)+std::string(".tif");
  GDALDataset *fout       = poDriver->Create(output_name.c_str(), width, segment_height, 1, GDT_Int32, NULL);
  if(fout==NULL){
    std::cerr<<"could not create output file."<<std::endl;
    return;
  }

  //The geotransform maps each grid cell to a point in an affine-transformed
  //projection of the actual terrain. The geostransform is specified as follows:
  //    Xgeo = GT(0) + Xpixel*GT(1) + Yline*GT(2)
  //    Ygeo = GT(3) + Xpixel*GT(4) + Yline*GT(5)
  //In case of north up images, the GT(2) and GT(4) coefficients are zero, and
  //the GT(1) is pixel width, and GT(5) is pixel height. The (GT(0),GT(3))
  //position is the top left corner of the top left pixel of the raster.
  double geotrans[6];
  fin->GetGeoTransform(geotrans);

  //We shift the top-left pixel of the image southward to the appropriate
  //coordinate
  geotrans[3] += my_node_number*segment_height*geotrans[5];
  fout->SetGeoTransform(geotrans);

  //Copy the projection from the input file
  const char* projection_string=fin->GetProjectionRef();
  fout->SetProjection(projection_string);

  //We will write all of the data to the first raster band
  GDALRasterBand *oband = fout->GetRasterBand(1);

  //Since each real cell will have an accumulation of at least one, a no_data
  //value of zero is perfectly appropriate.
  oband->SetNoDataValue(0);

  std::cerr<<"Writing out."<<std::endl;
  //Debugging mode prints an ASCII file suitable for reviewing by hand if the
  //input is small
  #ifdef DEBUG
    std::ofstream foutasc( std::string("output") + std::to_string(my_node_number) + std::string(".asc") );
  #endif
  for(int y=0;y<segment_height;y++){
    oband->RasterIO(GF_Write, 0, y, width, 1, accum[y].data(), width, 1, GDT_Int32, 0, 0);
    #ifdef DEBUG
      for(int x=0;x<width;x++){
        int temp = accum[y][x];
        foutasc<<accum[y][x]<<" ";
      }
      foutasc<<std::endl;
    #endif
  }

  GDALClose(fin);
  GDALClose(fout);
}











//The master node receives the flow accumulations of the top and bottom rows of
//each segment, along with the flow directions of these rows and an array
//indicating whether those rows pass their flow between each other or if the
//flow terminates somewhere within the segment. From this information, the
//master node generates a return to each non-master node which allows those
//nodes to finish their computation as though everything had been done on one
//node.
void DoMaster(int total_number_of_nodes, char *flowdir_fname){
  GDALAllRegister();              //Load GDAL drivers so we can read and write
  boost::mpi::communicator world; //Open MPI communication line

  //Open D8 flow directions raster so that we can get its essential properties
  GDALDataset *fin = (GDALDataset*)GDALOpen(flowdir_fname, GA_ReadOnly);
  if(fin==NULL){
    std::cerr<<"Could not open file: "<<flowdir_fname<<std::endl;
    return;
  }

  //Get the essential properties
  GDALRasterBand *flowband = fin->GetRasterBand(1);
  int no_data              = flowband->GetNoDataValue();
  int width                = flowband->GetXSize();
  GDALClose(fin);

  Links         links       (total_number_of_nodes*2,LinksRow       (width  ));
  Accum         accum       (total_number_of_nodes*2,AccumRow       (width  ));
  Accum         accum_orig;
  Flowdirs      flowdirs    (total_number_of_nodes*2,FlowdirsRow    (width  ));
  Dependencies  dependencies(total_number_of_nodes*2,DependenciesRow(width,0));

  //Receive data from each non-master node
  for(int noden=1;noden<=total_number_of_nodes;noden++){
    //Since the master node is #0, we need to shift each incoming node's data so
    //it fits in our zero-based internal storage arrays
    int datan = noden-1;
    world.recv(noden, TOP_LINKS_TAG,        links   [2*datan]);
    world.recv(noden, BOT_LINKS_TAG,        links   [2*datan+1]);
    world.recv(noden, TOP_ACCUMULATION_TAG, accum   [2*datan]);
    world.recv(noden, BOT_ACCUMULATION_TAG, accum   [2*datan+1]);
    world.recv(noden, TOP_FLOWDIRS_TAG,     flowdirs[2*datan]);
    world.recv(noden, BOT_FLOWDIRS_TAG,     flowdirs[2*datan+1]);
  }


  //Each cell that flows points to a neighbouring cell. But not every cell is
  //pointed at. Cells which are not pointed at are the peaks from which flow
  //originates. In order to calculate the flow accumulation we begin at peaks
  //and pass flow downwards. Once flow has been passed downwards, the cell
  //receiving the flow no longer needs to wait for the cell which passed the
  //flow. When the receiving cell has no cells on which it is waiting, it then
  //becomes a peak itself. The number of cells pointing at a cell is its
  //"dependency count". In this section of the code we find each cell's
  //dependency count. In contradistinction to the homologous code in DoNode(),
  //our neighbours here may be at the other ends of flow paths.
  for(int y=0;y<links.size();y++)
  for(int x=0;x<width;x++){
    int n = flowdirs[y][x];  //Which direction does this cell flow?

    if(n<=0 || n==no_data)   //This cell doesn't flow, or the cell is a no_data.
      continue;              //Move on to the next cell.

    int nx = x+dx[n];        //The cell does flow. Make a note of the
    int ny = y+dy[n];        //coordinates of its neighbour.

    //The cell flows off the edge of the grid. Move on to next cell.
    if(nx<0 || ny<0 || nx==width || ny==flowdirs.size())
      continue;

    //In the master node, each set of two rows represents an entire segment's
    //worth of cells. Moving up or down within this representation of a segment
    //means that we are not going to a neighbouring cell, but, rather, we are
    //going to the terminal cell of a flow path originating at the cell we are
    //currently considering.
    if(y/2==ny/2){      //These cells are part of the same segment
      nx = links[y][x]; //Find the x-coordinate of the flow path's terminal cell
      //This cell is the initial cell of a flow path which terminates internally
      //within the segment. Therefore, move on to the next cell.
      if(nx==FLOW_TERMINATES)
        continue;
      //In this case, we have used `y` to find out what segment we are on. To
      //locate the flow path's terminal cell, begin by setting `ny` to the top
      //y-coordinate of the segment.
      ny = (ny%2==0)?ny:ny-1;
      if(nx<=0){       //This path ends on the top of the segment
        nx = -nx;      //Map (-Inf,0] to [0,Inf]
      } else {         //This path ends on the bottom of the segment
        ny++;          //ny will now point to the bottom of the segment
        nx--;          //Map [1,Inf) to [0,Inf)
      }
    }

    //Decrement the dependencies of the "neighbouring" cell
    dependencies[ny][nx]++;
  }


  //Cells which receive flow from another segment have already propagated their
  //flow to the terminal cells of their flow paths in DoNode(). Therefore, we
  //should set the accumulation of these initial cells to zero to avoid having
  //it be propagated to the terminal cell twice.
  for(size_t y=0;y<links.size();y++)
  for(int x=0;x<width;x++){
    int fd = flowdirs[y][x]; //Flow direction of this cell
    //This cell goes nowhere or is on the bottom of a segment going up into the
    //segment
    if( y%2==1 && (fd==1 || fd==2 || fd==3 || fd==4 || fd==5 || fd<=0) )
      accum[y][x] = 0;
    //This cell gos nowhere or is on the top of a segment going down into the
    //segement
    else if( y%2==0 && (fd==1 || fd==5 || fd==8 || fd==7 || fd==6 || fd<=0) )
      accum[y][x] = 0;
  }

  //We save a copy of the accumulation grid as it stands now. All cells which
  //receive flow and pass it into their segments have accumulation values of
  //zero. All cells which pass flow out of their segments have positive values.
  //Later, we will want an output where only the cells passing flow into the
  //segments have positive values and all other cells are zero. This is used to
  //help generate that output.
  accum_orig = accum;

  //Cells which have no dependencies are the sources of flow.
  std::queue<GridCell> sources;
  //The top row of the top segment and the bottom row of the bottom segment
  //cannot be the source of flow, so ignore them
  for(size_t y=1;y<links.size()-1;y++)
  for(int x=0;x<width;x++)
    //A real cell with no dependencies is a source
    if(dependencies[y][x]==0 && flowdirs[y][x]!=no_data)
      sources.emplace(x,y);

  //Now that we have identified sources of flow, visit each of them and pass
  //their flow on to their neighbours.
  while(!sources.empty()){        //Are their sources left?
    GridCell c = sources.front(); //Yes! Get the next source, order unimportant.
    sources.pop();                //We're done with this source now.

    int n = flowdirs[c.y][c.x];   //Where does this cell flow to?
    if(n<=0 || n==no_data)        //This cell doesn't flow, or is not a real
      continue;                   //cell. Move on to the next source.

    int nx = c.x+dx[n];           //This cell does flow somewhere and is real,
    int ny = c.y+dy[n];           //make a note of where it flows to

    //This cell flows off the edge of the grid. Move on to the next source.
    if(nx<0 || ny<0 || nx==width || ny==flowdirs.size())
      continue;

    //As before, when we were calculating the dependencies, if cells are part of
    //the same segement, we need to move accumulation to the terminal cells of
    //the flow paths their originate.
    if(c.y/2==ny/2){          //Cells part of the same segment
      nx = links[c.y][c.x];   //Get x-coordinate of terminal cell
      if(nx==FLOW_TERMINATES) //Cell terminates internally in segment, move on
        continue;             //to the next source.

      //In this case, we have used `y` to find out what segment we are on. To
      //locate the flow path's terminal cell, begin by setting `ny` to the top
      //y-coordinate of the segment.
      ny = (ny%2==0)?ny:ny-1;
      if(nx<=0){              //This path ends on the top of the segment
        nx = -nx;             //Map (-Inf,0] to [0,Inf]
      } else {                //This path ends on the bottom of the segment
        ny++;                 //ny will now point to the bottom of the segment
        nx--;                 //Map [1,Inf) to [0,Inf)
      }
    }

    //Add accumulation to the "neighbouring cell"
    accum[ny][nx] += accum[c.y][c.x];
    //Decrement the dependencies of the "neighbouring" cell
    dependencies[ny][nx]--;

    //If the "neigbouring" cell has no dependencies, it becomes a source
    if(dependencies[ny][nx]==0)
      sources.emplace(nx,ny);
  }


  //We now need to calculate how much accumulation should be added to the inputs
  //of each flow path. Therefore, we need the terminal cells of flow paths to
  //have zero accumulation. But, a segment may pass flow directly to a terminal
  //cell without going through a flow path. Therefore, for each cell, we
  //subtract the flow passed to that cell from any flow paths which terminate at
  //it. We also subtract any flow the cell had originally because it was the
  //terminal cell of a flow path.
  for(int y=0;y<links.size();y++)
  for(int x=0;x<width;x++){
    //This cell had no accumulation, so it cannot need correcting
    if(accum[y][x]==0)
      continue;

    //Which cell does this cell flow to?
    int n = flowdirs[y][x];

    //This cell has no flow or is a no_data cell. Ignore it
    if(n<=0 || n==no_data)
      continue;

    //TODO: Should this maybe go above the `n<=0||n==no_data` check?
    accum[y][x] -= accum_orig[y][x];

    //This is a real cell and its pointing somewhere. Make a note of where.
    int nx = x+dx[n];
    int ny = y+dy[n];

    //The flow goes off the edge of the grid. Move on to the next cell.
    if(nx<0 || ny<0 || nx==width || ny==flowdirs.size())
      continue;

    //As before, when we were calculating the dependencies, if cells are part of
    //the same segement, we need to move accumulation to the terminal cells of
    //the flow paths their originate.
    if(y/2==ny/2){
      nx = links[y][x];
      if(nx==FLOW_TERMINATES) //Path we are going into ends on a side edge or internally
        continue;
      ny = (ny%2==0)?ny:ny-1; //Map ny to the nearest strip top
      if(nx<=0){       //This path ends on the top of the strip
        nx = -nx;      //Map (-Inf,0] to [0,Inf]
      } else {         //This path ends on the bottom of the strip
        ny++;          //ny will not point to the bottom of the strip
        nx--;          //Map [1,Inf) to [0,Inf)
      }

      //Subtract this cell's input from the cell at the end of its flow path.
      //The non-master nodes will propagate this cell's accumulation to that
      //cell anyway.
      accum[ny][nx] -= accum[y][x];
    }
  }

  //Send information to each of the non-master nodes so they can complete their
  //calculations.
  for(int i=1;i<=total_number_of_nodes;i++){
    int n=i-1;
    world.send(i, TOP_ACCUMULATION_TAG, accum[2*n]);
    world.send(i, BOT_ACCUMULATION_TAG, accum[2*n+1]);
  }
}

int main(int argc, char **argv){
  boost::mpi::environment  env;
  boost::mpi::communicator world;

  if(argc!=2){
    std::cerr<<"Syntax: "<<argv[0]<<" <DEM>"<<std::endl;
    MPI_Finalize();
    return -1;
  }

  if(world.rank()>0)
    doNode(world.rank()-1,world.size()-1,argv[1]);
  else
    DoMaster(world.size()-1,argv[1]);

  return 0;
}