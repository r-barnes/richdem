#include "gdal_priv.h"
#include <vector>
#include <deque>
#include <iostream>
#include <limits>
#include <cstdint>

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

#define FLOWDIR_NO_DATA 255


template<class T, GDALDataType gdt_t>
void ProcessFile(GDALRasterBand *inp, GDALRasterBand *out){
  T no_data  = inp->GetNoDataValue();
  int width  = inp->GetXSize();
  int height = inp->GetYSize();

  //Rotating queue to hold the elevation data
  std::deque< std::vector<T> > elev;

  //Since each real cell will have an accumulation of at least one, a no_data
  //value of zero is perfectly appropriate.
  out->SetNoDataValue(FLOWDIR_NO_DATA);

  //Top row points always up when there are data cells
  {
    //Read top row in
    std::vector<T>    temp    (width);
    std::vector<char> flowdirs(width,FLOWDIR_NO_DATA);
    inp -> RasterIO(GF_Read, 0, 0, width, 1, temp.data(), width, 1, gdt_t, 0, 0);
    elev.push_back(temp);

    for(int x=0;x<width;x++)
      if(temp[x]!=no_data)
        flowdirs[x] = 3;
    out->RasterIO(GF_Write, 0, 0, width, 1, flowdirs.data(), width, 1, GDT_Byte, 0, 0);
  }

  //Read in the elevation data, one row at a time
  std::vector<char> flowdirs(width);
  std::vector<T>    inp_elev(width);
  for(int y=1;y<height-1;y++){
    if(elev.size()==3){
      //Leftmost cell points left out of raster
      if(elev[1][0]==no_data)
        flowdirs[0] = FLOWDIR_NO_DATA;
      else
        flowdirs[0] = 1;

      for(int x=1;x<width-1;x++){
        if(elev[1][x]==no_data){
          flowdirs[x] = FLOWDIR_NO_DATA;
          continue;
        }

        char lowest_n = 0;
        T lowest_elev = std::numeric_limits<T>::max();
        for(char n=1;n<=8;n++){
          int nx = x+dx[n];
          int ny = 1+dy[n];
          if(elev[ny][nx]==no_data){
            lowest_n = n;
            break;
          } else if(elev[ny][nx]<lowest_elev){
            lowest_n    = n;
            lowest_elev = elev[ny][nx];
          }
        }
        flowdirs[x] = lowest_n;
        // std::cerr<<lowest_n<<std::endl;
      }

      //Rightmost cell points right out of raster
      if(elev[1][width-1]==no_data)
        flowdirs[width-1] = FLOWDIR_NO_DATA;
      else
        flowdirs[width-1] = 5;

      out->RasterIO(GF_Write, 0, y, width, 1, flowdirs.data(), width, 1, GDT_Byte, 0, 0);

      elev.pop_front();
    }

    //Read data using its native format

    inp->RasterIO(GF_Read, 0, y, width, 1, inp_elev.data(), width, 1, gdt_t, 0, 0);
    elev.push_back(inp_elev);
  }

  //Top row points always points down when there are data cells
  {
    //Read top row in
    std::vector<T>    temp    (width);
    std::vector<char> flowdirs(width,FLOWDIR_NO_DATA);
    inp -> RasterIO(GF_Read, 0, height-1, width, 1, temp.data(), width, 1, gdt_t, 0, 0);
    elev.push_back(temp);

    for(int x=0;x<width;x++)
      if(temp[x]!=no_data)
        flowdirs[x] = 7;
    out->RasterIO(GF_Write, 0, height-1, width, 1, flowdirs.data(), width, 1, GDT_Byte, 0, 0);
  }

}


int main(int argc, char **argv){
  if(argc!=4){
    std::cerr<<"Syntax: "<<argv[0]<<" <Input> <Output> <Format>"<<std::endl;
    return -1;
  }

  std::cerr<<"Running "<<GDALVersionInfo("--version")<<std::endl;

  GDALAllRegister();              //Load GDAL drivers so we can read and write

  //Open dataset
  GDALDataset *fin = (GDALDataset*)GDALOpen(argv[1], GA_ReadOnly);
  if(fin==NULL){
    std::cerr<<"Could not open input file: "<<argv[1]<<std::endl;
    return -1;
  }

  //The first band is assumed to contain elevations
  GDALRasterBand *inp = fin->GetRasterBand(1);
  int width           = inp->GetXSize();
  int height          = inp->GetYSize();

  //Load an output driver
  GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(argv[3]);
  if(poDriver==NULL){
    std::cerr<<"Could not open GDAL driver."<<std::endl;
    return -1;
  }

  //Open an output file
  GDALDataset *fout = poDriver->Create(argv[2], width, height, 1, GDT_Byte, NULL);
  if(fout==NULL){
    std::cerr<<"Could not create output file: "<<argv[2]<<std::endl;
    return -1;
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
  fout->SetGeoTransform(geotrans);

  //Copy the projection from the input file
  const char* projection_string=fin->GetProjectionRef();
  fout->SetProjection(projection_string);


  //We will write all of the data to the first raster band
  GDALRasterBand *oband = fout->GetRasterBand(1);

  switch(inp->GetRasterDataType()){
    case GDT_Int16:
      ProcessFile<int16_t, GDT_Int16>(inp,oband);
      break;
    case GDT_Int32:
      ProcessFile<int32_t, GDT_Int32>(inp,oband);
      break;
    case GDT_Float32:
      ProcessFile<float, GDT_Float32>(inp,oband);
      break;
    default:
      std::cerr<<"Unrecognised data type: "<<GDALGetDataTypeName(inp->GetRasterDataType())<<std::endl;
      return -1;
  }

  GDALClose(fin);
  GDALClose(fout);

  return 0;
}