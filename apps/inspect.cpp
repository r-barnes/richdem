#include <iostream>
#include <iomanip>
#include "../libs/common/Array2D.hpp"

template<class T>
void Master(std::string input, int minx, int maxx, int miny, int maxy, bool flipV, bool flipH){
  Array2D<T> inp(input, 0, 0, 0, 0, false);

  const int fwidth = 10;

  int cx = -1;
  int cy = -1;

  if(inp.geotransform[1]<0)
    flipH = !flipH;
  if(inp.geotransform[5]>0)
    flipV = !flipV;

  if(flipV)
    inp.flipVert();
  if(flipH)
    inp.flipHorz();

  if(maxy==-1){
    cx         = minx;
    cy         = maxx;
    int radius = miny;
    minx       = std::max(cx-radius,0);
    maxx       = std::min(cx+radius,inp.width()-1);
    miny       = std::max(cy-radius,0);
    maxy       = std::min(cy+radius,inp.height()-1);
  }

  if(minx>maxx){
    std::cerr<<"MinX must be smaller than MaxX!"<<std::endl;
    return;
  }

  if(miny>maxy){
    std::cerr<<"MinY must be smaller than MaxY!"<<std::endl;
    return;
  }

  std::cout<<std::setw(fwidth)<<" ";
  for(int x=minx;x<=maxx;x++)
    std::cout<<std::setw(fwidth)<<x;
  std::cout<<std::endl;

  for(int y=miny;y<=maxy;y++){
    std::cout<<std::setw(fwidth)<<y;
    for(int x=minx;x<=maxx;x++){
      if(cx==x && cy==y)
        std::cout<<"\033[93m";

      if(NativeTypeToGDAL<T>()==GDT_Byte)
        std::cout<<std::setw(10)<<(int)inp(x,y);
      else
        std::cout<<std::setw(10)<<inp(x,y);

      if(cx==x && cy==y)
        std::cout<<"\033[39m";
    }
    std::cout<<std::endl;
  }
}

int main(int argc, char **argv){
  int32_t total_height;
  int32_t total_width;
  GDALDataType file_type;

  bool flipV  = false;
  bool flipH  = false;
  std::string inputfile;
  int minx    = -1;
  int maxx    = -1;
  int miny    = -1;
  int maxy    = -1;

  try{
    for(int i=1;i<argc;i++){
      if(strcmp(argv[i],"--flipH")==0 || strcmp(argv[i],"-H")==0){
        flipH = true;
      } else if(strcmp(argv[i],"--flipV")==0 || strcmp(argv[i],"-V")==0){
        flipV = true;
      } else if(argv[i][0]=='-'){
        throw std::invalid_argument("Unrecognised flag: "+std::string(argv[i]));
      } else if(inputfile.size()==0){
        inputfile = argv[i];
      } else if(minx==-1){
        minx = std::stoi(argv[i]);
      } else if(maxx==-1){
        maxx = std::stoi(argv[i]);
      } else if(miny==-1){
        miny = std::stoi(argv[i]);
      } else if(maxy==-1){
        maxy = std::stoi(argv[i]);
      } else {
        throw std::invalid_argument("Too many arguments.");
      }
    }
    if(inputfile.size()==0 || minx==-1 || maxx==-1 || miny==-1)
      throw std::invalid_argument("Too few arguments.");
    if(maxy!=-1 && minx>maxx)
      throw std::invalid_argument("MinX must be smaller than MaxX!");
    if(maxy!=-1 && miny>maxy)
      throw std::invalid_argument("MinY must be smaller than MaxY!");
    if(minx<0)
      throw std::invalid_argument("MinX must be >=0!");
    if(maxx<0)
      throw std::invalid_argument("MaxX must be >=0!");
    if(miny<0)
      throw std::invalid_argument("MinY must be >=0!");
  } catch (const std::invalid_argument &ia){
    std::cerr<<argv[0]<<" [--flipV] [--flipH] <Center X> <Center Y> <Radius>"<<std::endl;
    std::cerr<<argv[0]<<" [--flipV] [--flipH] <Min X> <Max X> <Min Y> <Max Y>"<<std::endl;

    std::cerr<<"###Error: "<<ia.what()<<std::endl;

    return -1;
  }

  //Get the total dimensions of the input file
  getGDALDimensions(inputfile, total_height, total_width, file_type, NULL);

  switch(file_type){
    case GDT_Byte:
      Master<uint8_t >(inputfile,minx,maxx,miny,maxy, flipV, flipH);break;
    case GDT_UInt16:
      Master<uint16_t>(inputfile,minx,maxx,miny,maxy, flipV, flipH);break;
    case GDT_Int16:
      Master<int16_t >(inputfile,minx,maxx,miny,maxy, flipV, flipH);break;
    case GDT_UInt32:
      Master<uint32_t>(inputfile,minx,maxx,miny,maxy, flipV, flipH);break;
    case GDT_Int32:
      Master<int32_t >(inputfile,minx,maxx,miny,maxy, flipV, flipH);break;
    case GDT_Float32:
      Master<float   >(inputfile,minx,maxx,miny,maxy, flipV, flipH);break;
    case GDT_Float64:
      Master<double  >(inputfile,minx,maxx,miny,maxy, flipV, flipH);break;
    default:
      std::cerr<<"Unrecognised data type!"<<std::endl;
      return -1;
  }
}