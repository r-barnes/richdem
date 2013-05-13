/**
  @file
  Header and templated code for reading and writing data, primarily in the
  ArcGrid ASCII format.

  Richard Barnes (rbarnes@umn.edu), 2012
*/
#ifndef _data_io_included
#define _data_io_included

#include "interface.hpp"
#include <fstream>
#include <string>
#include <iostream>
#include <typeinfo>
//#include <fcntl.h> //Used for posix_fallocate

int write_arrows(const char filename[], const char_2d &flowdirs);

#define OUTPUT_DEM  1
#define OUTPUT_OMG  2

class must_be {
  public:
    std::string match;
    must_be(const std::string &match) : match(match) {}
};

std::istream& operator>>( std::istream &is, const must_be &a ){
  std::string inp;
  is >> inp;
  if(inp!=a.match){
    std::cerr<<"Failed to match required input string '" << a.match << "'. Found '" <<inp<<"'."<<std::endl;
    exit(-1); //TODO: Can we fail gracefully?
  }
  return is;
}



/**
  @brief  Reads an ArcGrid ASCII file
  @author Richard Barnes

  @param[in]  &filename     Name of ArcGrid ASCII file to read
  @param[out] &elevations   DEM object containing contents of file

  @todo Won't handle char data input correctly

  @returns 0 upon success
*/
template<class T>
int load_ascii_data(std::string filename, array2d<T> &elevations){
  std::ifstream fin;
  size_t file_size;
  int rows,columns;
  Timer load_time;
  ProgressBar progress;

  load_time.start();

  diagnostic_arg("Opening input ASCII-DEM file \"%s\"...",filename.c_str());
  fin.open(filename.c_str());
  if(!fin.good()){
    diagnostic("failed!\n");
    exit(-1);
  }
  diagnostic("succeeded.\n");

  diagnostic("Calculating file size...");
  fin.seekg(0, fin.end);
  file_size=fin.tellg();
  fin.seekg(0, fin.beg);
  diagnostic("succeeded.\n");

//  posix_fadvise(fileno(fin),0,0,POSIX_FADV_SEQUENTIAL);

  diagnostic("Reading DEM header...");
  fin>>must_be("ncols")         >>columns;
  fin>>must_be("nrows")         >>rows;
  fin>>must_be("xllcorner")     >>elevations.xllcorner;
  fin>>must_be("yllcorner")     >>elevations.yllcorner;
  fin>>must_be("cellsize")      >>elevations.cellsize;
  fin>>must_be("NODATA_value")  >>elevations.no_data;
  diagnostic("succeeded.\n");

  diagnostic_arg("The loaded DEM will require approximately %ldMB of RAM.\n",columns*rows*((long)sizeof(float))/1024/1024);

  diagnostic("Resizing elevation matrix...");  //TODO: Consider abstracting this block
  elevations.resize(columns,rows);
  diagnostic("succeeded.\n");

  diagnostic("%%Reading elevation matrix...\n");
  progress.start(file_size);

  elevations.data_cells=0;
  for(int y=0;y<rows;y++){
    progress.update(fin.tellg()); //Todo: Check to see if ftell fails here?
    for(int x=0;x<columns;x++){
      fin>>elevations(x,y);
      if(elevations(x,y)!=elevations.no_data)
        elevations.data_cells++;
    }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());

  fin.close();

  diagnostic_arg(
    "Read %ld cells, of which %ld contained data (%ld%%).\n",
    elevations.width()*elevations.height(), elevations.data_cells,
    elevations.data_cells*100/elevations.width()/elevations.height()
  );

  load_time.stop();
  diagnostic_arg("Read time was: %lfs\n", load_time.accumulated());

  return 0;
}


/**
  @brief  Writes an ArcGrid ASCII file or OmniGlyph file
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  &filename     Name of ArcGrid ASCII file to write
  @param[in]  &elevations   DEM object to write
  @param[in]  precision     Floating-point precision of output

  @returns 0 upon success
*/
template <class T>
int output_ascii_data(
  const std::string filename,
  const array2d<T> &output_grid,
  int precision=8
){
  std::ofstream fout;
  std::string outputsep=" ";
  int output_type=OUTPUT_DEM;
  Timer write_time;
  ProgressBar progress;

  write_time.start();

  diagnostic_arg("Opening ASCII output file \"%s\"...",filename.c_str());
  fout.open(filename.c_str());
  if(!fout.is_open()){
    diagnostic("failed!\n");
    exit(-1);  //TODO: Need to make this safer! Don't just close after all that work!
  }
  diagnostic("succeeded.\n");

  //OmniGlyph output
  if(filename.substr(filename.length()-4)==".omg"){
    outputsep="|";
    output_type=OUTPUT_OMG;
    diagnostic("Writing OmniGlyph file header...");
    fout<<"Contents: Pixel array"<<std::endl;
    fout<<std::endl;
    fout<<"Width:    "<<output_grid.width()<<std::endl;
    fout<<"Height:   "<<output_grid.height()<<std::endl;
    fout<<std::endl;
    fout<<"Spectral bands:   1"<<std::endl;
    fout<<"Bits per band:   32"<<std::endl;
    fout<<"Range of values:   "<<output_grid.min()<<","<<output_grid.max()<<std::endl;
    fout<<"Actual range:   "<<output_grid.no_data<<","<<output_grid.max()<<std::endl;  //TODO: Assumes no_data is a small negative value
    fout<<"Gamma exponent:   0."<<std::endl;
    fout<<"Resolution:   100 pixels per inch"<<std::endl;
    fout<<std::endl;
    fout<<"|"<<std::endl;
  } else {
    diagnostic("Writing ArcGrid ASCII file header...");
    fout<<"ncols\t\t"<<output_grid.width()<<std::endl;
    fout<<"nrows\t\t"<<output_grid.height()<<std::endl;
    fout<<"xllcorner\t"<<std::fixed<<std::setprecision(precision)<<output_grid.xllcorner<<std::endl;
    fout<<"yllcorner\t"<<std::fixed<<std::setprecision(precision)<<output_grid.yllcorner<<std::endl;
    fout<<"cellsize\t"<<std::fixed<<std::setprecision(precision)<<output_grid.cellsize<<std::endl;
    fout<<"NODATA_value\t"<<std::fixed<<std::setprecision(precision);
    if(sizeof(T)==1)  //TODO: Crude way of detecting chars and bools
      fout<<(int)output_grid.no_data<<std::endl;
    else
      fout<<output_grid.no_data<<std::endl;
  }
  diagnostic("succeeded.\n");

  diagnostic("%%Writing ArcGrid ASCII file data...\n");
  progress.start( output_grid.width()*output_grid.height() );
  fout.precision(precision);
  fout.setf(std::ios::fixed);
  for(int y=0;y<output_grid.height();y++){
    progress.update( y*output_grid.width() );
    if(output_type==OUTPUT_OMG)
      fout<<"|";
    for(int x=0;x<output_grid.width();x++)
      if(sizeof(T)==1)  //TODO: Crude way of detecting chars and bools
        fout<<(int)output_grid(x,y)<<outputsep;
      else
        fout<<output_grid(x,y)<<outputsep;
    fout<<std::endl;
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());

  fout.close();

  write_time.stop();
  diagnostic_arg("Write time was: %lf\n", write_time.accumulated());

  return 0;
}








/**
  @brief  Writes a floating-point grid file
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  &basename     Name, without extension, of output file
  @param[in]  &output_grid  DEM object to write

  @todo Does not check byte order (big-endian, little-endian)

  @returns 0 upon success
*/
template <class T>
int write_floating_data(
  const std::string basename,
  const array2d<T> &output_grid
){
  Timer write_time;
  ProgressBar progress;
  std::string fn_header(basename), fn_data(basename);

  //TODO: The section below should work, but is something of an abomination
  if(typeid(T)==typeid(float)){
    fn_header+=".hdr";
    fn_data+=".flt";
  } else if (typeid(T)==typeid(double)){
    fn_header+=".hdr";
    fn_data+=".dflt";
  } else {
    std::cerr<<"Cannot read floating type data into this format!"<<std::endl;
    exit(-1);
  }

  write_time.start();


  {
    diagnostic_arg("Opening floating-point header file \"%s\" for writing...",fn_header.c_str());
    std::ofstream fout;
    fout.open(fn_header.c_str());
    if(!fout.is_open()){
      diagnostic("failed!\n");
      exit(-1);  //TODO: Need to make this safer! Don't just close after all that work!
    }
    diagnostic("succeeded.\n");

    diagnostic("Writing floating-point header file...");
    fout<<"ncols\t\t"<<output_grid.width()<<std::endl;
    fout<<"nrows\t\t"<<output_grid.height()<<std::endl;
    fout<<"xllcorner\t"<<std::fixed<<std::setprecision(10)<<output_grid.xllcorner<<std::endl;
    fout<<"yllcorner\t"<<std::fixed<<std::setprecision(10)<<output_grid.yllcorner<<std::endl;
    fout<<"cellsize\t"<<std::fixed<<std::setprecision(10)<<output_grid.cellsize<<std::endl;
    fout<<"NODATA_value\t"<<std::fixed<<std::setprecision(10)<<output_grid.no_data<<std::endl;
    fout<<"BYTEORDER\tLSBFIRST"<<std::endl; //TODO
    fout.close();
    diagnostic("succeeded.\n");
  }


  diagnostic_arg("Opening floating-point data file \"%s\" for writing...",fn_data.c_str());

  {
    std::ofstream fout(fn_data.c_str(), std::ios::binary | std::ios::out);
    if(!fout.is_open()){
      diagnostic("failed!\n");
      exit(-1);  //TODO: Need to make this safer! Don't just close after all that work!
    }
    diagnostic("succeeded.\n");

    diagnostic("%%Writing floating-point data file...\n");
    progress.start( output_grid.width()*output_grid.height() );
    for(int y=0;y<output_grid.height();++y){
      progress.update( y*output_grid.width() );
      for(int x=0;x<output_grid.width();++x)
        fout.write(reinterpret_cast<const char*>(&output_grid(x,y)), std::streamsize(sizeof(T)));
    }
    fout.close();
    write_time.stop();
    diagnostic_arg(SUCCEEDED_IN,progress.stop());
  }

  diagnostic_arg("Write time was: %lf\n", write_time.accumulated());

  return 0;
}






/**
  @brief  Reads in a floating-point grid file
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  &basename     Name, without extension, of input file
  @param[in]  &grid         DEM object in which to store data

  @todo Does not check byte order (big-endian, little-endian)

  @returns 0 upon success
*/
template <class T>
int read_floating_data(
  const std::string basename,
  array2d<T> &grid
){
  Timer io_time;
  ProgressBar progress;
  std::string fn_header(basename), fn_data(basename);

  //TODO: The section below should work, but is something of an abomination
  if(typeid(T)==typeid(float)){
    fn_header+=".hdr";
    fn_data+=".flt";
  } else if (typeid(T)==typeid(double)){
    fn_header+=".hdr";
    fn_data+=".dflt";
  } else {
    std::cerr<<"Cannot read floating type data into this format!"<<std::endl;
    exit(-1);
  }

  int columns, rows;
  std::string byteorder;

  io_time.start();


  {
    std::ifstream fin;
    diagnostic_arg("Opening floating-point header file \"%s\" for reading...",fn_header.c_str());
    fin.open(fn_header.c_str());
    if(fin==NULL){
      diagnostic("failed!\n");
      exit(-1);
    }
    diagnostic("succeeded.\n");


    diagnostic("Reading DEM header...");
    fin>>must_be("ncols")         >>columns;
    fin>>must_be("nrows")         >>rows;
    fin>>must_be("xllcorner")     >>grid.xllcorner;
    fin>>must_be("yllcorner")     >>grid.yllcorner;
    fin>>must_be("cellsize")      >>grid.cellsize;
    fin>>must_be("NODATA_value")  >>grid.no_data;
    fin>>must_be("BYTEORDER")     >>byteorder;
    diagnostic("succeeded.\n");
    fin.close();
  }

  diagnostic_arg("The loaded DEM will require approximately %ldMB of RAM.\n",columns*rows*((long)sizeof(float))/1024/1024);

  diagnostic("Resizing grid...");  //TODO: Consider abstracting this block
  grid.resize(columns,rows);
  diagnostic("succeeded.\n");



  diagnostic_arg("Opening floating-point data file \"%s\" for reading...",fn_data.c_str());

  {
    std::ifstream fin(fn_data.c_str(), std::ios::binary | std::ios::in);
    if(!fin.is_open()){
      diagnostic("failed!\n");
      exit(-1);  //TODO: Need to make this safer! Don't just close after all that work!
    }
    diagnostic("succeeded.\n");


    diagnostic("%%Reading data...\n");
    progress.start(columns*rows);
    grid.data_cells=0;
    for(int y=0;y<rows;++y){
      progress.update(y*columns); //Todo: Check to see if ftell fails here?
      for(int x=0;x<columns;++x){
        fin.read(reinterpret_cast<char*>(&grid(x,y)), std::streamsize(sizeof(T)));
        if(grid(x,y)!=grid.no_data)
          grid.data_cells++;
      }
    }
    io_time.stop();
    diagnostic_arg(SUCCEEDED_IN,progress.stop());

  }

  diagnostic_arg("Write time was: %lf\n", io_time.accumulated());

  return 0;
}


/**
  @brief  Universal read function. Calls everything else.
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]   &filename    Name of file to read in
  @param[out]  &grid        DEM object in which to store data

  @returns 0 upon success
*/
template<class T>
int read_data(std::string filename, array2d<T> &grid){
  if( filename.substr(filename.size()-3)=="flt" )
    return read_floating_data(filename.substr(0,filename.size()-4), grid);
  else
    return load_ascii_data(filename, grid);
}

#endif
