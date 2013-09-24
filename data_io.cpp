/**
  @file
  Implementation code for reading and writing data, primarily in the
  ArcGrid ASCII format.

  Richard Barnes (rbarnes@umn.edu), 2012
*/
#include <stdio.h>
#include <stdlib.h>
#include "data_io.hpp"
#include "interface.hpp"
#include "data_structures.hpp"
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include "utility.hpp"
#include <fstream>
#include <algorithm>

//write_arrows
/**
  @brief  Writes an arrow representation of a 2D D8 flow direction array
  @author Richard Barnes

  @param[in] &filename     Name of file to write to
  @param[in] &elevations   2D array of D8 flow directions

  @returns 0 upon success
*/
int write_arrows(const char filename[], const char_2d &flowdirs){
  std::locale::global(std::locale(""));
  std::wofstream fout;
  Timer write_time;
  ProgressBar progress;

  write_time.start();

  diagnostic_arg("Opening arrow output file \"%s\"...",filename);
  fout.open(filename);
  if(!fout.is_open()){
    diagnostic("failed!\n");
    exit(-1);  //TODO: Need to make this safer! Don't just close after all that work!
  }
  diagnostic("succeeded.\n");

  diagnostic("%%Writing arrows...\n");
  progress.start( flowdirs.width()*flowdirs.height() );
  for(int y=0;y<flowdirs.height();++y){
    progress.update( y*flowdirs.width() );
    for(int x=0;x<flowdirs.width();++x){
      if(flowdirs(x,y)==flowdirs.no_data)  //TODO: Crude way of detecting chars and bools
        fout<<L" ";
      else if (flowdirs(x,y)==NO_FLOW)
        fout<<fd[0];
      else
        fout<<fd[flowdirs(x,y)];
      fout<<L" ";
    }
    fout<<std::endl;
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());

  fout.close();

  write_time.stop();
  diagnostic_arg("Write time was: %lf\n", write_time.accumulated());

  return 0;
}


/**
Provides an input stream operator to the must_be class, which enables it to
assert that the next token of the input must match the argument provided to
the constructor of the class.
*/
std::istream& operator>>( std::istream &is, const must_be &a ){
  std::string inp, mstr;
  size_t cpos=is.tellg();
  is >> inp;
  mstr=a.match;
  std::transform(inp.begin(), inp.end(), inp.begin(), ::tolower);
  std::transform(mstr.begin(), mstr.end(), mstr.begin(), ::tolower);
  if(inp!=mstr){
    is.seekg(cpos);
    std::cerr<<"Failed to match required input string '"<<a.match<<"'. Found '"<<inp<< "'."<<std::endl;
    throw std::string("Failed to match!");
  }
  return is;
}
