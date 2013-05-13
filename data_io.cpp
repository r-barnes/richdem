/**
  @file
  Implementation code for reading and writing data, primarily in the
  ArcGrid ASCII format.

  Richard Barnes (rbarnes@umn.edu), 2012
*/
#include <stdio.h>
#include <stdlib.h>
#include "interface.hpp"
#include "data_structures.hpp"
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include "utility.hpp"
#include <fstream>

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
