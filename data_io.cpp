/**
  @file
  Implementation code for reading and writing data, primarily in the
  ArcGrid ASCII format.

  Richard Barnes (rbarnes@umn.edu), 2012
*/
#include <stdio.h>
#include <stdlib.h>
#include "interface.h"
#include "data_structures.h"
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include "utility.h"
#include <fstream>

//load_ascii_data
/**
  @brief  Reads an ArcGrid ASCII file
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  &filename     Name of ArcGrid ASCII file to read
  @param[out] &elevations   DEM object containing contents of file

  @todo Has minor bug where progress becomes negative at some point.

  @returns 0 upon success
*/
int load_ascii_data(const char filename[], float_2d &elevations){
  FILE *fin;
  long long file_size;
  int rows,columns;
  Timer load_time;
  ProgressBar progress;

  load_time.start();

  errno=0;
  diagnostic_arg("Opening input ASCII-DEM file \"%s\"...",filename);
  fin=fopen(filename,"r");
  if(fin==NULL){
    diagnostic_arg("failed with error %d: \"%s\"!\n",errno,strerror(errno));
    exit(-1);
  }
  diagnostic("succeeded.\n");

  diagnostic("Calculating file size...");
  if(fseek(fin,-1,SEEK_END)!=0){
    diagnostic("failed! (Couldn't jump to end of file.)\n");
    exit(-1);
  }
  if((file_size=ftell(fin))==-1){
    diagnostic("failed! (Couldn't determine file size.)\n");
    exit(-1);
  }
  if(fseek(fin,0,SEEK_SET)!=0){
    diagnostic("failed! (Couldn't jump back to beginning of file.)\n");
    exit(-1);
  }
  diagnostic("succeeded.\n");

//  posix_fadvise(fileno(fin),0,0,POSIX_FADV_SEQUENTIAL);

  diagnostic("Reading DEM header...");
  if(fscanf(fin,"ncols %d nrows %d xllcorner %lf yllcorner %lf cellsize %lf NODATA_value %f",&columns, &rows, &elevations.xllcorner, &elevations.yllcorner, &elevations.cellsize, &elevations.no_data)!=6){
    diagnostic("failed!\n");
    exit(-1);
  }
  diagnostic("succeeded.\n");

  diagnostic_arg("The loaded DEM will require approximately %ldMB of RAM.\n",columns*rows*((long)sizeof(float))/1024/1024);

  diagnostic("Resizing elevation matrix...");  //TODO: Consider abstracting this block
  elevations.resize(columns,rows);
  diagnostic("succeeded.\n");

  diagnostic("%%Reading elevation matrix...\n");
  progress.start(file_size);
  float temp;
  elevations.data_cells=0;
  for(int y=0;y<rows;y++){
    progress.update(ftell(fin)); //Todo: Check to see if ftell fails here?
    for(int x=0;x<columns;x++){
      if (fscanf(fin,"%f", &temp)!=1){
        diagnostic("\n\tFailed! (Couldn't read or convert a value!)\n");
        exit(-1);
      }
      elevations(x,y)=temp;
      if(temp!=elevations.no_data)
        elevations.data_cells++;
    }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());

  fclose(fin);

  diagnostic_arg(
    "Read %ld cells, of which %ld contained data (%ld%%).\n",
    elevations.width()*elevations.height(), elevations.data_cells,
    elevations.data_cells*100/elevations.width()/elevations.height()
  );

  load_time.stop();
  diagnostic_arg("Read time was: %lfs\n", load_time.accumulated());

  return 0;
}








int load_ascii_data(const char filename[], char_2d &data){
  FILE *fin;
  long long file_size;
  int rows,columns;
  Timer load_time;
  ProgressBar progress;
  int temp;

  load_time.start();

  errno=0;
  diagnostic_arg("Opening input ASCII-DEM file \"%s\"...",filename);
  fin=fopen(filename,"r");
  if(fin==NULL){
    diagnostic_arg("failed with error %d: \"%s\"!\n",errno,strerror(errno));
    exit(-1);
  }
  diagnostic("succeeded.\n");

  diagnostic("Calculating file size...");
  if(fseek(fin,-1,SEEK_END)!=0){
    diagnostic("failed! (Couldn't jump to end of file.)\n");
    exit(-1);
  }
  if((file_size=ftell(fin))==-1){
    diagnostic("failed! (Couldn't determine file size.)\n");
    exit(-1);
  }
  if(fseek(fin,0,SEEK_SET)!=0){
    diagnostic("failed! (Couldn't jump back to beginning of file.)\n");
    exit(-1);
  }
  diagnostic("succeeded.\n");

//  posix_fadvise(fileno(fin),0,0,POSIX_FADV_SEQUENTIAL);

  diagnostic("Reading DEM header...");
//  if(fscanf(fin,"ncols %d nrows %d xllcorner %lf yllcorner %lf cellsize %lf NODATA_value %d",&columns, &rows, &data.xllcorner, &data.yllcorner, &data.cellsize, &temp)!=6){
  if(fscanf(fin,"ncols %d nrows %d xllcorner %lf yllcorner %lf cellsize %lf",&columns, &rows, &data.xllcorner, &data.yllcorner, &data.cellsize)!=5){
    diagnostic("failed!\n");
    exit(-1);
  }
  diagnostic("succeeded.\n");

//  data.no_data=(char)temp;

  diagnostic_arg(
    "The loaded DEM will require approximately %ldMB of RAM.\n",
    columns*rows*((long)sizeof(float))/1024/1024
  );

  diagnostic("Resizing data matrix...");
  data.resize(columns,rows);
  diagnostic("succeeded.\n");

  diagnostic("%%Reading data matrix...\n");
  progress.start(file_size);

  data.data_cells=0;
  for(int y=0;y<rows;y++){
    progress.update(ftell(fin)); //Todo: Check to see if ftell fails here?
    for(int x=0;x<columns;x++){
      if (fscanf(fin,"%d", &temp)!=1){
        diagnostic("\n\tFailed! (Couldn't read or convert a value!)\n");
        exit(-1);
      }
      data(x,y)=(char)temp;
      if(temp!=data.no_data)
        data.data_cells++;
    }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());

  fclose(fin);

  diagnostic_arg(
    "Read %ld cells, of which %ld contained data (%ld%%).\n",
    data.width()*data.height(), data.data_cells,
    data.data_cells*100/data.width()/data.height()
  );

  load_time.stop();
  diagnostic_arg("Read time was: %lfs\n", load_time.accumulated());

  return 0;
}


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
