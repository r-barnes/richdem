#include "utility.h"
#include "data_structures.h"
#include "data_io.h"
#include "d8_methods.h"
#include "dinf_methods.h"
#include "pit_fill.h"
#include "interface.h"
#include "flat_resolution.h"
#include "debug.h"
#include <string>
#include <sys/time.h>

int main(int argc, char **argv){
  float_2d elevations;
  char_2d flowdirs;
  int_2d is_upslope;
  uint_2d fmask;

  load_ascii_data(argv[1],elevations);
  elevations.low_pass_filter();
  barnes_flood(elevations);

  output_ascii_data(argv[2], elevations);

  return 0;

  d8_flow_directions(elevations,flowdirs);
  flat_mask(flowdirs,fmask);

  output_ascii_data(argv[2], fmask);

  int tflats=0;
  for(int x=0;x<fmask.width();++x)
  for(int y=0;y<fmask.height();++y)
    if(fmask(x,y)==1)
      tflats++;
  printf("%d in flats\n",tflats);

  return 0;
}
