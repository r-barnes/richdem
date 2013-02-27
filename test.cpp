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

  load_ascii_data(argv[1],elevations);
  elevations.low_pass_filter();
  barnes_flood_flowdirs(elevations, flowdirs);
  d8_upslope_cells(atoi(argv[3])-elevations.xllcorner, (elevations.yllcorner+elevations.height())-atoi(argv[4]), atoi(argv[5])-elevations.xllcorner, (elevations.yllcorner+elevations.height())-atoi(argv[6]), flowdirs, is_upslope);

  output_ascii_data(argv[2], is_upslope);

  return 0;
}
