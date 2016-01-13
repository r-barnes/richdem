#include "utility.hpp"
#include "data_structures.hpp"
#include "data_io.hpp"
#include "d8_methods.hpp"
#include "dinf_methods.hpp"
#include "interface.hpp"
#include "flat_resolution.hpp"
#include <string>
#include <sys/time.h>

int main(int argc, char **argv){
  float_2d elevations;

  //Read in ArcGrid ASCII for Floating Raster Grid formatted files
  read_data(argv[1], elevations);

  char_2d flowdirs(elevations);
  d8_flow_directions(elevations,flowdirs);

  int_2d flat_resolution_mask(elevations), groups;
  resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);
  d8_flow_flats(flat_resolution_mask,groups,flowdirs);


  return 0;
}
