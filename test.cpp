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
  float_2d elevations,elevations2;

  load_ascii_data(argv[1],elevations);
  write_floating_data(argv[2],elevations);
  read_floating_data(argv[2],elevations);
  output_ascii_data(argv[3],elevations,8);

  return 0;
}
