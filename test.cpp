#include "utility.hpp"
#include "data_structures.hpp"
#include "data_io.hpp"
#include "d8_methods.hpp"
#include "dinf_methods.hpp"
#include "pit_fill.hpp"
#include "interface.hpp"
#include "flat_resolution.hpp"
#include "debug.hpp"
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
