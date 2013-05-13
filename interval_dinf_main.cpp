#define  RICHDEM_VERSION  "0r~370"
#include "utility.hpp"
#include "data_structures.hpp"
#include "data_io.hpp"
#include "dinf_methods.hpp"
#include "pit_fill.hpp"
#include "interface.hpp"
#include "flat_resolution.hpp"
#include "interval_dinf.hpp"
#include "debug.hpp"
#include <string>
#include <boost/numeric/interval.hpp>

int main(int argc, char **argv){
  Timer running_calc_time,running_io_time,total_time;

  total_time.start();

  float_2d elevations;
  running_io_time.start();
  load_ascii_data(argv[1],elevations);
  running_io_time.stop();

  running_calc_time.start();
  barnes_flood(elevations);
  running_calc_time.stop();

  float_2d flowdirs(elevations);
  dinf_flow_directions(elevations,flowdirs);

  int_2d flat_resolution_mask(elevations), groups;
  resolve_flats_barnes(elevations,flowdirs,flat_resolution_mask,groups);
  dinf_flow_flats(flat_resolution_mask,groups,flowdirs);
  flat_resolution_mask.clear();
  groups.clear();

  grid_engine< boost::numeric::interval<double> > area;

  running_calc_time.start();
  dinf_upslope_area_interval(flowdirs, area);
  running_calc_time.stop();

  running_calc_time.start();
  float max_width=0,upper_of_max=0,lower_of_max=0;
  //#pragma omp parallel for collapse(2) reduction(max:max_width)
  for(int x=0;x<area.width();x++)
  for(int y=0;y<area.height();y++)
    if(boost::numeric::width(area(x,y))>max_width){
      max_width=boost::numeric::width(area(x,y));
      upper_of_max=boost::numeric::upper(area(x,y));
      lower_of_max=boost::numeric::lower(area(x,y));
    }
  diagnostic_arg("Maximum interval width on %s: %.12lf [%.12lf,%.12lf]\n",argv[1], max_width, lower_of_max, upper_of_max);
  running_calc_time.stop();

  total_time.stop();
  diagnostic_arg("Total time was: %lfs\n", total_time.accumulated());
  diagnostic_arg("Time spent in computation: %lfs\n",running_calc_time.accumulated());
  diagnostic_arg("Time spent in I/O: %lfs\n",running_io_time.accumulated());

  return 0;
}
