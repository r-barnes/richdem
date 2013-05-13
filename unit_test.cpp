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
#include "unit_test.hpp"

int main(){
  float_2d elevations;
  load_ascii_data("unit_test/bf03.dem", elevations);

  {
    float_2d slope_percent, unit_slope_percent;
    d8_slope(elevations, slope_percent, TATTRIB_SLOPE_PERCENT);
    load_ascii_data("unit_test/bf03_slope_percent.txt", unit_slope_percent);
    printf("Average SLOPE PERCENT difference: %lf\n",unit_avg_diff(slope_percent, unit_slope_percent));
    slope_percent.clear();
    unit_slope_percent.clear();
  }

  {
    float_2d aspect, unit_aspect;
    d8_aspect(elevations, aspect);
    load_ascii_data("unit_test/bf03_aspect.txt", unit_aspect);
    printf("Average ASPECT difference: %lf\n",unit_ang_avg_diff(aspect, unit_aspect));
    aspect.clear();
    unit_aspect.clear();
  }

  {
    float_2d planform_curvature, unit_planform_curvature;
    d8_planform_curvature(elevations, planform_curvature);
    load_ascii_data("unit_test/bf03_planform_curvature.txt", unit_planform_curvature);
    printf("Average PLANFORM CURVATURE difference: %lf\n",unit_avg_diff(planform_curvature, unit_planform_curvature));
    planform_curvature.clear();
    unit_planform_curvature.clear();
  }

  {
    float_2d profile_curvature, unit_profile_curvature;
    d8_profile_curvature(elevations, profile_curvature);
    load_ascii_data("unit_test/bf03_profile_curvature.txt", unit_profile_curvature);
    printf("Average PROFILE CURVATURE difference: %lf\n",unit_avg_diff(profile_curvature, unit_profile_curvature));
    profile_curvature.clear();
    unit_profile_curvature.clear();
  }

}
