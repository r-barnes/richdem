#ifndef _d8_methods_included
#define _d8_methods_included

#include "data_structures.h"

#define SLOPE_RISERUN	1
#define SLOPE_PERCENT	2
#define SLOPE_RADIAN	3
#define SLOPE_DEGREE	4

void d8_flow_directions(const float_2d &elevations, char_2d &flowdirs, bool init=true);
void d8_upslope_area(const char_2d &flowdirs, uint_2d &area);
void d8_slope(const float_2d &elevations, float_2d &slopes, int slope_type=SLOPE_RISERUN);

#endif
