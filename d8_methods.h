#ifndef _d8_methods_included
#define _d8_methods_included

#include "data_structures.h"

int d8_flow_directions(const float_2d &elevations, char_2d &flowdirs, bool init=true);
int d8_upslope_area(const char_2d &flowdirs, uint_2d &area);

#endif
