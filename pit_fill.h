#ifndef pit_fill_include
#define pit_fill_include
#include "data_structures.h"

void barnes_flood(float_2d &elevations);
void barnes_flood_flowdirs(const float_2d &elevations, char_2d &flowdirs);

#endif
