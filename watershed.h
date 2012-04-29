#ifndef _include_watersheds
#define _include_watersheds

void find_watersheds(float_2d &elevations, int_2d &labels);
void watershed_area(const int_2d &labels);

#endif
