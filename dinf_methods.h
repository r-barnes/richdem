#ifndef _dinf_methods_included
#define _dinf_methods_included

void dinf_flow_directions(const float_2d &elevations, float_2d &flowdirs, bool init=true);
void dinf_upslope_area(const float_2d &flowdirs, float_2d &area);

#endif
