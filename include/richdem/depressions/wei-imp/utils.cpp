#include "utils.h"
#include "dem.h"
#include <string>
#include <math.h>
#include <time.h>


/*
*	neighbor index
*	5  6  7
*	4     0
*	3  2  1
*/
int	ix[8]	= { 0,1, 1, 1, 0,-1,-1,-1 };
int	iy[8]	= { 1,1, 0,-1,-1,-1, 0, 1 };

