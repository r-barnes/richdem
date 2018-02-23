#ifndef UTILS_HEAD_H
#define UTILS_HEAD_H

#include <queue>
#include <algorithm>
#include "dem.h"


extern int	ix[8];
extern int	iy[8];

inline int Get_rowTo(int dir, int row)
{
	return( row + ix[dir] );
}
inline int Get_colTo(int dir, int col){
	return( col + iy[dir] );
}



int randomi(int m,int n);

#endif
