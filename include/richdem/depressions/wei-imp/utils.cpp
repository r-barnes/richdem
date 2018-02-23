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

int randomi(int m,int n){
	int pos, dis;
	srand((unsigned)time(NULL));
	if (m == n){
		return m;
	}
	else if (m > n){
		pos = n;
		dis = m - n + 1;
		return rand() % dis + pos;
	}
	else{
		pos = m;
		dis = n - m + 1;
		return rand() % dis + pos;
	}
}
