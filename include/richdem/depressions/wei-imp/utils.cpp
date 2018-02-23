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


void setNoData(unsigned char* data, int length, unsigned char noDataValue)
{
	if (data == NULL || length == 0)
	{
		return;
	}

	for (int i = 0; i < length; i++)
	{
		data[i] = noDataValue;
	}
}
void setNoData(float* data, int length, float noDataValue)
{
	for (int i = 0; i < length; i++)
	{
		data[i] = noDataValue;
	}
}


const unsigned char value[8] = {128, 64, 32, 16, 8, 4, 2, 1};


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
