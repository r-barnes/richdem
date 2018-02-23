#ifndef UTILS_HEAD_H
#define UTILS_HEAD_H

#include <queue>
#include <algorithm>
#include "dem.h"


void calculateStatistics(const CDEM& dem, double* min, double* max, double* mean, double* stdDev);

extern int	ix[8];
extern int	iy[8];

inline int Get_rowTo(int dir, int row)
{
	return( row + ix[dir] );
}
inline int Get_colTo(int dir, int col){
	return( col + iy[dir] );
}



void setNoData(unsigned char* data, int length, unsigned char noDataValue);
void setNoData(float* data, int length, float noDataValue);
void setFlag(int index, unsigned char* flagArray);
bool isProcessed(int index, const unsigned char* flagArray);
CDEM* diff(CDEM& demA, CDEM& demB);
float randomf(float m,float n);
int randomi(int m,int n);
extern const unsigned char value[8];

#endif
