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

class Flag
{
public:
	int width, height;
	unsigned char* flagArray;
public:
	~Flag()
	{
		Free();
	}
	bool Init(int width,int height)
	{
		flagArray=NULL;
		this->width=width;
		this->height=height;
		flagArray = new unsigned char[(width * height + 7) / 8]();
		return flagArray!=NULL;
	}
	void Free()
	{
		delete[] flagArray;
		flagArray=NULL;
	}
	void SetFlag(int row, int col)
	{
		int index=row*width+col;
		flagArray[index / 8] |= value[index % 8];
	}
	void SetFlags(int row, int col, Flag& flag)
	{
		int index=row*width+col;
		int bIndex=index / 8;
		int bShift=index % 8;
		flagArray[bIndex] |= value[bShift];
		flag.flagArray[bIndex] |= value[bShift];
	}
	int IsProcessed(int row, int col)
	{
		//if the cell is outside the DEM, is is regared as processed
		if (row<0 ||row>=height || col<0|| col>=width) return true;
		int index=row*width+col;
		return flagArray[index / 8]&value[index % 8];
	}
	int IsProcessedDirect(int row, int col)
	{
		int index=row*width+col;
		return flagArray[index / 8]&value[index % 8];
	}
};

#endif
