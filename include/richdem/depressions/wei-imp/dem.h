#ifndef CDEM_HEADER_H
#define CDEM_HEADER_H

#include <string>
#include <algorithm>
#include <fstream>
#include <queue>
#include <functional>

#define NO_DATA_VALUE -9999.0f


/*
*	reverse of flow directions
*	2	4	8
*	1	0	16
*	128	64	32
*/
static unsigned char inverse[8] = {16, 32, 64, 128, 1, 2, 4, 8};
/*
*	flow direction		
*	32	64	128		
*	16	0	1		
*	8	4	2		
*/
static unsigned char	dir[8] = {1, 2, 4, 8, 16, 32, 64, 128};
class CDEM
{
protected:
	float* pDem;
	int width, height;
public:
	CDEM()
	{
		pDem=NULL;
	}
	~CDEM()
	{
		delete[] pDem;
	}
	bool Allocate();

	void freeMem();

	void initialElementsNodata();
	float asFloat(int row,int col) const;
	void Set_Value(int row,int col, float z);
	bool is_NoData(int row, int col) const;
	void Assign_NoData();
	int Get_NY() const;
	int Get_NX() const;
	float* getDEMdata() const;
	void SetHeight(int height);
	void SetWidth(int width);
	void readDEM(const std::string& filePath);
	bool is_InGrid(int row, int col) const;
	float getLength(unsigned int dir);
	unsigned char getDirction(int row, int col, float spill);
};
#endif
