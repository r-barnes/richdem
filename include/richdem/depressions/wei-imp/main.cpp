#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <queue>
#include <algorithm>
#include "dem.h"
#include "Node.h"
#include "utils.h"
#include <time.h>
#include <list>
#include <unordered_map>
#include "perlin.h"



using namespace std;
using std::cout;
using std::endl;
using std::string;
using std::getline;
using std::fstream;
using std::ifstream;
using std::priority_queue;
using std::binary_function;


typedef std::vector<Node> NodeVector;
typedef std::priority_queue<Node, NodeVector, Node::Greater> PriorityQueue;


void fillDEM(char* inputFile, char* outputFilledPath);
//compute stats for a DEM
void calculateStatistics(const CDEM& dem, double* min, double* max, double* mean, double* stdDev)
{
	int width = dem.Get_NX();
	int height = dem.Get_NY();

	int validElements = 0;
	double minValue, maxValue;
	double sum = 0.0;
	double sumSqurVal = 0.0;
	for (int row = 0; row < height; row++)
	{
		for (int col = 0; col < width; col++)
		{
			if (!dem.is_NoData(row, col))
			{
				double value = dem.asFloat(row, col);
				
				if (validElements == 0)
				{
					minValue = maxValue = value;
				}
				validElements++;
				if (minValue > value)
				{
					minValue = value;
				}
				if (maxValue < value)
				{
					maxValue = value;
				}

				sum += value;
				sumSqurVal += (value * value);
			}
		}
	}

	double meanValue = sum / validElements;
	double stdDevValue = sqrt((sumSqurVal / validElements) - (meanValue * meanValue));
	*min = minValue;
	*max = maxValue;
	*mean = meanValue;
	*stdDev = stdDevValue;
}

//generate DEM with Perlin Noise 
void createPerlinNoiseDEM(char* outputFilePath,int rows,int cols){
	CDEM dem;
	float fre=(float)randomi(0,512)/10;
	dem.SetHeight(rows);
	dem.SetWidth(cols);
	if (!dem.Allocate()) {
		std::cout << "Failed to allocate memory correctly!" << std::endl;
		return;
	}
	for (int x = 0; x < rows; ++x){
		for (int y = 0; y < cols; ++y){
			float vec2[] = { x / fre , y/ fre };
			dem.Set_Value(x,y,noise2(vec2)*1000);
		}
	}
	double min, max, mean, stdDev;
	calculateStatistics(dem, &min, &max, &mean, &stdDev);
	CreateGeoTIFF(outputFilePath, dem.Get_NY(), dem.Get_NX(), 
		(void *)dem.getDEMdata(),GDALDataType::GDT_Float32, NULL,
		&min, &max, &mean, &stdDev, -9999);
	std::cout << "create dem successfully!" << std::endl;
    return;
}

int main(int argc, char* argv[])
{
	if (argc<4){
		cout<<"\nThis program fills depressions in DEMs using our proposed algorithm. It can also generate DEMs with large depressions using Perlin noise method."<<endl;
		cout<<"\n\nTo create a DEM with Perlin noise method, use the following command:\nfill perlinDEM rows cols dem.tif \nThis command generates a Perlin noise DEM with given rows and cols. "<<endl;
		cout<<"\nExample usage:, ./fill perlinDEM 100 100 dem.tif"<<endl;
		cout<<"\n\nTo fill a DEM using our proposed method, use the following command:\nfill fillDEM inputFilePath outputfilePath.\nThis command fills depressions in the given DEM"<<endl;
		cout<<"\nExample usage: ./fill fillDEM dem.tif demfilled.tif"<<endl;
		return 1;
	}
	char* function=argv[1];
	if (strcmp(function,"perlinDEM")==0){
		string sr(argv[2]);
		string sc(argv[3]);
		if (argc<5||atoi(&sr[0])<0||atoi(&sc[0])<0){
			cout<<"Unknown parameters,please input correct parameters!"<<endl;
			return 1;
		}
		string outputFilePath(argv[4]);
		createPerlinNoiseDEM(&outputFilePath[0],atoi(&sr[0]),atoi(&sc[0]));
	}
	else if(strcmp(function,"fillDEM")==0){
		string inputFilePath(argv[2]);
		string outputFilePath(argv[3]);
		fillDEM(&inputFilePath[0], &outputFilePath[0]);
	}
	else {
		cout<<"Unknown parameters,please input correct parameters!"<<endl;
		return 1;
	}
	return 0;
}


