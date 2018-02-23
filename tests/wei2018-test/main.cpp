#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <queue>
#include <algorithm>
#include <time.h>
#include <list>
#include <unordered_map>
#include "perlin.h"

#include <richdem/depressions/depressions.hpp>
#include <richdem/common/Array2D.hpp>
#include <stdexcept>


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



//generate DEM with Perlin Noise 
richdem::Array2D<double> createPerlinNoiseDEM(int width,int height){
	richdem::Array2D<double> dem(width,height);
	float fre = (float)randomi(0,512)/10;
	for (int y = 0; y < height; ++y)
	for (int x = 0; x < width; ++x){
		float vec2[] = { x / fre , y/ fre };
		dem(x,y)     = 1000*noise2(vec2);
	}

	return dem;
}

int main(int argc, char* argv[]){
	const int SIZE = 5000;

	auto dem = createPerlinNoiseDEM(SIZE,SIZE);

	auto dem2 = dem;

  dem.saveGDAL("/z/before.tif");

  richdem::PriorityFlood_Wei2018(dem);

  richdem::FillDepressions(dem2);

	dem.saveGDAL("/z/after_wei.tif");
	dem2.saveGDAL("/z/after_richdem.tif");

	dem2.setNoData(dem.noData());

	double diff = 0;
	for(int i=0;i<dem.size();i++)
		diff += std::abs(dem(i)-dem2(i));

	std::cout<<"Diff: "<<std::setprecision(20)<<diff<<std::endl;

  std::cout<<"Comparing "<<(int)(dem==dem2)<<std::endl;

	return 0;
}


