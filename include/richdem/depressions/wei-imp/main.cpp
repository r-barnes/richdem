#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <queue>
#include <algorithm>
#include "dem.h"
#include "utils.h"
#include <time.h>
#include <list>
#include <unordered_map>
#include "perlin.h"
#include "filldem.hpp"

#include <richdem/richdem.hpp>
#include <richdem/common/Array2D.hpp>
#include <stdexcept>


using namespace std;
using std::cout;
using std::endl;
using std::string;
using std::getline;
using std::fstream;
using std::ifstream;
using std::priority_queue;
using std::binary_function;




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
	const int SIZE=1000;

	auto dem = createPerlinNoiseDEM(SIZE,SIZE);

	auto dem2 = dem;



	// richdem::Array2D<float> rddem(SIZE,SIZE);
	// for(int y=0;y<SIZE;y++)
	// for(int x=0;x<SIZE;x++)
	// 	rddem(x,y) = dem.asFloat(y,x);

  dem.saveGDAL("/z/before.tif");

  richdem::fillDEM(dem);

	// richdem::Array2D<float> outdem(SIZE,SIZE);
	// for(int y=0;y<SIZE;y++)
	// for(int x=0;x<SIZE;x++)
	// 	outdem(x,y) = dem.asFloat(y,x);




  std::cout<<"Rich"<<std::endl;


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


