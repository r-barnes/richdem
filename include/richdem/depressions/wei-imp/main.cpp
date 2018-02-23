#include <stdio.h>
#include <iomanip>
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

#include <richdem/richdem.hpp>
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


typedef std::vector<Node> NodeVector;
typedef std::priority_queue<Node, NodeVector, Node::Greater> PriorityQueue;


void fillDEM(CDEM &dem);

//generate DEM with Perlin Noise 
CDEM createPerlinNoiseDEM(int rows,int cols){
	CDEM dem;
	float fre=(float)randomi(0,512)/10;
	dem.SetHeight(rows);
	dem.SetWidth(cols);
	if (!dem.Allocate()) {
		throw std::runtime_error("Failed to allocate memory correctly!");
	}
	for (int x = 0; x < rows; ++x){
		for (int y = 0; y < cols; ++y){
			float vec2[] = { x / fre , y/ fre };
			dem.Set_Value(x,y,noise2(vec2)*1000);
		}
	}

	return dem;
}

int main(int argc, char* argv[]){
	const int SIZE=1000;

	auto dem = createPerlinNoiseDEM(SIZE,SIZE);

	richdem::Array2D<float> rddem(SIZE,SIZE);
	for(int y=0;y<SIZE;y++)
	for(int x=0;x<SIZE;x++)
		rddem(x,y) = dem.asFloat(y,x);


  fillDEM(dem);

	richdem::Array2D<float> outdem(SIZE,SIZE);
	for(int y=0;y<SIZE;y++)
	for(int x=0;x<SIZE;x++)
		outdem(x,y) = dem.asFloat(y,x);




  std::cout<<"Rich"<<std::endl;

  rddem.saveGDAL("/z/before.tif");

  richdem::FillDepressions(rddem);

	rddem.saveGDAL("/z/after_richdem.tif");
	outdem.saveGDAL("/z/after_wei.tif");

	outdem.setNoData(rddem.noData());

	double diff = 0;
	for(int i=0;i<rddem.size();i++)
		diff += std::abs(outdem(i)-rddem(i));

	std::cout<<"Diff: "<<std::setprecision(20)<<diff<<std::endl;

  std::cout<<"Comparing "<<(int)(outdem==rddem)<<std::endl;

	return 0;
}


