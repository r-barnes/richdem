#include <iostream>
#include <fstream>
#include <cstdlib>

#define XBIG 10000
#define YBIG 10000
#define IN_GRID(x,y)   (x>=0 && y>=0 && x<x_max && y<y_max)
#define EDGE_GRID(x,y) (x==0 || y==0 || x==x_max-1 || y==y_max-1)
using namespace std;
char elevations[XBIG][YBIG];
int x_max,y_max;

int PrintDEM(){
	ofstream fout("square_grid.dem");
	fout<<"ncols\t\t"<<x_max<<endl;
	fout<<"nrows\t\t"<<y_max<<endl;
	fout<<"xllcorner\t\t0"<<endl;
	fout<<"yllcorner\t\t0"<<endl;
	fout<<"cellsize\t\t3"<<endl;
	fout<<"NODATA_value\t\t-1"<<endl;
	for(int y=0;y<y_max;y++){
		for(int x=0;x<x_max;x++)
			fout<<(int)elevations[x][y]<<" ";
		fout<<endl;
	}
	fout.close();
}

int GenerateDEM(){
	for(int x=0;x<x_max;x++)
	for(int y=0;y<y_max;y++)
		if(EDGE_GRID(x,y))
			elevations[x][y]=9;
		else
			elevations[x][y]=6;

	elevations[2][y_max-1]=5;
}

int main(int argc, char **argv){
	if(argc!=2){
		printf("%s <LENGTH OF ONE SIDE OF GRID>\n",argv[0]);
		return -1;
	}
	x_max=y_max=atoi(argv[1]);
	GenerateDEM();	
	PrintDEM();
}
