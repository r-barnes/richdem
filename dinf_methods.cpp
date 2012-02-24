#include "data_structures.h"
#include "utility.h"
#include "d8_methods.h"
#include "interface.h"
#include <omp.h>
#include <cmath>
#include <queue>

//D8 Directions
static int const dx[9]={0,-1,-1,0,1,1,1,0,-1};
static int const dy[9]={0,0,-1,-1,-1,0,1,1,1};
static int const inverse_flow[9]={0,5,6,7,8,1,2,3,4};

/*
814 citations
@article{tarboton1997new,
  title={A new method for the determination of flow directions and upslope areas in grid digital elevation models},
  author={Tarboton, D.G.},
  journal={Water resources research},
  volume={33},
  pages={309--319},
  year={1997},
  publisher={Citeseer}
}
*/


//Table 1 of Tarboton TODO (i=y-coordinate, j=x-coordinate)
int dy_e1[8]={0,-1,-1,0,0,1,1,0};
int dx_e1[8]={1,0,0,-1,-1,0,0,1};
int dy_e2[8]={-1,-1,-1,-1,1,1,1,1};
int dx_e2[8]={1,1,-1,-1,-1,-1,1,1};
double ac[8]={0,1,1,2,2,3,3,4};
double af[8]={1,-1,1,-1,1,-1,1,-1};

float dinf_FlowDir(const float_2d &elevations, const int x, const int y){
	double smax=0;
	int nmax=-1;
	double rmax=0;

	double e0,e1,e2,d1,d2,s1,s2,r,s;

	if (EDGE_GRID(x,y,elevations.size1(),elevations.size2())) return -1; //Edge cells do not have a defined flow direction
	if (elevations(x,y)==elev_NO_DATA) return dinf_NO_DATA; //Missing data

	for(int n=0;n<8;n++){
		if(!IN_GRID(x+dx_e1[n],y+dy_e1[n],elevations.size1(),elevations.size2())) continue;
		if(!IN_GRID(x+dx_e2[n],y+dy_e2[n],elevations.size1(),elevations.size2())) continue;
		if(elevations(x+dx_e1[n],y+dy_e1[n])==elev_NO_DATA) continue;
		if(elevations(x+dx_e2[n],y+dy_e2[n])==elev_NO_DATA) continue;

		e0=elevations(x,y);
		e1=elevations(x+dx_e1[n],y+dy_e1[n]);
		e2=elevations(x+dx_e2[n],y+dy_e2[n]);
		d1=1;
		d2=1;
		s1=(e0-e1)/d1;
		s2=(e1-e2)/d2;
		r=atan2(s2,s1);
		s=sqrt(s1*s1+s2*s2);
		if(r<0){
			r=0;
			s=s1;
		} else if(r>atan2(d2,d1)){
			r=atan2(d2,d1);
			s=(e0-e2)/sqrt(d1*d1+d2*d2);
		}
		if(s>smax){
			smax=s;
			nmax=n;
			rmax=r;
		}
	}

	double rg=-1;
	if(nmax!=-1)
		rg=(af[nmax]*rmax+ac[nmax]*M_PI/2);//*(180.0/M_PI); //TODO: Remove degree conversion

	return rg;
}

int dinf_flow_directions(const float_2d &elevations, float_2d &flowdirs){
	diagnostic_arg("The Dinf flow directions will require approximately %ldMB of RAM.\n",elevations.size1()*elevations.size2()*sizeof(float)/1024/1024);
	diagnostic("Resizing flow directions matrix...");
	try{
		flowdirs.resize(elevations.size1(),elevations.size2(),false);
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Calculating Dinf flow directions...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.size1();x++){
		progress_bar(x*omp_get_num_threads()*elevations.size2()*100/(elevations.size1()*elevations.size2()));
		for(int y=0;y<elevations.size2();y++)
			flowdirs(x,y)=dinf_FlowDir(elevations,x,y);
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}





/*
@inproceedings{Wallis2009,
address = {Las Vegas, Nevada, USA},
author = {Wallis, Chase and Watson, Dan and Tarboton, David and Wallace, Robert},
booktitle = {International Conference on Parallel and Distributed Processing Techniques and Applications},
file = {:home/rick/projects/watershed/papers/dinf\_refs/10.1.1.158.2864.pdf:pdf},
pages = {1--5},
title = {{Parallel Flow-Direction and Contributing Area Calculation for Hydrology Analysis in Digital Elevation Models}},
year = {2009}
}

Also:
@article{wallaceparallel,
  title={Parallel Algorithms for Processing Hydrologic Properties from Digital Terrain},
  author={Wallace, RM and Tarboton, DG and Watson, DW and Schreuders, KAT and Tesfa, TK}
}
*/
//234
//105
//876
bool does_cell_flow_into_me(const int x, const int y, int n, const float_2d &flowdirs){ //TODO: Maybe should have some kind of tolerance value - there is a little leakage happening due to rounding when water should be flowing directly into a single cell
	if(n==-1) return false; //TODO: This should not happen.
	if(n==0) { //TODO
		diagnostic("A terrible mistake has occurred in the procedure claling the 'does_cell_flow_into_me' procedure.\n");
		return false;
	}
	n=inverse_flow[n];
	if(     n==1 && flowdirs(x,y)>3*M_PI/4 && flowdirs(x,y)<5*M_PI/4)
		return true;
	else if(n==2 && flowdirs(x,y)>2*M_PI/4 && flowdirs(x,y)<4*M_PI/4)
		return true;
	else if(n==3 && flowdirs(x,y)>1*M_PI/4 && flowdirs(x,y)<3*M_PI/4)
		return true;
	else if(n==4 && flowdirs(x,y)>0*M_PI/4 && flowdirs(x,y)<2*M_PI/4)
		return true;
	else if(n==5 && flowdirs(x,y)>7*M_PI/4 && flowdirs(x,y)<1*M_PI/4)
		return true;
	else if(n==6 && flowdirs(x,y)>6*M_PI/4 && flowdirs(x,y)<8*M_PI/4)
		return true;
	else if(n==7 && flowdirs(x,y)>5*M_PI/4 && flowdirs(x,y)<7*M_PI/4)
		return true;
	else if(n==8 && flowdirs(x,y)>4*M_PI/4 && flowdirs(x,y)<6*M_PI/4)
		return true;
	else
		return false;
}

int dinf_upslope_area(const float_2d &flowdirs){
	char_2d dependency;
	std::queue<grid_cell*> sources;

	diagnostic_arg("The sources queue will require at most approximately %ldMB of RAM.\n",flowdirs.size1()*flowdirs.size2()*sizeof(grid_cell)/1024/1024);

	diagnostic_arg("The dependency matrix will require approximately %ldMB of RAM.\n",flowdirs.size1()*flowdirs.size2()*sizeof(char)/1024/1024);
	diagnostic("Resizing dependency matrix...");
	try{
		dependency.resize(flowdirs.size1(),flowdirs.size2(),false);
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Calculating dependency matrix...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<flowdirs.size1();x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.size2()*100/(flowdirs.size1()*flowdirs.size2()));
		for(int y=0;y<flowdirs.size2();y++){
			dependency(x,y)=0;
			for(int n=1;n<=8;n++){
				if(!IN_GRID(x+dx[n],y+dy[n],flowdirs.size1(),flowdirs.size2())) continue;
				if(does_cell_flow_into_me(x+dx[n],y+dy[n],n,flowdirs))
					dependency(x,y)++;
			}
		}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	diagnostic("Locating source cells...\n");
	progress_bar(-1);
	for(int x=0;x<flowdirs.size1();x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.size2()*100/(flowdirs.size1()*flowdirs.size2()));
		for(int y=0;y<flowdirs.size2();y++)
			if(dependency(x,y)==0)
				sources.push(new grid_cell(x,y));
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}
