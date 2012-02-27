#include "data_structures.h"
#include "utility.h"
#include "d8_methods.h"
#include "interface.h"
#include <omp.h>
#include <cmath>
#include <queue>

#define NO_DINF_FLOW -1

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


//Table 1 of Tarboton
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
	if (elevations(x,y)==elevations.no_data) return dinf_NO_DATA; //Missing data

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

	double rg=NO_DINF_FLOW;
	if(nmax!=-1)
		rg=(af[nmax]*rmax+ac[nmax]*M_PI/2);

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


/*
We must convert the Dinf angle system to cells within the D8 system
I use the following grid for the D8 system
//234
//105
//876
To convert Dinf to this, take
(int)(flowdir/45)      D8
      0                4,5
      1                3,4
      2                2,3
      3                1,2
      4                1,8
      5                7,8
      6                6,7
      7                5,6
*/
//These arrays have a 9th element which repeats the 8th element because floating point rounding errors occassionally result in the 9th element being accessed.
int dinf_to_d8_low[9]= {4,3,2,1,1,7,6,5,5};
int dinf_to_d8_high[9]={5,4,3,2,8,8,7,6,6};

bool does_cell_flow_into_me(float n_flowdir, int n){
	if(n_flowdir==NO_DINF_FLOW || n_flowdir==dinf_NO_DATA) return false;
	n=inverse_flow[n];
	int flown=(int)(n_flowdir/(M_PI/4));
	return (dinf_to_d8_low[flown]==n || dinf_to_d8_high[flown]==n);
}

//This reacts correctly if the flow direction wedge number exceeds 7.
float proportion_i_get(float flowdir, int n){
//	if(!does_cell_flow_into_me(flowdir,n)) return 0;
	int dinf_wedge=(int)(flowdir/(M_PI/4));
	float normalized_angle=flowdir-(M_PI/4)*dinf_wedge;

	if(n%2==dinf_wedge%2)
		return normalized_angle/(M_PI/4);
	else
		return 1-normalized_angle/(M_PI/4);
}

int dinf_upslope_area(const float_2d &flowdirs){
	char_2d dependency;
	float_2d area;
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

	diagnostic_arg("The area matrix will require approximately %ldMB of RAM.\n",flowdirs.size1()*flowdirs.size2()*sizeof(float)/1024/1024);
	diagnostic("Resizing the area matrix...");
	try{
		area.resize(flowdirs.size1(),flowdirs.size2(),false);
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
				if(does_cell_flow_into_me(flowdirs(x+dx[n],y+dy[n]),n))
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
			if(flowdirs(x,y)!=dinf_NO_DATA && dependency(x,y)==0)
				sources.push(new grid_cell(x,y));
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	diagnostic("Calculating up-slope areas...\n");
	progress_bar(-1);
	long int ccount=0;
	while(sources.size()>0){
		grid_cell *c=sources.front();
		sources.pop();

		ccount++;
		progress_bar(ccount*100/flowdirs.data_cells);

		area(c->x,c->y)=1;
		for(int n=1;n<=8;n++){
			if(!IN_GRID(c->x+dx[n],c->y+dy[n],flowdirs.size1(),flowdirs.size2())) continue;
			if(does_cell_flow_into_me(flowdirs(c->x+dx[n],c->y+dy[n]),n))
				area(c->x,c->y)+=proportion_i_get(flowdirs(c->x+dx[n],c->y+dy[n]),n)*area(c->x+dx[n],c->y+dy[n]);
		}
		if(flowdirs(c->x,c->y)!=dinf_NO_DATA && flowdirs(c->x,c->y)!=NO_DINF_FLOW){
			int n_low= dinf_to_d8_low [(int)(flowdirs(c->x,c->y)/(M_PI/4))];
			int n_high=dinf_to_d8_high[(int)(flowdirs(c->x,c->y)/(M_PI/4))];
			if( (--dependency(c->x+dx[n_low],c->y+dy[n_low]))==0)
				sources.push(new grid_cell(c->x+dx[n_low],c->y+dy[n_low]));
			if( (--dependency(c->x+dx[n_high],c->y+dy[n_high]))==0)
				sources.push(new grid_cell(c->x+dx[n_high],c->y+dy[n_high]));
		}

		delete c;
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}
