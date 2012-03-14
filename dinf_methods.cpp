#include "data_structures.h"
#include "utility.h"
#include "dinf_methods.h"
#include "interface.h"
#include <omp.h>
#include <cmath>
#include <queue>

#define NO_DINF_FLOW -1

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

	if (elevations(x,y)==elevations.no_data) return dinf_NO_DATA; //Missing data
	if (EDGE_GRID(x,y,elevations.width(),elevations.height())) return -1; //Edge cells do not have a defined flow direction

	for(int n=0;n<8;n++){
		if(!IN_GRID(x+dx_e1[n],y+dy_e1[n],elevations.width(),elevations.height())) continue;
		if(!IN_GRID(x+dx_e2[n],y+dy_e2[n],elevations.width(),elevations.height())) continue;
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

void dinf_flow_directions(const float_2d &elevations, float_2d &flowdirs, bool init){
	diagnostic_arg("The Dinf flow directions will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(float)/1024/1024);
	diagnostic("Resizing flow directions matrix...");
	flowdirs.resize(elevations.width(),elevations.height(),!init);
	diagnostic("succeeded.\n");

	if(init){
		diagnostic("Setting no_data value on flowdirs matrix...");
		flowdirs.no_data=dinf_NO_DATA;
		diagnostic("succeeded.\n");

		diagnostic("Initializing Dinf flow directions...\n");
		flowdirs.init(NO_FLOW);
		diagnostic("\tsucceeded.\n");
	}

	diagnostic("Calculating Dinf flow directions...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++){
		progress_bar(x*omp_get_num_threads()*elevations.height()*100/(elevations.width()*elevations.height()));
		for(int y=0;y<elevations.height();y++)
			if(flowdirs(x,y)==NO_FLOW)
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

void dinf_upslope_area(const float_2d &flowdirs, float_2d &area){
	char_2d dependency;
	std::queue<grid_cell> sources;

	diagnostic_arg("The sources queue will require at most approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(grid_cell)/1024/1024);

	diagnostic_arg("The dependency matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(char)/1024/1024);
	diagnostic("Resizing dependency matrix...");
	dependency.resize(flowdirs.width(),flowdirs.height(),false);
	diagnostic("succeeded.\n");

	diagnostic_arg("The area matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(float)/1024/1024);
	diagnostic("Resizing the area matrix...");
	area.resize(flowdirs.width(),flowdirs.height(),false);
	diagnostic("succeeded.\n");

	diagnostic("Setting no_data value on area matrix...");
	area.no_data=dinf_NO_DATA;
	diagnostic("succeeded.\n");

	diagnostic("Calculating dependency matrix & setting no_data cells...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<flowdirs.width();x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.height()*100/(flowdirs.width()*flowdirs.height()));
		for(int y=0;y<flowdirs.height();y++){
			dependency(x,y)=0;
			if(flowdirs(x,y)==dinf_NO_DATA){
				area(x,y)=dinf_NO_DATA;
				continue;
			}
			for(int n=1;n<=8;n++){
				if(!IN_GRID(x+dx[n],y+dy[n],flowdirs.width(),flowdirs.height())) continue;
				if(does_cell_flow_into_me(flowdirs(x+dx[n],y+dy[n]),n))
					dependency(x,y)++;
			}
		}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	diagnostic("Locating source cells...\n");
	progress_bar(-1);
	for(int x=0;x<flowdirs.width();x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.height()*100/(flowdirs.width()*flowdirs.height()));
		for(int y=0;y<flowdirs.height();y++)
			if(flowdirs(x,y)!=dinf_NO_DATA && dependency(x,y)==0)
				sources.push(grid_cell(x,y));
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	diagnostic("Calculating up-slope areas...\n");
	progress_bar(-1);
	long int ccount=0;
	while(sources.size()>0){
		grid_cell *c=&sources.front();

		ccount++;
		progress_bar(ccount*100/flowdirs.data_cells);

		if(flowdirs(c->x,c->y)==dinf_NO_DATA){ //TODO: May not be necessary due to code elsewhere
			area(c->x,c->y)=dinf_NO_DATA;
			sources.pop();
			continue;
		}

		area(c->x,c->y)=1;
		for(int n=1;n<=8;n++){
			if(!IN_GRID(c->x+dx[n],c->y+dy[n],flowdirs.width(),flowdirs.height())) continue;
			if(does_cell_flow_into_me(flowdirs(c->x+dx[n],c->y+dy[n]),n))
				area(c->x,c->y)+=proportion_i_get(flowdirs(c->x+dx[n],c->y+dy[n]),n)*area(c->x+dx[n],c->y+dy[n]);
		}
		if(flowdirs(c->x,c->y)!=dinf_NO_DATA && flowdirs(c->x,c->y)!=NO_DINF_FLOW){
			int n_low= dinf_to_d8_low [(int)(flowdirs(c->x,c->y)/(M_PI/4))];
			int n_high=dinf_to_d8_high[(int)(flowdirs(c->x,c->y)/(M_PI/4))];
			if( (--dependency(c->x+dx[n_low],c->y+dy[n_low]))==0)
				sources.push(grid_cell(c->x+dx[n_low],c->y+dy[n_low]));
			if( (--dependency(c->x+dx[n_high],c->y+dy[n_high]))==0)
				sources.push(grid_cell(c->x+dx[n_high],c->y+dy[n_high]));
		}

		sources.pop();
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}
