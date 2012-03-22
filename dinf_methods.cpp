#include "data_structures.h"
#include "utility.h"
#include "dinf_methods.h"
#include "interface.h"
#include <omp.h>
#include <cmath>
#include <queue>

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
const int dinf_to_d8_low[9]= {4,3,2,1,1,7,6,5,5};
const int dinf_to_d8_high[9]={5,4,3,2,8,8,7,6,6};

bool does_cell_flow_into_me(float n_flowdir, int n){
	if(n_flowdir==NO_FLOW || n_flowdir==dinf_NO_DATA) return false;
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
		if(flowdirs(c->x,c->y)!=dinf_NO_DATA && flowdirs(c->x,c->y)!=NO_FLOW){
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
