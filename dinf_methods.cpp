#include "data_structures.h"
#include "utility.h"
#include "dinf_methods.h"
#include "interface.h"
#include <omp.h>
#include <cmath>
#include <queue>
#include "visualize.h"

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
	area.init(-1);
	area.no_data=dinf_NO_DATA;
	diagnostic("succeeded.\n");

	diagnostic("Calculating dependency matrix & setting no_data cells...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<flowdirs.width();x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.height()*100/(flowdirs.width()*flowdirs.height()));
		for(int y=0;y<flowdirs.height();y++){
			dependency(x,y)=0;
			if(flowdirs(x,y)==flowdirs.no_data){
				area(x,y)=area.no_data;
				continue;
			}
			for(int n=1;n<=8;n++){
				if(!IN_GRID(x+dx[n],y+dy[n],flowdirs.width(),flowdirs.height()))
					continue;
				else if (flowdirs(x+dx[n],y+dy[n])==NO_FLOW)
					continue;
				else if(does_cell_flow_into_me(flowdirs(x+dx[n],y+dy[n]),n))
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
			if(flowdirs(x,y)==flowdirs.no_data)
				continue;
			else if(flowdirs(x,y)==NO_FLOW)
				continue;
			else if(dependency(x,y)==0)
				sources.push(grid_cell(x,y));
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");

	diagnostic("Calculating up-slope areas...\n");
	progress_bar(-1);
	long int ccount=0;
	while(sources.size()>0){
		grid_cell c=sources.front();
		sources.pop();

		ccount++;
		progress_bar(ccount*100/flowdirs.data_cells);

		area(c.x,c.y)=1;
		for(int n=1;n<=8;n++){
			if(!IN_GRID(c.x+dx[n],c.y+dy[n],flowdirs.width(),flowdirs.height()))
				continue;
			if(does_cell_flow_into_me(flowdirs(c.x+dx[n],c.y+dy[n]),n))
				area(c.x,c.y)+=proportion_i_get(flowdirs(c.x+dx[n],c.y+dy[n]),n)*area(c.x+dx[n],c.y+dy[n]);
		}
		if(flowdirs(c.x,c.y)!=flowdirs.no_data && flowdirs(c.x,c.y)!=NO_FLOW){
			int n_low= dinf_to_d8_low [(int)(flowdirs(c.x,c.y)/(M_PI/4))];
			int n_high=dinf_to_d8_high[(int)(flowdirs(c.x,c.y)/(M_PI/4))];
			int nlx=c.x+dx[n_low],nly=c.y+dy[n_low];
			int nhx=c.x+dx[n_high],nhy=c.y+dy[n_high];
			if( IN_GRID(nlx,nly,flowdirs.width(),flowdirs.height()) && area(nlx,nly)!=area.no_data && (--dependency(nlx,nly))==0)
				sources.push(grid_cell(nlx,nly));
			if( IN_GRID(nhx,nhy,flowdirs.width(),flowdirs.height()) && area(nhx,nhy)!=area.no_data && (--dependency(nhx,nhy))==0)
				sources.push(grid_cell(nhx,nhy));
		}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}








float dinf_masked_FlowDir(const int_2d &flat_resolution_mask, const int_2d &groups, const int x, const int y){
	double smax=0;
	int nmax=-1;
	double rmax=0;

	double e0,e1,e2,d1,d2,s1,s2,r,s;

	//Yes, this should be 0-8, this is the Tarboton neighbour system
	for(int n=0;n<8;n++){
		if(groups(x+dx_e1[n],y+dy_e1[n])!=groups(x,y)) continue;
		if(groups(x+dx_e2[n],y+dy_e2[n])!=groups(x,y)) continue;
		//if(elevations(x+dx_e1[n],y+dy_e1[n])==elevations.no_data) continue;
		//if(elevations(x+dx_e2[n],y+dy_e2[n])==elevations.no_data) continue;

		e0=flat_resolution_mask(x,y);
		e1=flat_resolution_mask(x+dx_e1[n],y+dy_e1[n]);
		e2=flat_resolution_mask(x+dx_e2[n],y+dy_e2[n]);
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

	double rg=NO_FLOW;
	if(nmax!=-1)
		rg=(af[nmax]*rmax+ac[nmax]*M_PI/2);

	return rg;
}

void dinf_flow_flats(const int_2d &flat_resolution_mask, const int_2d &groups, float_2d &flowdirs){
	diagnostic("Calculating Dinf flow directions using flat mask...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=1;x<flat_resolution_mask.width()-1;x++){
		progress_bar(x*omp_get_num_threads()*flat_resolution_mask.height()*100/(flat_resolution_mask.width()*flat_resolution_mask.height()));
		for(int y=1;y<flat_resolution_mask.height()-1;y++)
			if(flat_resolution_mask(x,y)==flat_resolution_mask.no_data)
				continue;
			else if(flowdirs(x,y)==NO_FLOW)
				flowdirs(x,y)=dinf_masked_FlowDir(flat_resolution_mask,groups,x,y);
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}
