#ifndef _dinf_methods_included
#define _dinf_methods_included

#include "interface.h"
#include "data_structures.h"
#include <omp.h>
#include "utility.h"

void dinf_upslope_area(const float_2d &flowdirs, float_2d &area);

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
static const int dy_e1[8]={0,-1,-1,0,0,1,1,0};
static const int dx_e1[8]={1,0,0,-1,-1,0,0,1};
static const int dy_e2[8]={-1,-1,-1,-1,1,1,1,1};
static const int dx_e2[8]={1,1,-1,-1,-1,-1,1,1};
static const double ac[8]={0,1,1,2,2,3,3,4};
static const double af[8]={1,-1,1,-1,1,-1,1,-1};

template <class T>
float dinf_FlowDir(const array2d<T> &elevations, const int x, const int y){
	double smax=0;
	int nmax=-1;
	double rmax=0;

	double e0,e1,e2,d1,d2,s1,s2,r,s;

	if (EDGE_GRID(x,y,elevations.width(),elevations.height())){
		if(x==0 && y==0)
			return 3*M_PI/4;	//D8: 2
		else if(x==0 && y==elevations.height()-1)
			return 5*M_PI/4;	//D8: 8
		else if(x==elevations.width()-1 && y==0)
			return 1*M_PI/4;	//D8: 4
		else if(x==elevations.width()-1 && y==elevations.height()-1)
			return 7*M_PI/4;	//D8: 6
		else if(x==0)
			return 4*M_PI/4;	//D8: 1
		else if(x==elevations.width()-1)
			return 0*M_PI/4;	//D8: 5
		else if(y==0)
			return 2*M_PI/4;	//D8: 3
		else if(y==elevations.height()-1)
			return 6*M_PI/4;	//D8: 7
	}
	
	//Since I am not on the edge of the grid if I've made it this far, may neighbours cannot be off the grid
	//Very negative no_data's should be acceptable, and suck water of the grid.
	for(int n=0;n<8;n++){
		//if(elevations(x+dx_e1[n],y+dy_e1[n])==elevations.no_data) continue;
		//if(elevations(x+dx_e2[n],y+dy_e2[n])==elevations.no_data) continue;

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

	double rg=NO_FLOW;
	if(nmax!=-1)
		rg=(af[nmax]*rmax+ac[nmax]*M_PI/2);

	return rg;
}

template <class T>
void dinf_flow_directions(const array2d<T> &elevations, float_2d &flowdirs){
	diagnostic_arg("The Dinf flow directions will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(float)/1024/1024);
	diagnostic("Resizing flow directions matrix...");
	flowdirs.resize(elevations.width(),elevations.height(),!init);
	diagnostic("succeeded.\n");

	diagnostic("Setting no_data value on flowdirs matrix...");
	flowdirs.no_data=dinf_NO_DATA;
	diagnostic("succeeded.\n");

	diagnostic("Initializing Dinf flow directions...\n");
	flowdirs.init(NO_FLOW);
	diagnostic("\tsucceeded.\n");

	diagnostic("Calculating Dinf flow directions...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++){
		progress_bar(x*omp_get_num_threads()*elevations.height()*100/(elevations.width()*elevations.height()));
		for(int y=0;y<elevations.height();y++)
			if(elevations(x,y)==elevations.no_data)
				flowdirs(x,y)=dinf_NO_DATA;
			else if(flowdirs(x,y)==NO_FLOW)
				flowdirs(x,y)=dinf_FlowDir(elevations,x,y);
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}

#endif
