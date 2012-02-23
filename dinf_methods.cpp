#include "data_structures.h"
#include "utility.h"
#include "d8_methods.h"
#include "interface.h"
#include <omp.h>
#include <cmath>

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

	if (EDGE_GRID(x,y,elevations.size1(),elevations.size2())) return 0;
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
		rg=af[nmax]*rmax+ac[nmax]*M_PI/2;

	return rg*(180.0/M_PI); //TODO: Converts to degrees
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

	print_flow(flowdirs);
}
