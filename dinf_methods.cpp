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
int ac[8]={0,1,1,2,2,3,3,4};
int af[8]={1,-1,1,-1,1,-1,1,-1};

float dinf_FlowDir(float_2d &elevations, int x, int y){
	float smax=0;
	int nmax=0;
	float rmax=0;
	float flowdir=NO_FLOW;
	double e0,e1,e2,d1,d2,s1,s2,r,s;
//	flat tflat={x,y};
//	bool has_different=false,has_same=false;

	if (elevations(x,y)==elev_NO_DATA) return dinf_NO_DATA; //Missing data

	for(int n=0;n<8;n++){
		if(elevations(x+dx_e1[n],y+dy_e1[n])==elev_NO_DATA) continue;
		if(elevations(x+dx_e2[n],y+dy_e2[n])==elev_NO_DATA) continue;
		if(!IN_GRID(x+dx_e1[n],y+dy_e1[n],elevations.size1(),elevations.size2())) continue;
		if(!IN_GRID(x+dx_e2[n],y+dy_e2[n],elevations.size1(),elevations.size2())) continue;
//		if(elevations(x+dx[n],y+dy[n])!=elevations(x,y))
//			has_different=true;
//		else
//			has_same=true;

		//TODO: Need to detect if slope is actually going up!
		e0=elevations(x,y);
		e1=elevations(x+dx_e1[n],y+dy_e1[n]);
		e2=elevations(x+dx_e2[n],y+dy_e2[n]);
		d1=1;
		d2=1;
		s1=(e0-e1)/d1;
		s2=(e1-e2)/d2;
		r=atan(s2/s1);
		s=sqrt(s1*s1+s2*s2);
		if(r<0){
			r=0;s=s1;
		} else if(r>atan(d2/d1)){
			r=atan(d2/d1);
			s=(e0-e2)/sqrt(d1*d1+d2*d2);
		}
		if(s>smax){
			smax=s;
			nmax=n;
			rmax=r;
		}
	}

	double rg=-1;
	if(nmax!=0)
		rg=af[nmax]*rmax+ac[nmax]*M_PI/2;

//	if(has_same && has_different){
//		tflat.x=x;
//		tflat.y=y;
//		flats.push(tflat);
//	}

	return rg;
}

int dinf_flow_directions(float_2d elevations, float_2d flowdirs){
	diagnostic_arg("The Dinf flow directions will require approximately %ldMB of RAM.\n",elevations.size1()*elevations.size2()*sizeof(float)/1024/1024);
	diagnostic("Resizing flow directions matrix...");
	try{
		flowdirs.resize(elevations.size1(),elevations.size2());
	} catch (std::exception &e){
		diagnostic("failed!\n");
		return -1;
	}
	diagnostic("succeeded.\n");

	diagnostic("Calculating Dinf flow directions...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.size1();x++){
		progress_bar(x*omp_get_num_threads()*elevations.size2()*100/(elevations.size1()*elevations.size2())); //Todo: Should I check to see if ftell fails here?
		for(int y=0;y<elevations.size2();y++){
			if(EDGE_GRID(x,y,elevations.size1(),elevations.size2())){ //Edges of grid undefined
				flowdirs(x,y)=0;
				continue;
			}
			if(flowdirs(x,y)!=0) continue; //Flow direction already set
				flowdirs(x,y)=dinf_FlowDir(elevations,x,y);
		}
	}
	progress_bar(-1);
	diagnostic("\tsucceeded.\n");
}
