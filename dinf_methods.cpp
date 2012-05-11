#include "data_structures.h"
#include "utility.h"
#include "dinf_methods.h"
#include "interface.h"
#include "debug.h"
#include <cmath>
#include <queue>
#ifdef _OPENMP
	#include <omp.h>
#endif

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

//321
//4 0
//567
void where_do_i_flow(float flowdir, int &nhigh, int &nlow){
	//If it is very close to being directed into just one cell
	//then we direct it into just one cell. If we mistakenly direct
	//it into 2 cells, then we may create unresolvable loops in the
	//flow accumulation algorithm, whereas nothing is lost if we
	//leave out one of the two cells (provided a negligible amount
	//of flow is directed to the one we leave out).
	assert(flowdir>=0 && flowdir<=2*M_PI+1e-6);

	float temp=flowdir/(M_PI/4.);

	if(fabs(temp-(int)temp)<1e-6){
		nlow=-1;
		nhigh=(int)ROUND(temp);
	} else {
		nlow=(int)temp;
		nhigh=nlow+1;
	}

	//8 is not technically a direction, but, since things move in a circle,
	//it overlaps with 0. It should _never_ be greater than 8.
	assert(nhigh>=0 && nhigh<=8);
}

//This reacts correctly if the flow direction wedge number exceeds 7.
void area_proportion(float flowdir, int nhigh, int nlow, float &phigh, float &plow){
	if(nlow==-1){
		phigh=1;
		plow=0;
	} else {
		phigh=(nhigh*(M_PI/4.0)-flowdir)/(M_PI/4.0);
		plow=1-phigh;
	}

	assert(phigh+plow==1);	//TODO: This isn't necessarily so in floating-point... or is it?
}

/*//TODO: Debugging code used for checking for loops. Since loops should not occur in the output of the production code, this is not needed.
bool is_loop(const float_2d &flowdirs, int n, int x, int y, int c2x, int c2y){
	int nh,nl;
	if(! flowdirs.in_grid(c2x, c2y) || flowdirs(c2x,c2y)==flowdirs.no_data || flowdirs(c2x,c2y)==NO_FLOW)
		return false;
	where_do_i_flow(flowdirs(c2x,c2y),nh,nl);
	if(n==dinf_d8_inverse[nh] || (nl!=-1 && n==dinf_d8_inverse[nl])){
		printf("Beware dir %d (%d and %d).\n",n,nh,nl);
		flowdirs.surroundings(x,y,8);
		return true;
	}
	return false;
}*/

void dinf_upslope_area(const float_2d &flowdirs, float_2d &area){
	char_2d dependency;
	std::queue<grid_cell> sources;

	diagnostic_arg("The sources queue will require at most approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*((long)sizeof(grid_cell))/1024/1024);

	diagnostic_arg("Setting up the dependency matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*((long)sizeof(char))/1024/1024);
	diagnostic("Resizing dependency matrix...");
	dependency.resize(flowdirs.width(),flowdirs.height(),false);
	dependency.init(0);
	diagnostic("succeeded.\n");

	diagnostic_arg("Setting up the area matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*((long)sizeof(float))/1024/1024);
	diagnostic("Resizing the area matrix...");
	area.resize(flowdirs.width(),flowdirs.height(),false);
	area.init(0);
	area.no_data=dinf_NO_DATA;
	diagnostic("succeeded.\n");

	diagnostic("%%Calculating dependency matrix & setting no_data cells...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<flowdirs.width();x++){
		progress_bar(x*flowdirs.height()*100/(flowdirs.width()*flowdirs.height()));
		for(int y=0;y<flowdirs.height();y++){
			if(flowdirs(x,y)==flowdirs.no_data){
				area(x,y)=area.no_data;
				dependency(x,y)=9;	//Note: This is an unnecessary safety precaution
				continue;
			}
			if(flowdirs(x,y)==NO_FLOW)
				continue;
			int n_high,n_low;
			int nhx,nhy,nlx,nly;
			where_do_i_flow(flowdirs(x,y),n_high,n_low);
			nhx=x+dinf_dx[n_high],nhy=y+dinf_dy[n_high];
			if(n_low!=-1){
				nlx=x+dinf_dx[n_low];
				nly=y+dinf_dy[n_low];
			}
			if( n_low!=-1 && flowdirs.in_grid(nlx,nly) && flowdirs(nlx,nly)!=flowdirs.no_data )
				dependency(nlx,nly)++;
			if( flowdirs.in_grid(nhx,nhy) && flowdirs(nhx,nhy)!=flowdirs.no_data )
				dependency(nhx,nhy)++;
		}
	}
	diagnostic_arg(SUCCEEDED_IN,progress_bar(-1));

	diagnostic("%%Locating source cells...\n");
	progress_bar(-1);
	for(int x=0;x<flowdirs.width();x++){
		progress_bar(x*flowdirs.height()*100/(flowdirs.width()*flowdirs.height()));
		for(int y=0;y<flowdirs.height();y++)
			if(flowdirs(x,y)==flowdirs.no_data)
				continue;
			else if(flowdirs(x,y)==NO_FLOW)
				continue;
			else if(dependency(x,y)==0)
				sources.push(grid_cell(x,y));
	}
	diagnostic_arg(SUCCEEDED_IN,progress_bar(-1));

	diagnostic("%%Calculating up-slope areas...\n");
	progress_bar(-1);
	long int ccount=0;
	while(sources.size()>0){
		grid_cell c=sources.front();
		sources.pop();

		ccount++;
		progress_bar(ccount*100/flowdirs.data_cells);
		area(c.x,c.y)+=1;

		if(flowdirs(c.x,c.y)==flowdirs.no_data || flowdirs(c.x,c.y)==NO_FLOW)
			continue;

		int n_high,n_low,nhx,nhy,nlx,nly;
		where_do_i_flow(flowdirs(c.x,c.y),n_high,n_low);
		nhx=c.x+dinf_dx[n_high],nhy=c.y+dinf_dy[n_high];

		float phigh,plow;
		area_proportion(flowdirs(c.x,c.y), n_high, n_low, phigh, plow);
		if(flowdirs.in_grid(nhx,nhy))
			area(nhx,nhy)+=area(c.x,c.y)*phigh;

		if(n_low!=-1){
			nlx=c.x+dinf_dx[n_low];
			nly=c.y+dinf_dy[n_low];
			if(flowdirs.in_grid(nlx,nly)){
				area(nlx,nly)+=area(c.x,c.y)*plow;
				if(flowdirs(nlx,nly)!=flowdirs.no_data && (--dependency(nlx,nly))==0){
					sources.push(grid_cell(nlx,nly));
				}
			}
		}

		if( flowdirs.in_grid(nhx,nhy) && flowdirs(nhx,nhy)!=flowdirs.no_data && (--dependency(nhx,nhy))==0)
			sources.push(grid_cell(nhx,nhy));
	}
	diagnostic_arg(SUCCEEDED_IN,progress_bar(-1));
}




static const float d8_to_dinf[9]={-1, 4*M_PI/4, 3*M_PI/4, 2*M_PI/4, 1*M_PI/4, 0, 7*M_PI/4, 6*M_PI/4, 5*M_PI/4};

float dinf_masked_FlowDir(const int_2d &flat_resolution_mask, const int_2d &groups, const int x, const int y){
	double smax=0;
	int nmax=-1;
	double rmax=0;

	double e0,e1,e2,d1,d2,s1,s2,r,s;

	//Yes, this should be 0-8, this is the Tarboton neighbour system
	for(int n=0;n<8;n++){
		//TODO: Can these ever give !IN_GRID errors?
		if(groups(x+dx_e1[n],y+dy_e1[n])!=groups(x,y)) continue;
		if(groups(x+dx_e2[n],y+dy_e2[n])!=groups(x,y)) continue;

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
	else
		for(int n=1;n<=8;n++)	//TODO: I have a feeling this is potentially unsafe as it may create dependency loops. Does it?
			if(groups(x+dx[n],y+dy[n])==groups(x,y) && flat_resolution_mask(x+dx[n],y+dy[n])<flat_resolution_mask(x,y)){
				rg=d8_to_dinf[n];
				break;
			}

	return rg;
}

void dinf_flow_flats(const int_2d &flat_resolution_mask, const int_2d &groups, float_2d &flowdirs){
	diagnostic("%%Calculating Dinf flow directions using flat mask...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=1;x<flat_resolution_mask.width()-1;x++){
		progress_bar(x*flat_resolution_mask.height()*100/(flat_resolution_mask.width()*flat_resolution_mask.height()));
		for(int y=1;y<flat_resolution_mask.height()-1;y++)
			if(flat_resolution_mask(x,y)==flat_resolution_mask.no_data)
				continue;
			else if(flowdirs(x,y)==NO_FLOW)
				flowdirs(x,y)=dinf_masked_FlowDir(flat_resolution_mask,groups,x,y);
	}
	diagnostic_arg(SUCCEEDED_IN,progress_bar(-1));
}
