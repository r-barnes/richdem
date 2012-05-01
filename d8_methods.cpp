#include "data_structures.h"
#include "utility.h"
#include "d8_methods.h"
#include "interface.h"
#include <omp.h>
#include <queue>
#include <limits>
#include <cmath>



//234
//105
//876
int d8_masked_FlowDir(const int_2d &flat_resolution_mask, const int_2d &groups, const int x, const int y){
	int minimum_elevation=flat_resolution_mask(x,y);
	int flowdir=NO_FLOW;

	for(int n=1;n<=8;n++){
		int nx=x+dx[n];
		int ny=y+dy[n];
		if(	groups(nx,ny)==groups(x,y) && (flat_resolution_mask(nx,ny)<minimum_elevation || (flat_resolution_mask(nx,ny)==minimum_elevation && flowdir>0 && flowdir%2==0 && n%2==1)) ){
			minimum_elevation=flat_resolution_mask(nx,ny);
			flowdir=n;
		}
	}

	return flowdir;
}

void d8_flow_flats(const int_2d &flat_resolution_mask, const int_2d &groups, char_2d &flowdirs){
	diagnostic("Calculating D8 flow directions using flat mask...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=1;x<flat_resolution_mask.width()-1;x++){
		progress_bar(x*omp_get_num_threads()*flat_resolution_mask.height()*100/(flat_resolution_mask.width()*flat_resolution_mask.height()));
		for(int y=1;y<flat_resolution_mask.height()-1;y++)
			if(flat_resolution_mask(x,y)==flat_resolution_mask.no_data)
				continue;
			else if (flowdirs(x,y)==NO_FLOW)
				flowdirs(x,y)=d8_masked_FlowDir(flat_resolution_mask,groups,x,y);
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
}





void d8_upslope_area(const char_2d &flowdirs, int_2d &area){
	char_2d dependency;
	std::queue<grid_cell> sources;

	diagnostic_arg("The sources queue will require at most approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(grid_cell)/1024/1024);

	diagnostic_arg("The dependency matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(char)/1024/1024);
	diagnostic("Resizing dependency matrix...");
	dependency.resize(flowdirs.width(),flowdirs.height(),false);
	diagnostic("succeeded.\n");

	diagnostic_arg("The area matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*sizeof(unsigned int)/1024/1024);
	diagnostic("Resizing the area matrix...");
	area.resize(flowdirs.width(),flowdirs.height(),false);
	diagnostic("succeeded.\n");
	diagnostic("Initializing the area matrix...");
	area.init(-1);
	area.no_data=d8_NO_DATA;
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
			for(int n=1;n<=8;n++)
				if(!flowdirs.in_grid(x+dx[n],y+dy[n]))
					continue;
				else if(flowdirs(x+dx[n],y+dy[n])==NO_FLOW)
					continue;
				else if(flowdirs(x+dx[n],y+dy[n])==flowdirs.no_data)
					continue;
				else if(n==inverse_flow[flowdirs(x+dx[n],y+dy[n])])
					dependency(x,y)++;
		}
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));

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
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));

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
			if(!flowdirs.in_grid(c.x+dx[n],c.y+dy[n]))
				continue;
			if(flowdirs(c.x+dx[n],c.y+dy[n])!=NO_FLOW && n==inverse_flow[flowdirs(c.x+dx[n],c.y+dy[n])])
				area(c.x,c.y)+=area(c.x+dx[n],c.y+dy[n]);
		}
		if(flowdirs(c.x,c.y)!=NO_FLOW){
			int nx=c.x+dx[flowdirs(c.x,c.y)];
			int ny=c.y+dy[flowdirs(c.x,c.y)];
			if( flowdirs.in_grid(nx,ny) && area(nx,ny)!=area.no_data && (--dependency(nx,ny))==0)
				sources.push(grid_cell(nx,ny));
		}
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
}







inline void d8_slope_and_aspect(const float_2d &elevations, int x0, int y0, float &slope, float &aspect, float &curvature){
//Slope derived from ArcGIS help at:
//http://webhelp.esri.com/arcgiSDEsktop/9.3/index.cfm?TopicName=How%20Slope%20works
//Cells are identified as
//Aspect derived from AcrGIS help at:
//http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#/How_Aspect_works/00q900000023000000/
//Curvature dervied from ArcGIS help at:
//http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//00q90000000t000000
//http://blogs.esri.com/esri/arcgis/2010/10/27/understanding-curvature-rasters/

//a b c
//d e f
//g h i
	float a,b,c,d,e,f,g,h,i;
	a=b=c=d=e=f=g=h=i=elevations(x0,y0);
	if(x0-1>0 && y0-1>0)                                    a=elevations(x0-1,y0-1);
	if(x0-1>0)                                              d=elevations(x0-1,y0);
	if(x0-1>0 && y0+1<elevations.height())                  g=elevations(x0-1,y0+1);
	if(y0-1>0)                                              b=elevations(x0,y0-1);
	if(y0+1<elevations.height())                            h=elevations(x0,y0+1);
	if(x0+1<elevations.width() && y0-1>0)                   c=elevations(x0+1,y0-1);
	if(x0+1<elevations.width())                             f=elevations(x0+1,y0);
	if(x0+1<elevations.width() && y0+1<elevations.height()) i=elevations(x0+1,y0+1);
	if(a==elevations.no_data) a=e;
	if(b==elevations.no_data) b=e;
	if(c==elevations.no_data) c=e;
	if(d==elevations.no_data) d=e;
	if(f==elevations.no_data) f=e;
	if(g==elevations.no_data) g=e;
	if(h==elevations.no_data) h=e;
	if(i==elevations.no_data) i=e;

	//TODO: Could use actual cell sizes below
	float dzdx=( (c+2*f+i) - (a+2*d+g)) / 8; //(8*x_cell_size);
	float dzdy=( (g+2*h+i) - (a+2*b+c)) / 8; //(8*y_cell_size);
	slope=sqrt(dzdx*dzdx+dzdy*dzdy);


	aspect=180.0/M_PI*atan2(dzdy,-dzdx);
	if(aspect<0)
		aspect=90-aspect;
	else if(aspect>90.0)
		aspect=360.0-aspect+90.0;
	else
		aspect=90.0-aspect;
	if(slope==0)
		aspect=-1;	//Special value denoting a flat

//Z1 Z2 Z3
//Z4 Z5 Z6
//Z7 Z8 Z9
	float D=( (d+f)/2 - e) / elevations.cellsize;	//D = [(Z4 + Z6) /2 - Z5] / L2
	float E=( (b+h)/2 - e) / elevations.cellsize;	//E = [(Z2 + Z8) /2 - Z5] / L2
	curvature=-2*(D+E)*100;
}


void d8_slope(const float_2d &elevations, float_2d &slopes, int slope_type){
	diagnostic_arg("The slope matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(float)/1024/1024);
	diagnostic("Setting up the slope matrix...");
	slopes.resize(elevations.width(),elevations.height());
	slopes.no_data=-1;
	diagnostic("succeeded.\n");

	diagnostic("Calculating slopes...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++){
		progress_bar(x*omp_get_num_threads()*elevations.height()*100/(elevations.width()*elevations.height()));
		for(int y=0;y<elevations.height();y++){
			if(elevations(x,y)==elevations.no_data){
				slopes(x,y)=slopes.no_data;
				continue;
			}
			float rise_over_run,aspect,curvature;
			d8_slope_and_aspect(elevations,x,y,rise_over_run,aspect,curvature);
			switch(slope_type){
				case SLOPE_RISERUN:
					slopes(x,y)=rise_over_run;
					break;
				case SLOPE_PERCENT:
					slopes(x,y)=rise_over_run*100;
					break;
				case SLOPE_RADIAN:
					slopes(x,y)=atan(rise_over_run);
					break;
				case SLOPE_DEGREE:
					slopes(x,y)=atan(rise_over_run)*180/M_PI;
					break;
			}
		}
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
}



void d8_aspect(const float_2d &elevations, float_2d &aspects){
	diagnostic_arg("The aspects matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(float)/1024/1024);
	diagnostic("Setting up the aspects matrix...");
	aspects.resize(elevations.width(),elevations.height());
	aspects.no_data=-2;
	diagnostic("succeeded.\n");

	diagnostic("Calculating aspects...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++){
		progress_bar(x*omp_get_num_threads()*elevations.height()*100/(elevations.width()*elevations.height()));
		for(int y=0;y<elevations.height();y++){
			if(elevations(x,y)==elevations.no_data){
				aspects(x,y)=aspects.no_data;
				continue;
			}
			float rise_over_run,aspect,curvature;
			d8_slope_and_aspect(elevations,x,y,rise_over_run,aspect,curvature);
			aspects(x,y)=aspect;
		}
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
}





void d8_curvature(const float_2d &elevations, float_2d &curvatures){
	diagnostic_arg("The curvatures matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*sizeof(float)/1024/1024);
	diagnostic("Setting up the curvatures matrix...");
	curvatures.resize(elevations.width(),elevations.height());
	curvatures.no_data=-2;
	diagnostic("succeeded.\n");

	diagnostic("Calculating curvatures...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++){
		progress_bar(x*omp_get_num_threads()*elevations.height()*100/(elevations.width()*elevations.height()));
		for(int y=0;y<elevations.height();y++){
			if(elevations(x,y)==elevations.no_data){
				curvatures(x,y)=curvatures.no_data;
				continue;
			}
			float rise_over_run,aspect,curvature;
			d8_slope_and_aspect(elevations,x,y,rise_over_run,aspect,curvature);
			curvatures(x,y)=curvature;
		}
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
}







//TODO: ArcGIS calculates flow directions differently, so the procedure below is not valid
//http://webhelp.esri.com/arcgiSDEsktop/9.3/index.cfm?TopicName=Determining_flow_direction
//http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=flow_direction
void d8_arcgis_convert(char_2d &flowdirs){
//234   32 64 128
//105   16     1
//876    8  4  2
	const char arcgis_flowdirs[9]={0,16,32,64,128,1,2,4,8};
	diagnostic("Converting flow directions to ArcGIS format...\n");
	progress_bar(-1);
	#pragma omp parallel for
	for(int x=0;x<flowdirs.width();x++){
		progress_bar(x*omp_get_num_threads()*flowdirs.height()*100/(flowdirs.width()*flowdirs.height()));
		for(int y=0;y<flowdirs.height();y++){
			if(flowdirs(x,y)==flowdirs.no_data) continue;
			flowdirs(x,y)=arcgis_flowdirs[flowdirs(x,y)];
		}
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs.\033[39m\n",progress_bar(-1));
}
