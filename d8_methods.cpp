#include "data_structures.h"
#include "utility.h"
#include "d8_methods.h"
#include "interface.h"
#include <queue>
#include <limits>
#include <cmath>
#include <stack>
#include <map>
#ifdef _OPENMP
	#include <omp.h>
#endif

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
	ProgressBar progress;

	diagnostic("%%Calculating D8 flow directions using flat mask...\n");
	progress.start( flat_resolution_mask.width()*flat_resolution_mask.height() );
	#pragma omp parallel for
	for(int x=1;x<flat_resolution_mask.width()-1;x++){
		progress.update( x*flat_resolution_mask.height() );
		for(int y=1;y<flat_resolution_mask.height()-1;y++)
			if(flat_resolution_mask(x,y)==flat_resolution_mask.no_data)
				continue;
			else if (flowdirs(x,y)==NO_FLOW)
				flowdirs(x,y)=d8_masked_FlowDir(flat_resolution_mask,groups,x,y);
	}
	diagnostic_arg(SUCCEEDED_IN,progress.stop());
}





void d8_upslope_area(const char_2d &flowdirs, int_2d &area){
	char_2d dependency;
	std::queue<grid_cell> sources;
	ProgressBar progress;

	diagnostic_arg("The sources queue will require at most approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*((long)sizeof(grid_cell))/1024/1024);

	diagnostic_arg("The dependency matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*((long)sizeof(char))/1024/1024);
	diagnostic("Resizing dependency matrix...");
	dependency.resize(flowdirs.width(),flowdirs.height(),false);
	diagnostic("succeeded.\n");

	diagnostic_arg("The area matrix will require approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*((long)sizeof(unsigned int))/1024/1024);
	diagnostic("Resizing the area matrix...");
	area.resize(flowdirs.width(),flowdirs.height(),false);
	diagnostic("succeeded.\n");
	diagnostic("Initializing the area matrix...");
	area.init(0);
	area.no_data=d8_NO_DATA;
	diagnostic("succeeded.\n");

	diagnostic("%%Calculating dependency matrix & setting no_data cells...\n");
	progress.start( flowdirs.width()*flowdirs.height() );
	#pragma omp parallel for
	for(int x=0;x<flowdirs.width();x++){
		progress.update( x*flowdirs.height() );
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
					++dependency(x,y);
		}
	}
	diagnostic_arg(SUCCEEDED_IN,progress.stop());

	diagnostic("%%Locating source cells...\n");
	progress.start( flowdirs.width()*flowdirs.height() );
	for(int x=0;x<flowdirs.width();x++){
		progress.update( x*flowdirs.height() );
		for(int y=0;y<flowdirs.height();y++)
			if(flowdirs(x,y)==flowdirs.no_data)
				continue;
			else if(dependency(x,y)==0)
				sources.push(grid_cell(x,y));
	}
	diagnostic_arg(SUCCEEDED_IN,progress.stop());

	diagnostic("%%Calculating up-slope areas...\n");
	progress.start(flowdirs.data_cells);
	long int ccount=0;
	while(sources.size()>0){
		grid_cell c=sources.front();
		sources.pop();

		ccount++;
		progress.update(ccount);

		area(c.x,c.y)+=1;

		if(flowdirs(c.x,c.y)==NO_FLOW)
			continue;

		int nx=c.x+dx[flowdirs(c.x,c.y)];
		int ny=c.y+dy[flowdirs(c.x,c.y)];
		if(flowdirs.in_grid(nx,ny) && area(nx,ny)!=area.no_data){
			area(nx,ny)+=area(c.x,c.y);
			if((--dependency(nx,ny))==0)
				sources.push(grid_cell(nx,ny));
		}
	}
	diagnostic_arg(SUCCEEDED_IN,progress.stop());
}






//Procedure:	d8_terrain_attrib_helper
//Description:
//		This calculates a variety of terrain attributes according
//		to the work of Burrough 1998's "Principles of Geographical
//		Information Systems" (p. 190). Algorithms and helpful ArcGIS
//		pages are noted in comments in the code.
//Restrictions:
//		This function should never be called on a NoData cell
//Inputs:
//		elevations		A 2D array of cell elevations
//		x0				X coordinate of cell to perform calculation on
//		y0				Y coordinate of cell to perform calculation on
//		rise_over_run	Returns rise-over-run slope
//		aspect			Returns aspect in degrees [0,360). 0 is north,
//						degrees increase in a clockwise fashion.
//		curvature		Returns the difference of profile and planform
//						curvatures
//		profile_curvature	Returns the profile curvature
//		planform_curvature	Returns the planform curvature
//Requirements:
//		None
//Effects:
//		Variables noted as "Returns" above are altered
//Returns:
//		None
inline static void d8_terrain_attrib_helper(const float_2d &elevations, int x0, int y0, float &rise_over_run, float &aspect, float &curvature, float &profile_curvature, float &planform_curvature){
/*
Slope derived from ArcGIS help at:
http://webhelp.esri.com/arcgiSDEsktop/9.3/index.cfm?TopicName=How%20Slope%20works
http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#/How_Slope_works/009z000000vz000000/

Cells are identified as
Aspect derived from AcrGIS help at:
http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#/How_Aspect_works/00q900000023000000/
http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#/How_Aspect_works/009z000000vp000000/

Curvature dervied from ArcGIS help at:
http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//00q90000000t000000
http://blogs.esri.com/esri/arcgis/2010/10/27/understanding-curvature-rasters/

Expanded ArcGIS curvature info
http://support.esri.com/en/knowledgebase/techarticles/detail/21942

Burrough 1998's "Principles of Geographical Information Systems" explains all these algorithms
*/

//a b c
//d e f
//g h i
	double a,b,c,d,e,f,g,h,i;
	//Deal with grid edges and NoData values in the manner suggested by
	//ArcGIS. Note that this function should never be called on a NoData cell
	a=b=c=d=e=f=g=h=i=elevations(x0,y0);
	if(elevations.in_grid(x0-1,y0-1))			a=elevations(x0-1,y0-1);
	if(elevations.in_grid(x0-1,y0))             d=elevations(x0-1,y0);
	if(elevations.in_grid(x0-1,y0+1))           g=elevations(x0-1,y0+1);
	if(elevations.in_grid(x0  ,y0-1))           b=elevations(x0,y0-1);
	if(elevations.in_grid(x0  ,y0+1))           h=elevations(x0,y0+1);
	if(elevations.in_grid(x0+1,y0-1))           c=elevations(x0+1,y0-1);
	if(elevations.in_grid(x0+1,y0))             f=elevations(x0+1,y0);
	if(elevations.in_grid(x0+1,y0+1))           i=elevations(x0+1,y0+1);
	if(a==elevations.no_data) a=e;
	if(b==elevations.no_data) b=e;
	if(c==elevations.no_data) c=e;
	if(d==elevations.no_data) d=e;
	if(f==elevations.no_data) f=e;
	if(g==elevations.no_data) g=e;
	if(h==elevations.no_data) h=e;
	if(i==elevations.no_data) i=e;

	//TODO: Convert elevations to meters? (1ft=0.3048m)
/*	a*=0.3048;
	b*=0.3048;
	c*=0.3048;
	d*=0.3048;
	e*=0.3048;
	f*=0.3048;
	g*=0.3048;
	h*=0.3048;
	i*=0.3048;*/

	double dzdx,dzdy;
	//Aspect calculation in the manner of Horn 1981
	//ArcGIS doesn't use cell size fo aspect calculations.
	dzdx=( (c+2*f+i) - (a+2*d+g)) / 8;
	dzdy=( (g+2*h+i) - (a+2*b+c)) / 8;
	aspect=180.0/M_PI*atan2(dzdy,-dzdx);
	if(aspect<0)
		aspect=90-aspect;
	else if(aspect>90.0)
		aspect=360.0-aspect+90.0;
	else
		aspect=90.0-aspect;

	//Slope calculation in the manner of Horn 1981
	//But cellsize is accounted for in slope
	dzdx/=elevations.cellsize;
	dzdy/=elevations.cellsize;
	rise_over_run=sqrt(dzdx*dzdx+dzdy*dzdy);

	if(rise_over_run==0)
		aspect=-1;	//Special value denoting a flat



//Z1 Z2 Z3   a b c
//Z4 Z5 Z6   d e f
//Z7 Z8 Z9   g h i
	//Curvatures in the manner of Zevenbergen and Thorne 1987
	double L=elevations.cellsize;		//TODO: Should be in the same units as z
	double D=( (d+f)/2 - e) / L / L;	//D = [(Z4 + Z6) /2 - Z5] / L^2
	double E=( (b+h)/2 - e) / L / L;	//E = [(Z2 + Z8) /2 - Z5] / L^2
	double F=(-a+c+g-i)/4/L/L;			//F=(-Z1+Z3+Z7-Z9)/(4L^2)
	double G=(-d+f)/2/L;				//G=(-Z4+Z6)/(2L)
	double H=(b-h)/2/L;					//H=(Z2-Z8)/(2L)
	curvature=(float)(-2*(D+E)*100);

	profile_curvature=(float)(2*(D*G*G+E*H*H+F*G*H)/(G*G+H*H)*100);
	planform_curvature=(float)(-2*(D*H*H+E*G*G-F*G*H)/(G*G+H*H)*100);
}




//Procedure:	d8_terrain_attribute
//Description:
//		This calculates a variety of terrain attributes according
//		to the work of Burrough 1998's "Principles of Geographical
//		Information Systems" (p. 190). This function calls
//		d8_terrain_attrib_helper to calculate the actual attributes.
//		This function may perform some post-processing (such as on
//		slope), but it's purpose is essentially to just scan the grid
//		and pass off the work.
//Restrictions:
//		None
//Inputs:
//		elevations		A 2D array of cell elevations
//		attribs			A 2D array to hold the calculated attribute
//		attrib			The attribute to be calculated
//Attributes:
//		Values "attrib" can take:
//			TATTRIB_SLOPE
//			TATTRIB_CURVATURE
//			TATTRIB_PLANFORM_CURVATURE
//			TATTRIB_PROFILE_CURVATURE
//			TATTRIB_ASPECT
//			TATTRIB_SLOPE_RISERUN
//			TATTRIB_SLOPE_PERCENT
//			TATTRIB_SLOPE_RADIAN
//			TATTRIB_SLOPE_DEGREE
//Requirements:
//		None
//Effects:
//		"attribs" is altered to contain the calculated attribute
//Returns:
//		None
void d8_terrain_attribute(const float_2d &elevations, float_2d &attribs, int attrib){
	ProgressBar progress;

	diagnostic_arg("The attribute #%d matrix will require approximately %ldMB of RAM.\n", attrib, elevations.width()*elevations.height()*((long)sizeof(float))/1024/1024);
	diagnostic("Setting up the attribute matrix...");
	attribs.copyprops(elevations);
	attribs.no_data=-2;	//TODO: Should push this out to the calling helper functions
	diagnostic("succeeded.\n");

	diagnostic_arg("%%Calculating terrain attribute #%d...\n",attrib);
	progress.start( elevations.width()*elevations.height() );
	#pragma omp parallel for
	for(int x=0;x<elevations.width();x++){
		progress.update( x*elevations.height() );
		for(int y=0;y<elevations.height();y++){
			if(elevations(x,y)==elevations.no_data){
				attribs(x,y)=attribs.no_data;
				continue;
			}
			float rise_over_run, aspect, curvature, profile_curvature, planform_curvature;
			d8_terrain_attrib_helper(elevations, x, y, rise_over_run, aspect, curvature, profile_curvature, planform_curvature);
			switch(attrib){
				case TATTRIB_SLOPE:
					attribs(x,y)=rise_over_run;break;
				case TATTRIB_CURVATURE:
					attribs(x,y)=curvature;break;
				case TATTRIB_PLANFORM_CURVATURE:
					attribs(x,y)=planform_curvature;break;
				case TATTRIB_PROFILE_CURVATURE:
					attribs(x,y)=profile_curvature;break;
				case TATTRIB_ASPECT:
					attribs(x,y)=aspect;break;
				case TATTRIB_SLOPE_RISERUN:
					attribs(x,y)=rise_over_run;
					break;
				case TATTRIB_SLOPE_PERCENT:
					attribs(x,y)=rise_over_run*100;
					break;
				case TATTRIB_SLOPE_RADIAN:
					attribs(x,y)=atan(rise_over_run);
					break;
				case TATTRIB_SLOPE_DEGREE:
					attribs(x,y)=atan(rise_over_run)*180/M_PI;
					break;
			}
		}
	}
	diagnostic_arg(SUCCEEDED_IN,progress.stop());
}

void d8_slope(const float_2d &elevations, float_2d &slopes, int slope_type){
	diagnostic("\n###Slope attribute calculation\n");
	d8_terrain_attribute(elevations, slopes, slope_type);
}


void d8_aspect(const float_2d &elevations, float_2d &aspects){
	diagnostic("\n###Aspect attribute calculation\n");
	d8_terrain_attribute(elevations, aspects, TATTRIB_ASPECT);
}

void d8_curvature(const float_2d &elevations, float_2d &curvatures){
	diagnostic("\n###Curvature attribute calculation\n");
	d8_terrain_attribute(elevations, curvatures, TATTRIB_CURVATURE);
}

void d8_planform_curvature(const float_2d &elevations, float_2d &planform_curvatures){
	diagnostic("\n###Planform curvature attribute calculation\n");
	d8_terrain_attribute(elevations, planform_curvatures,  TATTRIB_PLANFORM_CURVATURE);
}

void d8_profile_curvature(const float_2d &elevations, float_2d &profile_curvatures){
	diagnostic("\n###Profile curvature attribute calculation\n");
	d8_terrain_attribute(elevations, profile_curvatures,  TATTRIB_PROFILE_CURVATURE);
}











//FindWatersheds works in the same way as BarnesFlood, save that it labels watersheds, working inwards from the edges of the DEM. This is helpful for checking flow accumulation
//Procedure:	FindWatersheds
//Description:
//		Same as BarnesFlood. Labels starts out as no_data. If it is found that a
//		no_data labels cell coincides with a non-no_data elevations cell, then this
//		is the beginning of a new watershed. Cells which are flooded from a labeled
//		cell take on that cell's label
//Inputs:
//		elevations		A 2D array of cell elevations
//		labels			A 2D array to hold the watershed labels
//Requirements:
//		None
//Effects:
//		"labels" will be altered to contain the labels of the watersheds
//Returns:
//		None
void find_watersheds(float_2d &elevations, int_2d &labels, bool alter_elevations){
	std::priority_queue<grid_cellz, std::vector<grid_cellz>, grid_cellz_compare> open;
	std::stack<grid_cellz, std::vector<grid_cellz> > meander;
	bool_2d closed;
	unsigned long processed_cells=0;
	unsigned long pitc=0,openc=0;
	int clabel=1;	//TODO: Thought this was more clear than zero in the results.
	ProgressBar progress;

	diagnostic("\n###Barnes Flood+Watershed Labels\n");
	diagnostic_arg("The closed matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*((long)sizeof(bool))/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	closed.resize(elevations.width(),elevations.height());
	closed.init(false);
	diagnostic("succeeded.\n");

	diagnostic_arg("The labels matrix will require approximately %ldMB of RAM.\n",elevations.width()*elevations.height()*((long)sizeof(bool))/1024/1024);
	diagnostic("Setting up boolean flood array matrix...");
	labels.copyprops(elevations);
	labels.no_data=-1;
	labels.init(labels.no_data);
	diagnostic("succeeded.\n");

	diagnostic_arg("The open priority queue will require approximately %ldMB of RAM.\n",(elevations.width()*2+elevations.height()*2)*((long)sizeof(grid_cellz))/1024/1024);
	diagnostic("Adding cells to the open priority queue...");
	for(int x=0;x<elevations.width();x++){
		open.push(grid_cellz(x,0,elevations(x,0) ));
		open.push(grid_cellz(x,elevations.height()-1,elevations(x,elevations.height()-1) ));
		closed(x,0)=true;
		closed(x,elevations.height()-1)=true;
	}
	for(int y=1;y<elevations.height()-1;y++){
		open.push(grid_cellz(0,y,elevations(0,y)	));
		open.push(grid_cellz(elevations.width()-1,y,elevations(elevations.width()-1,y) ));
		closed(0,y)=true;
		closed(elevations.width()-1,y)=true;
	}
	diagnostic("succeeded.\n");

	diagnostic("%%Performing the Barnes Flood+Watershed Labels...\n");
	progress.start( elevations.width()*elevations.height() );
	while(open.size()>0 || meander.size()>0){
		grid_cellz c;
		if(meander.size()>0){
			c=meander.top();
			meander.pop();
			pitc++;
		} else {
			c=open.top();
			open.pop();
			openc++;
		}
		processed_cells++;

		if(labels(c.x,c.y)==labels.no_data && elevations(c.x,c.y)!=elevations.no_data)	//Implies a cell without a label which borders the edge of the DEM or a region of no_data
			labels(c.x,c.y)=clabel++;

		for(int n=1;n<=8;n++){
			int nx=c.x+dx[n];
			int ny=c.y+dy[n];
			if(!elevations.in_grid(nx,ny)) continue;
			if(closed(nx,ny)) 
				continue;

			labels(nx,ny)=labels(c.x,c.y);

			closed(nx,ny)=true;
			if(elevations(nx,ny)<=c.z){
				if(alter_elevations)
					elevations(nx,ny)=c.z;
				meander.push(grid_cellz(nx,ny,c.z));
			} else
				open.push(grid_cellz(nx,ny,elevations(nx,ny)));
		}
		progress.update(processed_cells);
	}
	diagnostic_arg("\t\033[96msucceeded in %.2lfs\033[39m\n",progress.stop());
	diagnostic_arg("%ld cells processed. %ld in pits, %ld not in pits.\n",processed_cells,pitc,openc);
}


void watershed_area(const int_2d &labels){
	std::map<int, int> wsheds;
	for(int x=0;x<labels.width();x++)
	for(int y=0;y<labels.height();y++)
		if(labels(x,y)!=labels.no_data)
			wsheds[labels(x,y)]++;

	for(std::map<int, int>::iterator i=wsheds.begin();i!=wsheds.end();i++)
		printf("Watershed %d has area %d\n",(*i).first,(*i).second);
}














/**
	@brief  Calculates the SPI terrain attribute
	@author Richard Barnes

	@param[in]	&flow_accumulation	A flow accumulation grid (#dinf_upslope_area)
	@param[in]	&percent_slope		A percent_slope grid (#d8_slope)
	@param[out]	&result				Altered to return the calculated SPI

	@pre flow_accumulation and percent_slope must be the same size

	@post result takes the properties and dimensions of flow_accumulation

	@todo Generalize for float and int grids
*/
void d8_SPI(const float_2d &flow_accumulation, const float_2d &percent_slope, float_2d &result){
	Timer timer;

	diagnostic("\n###d8_SPI\n");

	if(flow_accumulation.width()!=percent_slope.width() || flow_accumulation.height()!=percent_slope.height()){
		diagnostic("Couldn't calculate SPI! The input matricies were of unequal dimensions!\n");
		exit(-1);
	}

	diagnostic_arg("The SPI matrix will require approximately %ldMB of RAM.\n", flow_accumulation.width()*flow_accumulation.height()*((long)sizeof(float))/1024/1024);
	diagnostic("Setting up the SPI matrix...");
	result.copyprops(flow_accumulation);
	result.no_data=-1;	//Log(x) can't take this value of real inputs, so we're good
	diagnostic("succeeded.\n");

	diagnostic("Calculating SPI...\n");
	timer.start();
	#pragma omp parallel for collapse(2)
	for(int x=0;x<flow_accumulation.width();x++)
		for(int y=0;y<flow_accumulation.height();y++)
			if(flow_accumulation(x,y)==flow_accumulation.no_data || percent_slope(x,y)==percent_slope.no_data)
				result(x,y)=result.no_data;
			else
				result(x,y)=log( flow_accumulation.cellsize*(flow_accumulation(x,y)+0.001) * ((percent_slope(x,y)/100) + 0.001) );
	diagnostic_arg("succeeded in %lfs.\n",timer.lap());
}







/**
	@brief  Calculates the CTI terrain attribute
	@author Richard Barnes

	@param[in]	&flow_accumulation	A flow accumulation grid (#dinf_upslope_area)
	@param[in]	&percent_slope		A percent_slope grid (#d8_slope)
	@param[out]	&result				Altered to return the calculated SPI

	@pre flow_accumulation and percent_slope must be the same size

	@post result takes the properties and dimensions of flow_accumulation

	@todo Generalize for float and int grids
*/
void d8_CTI(const float_2d &flow_accumulation, const float_2d &percent_slope, float_2d &result){
	Timer timer;

	diagnostic("\n###d8_CTI\n");

	if(flow_accumulation.width()!=percent_slope.width() || flow_accumulation.height()!=percent_slope.height()){
		diagnostic("Couldn't calculate CTI! The input matricies were of unequal dimensions!\n");
		exit(-1);
	}

	diagnostic_arg("The CTI matrix will require approximately %ldMB of RAM.\n", flow_accumulation.width()*flow_accumulation.height()*((long)sizeof(float))/1024/1024);
	diagnostic("Setting up the CTI matrix...");
	result.copyprops(flow_accumulation);
	result.no_data=-1;	//Log(x) can't take this value of real inputs, so we're good
	diagnostic("succeeded.\n");

	diagnostic("Calculating CTI...");
	timer.start();
	#pragma omp parallel for collapse(2)
	for(int x=0;x<flow_accumulation.width();x++)
		for(int y=0;y<flow_accumulation.height();y++)
			if(flow_accumulation(x,y)==flow_accumulation.no_data || percent_slope(x,y)==percent_slope.no_data)
				result(x,y)=result.no_data;
			else
				result(x,y)=log( flow_accumulation.cellsize*(flow_accumulation(x,y)+0.001) / ((percent_slope(x,y)/100) + 0.001) );
	diagnostic_arg("succeeded in %lfs.\n",timer.lap());
}
