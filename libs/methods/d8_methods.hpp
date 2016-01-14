#ifndef _richdem_d8_methods_hpp_
#define _richdem_d8_methods_hpp_

#include "../common/Array2D.hpp"

/// Used with #d8_terrain_attribute to get an ill-defined curvature thing
/// (TODO)
#define TATTRIB_CURVATURE           2

/// Used with #d8_terrain_attribute to get planform curvature as per
/// Zevenbergen and Thorne 1987 */
#define TATTRIB_PLANFORM_CURVATURE  3

/// Used with #d8_terrain_attribute to get profile curvature as per
/// Zevenbergen and Thorne 1987 */
#define TATTRIB_PROFILE_CURVATURE   4

/// Used with #d8_terrain_attribute to get aspect as per Horn 1981
#define TATTRIB_ASPECT              5

/// Used with #d8_terrain_attribute to get slope as per Horn 1981
#define TATTRIB_SLOPE_RISERUN       6

/// Used with #d8_terrain_attribute to get slope as per Horn 1981, the value
/// is multiplied by 100
#define TATTRIB_SLOPE_PERCENT       7

/// Used with #d8_terrain_attribute to get slope as per Horn 1981, an arc
/// tangent of the value is taken
#define TATTRIB_SLOPE_RADIAN        8

/// Used with #d8_terrain_attribute to get slope as per Horn 1981, an arc 
/// tangent of the value is taken and converted to degrees
#define TATTRIB_SLOPE_DEGREE        9


//d8_upslope_area
/**
  @brief  Calculates the D8 up-slope area, given the D8 flow directions
  @author Richard Barnes (rbarnes@umn.edu)

  This calculates the D8 up-slope area of a grid of D8 flow directions using
  by calculating each cell's dependency on its neighbours and then using
  a priority-queue to process cells in a top-of-the-watershed-down fashion

  @param[in]  &flowdirs  A D8 flowdir grid from d8_flow_directions()
  @param[out] &area      Returns the up-slope area of each cell
*/
template<class T, class U>
void d8_upslope_area(const Array2D<T> &flowdirs, Array2D<U> &area){
  char_2d dependency;
  std::queue<grid_cell> sources;
  ProgressBar progress;

  diagnostic("\n###D8 Upslope Area\n");

  diagnostic_arg(
    "The sources queue will require at most approximately %ldMB of RAM.\n",
    flowdirs.width()*flowdirs.height()*((long)sizeof(grid_cell))/1024/1024
  );

  diagnostic("Resizing dependency matrix...");
  dependency.copyprops(flowdirs);
  diagnostic("succeeded.\n");

  diagnostic("Setting up the area matrix...");
  area.copyprops(flowdirs);
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
        else if(n==inverse_flow[(int)flowdirs(x+dx[n],y+dy[n])])
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

    int nx=c.x+dx[(int)flowdirs(c.x,c.y)];
    int ny=c.y+dy[(int)flowdirs(c.x,c.y)];
    if(flowdirs.in_grid(nx,ny) && area(nx,ny)!=area.no_data){
      area(nx,ny)+=area(c.x,c.y);
      if((--dependency(nx,ny))==0)
        sources.push(grid_cell(nx,ny));
    }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
}




//d8_upslope_cells
/**
  @brief  Calculates which cells ultimately D8-flow through a given cell
  @author Richard Barnes (rbarnes@umn.edu)

  Given the coordinates x, y of a cell, this returns a grid indicating
  which cells ultimately flow into the indicated cell.
  1=Upslope cell
  2=Member of initializing line
  All other cells have a no_data value

  @param[in]  &flowdirs  A D8 flowdir grid from d8_flow_directions()
  @param[out] &area      Returns the up-slope area of each cell
*/
template<class T, class U>
void d8_upslope_cells(
  int x0,
  int y0,
  int x1,
  int y1,
  const Array2D<T> &flowdirs,
  Array2D<U>       &upslope_cells
){
  diagnostic("Setting up the upslope_cells matrix...");
  upslope_cells.copyprops(flowdirs);
  upslope_cells.init(d8_NO_DATA);
  upslope_cells.no_data=d8_NO_DATA;
  diagnostic("succeeded.\n");
  ProgressBar progress;

  std::queue<grid_cell> expansion;

  if(x0>x1){
    std::swap(x0,x1);
    std::swap(y0,y1);
  }

  //Modified Bresenham Line-Drawing Algorithm
  int deltax     = x1-x0;
  int deltay     = y1-y0;
  float error    = 0;
  float deltaerr = (float)deltay/(float)deltax;

  if (deltaerr<0)
    deltaerr = -deltaerr;

  diagnostic_arg("Line slope is %f\n",deltaerr);
  int y=y0;
  for(int x=x0;x<=x1;x++){
    expansion.push(grid_cell(x,y));
    upslope_cells(x,y)=2;
    error+=deltaerr;
    if (error>=0.5) {
      expansion.push(grid_cell(x+1,y));
      upslope_cells(x+1,y) = 2;
      y                   += sgn(deltay);
      error               -= 1;
    }
  }

  progress.start(flowdirs.data_cells);
  long int ccount=0;
  while(expansion.size()>0){
    grid_cell c=expansion.front();
    expansion.pop();

    progress.update(ccount++);

    for(int n=1;n<=8;n++)
      if(!flowdirs.in_grid(c.x+dx[n],c.y+dy[n]))
        continue;
      else if(flowdirs(c.x+dx[n],c.y+dy[n])==NO_FLOW)
        continue;
      else if(flowdirs(c.x+dx[n],c.y+dy[n])==flowdirs.no_data)
        continue;
      else if(upslope_cells(c.x+dx[n],c.y+dy[n])==upslope_cells.no_data && n==inverse_flow[flowdirs(c.x+dx[n],c.y+dy[n])]){
        expansion.push(grid_cell(c.x+dx[n],c.y+dy[n]));
        upslope_cells(c.x+dx[n],c.y+dy[n])=1;
      }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
  diagnostic_arg("Found %ld up-slope cells.\n",ccount);
}




//d8_SPI
/**
  @brief  Calculates the SPI terrain attribute
  @author Richard Barnes (rbarnes@umn.edu)

  \f$(\textit{CellSize}\cdot\textit{FlowAccumulation}+0.001)\cdot(\frac{1}{100}\textit{PercentSlope}+0.001)\f$

  @param[in]   &flow_accumulation
    A flow accumulation grid (dinf_upslope_area())
  @param[in]   &percent_slope      
    A percent_slope grid (d8_slope())
  @param[out]  &result            
    Altered to return the calculated SPI

  @pre \pname{flow_accumulation} and \pname{percent_slope} must be the same size

  @post \pname{result} takes the properties and dimensions of \pname{flow_accumulation}

  @todo Generalize for float and int grids
*/

template<class T, class U, class V>
void d8_SPI(
  const Array2D<T> &flow_accumulation,
  const Array2D<U> &percent_slope,
        Array2D<V> &result
){
  Timer timer;

  diagnostic("\n###d8_SPI\n");

  if(flow_accumulation.width()!=percent_slope.width() || flow_accumulation.height()!=percent_slope.height()){
    diagnostic("Couldn't calculate SPI! The input matricies were of unequal dimensions!\n");
    exit(-1);
  }

  diagnostic("Setting up the SPI matrix...");
  result.copyprops(flow_accumulation);
  result.no_data=-1;  //Log(x) can't take this value of real inputs, so we're good
  diagnostic("succeeded.\n");

  diagnostic("Calculating SPI...\n");
  timer.start();
  #pragma omp parallel for collapse(2)
  for(int x=0;x<flow_accumulation.width();x++)
    for(int y=0;y<flow_accumulation.height();y++)
      if(flow_accumulation(x,y)==flow_accumulation.no_data || percent_slope(x,y)==percent_slope.no_data)
        result(x,y)=result.no_data;
      else
        result(x,y)=log( (flow_accumulation(x,y)/flow_accumulation.cellsize + 0.001) * (percent_slope(x,y)/100+0.001) );
  diagnostic_arg("succeeded in %lfs.\n",timer.lap());
}






//d8_CTI
/**
  @brief  Calculates the CTI terrain attribute
  @author Richard Barnes (rbarnes@umn.edu)

  \f$\log{\frac{\textit{CellSize}\cdot\textit{FlowAccumulation}+0.001}{\frac{1}{100}\textit{PercentSlope}+0.001}}\f$

  @param[in]  &flow_accumulation 
    A flow accumulation grid (dinf_upslope_area())
  @param[in]  &percent_slope     
    A percent_slope grid (d8_slope())
  @param[out] &result             
    Altered to return the calculated SPI

  @pre \pname{flow_accumulation} and \pname{percent_slope} must be the same size

  @post \pname{result} takes the properties and dimensions of \pname{flow_accumulation}

  @todo Generalize for float and int grids
*/
template<class T, class U, class V>
void d8_CTI(
  const Array2D<T> &flow_accumulation,
  const Array2D<U> &percent_slope,
        Array2D<V> &result
){
  Timer timer;

  diagnostic("\n###d8_CTI\n");

  if(flow_accumulation.width()!=percent_slope.width() || flow_accumulation.height()!=percent_slope.height()){
    diagnostic("Couldn't calculate CTI! The input matricies were of unequal dimensions!\n");
    exit(-1);
  }

  diagnostic("Setting up the CTI matrix...");
  result.copyprops(flow_accumulation);
  result.no_data=-1;  //Log(x) can't take this value of real inputs, so we're good
  diagnostic("succeeded.\n");

  diagnostic("Calculating CTI...");
  timer.start();
  #pragma omp parallel for collapse(2)
  for(int x=0;x<flow_accumulation.width();x++)
    for(int y=0;y<flow_accumulation.height();y++)
      if(flow_accumulation(x,y)==flow_accumulation.no_data || percent_slope(x,y)==percent_slope.no_data)
        result(x,y)=result.no_data;
      else
        result(x,y)=log( (flow_accumulation(x,y)/flow_accumulation.cellsize + 0.001) / (percent_slope(x,y)/100+0.001) );
  diagnostic_arg("succeeded in %lfs.\n",timer.lap());
}



//d8_terrain_attrib_helper
/**
  @brief  Calculate a variety of terrain attributes
  @author Richard Barnes (rbarnes@umn.edu), Burrough (1998)

  This calculates a variety of terrain attributes according
  to the work of Burrough 1998's "Principles of Geographical
  Information Systems" (p. 190). Algorithms and helpful ArcGIS
  pages are noted in comments in the code.

  @param[in]  &elevations      An elevation grid
  @param[in]  x0          X coordinate of cell to perform calculation on
  @param[in]  y0          Y coordinate of cell to perform calculation on
  @param[out]  &rise_over_run
    Returns rise-over-run slope as per Horn 1981
  @param[out]  &aspect
    Returns aspect as per Horn 1981 in degrees [0,360). Degrees increase
    in a clockwise fashion. 0 is north, -1 indicates a flat surface.
  @param[out]  &curvature
    Returns the difference of profile and planform curvatures
    (TODO: Clarify, this is from ArcGIS and was poorly described)
  @param[out]  &profile_curvature
    Returns the profile curvature as per Zevenbergen and Thorne 1987.
    0 indicates a flat surface.
  @param[out]  &planform_curvature  
    Returns the planform curvature as per Zevenbergen and Thorne 1987.
    0 indicates a flat surface.

  @pre This function should never be called on a NoData cell
*/
template<class T>
inline static void d8_terrain_attrib_helper(
  const Array2D<T> &elevations,
  int   x0,
  int   y0,
  float zscale,
  float &rise_over_run,
  float &aspect,
  float &curvature,
  float &profile_curvature,
  float &planform_curvature
){
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
  if(elevations.in_grid(x0-1,y0-1))  a = elevations(x0-1,y0-1);
  if(elevations.in_grid(x0-1,y0))    d = elevations(x0-1,y0);
  if(elevations.in_grid(x0-1,y0+1))  g = elevations(x0-1,y0+1);
  if(elevations.in_grid(x0  ,y0-1))  b = elevations(x0,y0-1);
  if(elevations.in_grid(x0  ,y0+1))  h = elevations(x0,y0+1);
  if(elevations.in_grid(x0+1,y0-1))  c = elevations(x0+1,y0-1);
  if(elevations.in_grid(x0+1,y0))    f = elevations(x0+1,y0);
  if(elevations.in_grid(x0+1,y0+1))  i = elevations(x0+1,y0+1);
  if(a==elevations.no_data)          a = e;
  if(b==elevations.no_data)          b = e;
  if(c==elevations.no_data)          c = e;
  if(d==elevations.no_data)          d = e;
  if(f==elevations.no_data)          f = e;
  if(g==elevations.no_data)          g = e;
  if(h==elevations.no_data)          h = e;
  if(i==elevations.no_data)          i = e;

  //TODO: Convert elevations to meters? (1ft=0.3048m)
  a *= zscale;
  b *= zscale;
  c *= zscale;
  d *= zscale;
  e *= zscale;
  f *= zscale;
  g *= zscale;
  h *= zscale;
  i *= zscale;

  double dzdx,dzdy;
  //Aspect calculation in the manner of Horn 1981
  //ArcGIS doesn't use cell size fo aspect calculations.
  dzdx   = ( (c+2*f+i) - (a+2*d+g) ) / 8;
  dzdy   = ( (g+2*h+i) - (a+2*b+c) ) / 8;
  aspect = 180.0/M_PI*atan2(dzdy,-dzdx);
  if(aspect<0)
    aspect = 90-aspect;
  else if(aspect>90.0)
    aspect = 360.0-aspect+90.0;
  else
    aspect = 90.0-aspect;

  //Slope calculation in the manner of Horn 1981
  //But cellsize is accounted for in slope
  dzdx /= elevations.cellsize;
  dzdy /= elevations.cellsize;

  rise_over_run = sqrt(dzdx*dzdx+dzdy*dzdy);

  if(rise_over_run==0){
    aspect             = -1;  //Special value denoting a flat
    curvature          = 0;
    profile_curvature  = 0;
    planform_curvature = 0;
    return;
  }






  //Z1 Z2 Z3   a b c
  //Z4 Z5 Z6   d e f
  //Z7 Z8 Z9   g h i
  //Curvatures in the manner of Zevenbergen and Thorne 1987
  double L  = elevations.cellsize;     //TODO: Should be in the same units as z
  double D  = ( (d+f)/2 - e) / L / L;  //D = [(Z4 + Z6) /2 - Z5] / L^2
  double E  = ( (b+h)/2 - e) / L / L;  //E = [(Z2 + Z8) /2 - Z5] / L^2
  double F  = (-a+c+g-i)/4/L/L;        //F=(-Z1+Z3+Z7-Z9)/(4L^2)
  double G  = (-d+f)/2/L;              //G=(-Z4+Z6)/(2L)
  double H  = (b-h)/2/L;               //H=(Z2-Z8)/(2L)
  curvature = (float)(-2*(D+E)*100);

  if(G==0 && H==0){
    profile_curvature  = 0;
    planform_curvature = 0;
  } else {
    profile_curvature  = (float)(2*(D*G*G+E*H*H+F*G*H)/(G*G+H*H)*100);
    planform_curvature = (float)(-2*(D*H*H+E*G*G-F*G*H)/(G*G+H*H)*100);
  }
}




//d8_terrain_attribute
/**
  @brief  Calculate a variety of terrain attributes
  @author Richard Barnes (rbarnes@umn.edu), Burrough (1998)

  This calculates a variety of terrain attributes according
  to the work of Burrough 1998's "Principles of Geographical
  Information Systems" (p. 190). This function calls
  d8_terrain_attrib_helper to calculate the actual attributes.
  This function may perform some post-processing (such as on
  slope), but it's purpose is essentially to just scan the grid
  and pass off the work to d8_terrain_attrib_helper().

  Possible attribute values are
  <ul>
    <li>#TATTRIB_CURVATURE</li>
    <li>#TATTRIB_PLANFORM_CURVATURE</li>
    <li>#TATTRIB_PROFILE_CURVATURE</li>
    <li>#TATTRIB_ASPECT</li>
    <li>#TATTRIB_SLOPE_RISERUN</li>
    <li>#TATTRIB_SLOPE_PERCENT</li>
    <li>#TATTRIB_SLOPE_RADIAN</li>
    <li>#TATTRIB_SLOPE_DEGREE</li>
  </ul>

  @param[in]  &elevations      An elevation grid
  @param[out]  &attribs      A grid to hold the results
  @param[in]  &attrib        The attribute to be calculated

  @post \pname{attribs} takes the properties and dimensions of \pname{elevations}
*/
template<class T>
void d8_terrain_attribute(
  const Array2D<T> &elevations,
  Array2D<float>   &attribs,
  int   attrib,
  float zscale
){
  ProgressBar progress;

  diagnostic_arg("Setting up the attribute #%d matrix...", attrib);
  attribs.copyprops(elevations);
  attribs.no_data=-99999;  //TODO: Should push this out to the calling helper functions
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
      float rise_over_run, aspect, curvature;
      float profile_curvature, planform_curvature;
      d8_terrain_attrib_helper(
        elevations, x, y, zscale, rise_over_run, aspect, curvature,
        profile_curvature, planform_curvature
      );
      switch(attrib){
        case TATTRIB_CURVATURE:
          attribs(x,y) = curvature;
          break;
        case TATTRIB_PLANFORM_CURVATURE:
          attribs(x,y) = planform_curvature;
          break;
        case TATTRIB_PROFILE_CURVATURE:
          attribs(x,y) = profile_curvature;
          break;
        case TATTRIB_ASPECT:
          attribs(x,y) = aspect;
          break;
        case TATTRIB_SLOPE_RISERUN:
          attribs(x,y) = rise_over_run;
          break;
        case TATTRIB_SLOPE_PERCENT:
          attribs(x,y) = rise_over_run*100;
          break;
        case TATTRIB_SLOPE_RADIAN:
          attribs(x,y) = atan(rise_over_run);
          break;
        case TATTRIB_SLOPE_DEGREE:
          attribs(x,y) = atan(rise_over_run)*180/M_PI;
          break;
      }
    }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
}


//d8_slope
/**
  @brief  Calculates the slope
  @author Richard Barnes (rbarnes@umn.edu), Horn (1981)

  Calculates the slope using Horn 1981, as per Burrough 1998's
  "Principles of Geographical Information Systems" (p. 190)

  Possible slope types are
  <ul>
    <li>#TATTRIB_SLOPE_RISERUN</li>
    <li>#TATTRIB_SLOPE_PERCENT</li>
    <li>#TATTRIB_SLOPE_RADIAN</li>
    <li>#TATTRIB_SLOPE_DEGREE</li>
  </ul>

  @param[in]  &elevations   An elevation grid
  @param[out] &slopes       A slope grid
  @param[in]  &slope_type   A type, as in the description
*/
template<class T>
void d8_slope(
  const Array2D<T> &elevations,
  Array2D<float>   &slopes,
  int slope_type,
  float zscale
){
  diagnostic("\n###Slope attribute calculation\n");
  d8_terrain_attribute(elevations, slopes, slope_type, zscale);
}

template<class T>
void d8_aspect(
  const Array2D<T> &elevations,
  Array2D<float>   &aspects,
  float zscale
){
  diagnostic("\n###Aspect attribute calculation\n");
  d8_terrain_attribute(elevations, aspects, TATTRIB_ASPECT, zscale);
}

void d8_curvature(
  const Array2D<T> &elevations, 
  Array2D<float>   &curvatures, 
  float zscale
){
  diagnostic("\n###Curvature attribute calculation\n");
  d8_terrain_attribute(elevations, curvatures, TATTRIB_CURVATURE, zscale);
}

template<class T>
void d8_planform_curvature(
  const Array2D<T> &elevations,
  Array2D<float>   &planform_curvatures,
  float zscale
){
  diagnostic("\n###Planform curvature attribute calculation\n");
  d8_terrain_attribute(elevations, planform_curvatures,  TATTRIB_PLANFORM_CURVATURE, zscale);
}

template<class T>
void d8_profile_curvature(
  const Array2D<T> &elevations,
  Array2D<float>   &profile_curvatures,
  float zscale
){
  diagnostic("\n###Profile curvature attribute calculation\n");
  d8_terrain_attribute(elevations, profile_curvatures,  TATTRIB_PROFILE_CURVATURE, zscale);
}

#endif
