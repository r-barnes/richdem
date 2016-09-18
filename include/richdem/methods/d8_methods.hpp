/**
  @file
  @brief Defines a number of functions for calculating terrain attributes

  Richard Barnes (rbarnes@umn.edu), 2015
*/
//TODO: Curvature may be "Least Squares Fitted Plane" per http://gis.stackexchange.com/questions/37066/how-to-calculate-terrain-curvature
#ifndef _richdem_d8_methods_hpp_
#define _richdem_d8_methods_hpp_

#include "richdem/common/Array2D.hpp"
#include "richdem/common/constants.hpp"
#include "richdem/common/grid_cell.hpp"
#include "richdem/common/ProgressBar.hpp"
#include <queue>

/**
  @brief  Returns the sign (+1, -1, 0) of a number. Branchless.
  @author Richard Barnes (rbarnes@umn.edu)

  @param[in]  val  Input value

  @return
    -1 for a negative input, +1 for a positive input, and 0 for a zero input
*/
template <class T>
static T sgn(T val){
  return (T(0) < val) - (val < T(0));
}






/**
  @brief  Calculates the D8 flow accumulation, given the D8 flow directions
  @author Richard Barnes (rbarnes@umn.edu)

  This calculates the D8 flow accumulation of a grid of D8 flow directions by
  calculating each cell's dependency on its neighbours and then using a
  priority-queue to process cells in a top-of-the-watershed-down fashion

  @param[in]  &flowdirs  A D8 flowdir grid from d8_flow_directions()
  @param[out] &area      Returns the up-slope area of each cell
*/
template<class T, class U>
void d8_flow_accum(const Array2D<T> &flowdirs, Array2D<U> &area){
  std::queue<GridCell> sources;
  ProgressBar progress;

  std::cerr<<"\nA D8 Flow Accumulation"<<std::endl;
  std::cerr<<"C TODO"<<std::endl;

  std::cerr<<"The sources queue will require at most approximately "
           <<(flowdirs.size()*((long)sizeof(GridCell))/1024/1024)
           <<"MB of RAM."<<std::endl;

  std::cerr<<"p Resizing dependency matrix..."<<std::endl;
  Array2D<int8_t> dependency(flowdirs,0);

  std::cerr<<"p Setting up the area matrix..."<<std::endl;
  area.resize(flowdirs,0);
  area.setNoData(-1);

  std::cerr<<"p Calculating dependency matrix & setting noData() cells..."<<std::endl;
  progress.start( flowdirs.size() );
  #pragma omp parallel for
  for(int y=0;y<flowdirs.height();y++){
    progress.update( y*flowdirs.width() );
    for(int x=0;x<flowdirs.width();x++){
      if(flowdirs.isNoData(x,y)){
        area(x,y) = area.noData();
        continue;
      }

      int n = flowdirs(x,y); //The neighbour this cell flows into
      if(n==NO_FLOW)         //This cell does not flow into a neighbour
        continue;

      int nx = x+dx[n];      //x-coordinate of the neighbour
      int ny = y+dy[n];      //y-coordinate of the neighbour

      //Neighbour is not on the grid
      if(!flowdirs.inGrid(nx,ny))
        continue;

      //Neighbour is valid and is part of the grid. The neighbour depends on this
      //cell, so increment its dependency count.
      ++dependency(nx,ny);
    }
  }
  std::cerr<<"t Dependency calculation time = "<<progress.stop()<<" s"<<std::endl;

  std::cerr<<"p Locating source cells..."<<std::endl;
  for(int y=0;y<flowdirs.height();y++)
  for(int x=0;x<flowdirs.width();x++)
    if(dependency(x,y)==0 && !flowdirs.isNoData(x,y))
      sources.emplace(x,y);

  std::cerr<<"p Calculating flow accumulation areas..."<<std::endl;
  progress.start(flowdirs.numDataCells());
  long int ccount=0;
  while(sources.size()>0){
    GridCell c=sources.front();
    sources.pop();

    ccount++;
    progress.update(ccount);

    area(c.x,c.y)++;

    int n = flowdirs(c.x,c.y);

    if(n==NO_FLOW)
      continue;

    int nx=c.x+dx[n];
    int ny=c.y+dy[n];

    if(!flowdirs.inGrid(nx,ny))
      continue;
    if(flowdirs.isNoData(nx,ny))
      continue;

    area(nx,ny)+=area(c.x,c.y);
    --dependency(nx,ny);

    if(dependency(nx,ny)==0)
      sources.emplace(nx,ny);
  }
  std::cerr<<"t Flow accumulation calculation time = "<<progress.stop()<<" s"<<std::endl;

  //TODO: Explain this better
  int loops=0;
  for(int i=-1;i>=-8;i--)
    loops+=dependency.countval(-1);
  std::cerr<<"m Input contained at least = "<<loops<<" loops"<<std::endl;
}




//d8_upslope_cells
/**
  @brief  Calculates which cells ultimately D8-flow through a given cell
  @author Richard Barnes (rbarnes@umn.edu)

  Given the coordinates x, y of a cell, this returns a grid indicating
  which cells ultimately flow into the indicated cell.
  1=Upslope cell
  2=Member of initializing line
  All other cells have a noData() value

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
  std::cerr<<"Setting up the upslope_cells matrix..."<<std::flush;
  upslope_cells.resize(flowdirs);
  upslope_cells.setAll(FLOWDIR_NO_DATA);
  upslope_cells.setNoData(FLOWDIR_NO_DATA);
  std::cerr<<"succeeded."<<std::endl;
  ProgressBar progress;

  std::queue<GridCell> expansion;

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

  std::cerr<<"Line slope is "<<deltaerr<<std::endl;
  int y=y0;
  for(int x=x0;x<=x1;x++){
    expansion.push(GridCell(x,y));
    upslope_cells(x,y)=2;
    error+=deltaerr;
    if (error>=0.5) {
      expansion.push(GridCell(x+1,y));
      upslope_cells(x+1,y) = 2;
      y                   += sgn(deltay);
      error               -= 1;
    }
  }

  progress.start(flowdirs.data_cells);
  long int ccount=0;
  while(expansion.size()>0){
    GridCell c=expansion.front();
    expansion.pop();

    progress.update(ccount++);

    for(int n=1;n<=8;n++)
      if(!flowdirs.inGrid(c.x+dx[n],c.y+dy[n]))
        continue;
      else if(flowdirs(c.x+dx[n],c.y+dy[n])==NO_FLOW)
        continue;
      else if(flowdirs(c.x+dx[n],c.y+dy[n])==flowdirs.noData())
        continue;
      else if(upslope_cells(c.x+dx[n],c.y+dy[n])==upslope_cells.noData() && n==d8_inverse[flowdirs(c.x+dx[n],c.y+dy[n])]){
        expansion.push(GridCell(c.x+dx[n],c.y+dy[n]));
        upslope_cells(c.x+dx[n],c.y+dy[n])=1;
      }
  }
  std::cerr<<"Succeeded in "<<progress.stop()<<std::endl;
  std::cerr<<"Found "<<ccount<<" up-slope cells."<<std::endl;
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
  const Array2D<U> &riserun_slope,
        Array2D<V> &result
){
  Timer timer;

  std::cerr<<"\n###d8_SPI"<<std::endl;

  if(flow_accumulation.width()!=riserun_slope.width() || flow_accumulation.height()!=riserun_slope.height()){
    std::cerr<<"Couldn't calculate SPI! The input matricies were of unequal dimensions!"<<std::endl;
    exit(-1);
  }

  std::cerr<<"Setting up the SPI matrix..."<<std::flush;
  result.resize(flow_accumulation);
  result.noData()=-1;  //Log(x) can't take this value of real inputs, so we're good
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"Calculating SPI..."<<std::endl;
  timer.start();
  #pragma omp parallel for collapse(2)
  for(int x=0;x<flow_accumulation.width();x++)
    for(int y=0;y<flow_accumulation.height();y++)
      if(flow_accumulation(x,y)==flow_accumulation.noData() || riserun_slope(x,y)==riserun_slope.noData())
        result(x,y)=result.noData();
      else
        result(x,y)=log( (flow_accumulation(x,y)/flow_accumulation.getCellArea()) * (riserun_slope(x,y)+0.001) );
  std::cerr<<"succeeded in "<<timer.stop()<<"s."<<std::endl;
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
  const Array2D<U> &riserun_slope,
        Array2D<V> &result
){
  Timer timer;

  std::cerr<<"\n###d8_CTI"<<std::endl;

  if(flow_accumulation.width()!=riserun_slope.width() || flow_accumulation.height()!=riserun_slope.height()){
    std::cerr<<"Couldn't calculate CTI! The input matricies were of unequal dimensions!"<<std::endl;
    exit(-1);
  }

  std::cerr<<"Setting up the CTI matrix..."<<std::flush;
  result.resize(flow_accumulation);
  result.setNoData(-1);  //Log(x) can't take this value of real inputs, so we're good
  std::cerr<<"succeeded."<<std::endl;

  std::cerr<<"Calculating CTI..."<<std::flush;
  timer.start();
  #pragma omp parallel for collapse(2)
  for(int x=0;x<flow_accumulation.width();x++)
    for(int y=0;y<flow_accumulation.height();y++)
      if(flow_accumulation(x,y)==flow_accumulation.noData() || riserun_slope(x,y)==riserun_slope.noData())
        result(x,y)=result.noData();
      else
        result(x,y)=log( (flow_accumulation(x,y)/flow_accumulation.getCellArea()) / (riserun_slope(x,y)+0.001) );
  std::cerr<<"succeeded in "<<timer.stop()<<"s."<<std::endl;
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
class TerrainAttributator {
 private:
  double a,b,c,d,e,f,g,h,i;
  double L,D,E,F,G,H;
  double zscale;

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
  */

  //a b c
  //d e f
  //g h i
  //Deal with grid edges and NoData values in the manner suggested by
  //ArcGIS. Note that this function should never be called on a NoData cell
  void setup(const Array2D<T> &elevations, const int x, const int y){
    a=b=c=d=e=f=g=h=i=elevations(x,y);
    if(elevations.inGrid(x-1,y-1))  a = elevations(x-1,y-1);
    if(elevations.inGrid(x-1,y  ))  d = elevations(x-1,y  );
    if(elevations.inGrid(x-1,y+1))  g = elevations(x-1,y+1);
    if(elevations.inGrid(x  ,y-1))  b = elevations(x,  y-1);
    if(elevations.inGrid(x  ,y+1))  h = elevations(x,  y+1);
    if(elevations.inGrid(x+1,y-1))  c = elevations(x+1,y-1);
    if(elevations.inGrid(x+1,y  ))  f = elevations(x+1,y  );
    if(elevations.inGrid(x+1,y+1))  i = elevations(x+1,y+1);
    if(a==elevations.noData())      a = e;
    if(b==elevations.noData())      b = e;
    if(c==elevations.noData())      c = e;
    if(d==elevations.noData())      d = e;
    if(f==elevations.noData())      f = e;
    if(g==elevations.noData())      g = e;
    if(h==elevations.noData())      h = e;
    if(i==elevations.noData())      i = e;

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
  }

  void curvatureSetup(const Array2D<T> &elevations, int x, int y){
    setup(elevations,x,y);

    //Z1 Z2 Z3   a b c
    //Z4 Z5 Z6   d e f
    //Z7 Z8 Z9   g h i
    //Curvatures in the manner of Zevenbergen and Thorne 1987
    double L  = elevations.getCellArea(); //TODO: Should be in the same units as z
    double D  = ( (d+f)/2 - e) / L / L;   //D = [(Z4 + Z6) /2 - Z5] / L^2
    double E  = ( (b+h)/2 - e) / L / L;   //E = [(Z2 + Z8) /2 - Z5] / L^2
    double F  = (-a+c+g-i)/4/L/L;         //F=(-Z1+Z3+Z7-Z9)/(4L^2)
    double G  = (-d+f)/2/L;               //G=(-Z4+Z6)/(2L)
    double H  = (b-h)/2/L;                //H=(Z2-Z8)/(2L)
  }

 public:
  typedef double (TerrainAttributator<T>::*FcnPtr)(const Array2D<T> &elevations, int x, int y);

  TerrainAttributator(double zscale){
    this->zscale = zscale;
  }

  ///@brief  Calculates aspect in degrees in the manner of Horn 1981
  ///@return Aspect in degrees in the manner of Horn 1981
  //ArcGIS doesn't use cell size for aspect calculations.
  double aspect(const Array2D<T> &elevations, int x, int y){
    setup(elevations,x,y);

    double dzdx = ( (c+2*f+i) - (a+2*d+g) ) / 8; //TODO? Divide by delta x, according to Horn
    double dzdy = ( (g+2*h+i) - (a+2*b+c) ) / 8; //TODO? Divide by delta y, according to Horn
    aspect      = 180.0/M_PI*atan2(dzdy,-dzdx);
    if(aspect<0)
      return 90-aspect;
    else if(aspect>90.0)
      return 360.0-aspect+90.0;
    else
      return 90.0-aspect;
  }

  ///@brief  Calculates the rise/run slope along the maximum gradient on a fitted surface over a 3x3 be neighbourhood in the manner of Horn 1981
  ///@return Rise/run slope
  double slope_riserun(const Array2D<T> &elevations, int x, int y){
    setup(elevations,x,y);

    //But cellsize is accounted for in slope
    double dzdx = ( (c+2*f+i) - (a+2*d+g) ) / 8; //TODO? Divide by delta x, according to Horn
    double dzdy = ( (g+2*h+i) - (a+2*b+c) ) / 8; //TODO? Divide by delta y, according to Horn

    //TODO: Incorporate zscale. The above should do it.
    //The above fits are surface to a 3x3 neighbour hood. This returns the slope
    //along the direction of maximum gradient.
    return sqrt(dzdx*dzdx+dzdy*dzdy);
  }

  double curvature(const Array2D<T> &elevations, int x, int y){
    curvatureSetup(elevations,x,y);

    return (-2*(D+E)*100);
  }

  double planform_curvature(const Array2D<T> &elevations, int x, int y){
    curvatureSetup(elevations,x,y);

    if(G==0 && H==0)
      return 0;
    else
      return (-2*(D*H*H+E*G*G-F*G*H)/(G*G+H*H)*100);
  }

  double profile_curvature(const Array2D<T> &elevations, int x, int y){
    curvatureSetup(elevations,x,y);

    if(G==0 && H==0)
      return 0;
    else
      return (2*(D*G*G+E*H*H+F*G*H)/(G*G+H*H)*100);
  }

  double slope_percent(const Array2D<T> &elevations, int x, int y){
    return slope_riserun(elevations,x,y)*100;
  }

  double slope_radian(const Array2D<T> &elevations, int x, int y){
    return atan(slope_riserun(elevations,x,y));
  }

  double slope_degree(const Array2D<T> &elevations, int x, int y){
    return atan(slope_riserun(elevations,x,y))*180/M_PI;
  }

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
  void process(const Array2D<T> &elevations, Array2D<float> &attribs, FcnPtr fcn){
    attribs.resize(elevations);
    attribs.setNoData(-9999);  //TODO: Should push this out to the calling helper functions

    for(int y=0;y<elevations.height();y++)
    for(int x=0;x<elevations.width();x++)
      if(elevations.isNoData(x,y))
        attribs(x,y) = attribs.noData();
      else
        attribs(x,y) = (this->*fcn)(elevations,x,y);
  }
};



/**
  @brief  Calculates the slope as rise/run
  @author Richard Barnes (rbarnes@umn.edu), Horn (1981)

  Calculates the slope using Horn 1981, as per Burrough 1998's
  "Principles of Geographical Information Systems" (p. 190)

  @param[in]  &elevations   An elevation grid
  @param[out] &slopes       A slope grid
*/
template<class T>
void d8_slope_riserun(
  const Array2D<T> &elevations,
  Array2D<float>   &slopes,
  float zscale
){
  std::cerr<<"\nA Slope calculation (rise/run)"<<std::endl;
  std::cerr<<"C Horn 1981 (TODO)"<<std::endl;
  TerrainAttributator<T> ta(zscale);
  ta.process(elevations, slopes, &TerrainAttributator<T>::slope_riserun);
}

/**
  @brief  Calculates the slope as percentage
  @author Richard Barnes (rbarnes@umn.edu), Horn (1981)

  Calculates the slope using Horn 1981, as per Burrough 1998's
  "Principles of Geographical Information Systems" (p. 190)

  @param[in]  &elevations   An elevation grid
  @param[out] &slopes       A slope grid
*/
template<class T>
void d8_slope_percentage(
  const Array2D<T> &elevations,
  Array2D<float>   &slopes,
  float zscale
){
  std::cerr<<"\nA Slope calculation (percenage)"<<std::endl;
  std::cerr<<"C Horn 1981 (TODO)"<<std::endl;
  TerrainAttributator<T> ta(zscale);
  ta.process(elevations, slopes, &TerrainAttributator<T>::slope_percent);
}

/**
  @brief  Calculates the slope as degrees
  @author Richard Barnes (rbarnes@umn.edu), Horn (1981)

  Calculates the slope using Horn 1981, as per Burrough 1998's
  "Principles of Geographical Information Systems" (p. 190)

  @param[in]  &elevations   An elevation grid
  @param[out] &slopes       A slope grid
*/
template<class T>
void d8_slope_degrees(
  const Array2D<T> &elevations,
  Array2D<float>   &slopes,
  float zscale
){
  std::cerr<<"\nA Slope calculation (degrees)"<<std::endl;
  std::cerr<<"C Horn, B.K.P., 1981. Hill shading and the reflectance map. Proceedings of the IEEE 69, 14–47. doi:10.1109/PROC.1981.11918"<<std::endl;
  TerrainAttributator<T> ta(zscale);
  ta.process(elevations, slopes, &TerrainAttributator<T>::slope_degree);
}

/**
  @brief  Calculates the slope as radians
  @author Richard Barnes (rbarnes@umn.edu), Horn (1981)

  Calculates the slope using Horn 1981, as per Burrough 1998's
  "Principles of Geographical Information Systems" (p. 190)

  @param[in]  &elevations   An elevation grid
  @param[out] &slopes       A slope grid
*/
template<class T>
void d8_slope_radians(
  const Array2D<T> &elevations,
  Array2D<float>   &slopes,
  float zscale
){
  std::cerr<<"\nA Slope calculation (radians)"<<std::endl;
  std::cerr<<"C Horn, B.K.P., 1981. Hill shading and the reflectance map. Proceedings of the IEEE 69, 14–47. doi:10.1109/PROC.1981.11918"<<std::endl;
  TerrainAttributator<T> ta(zscale);
  ta.process(elevations, slopes, &TerrainAttributator<T>::slope_radian);
}

/**
  @brief  Calculates the terrain aspect
  @author Richard Barnes (rbarnes@umn.edu), Horn (1981)

  Calculates the aspect per Horn 1981, as described by Burrough 1998's
  "Principles of Geographical Information Systems" (p. 190)

  @param[in]  &elevations   An elevation grid
  @param[out] &aspects      An aspect grid
*/
template<class T>
void d8_aspect(
  const Array2D<T> &elevations,
  Array2D<float>   &aspects,
  float zscale
){
  std::cerr<<"\nA Aspect attribute calculation"<<std::endl;
  std::cerr<<"C Horn, B.K.P., 1981. Hill shading and the reflectance map. Proceedings of the IEEE 69, 14–47. doi:10.1109/PROC.1981.11918"<<std::endl;
  TerrainAttributator<T> ta(zscale);
  ta.process(elevations, aspects, &TerrainAttributator<T>::aspect);
}

/**
  @brief  Calculates the terrain curvature per Zevenbergen and Thorne 1987
  @author Richard Barnes (rbarnes@umn.edu), Horn (1981)

  Calculates the curvature per Zevenbergen and Thorne 1987, as described by
  Burrough 1998's "Principles of Geographical Information Systems" (p. 190)

  @param[in]  &elevations   An elevation grid
  @param[out] &curvatures   A curvature grid
*/
template<class T>
void d8_curvature(
  const Array2D<T> &elevations, 
  Array2D<float>   &curvatures, 
  float zscale
){
  std::cerr<<"\nA Curvature attribute calculation"<<std::endl;
  std::cerr<<"C Zevenbergen, L.W., Thorne, C.R., 1987. Quantitative analysis of land surface topography. Earth surface processes and landforms 12, 47–56."<<std::endl;
  TerrainAttributator<T> ta(zscale);
  ta.process(elevations, curvatures, &TerrainAttributator<T>::curvature);
}


/**
  @brief  Calculates the terrain planform curvature per Zevenbergen and Thorne 1987
  @author Richard Barnes (rbarnes@umn.edu), Horn (1981)

  Calculates the curvature per Zevenbergen and Thorne 1987, as described by
  Burrough 1998's "Principles of Geographical Information Systems" (p. 190)

  @param[in]  &elevations          An elevation grid
  @param[out] &planform_curvatures A planform curvature grid
*/
template<class T>
void d8_planform_curvature(
  const Array2D<T> &elevations,
  Array2D<float>   &planform_curvatures,
  float zscale
){
  std::cerr<<"\nA Planform curvature attribute calculation"<<std::endl;
  std::cerr<<"C Zevenbergen, L.W., Thorne, C.R., 1987. Quantitative analysis of land surface topography. Earth surface processes and landforms 12, 47–56."<<std::endl;
  TerrainAttributator<T> ta(zscale);
  ta.process(elevations, planform_curvatures, &TerrainAttributator<T>::planform_curvature);
}

/**
  @brief  Calculates the terrain profile curvature per Zevenbergen and Thorne 1987
  @author Richard Barnes (rbarnes@umn.edu), Horn (1981)

  Calculates the curvature per Zevenbergen and Thorne 1987, as described by
  Burrough 1998's "Principles of Geographical Information Systems" (p. 190)

  @param[in]  &elevations         An elevation grid
  @param[out] &profile_curvatures A profile curvature grid
*/
template<class T>
void d8_profile_curvature(
  const Array2D<T> &elevations,
  Array2D<float>   &profile_curvatures,
  float zscale
){
  std::cerr<<"\nA Profile curvature attribute calculation"<<std::endl;
  std::cerr<<"C Zevenbergen, L.W., Thorne, C.R., 1987. Quantitative analysis of land surface topography. Earth surface processes and landforms 12, 47–56."<<std::endl;
  TerrainAttributator<T> ta(zscale);
  ta.process(elevations, profile_curvatures, &TerrainAttributator<T>::profile_curvature);
}

#endif
