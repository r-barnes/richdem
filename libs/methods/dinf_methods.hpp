#ifndef _richdem_dinf_methods_hpp_
#define _richdem_dinf_methods_hpp_

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
static void where_do_i_flow(float flowdir, int &nhigh, int &nlow){
  //If it is very close to being directed into just one cell
  //then we direct it into just one cell. If we mistakenly direct
  //it into 2 cells, then we may create unresolvable loops in the
  //flow accumulation algorithm, whereas nothing is lost if we
  //leave out one of the two cells (provided a negligible amount
  //of flow is directed to the one we leave out).
  assert(flowdir>=0 && flowdir<=2*M_PI+1e-6);

  flowdir /= (M_PI/4.);

  if(fabs(flowdir-(int)flowdir)<1e-6){
    nlow  = -1;
    nhigh = (int)round(flowdir);
  } else {
    nlow  = (int)flowdir;
    nhigh = nlow+1;
  }

  //8 is not technically a direction, but, since things move in a circle,
  //it overlaps with 0. It should _never_ be greater than 8.
  assert(nhigh>=0 && nhigh<=8);
}

//This reacts correctly if the flow direction wedge number exceeds 7.
static void area_proportion(float flowdir, int nhigh, int nlow, float &phigh, float &plow){
  if(nlow==-1){
    phigh = 1;
    plow  = 0;
  } else {
    phigh = (nhigh*(M_PI/4.0)-flowdir)/(M_PI/4.0);
    plow  = 1-phigh;
  }

  assert(phigh+plow==1);  //TODO: This isn't necessarily so in floating-point... or is it?
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
  ProgressBar progress;

  diagnostic("\n###Dinf Upslope Area\n");

  diagnostic_arg("The sources queue will require at most approximately %ldMB of RAM.\n",flowdirs.width()*flowdirs.height()*((long)sizeof(grid_cell))/1024/1024);

  diagnostic("Setting up the dependency matrix...");
  dependency.copyprops(flowdirs);
  dependency.init(0);
  diagnostic("succeeded.\n");

  diagnostic("Setting up the area matrix...");
  area.copyprops(flowdirs);
  area.init(0);
  area.no_data=dinf_NO_DATA;
  diagnostic("succeeded.\n");

  bool has_cells_without_flow_directions=false;
  diagnostic("%%Calculating dependency matrix & setting no_data cells...\n");
  progress.start( flowdirs.width()*flowdirs.height() );
  #pragma omp parallel for reduction(|:has_cells_without_flow_directions)
  for(int x=0;x<flowdirs.width();x++){
    progress.update( x*flowdirs.height() );
    for(int y=0;y<flowdirs.height();y++){
      if(flowdirs(x,y)==flowdirs.no_data){
        area(x,y)       = area.no_data;
        dependency(x,y) = 9;  //Note: This is an unnecessary safety precaution
        continue;
      }
      if(flowdirs(x,y)==NO_FLOW){
        has_cells_without_flow_directions=true;
        continue;
      }
      int n_high,n_low;
      int nhx,nhy,nlx,nly;
      where_do_i_flow(flowdirs(x,y),n_high,n_low);
      nhx=x+dinf_dx[n_high],nhy=y+dinf_dy[n_high];
      if(n_low!=-1){
        nlx = x+dinf_dx[n_low];
        nly = y+dinf_dy[n_low];
      }
      if( n_low!=-1 && flowdirs.in_grid(nlx,nly) && flowdirs(nlx,nly)!=flowdirs.no_data )
        dependency(nlx,nly)++;
      if( flowdirs.in_grid(nhx,nhy) && flowdirs(nhx,nhy)!=flowdirs.no_data )
        dependency(nhx,nhy)++;
    }
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
  if(has_cells_without_flow_directions)
    diagnostic("\033[91mNot all cells had defined flow directions! This implies that there will be digital dams!\033[39m\n");

  diagnostic("%%Locating source cells...\n");
  progress.start( flowdirs.width()*flowdirs.height() );
  for(int x=0;x<flowdirs.width();x++){
    progress.update( x*flowdirs.height() );
    for(int y=0;y<flowdirs.height();y++)
      if(flowdirs(x,y)==flowdirs.no_data)
        continue;
      else if(flowdirs(x,y)==NO_FLOW)
        continue;
      else if(dependency(x,y)==0)
        sources.push(grid_cell(x,y));
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());

  diagnostic("%%Calculating up-slope areas...\n");
  progress.start( flowdirs.data_cells );
  long int ccount=0;
  while(sources.size()>0){
    grid_cell c=sources.front();
    sources.pop();

    ccount++;
    progress.update(ccount);

    if(flowdirs(c.x,c.y)==flowdirs.no_data)  //TODO: This line shouldn't be necessary since NoData's do not get added below
      continue;

    area(c.x,c.y)+=1;

    if(flowdirs(c.x,c.y)==NO_FLOW)
      continue;

    int n_high,n_low,nhx,nhy,nlx,nly;
    where_do_i_flow(flowdirs(c.x,c.y),n_high,n_low);
    nhx = c.x+dinf_dx[n_high];
    nhy = c.y+dinf_dy[n_high];

    float phigh,plow;
    area_proportion(flowdirs(c.x,c.y), n_high, n_low, phigh, plow);
    if(flowdirs.in_grid(nhx,nhy) && flowdirs(nhx,nhy)!=flowdirs.no_data)
      area(nhx,nhy)+=area(c.x,c.y)*phigh;

    if(n_low!=-1){
      nlx = c.x+dinf_dx[n_low];
      nly = c.y+dinf_dy[n_low];
      if(flowdirs.in_grid(nlx,nly) && flowdirs(nlx,nly)!=flowdirs.no_data){
        area(nlx,nly)+=area(c.x,c.y)*plow;
        if((--dependency(nlx,nly))==0)
          sources.push(grid_cell(nlx,nly));
      }
    }

    if( flowdirs.in_grid(nhx,nhy) && flowdirs(nhx,nhy)!=flowdirs.no_data && (--dependency(nhx,nhy))==0)
      sources.push(grid_cell(nhx,nhy));
  }
  diagnostic_arg(SUCCEEDED_IN,progress.stop());
}




static const float d8_to_dinf[9]={-1, 4*M_PI/4, 3*M_PI/4, 2*M_PI/4, 1*M_PI/4, 0, 7*M_PI/4, 6*M_PI/4, 5*M_PI/4};



#endif
