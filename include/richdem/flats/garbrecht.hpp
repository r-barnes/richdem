#ifndef _richdem_garbrecht_hpp_
#define _richdem_garbrecht_hpp_

#include <deque>
#include <cstdint>
#include <iostream>
#include "richdem/common/Array2D.hpp"
#include "richdem/common/grid_cell.hpp"
#include "richdem/flowmet/d8_flowdirs.hpp"
#include "richdem/common/logger.hpp"

namespace richdem {

typedef std::deque<grid_cell> flat_type;

void FindFlats(
	const Array2D<uint8_t> &flowdirs,
	flat_type &flats
){
	for(int x=0;x<flowdirs.width();++x)
	for(int y=0;y<flowdirs.height();++y)
		if(flowdirs(x,y)==NO_FLOW)
			flats.push_back(grid_cell(x,y));
}

template<class T>
void GradientTowardsLower(
	const Array2D<T>       &elevations,
	const Array2D<uint8_t> &flowdirs,
	flat_type              &flats,
	Array2D<int32_t>       &inc1
){
	int loops              = 0;
	int number_incremented = -1;
	RDLOG_PROGRESS<<"Setting up the inc1 matrix...";
  inc1.resize(elevations);
  inc1.init(0);
	RDLOG_PROGRESS<<"Calculating inc1 matrix...";
	while(number_incremented!=0){
		number_incremented=0;
		for(int i=0;i<(int)flats.size();i++){
			bool increment_elevation=true;
			int x=flats[i].x;
			int y=flats[i].y;
			for(int n=1;n<=8;n++){
				if(	elevations(x+dx[n],y+dy[n])<elevations(x,y) && 
					flowdirs(x+dx[n],y+dy[n])!=NO_FLOW && flowdirs(x+dx[n],y+dy[n])!=flowdirs.noData()
				){
					increment_elevation=false;
					break;
				}
				else if(inc1(x+dx[n],y+dy[n])<loops && 
						elevations(x+dx[n],y+dy[n])==elevations(x,y)
				){
					increment_elevation=false;
					break;
				}	
			}
			if(increment_elevation){
				inc1(x,y)++;
				number_incremented++;
			}
		}
		loops++;
	}
}

template<class T>
void GradientAwayFromHigher(
	const Array2D<T>       &elevations,
	const Array2D<uint8_t> &flowdirs,
	flat_type              &flats,
	Array2D<int32_t>       &inc2
){
	int loops                       = 0;
	unsigned int number_incremented = 0;
  RDLOG_PROGRESS<<"Setting up the inc2 matrix...";
  inc2.resize(elevations);
  inc2.init(0);
  RDLOG_PROGRESS<<"Calculting inc2 matrix...";
	while(number_incremented<flats.size()){
		for(int i=0;i<(int)flats.size();i++){
			int x=flats[i].x;
			int y=flats[i].y;
			if(inc2(x,y)>0){
				inc2(x,y)++;
				continue;
			}
		}
		for(int i=0;i<(int)flats.size();i++){
			bool has_higher=false,has_lower=false;
			int x=flats[i].x;
			int y=flats[i].y;
			if(inc2(x,y)>0) continue;
			for(int n=1;n<=8;n++){
				if( !has_higher &&
					(elevations(x+dx[n],y+dy[n])>elevations(x,y) ||
					inc2(x+dx[n],y+dy[n])==2)
				)
					has_higher=true;
				else if( !has_lower &&
						elevations(x+dx[n],y+dy[n])<elevations(x,y)
				)
					has_lower=true;
			}
			if(has_higher && !has_lower){
				inc2(x,y)++;
				number_incremented++;
			}
		}
		loops++;
	}
}

template<class T>
void CombineGradients(
	Array2D<T>             &elevations,
	const Array2D<int32_t> &inc1,
	const Array2D<int32_t> &inc2, 
	float epsilon //TODO
){
  RDLOG_PROGRESS<<"Combining the gradients...";
	for(int x=0;x<elevations.width();++x)
	for(int y=0;y<elevations.height();++y)
		elevations(x,y)+=(inc1(x,y)+inc2(x,y))*epsilon;
}



template<class T>
void GarbrechtAlg(Array2D<T> &elevations, Array2D<uint8_t> &flowdirs){
  Timer flat_resolution_timer;

  flat_resolution_timer.start();
  flat_type flats;
  FindFlats(flowdirs,flats);

  Array2D<int32_t> inc1, inc2;
  GradientTowardsLower  (elevations, flowdirs, flats, inc1);
  GradientAwayFromHigher(elevations, flowdirs, flats, inc2);
  CombineGradients(elevations, inc1, inc2, 0.001);
  flat_resolution_timer.stop();

  d8_flow_directions(elevations,flowdirs);

  std::cout<<flat_resolution_timer.accumulated()<<" seconds were used to resolve flats."<<std::endl;
}

}

#endif
