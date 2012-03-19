#ifndef _debug_included
#define _debug_included

#include <deque>

#define PRINT(ARR,PREC,WIDTH) std::cout<<std::setprecision(PREC)<<std::setw(WIDTH)<<ARR<<std::endl;

void print_edges(float_2d &elevations, std::deque<grid_cell> &low_edges, std::deque<grid_cell> &high_edges);

template <class T> //TODO: Needs error checking for dimensions, et cetera
void array2ddiff(const array2d<T> &arr1, const array2d<T> &arr2, array2d<T> &result){
	for(int x=0;x<arr1.width();x++)
	for(int y=0;y<arr2.height();y++)
		if(arr1(x,y)==arr1.no_data || arr2(x,y)==arr2.no_data)
			result(x,y)=arr1.no_data;
		else
			result(x,y)=arr1(x,y)-arr2(x,y);
}

#endif
