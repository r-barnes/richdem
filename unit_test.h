#ifndef _unit_test_included
#define _unit_test_included

#include <cmath>
#include "interface.h"

template <class T>
double unit_avg_diff(const array2d<T> &arr1, const array2d<T> &arr2){
	diagnostic("Finding average difference between two arrays...");
	if(arr1.width()!=arr2.width() || arr1.height()!=arr2.height())
		diagnostic("failed! The arrays do not have the same dimensions!\n");

	double diff=0;
	int ccount=0;
	#pragma omp parallel for collapse(2) reduction(+:diff) reduction(+:ccount)
	for(int x=0;x<arr1.width();x++)
	for(int y=0;y<arr2.height();y++){
		if(arr1(x,y)==arr1.no_data || arr2(x,y)==arr2.no_data)
			continue;
//		diagnostic_arg("Percent slopes. RichDEM: %f, ArcGIS: %f\n",arr1(x,y),arr2(x,y));
		diff+=fabs(arr1(x,y)-arr2(x,y));
		ccount++;
	}
	diagnostic("success!\n");

	return (float)(diff/ccount);
}

#endif
