#include <stdio.h>
#include "interface.h"
#include <sys/time.h>
#include <omp.h>

#define PROGRESS_BAR "=================================================="

double timediff(const timeval &startTime){
	timeval endTime;
	long seconds, useconds;
	double duration;

	gettimeofday(&endTime, NULL);

	seconds  = endTime.tv_sec  - startTime.tv_sec;
	useconds = endTime.tv_usec - startTime.tv_usec;

	return seconds + useconds/1000000.0;
}

void progress_bar(int percent){
	static int old_percent=-2;
	static timeval startTime;

	if(omp_get_thread_num()!=0) return;

	if(percent==old_percent)
		return;
	if(percent==-1){
		gettimeofday(&startTime, NULL);
		fprintf(stderr,"\r\033[2K",0,PROGRESS_BAR); //The ANSI code used in this line clears the entire line
		old_percent=0;
		return;
	}
	if(percent>100)
		percent=100;
	fprintf(stderr,"\r\033[2K[%-50.*s] (%d%% - %.1lfs left - %d thread)", percent/2, PROGRESS_BAR, percent, timediff(startTime)/percent*(100-percent),omp_get_num_threads()); //The ANSI code used in this line clears the entire line
	old_percent=percent;
}
