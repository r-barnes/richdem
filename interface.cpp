#include <stdio.h>
#include "interface.h"

#define PROGRESS_BAR "=================================================="

void progress_bar(int percent){
	static int old_percent=-2;
	if(percent==old_percent)
		return;
	if(percent==-1){
		fprintf(stderr,"\r %-50.*s         \r",0,PROGRESS_BAR);
		old_percent=0;
		return;
	}
	if(percent>100)
		percent=100;
	fprintf(stderr,"\r[%-50.*s] (%d%%) ",percent/2,PROGRESS_BAR,percent);
	old_percent=percent;
}
