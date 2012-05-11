#ifndef _interface_included
#define _interface_included

#include <cstdio>
#include <sys/time.h>
#define diagnostic_arg(message,...) fprintf(stderr,message,__VA_ARGS__)
#define diagnostic(message) fprintf(stderr,message)
double progress_bar(int percent);
double timediff(const timeval &startTime);

#ifndef ARCGIS
	#define SUCCEEDED_IN "\t\033[96msucceeded in %.2lfs.\033[39m\n"
#else
	#define SUCCEEDED_IN "\tsucceeded in %.2lfs.\n"
#endif

#endif
