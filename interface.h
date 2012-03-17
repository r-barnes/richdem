#ifndef _interface_included
#define _interface_included

#include <cstdio>
#include <sys/time.h>
#define diagnostic_arg(message,...) fprintf(stderr,message,__VA_ARGS__)
#define diagnostic(message) fprintf(stderr,message)
double progress_bar(int percent);
double timediff(const timeval &startTime);

#endif
