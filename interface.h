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



class ProgressBar{
	private:
		long total_work;
		long old_percent;
		Timer timer;
	public:
		ProgressBar(int total_work0) : total_work(total_work0){}
		void start(){
			if(omp_get_thread_num()!=0)
				return;

			timer.start();
			old_percent=0;
			#ifndef ARCGIS
				//This ANSI code clears the entire line
				fprintf(stderr,"\r\033[2K");
			#else
				//The ArcGIS python script resets its progress bar when it sees this
				fprintf(stderr,"P%%c\n");
			#endif
		}
		void update(long work_done){
			if(omp_get_thread_num()!=0)
				return;

			long percent;
			percent=work_done*omp_get_num_threads()*100/work_total;
			if(percent>100)
				percent=100;
			if(percent==old_percent)
				return;
			old_percent=percent;

			#ifndef ARCGIS
				fprintf(stderr,"\r\033[2K[%-50.*s] (%d%% - %.1lfs left - %d thread)", percent/2, PROGRESS_BAR, percent, timer.lap()/percent*(100-percent),omp_get_num_threads()); //This ANSI code clears the entire line
			#else
				fprintf(stderr,"P%%%d\n", percent); //The ArcGIS python script adjusts its progress bar when it sees this
			#endif
		}
		double stop(){
			if(omp_get_thread_num()!=0)
				return;

			#ifndef ARCGIS
				//This ANSI code clears the entire line
				fprintf(stderr,"\r\033[2K");
			#else
				//The ArcGIS python script resets its progress bar when it sees this
				fprintf(stderr,"P%%c\n");
			#endif

			timer.stop();
			return timer.accumulated();
		}
		double time_it_took(){
			return timer.accumulated();
		}
};



#endif
