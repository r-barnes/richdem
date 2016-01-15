/**
  @file
  Defines timing functions and progress bars.

  Richard Barnes (rbarnes@umn.edu), 2016
*/
#ifndef _richdem_utility_hpp_
#define _richdem_utility_hpp_

#include <string>
#include <sys/time.h>

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num()  0
  #define omp_get_num_threads() 1
#endif

class Timer{
  private:
    timeval start_time; ///<Last time the timer was started
    double accumulated_time; ///<Accumulated running time since creation
    bool running; ///<True when the timer is running
    ///Number of seconds between two time objects
    double timediff(timeval beginning, timeval end){
      long seconds, useconds;
      seconds  = end.tv_sec  - beginning.tv_sec;
      useconds = end.tv_usec - beginning.tv_usec;
      return seconds + useconds/1000000.0;
    }
  public:
    Timer(){
      accumulated_time=0;
      running=false;
    }
    void start(){
      if(running)
        throw "Timer was already started!";
      running=true;
      gettimeofday(&start_time, NULL);
    }
    void stop(){
      if(!running)
        throw "Timer was already stopped!";
      running=false;
      timeval end_time;
      gettimeofday(&end_time, NULL);

      accumulated_time+=timediff(start_time,end_time);
    }
    double accumulated(){
      if(running)
        throw "Timer is still running!";
      return accumulated_time;
    }
    double lap(){
      if(!running)
        throw "Timer was not started!";
      timeval lap_time;
      gettimeofday(&lap_time, NULL);
      return timediff(start_time,lap_time);
    }
};

#define PROGRESS_BAR "=================================================="

class ProgressBar{
  private:
    long total_work;
    int old_percent;
    Timer timer;
  public:
    void start(long total_work0){
      timer.start();
      total_work  = total_work0;
      old_percent = 0;
      //This ANSI code clears the entire line
      std::cerr<<"\r\033[2K"<<std::flush;
    }
    void update(long work_done){
      if(omp_get_thread_num()!=0)
        return;

      int percent;
      percent=(int)(work_done*omp_get_num_threads()*100/total_work);
      if(percent>100)
        percent=100;
      if(percent==old_percent)
        return;
      old_percent=percent;

      std::cerr<<"\r\033[2K["
               <<std::string(percent/2, '=')
               <<"] ("
               <<percent<<"% - "<<timer.lap()/percent*(100-percent)
               <<omp_get_num_threads()<< " threads)"<<std::flush;
    }
    double stop(){
      //This ANSI code clears the entire line
      std::cerr<<"\r\033[2K"<<std::flush;

      timer.stop();
      return timer.accumulated();
    }
    double time_it_took(){
      return timer.accumulated();
    }
};

#endif
