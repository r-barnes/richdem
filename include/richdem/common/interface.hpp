/**
  @file
  Defines timing functions and progress bars.

  Richard Barnes (rbarnes@umn.edu), 2016
*/
#ifndef _richdem_utility_hpp_
#define _richdem_utility_hpp_

#include <string>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <stdexcept>

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num()  0
  #define omp_get_num_threads() 1
#endif

class Timer{
 private:
  timeval start_time;      ///<Last time the timer was started
  double accumulated_time; ///<Accumulated running time since creation
  bool running;            ///<True when the timer is running

  ///Number of seconds between two time objects
  double timediff(timeval beginning, timeval end){
    long seconds, useconds;
    seconds  = end.tv_sec  - beginning.tv_sec;
    useconds = end.tv_usec - beginning.tv_usec;
    return seconds + useconds/1000000.0;
  }
 public:
  Timer(){
    accumulated_time = 0;
    running          = false;
  }
  void start(){
    if(running)
      std::logic_error("Timer was already started!");
    running=true;
    gettimeofday(&start_time, NULL);
  }
  double stop(){
    if(!running)
      std::logic_error("Timer was already stopped!");
    running=false;
    timeval end_time;
    gettimeofday(&end_time, NULL);

    accumulated_time += timediff(start_time,end_time);

    return accumulated_time;
  }
  double accumulated(){
    if(running)
      std::logic_error("Timer is still running!");
    return accumulated_time;
  }
  double lap(){
    if(!running)
      std::logic_error("Timer was not started!");
    timeval lap_time;
    gettimeofday(&lap_time, NULL);
    return timediff(start_time,lap_time);
  }
};

#define PROGRESS_BAR "=================================================="

class ProgressBar{
  private:
    uint64_t total_work;
    uint64_t next_update;
    uint64_t call_diff;
    uint16_t old_percent;
    Timer    timer;

    void clearConsoleLine() const {
      std::cerr<<"\r\033[2K"<<std::flush;
    }

  public:
    void start(uint64_t total_work){
      timer = Timer();
      timer.start();
      this->total_work = total_work;
      next_update      = 0;
      call_diff        = total_work/200;
      old_percent      = 0;
      clearConsoleLine();
    }

    void update(uint64_t work_done){
      #ifdef NOPROGRESS
        return;
      #endif

      if(omp_get_thread_num()!=0)
        return;

      if(work_done<next_update)
        return;

      next_update += call_diff;

      uint16_t percent = (uint8_t)(work_done*omp_get_num_threads()*100/total_work);
      if(percent>100)
        percent=100;
      if(percent==old_percent)
        return;
      old_percent=percent;

      std::cerr<<"\r\033[2K["
               <<std::string(percent/2, '=')<<std::string(50-percent/2, ' ')
               <<"] ("
               <<percent<<"% - "
               <<std::fixed<<std::setprecision(1)<<timer.lap()/percent*(100-percent)
               <<"s - "
               <<omp_get_num_threads()<< " threads)"<<std::flush;
    }

    double stop(){
      clearConsoleLine();

      timer.stop();
      return timer.accumulated();
    }

    double time_it_took(){
      return timer.accumulated();
    }
};

#endif
