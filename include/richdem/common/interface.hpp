/**
  @file
  @brief Defines a handy progress bar so the user doesn't get bored or panicked.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_utility_hpp_
#define _richdem_utility_hpp_

#include <string>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <stdexcept>
#include "richdem/common/timer.hpp"

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num()  0
  #define omp_get_num_threads() 1
#endif

class ProgressBar{
  private:
    long total_work;
    int call_count;
    int old_call_count;
    int old_percent;
    int call_diff;
    Timer timer;
  public:
    void start(long total_work0){
      timer = Timer();
      timer.start();
      total_work     = total_work0;
      old_percent    = 0;
      call_count     = 0;
      old_call_count = 0;
      call_diff      = total_work/200;
      //This ANSI code clears the entire line
      std::cerr<<"\r\033[2K"<<std::flush;
    }
    void update(long work_done){
      #ifdef NOPROGRESS
        return;
      #endif
      if(omp_get_thread_num()!=0)
        return;
      ++call_count;
      if(call_count-old_call_count!=call_diff)
        return;

      old_call_count = call_count;

      int percent;
      percent=(int)(work_done*omp_get_num_threads()*100/total_work);
      if(percent>100)
        percent=100;
      if(percent==old_percent)
        return;
      old_percent=percent;

      if(total_work!=-1)
        std::cerr<<"\r\033[2K["
                 <<std::string(percent/2, '=')<<std::string(50-percent/2, ' ')
                 <<"] ("
                 <<percent<<"% - "
                 <<std::fixed<<std::setprecision(1)<<timer.lap()/percent*(100-percent)
                 <<"s - "
                 <<omp_get_num_threads()<< " threads)"<<std::flush;
      else
        std::cerr<<"\r\033[2K"<<std::setw(12)<<work_done<<" cells completed."<<std::endl;
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
