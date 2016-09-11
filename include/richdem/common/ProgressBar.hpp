/**
  @file
  @brief Defines a handy progress bar object so users don't get impatient.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_progress_bar_hpp_
#define _richdem_progress_bar_hpp_

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

///The progress bar object accepts
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
