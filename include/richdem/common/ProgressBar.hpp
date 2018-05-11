/**
  @file
  @brief Defines a handy progress bar object so users don't get impatient.

  The progress bar indicates to the user how much work has been completed, how
  much is left, and how long it is estimated to take. It accounts for
  multithreading by assuming uniform progress by all threads.

  Define the global macro `NOPROGRESS` disables the progress bar, which may 
  speed up the program.

  The progress bar looks like this:

      [===================================               ] (70% - 0.2s - 1 threads)

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_progress_bar_hpp_
#define _richdem_progress_bar_hpp_

#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "richdem/common/timer.hpp"

namespace richdem {

///Macros used to disguise the fact that we do not have multithreading enabled.
#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num()  0
  #define omp_get_num_threads() 1
#endif

///@brief Manages a console-based progress bar to keep the user entertained.
///
///Defining the global `NOPROGRESS` will
///disable all progress operations, potentially speeding up a program. The look
///of the progress bar is shown in ProgressBar.hpp.
class ProgressBar{
  private:
    uint32_t total_work;    ///< Total work to be accomplished
    uint32_t next_update;   ///< Next point to update the visible progress bar
    uint32_t call_diff;     ///< Interval between updates in work units
    uint32_t work_done;
    uint16_t old_percent;   ///< Old percentage value (aka: should we update the progress bar) TODO: Maybe that we do not need this
    Timer    timer;         ///< Used for generating ETA

    ///Clear current line on console so a new progress bar can be written
    void clearConsoleLine() const {
      std::cerr<<"\r\033[2K"<<std::flush;
    }

  public:
    ///@brief Start/reset the progress bar.
    ///@param total_work  The amount of work to be completed, usually specified in cells.
    void start(uint32_t total_work){
      timer = Timer();
      timer.start();
      this->total_work = total_work;
      next_update      = 0;
      call_diff        = total_work/200;
      old_percent      = 0;
      work_done        = 0;
      clearConsoleLine();
    }

    ///@brief Update the visible progress bar, but only if enough work has been done.
    ///
    ///Define the global `NOPROGRESS` flag to prevent this from having an
    ///effect. Doing so may speed up the program's execution.
    void update(uint32_t work_done0){
      #ifdef NOPROGRESS
        return;
      #endif

      if(omp_get_thread_num()!=0)
        return;

      work_done = work_done0;

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

    ///Increment by one the work done and update the progress bar
    ProgressBar& operator++(){
      work_done++;
      update(work_done);
      return *this;
    }

    ///Stop the progress bar. Throws an exception if it wasn't started.
    ///@return The number of seconds the progress bar was running.
    double stop(){
      clearConsoleLine();

      timer.stop();
      return timer.accumulated();
    }

    ///@return Return the time the progress bar ran for.
    double time_it_took(){
      return timer.accumulated();
    }

    uint32_t cellsProcessed() const {
      return work_done;
    }
};

}

#endif
