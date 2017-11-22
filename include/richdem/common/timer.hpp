/**
  @file
  @brief Defines the Timer class, which is used for timing code.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_timer_
#define _richdem_timer_

#include <sys/time.h>
#include <stdexcept>

namespace richdem {

///@brief Used to time how intervals in code.
///
///Such as how long it takes a given function to run, or how long I/O has taken.
class Timer{
  private:
    timeval start_time;      ///< Last time the timer was started
    double accumulated_time; ///< Accumulated running time since creation
    bool running;            ///< True when the timer is running

    ///Number of (fractional) seconds between two time objects
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

    ///Start the timer. Throws an exception if timer was already running.
    void start(){
      if(running)
        throw std::runtime_error("Timer was already started!");
      running=true;
      gettimeofday(&start_time, NULL);
    }

    ///Stop the timer. Throws an exception if timer was already stopped.
    ///Calling this adds to the timer's accumulated time.
    ///@return The accumulated time in seconds.
    double stop(){
      if(!running)
        throw std::runtime_error("Timer was already stopped!");
      running=false;
      timeval end_time;
      gettimeofday(&end_time, NULL);

      accumulated_time+=timediff(start_time,end_time);

      return accumulated_time;
    }

    ///Returns the timer's accumulated time. Throws an exception if the timer is
    ///running.
    double accumulated(){
      if(running)
        throw std::runtime_error("Timer is still running!");
      return accumulated_time;
    }

    ///Returns the time between when the timer was started and the current
    ///moment. Throws an exception if the timer is not running.
    double lap(){
      if(!running)
        throw std::runtime_error("Timer was not started!");
      timeval lap_time;
      gettimeofday(&lap_time, NULL);
      return timediff(start_time,lap_time);
    }

    ///Stops the timer and resets its accumulated time. No exceptions are thrown
    ///ever.
    void reset(){
      accumulated_time = 0;
      running          = false;
    }
};

}

#endif
