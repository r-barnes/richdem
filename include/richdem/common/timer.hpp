/**
  @file
  @brief Defines the Timer class, which is used for timing code.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_timer_
#define _richdem_timer_

#include <chrono>
#include <stdexcept>

namespace richdem {

///@brief Used to time how intervals in code.
///
///Such as how long it takes a given function to run, or how long I/O has taken.
class Timer{
  private:
    typedef std::chrono::high_resolution_clock clock;
    typedef std::chrono::duration<double, std::ratio<1> > second;

    std::chrono::time_point<clock> start_time;  ///< Last time the timer was started
    double accumulated_time = 0;                ///< Accumulated running time since creation
    bool   running          = false;            ///< True when the timer is running

    ///Number of (fractional) seconds between two time objects
    double timediff(const std::chrono::time_point<clock> &start, const std::chrono::time_point<clock> &end){
      return static_cast<unsigned long int>(std::chrono::duration_cast<std::chrono::seconds>(end - start).count());
    }
  public:
    ///Creates a Timer which is not running and has no accumulated time
    Timer() = default;

    ///Start the timers. Throws an exception if timer was already running.
    void start(){
      if(running)
        throw std::runtime_error("Timer was already started!");
      running=true;
      start_time = clock::now();
    }

    ///Stop the timer. Throws an exception if timer was already stopped.
    ///Calling this adds to the timer's accumulated time.
    ///
    ///@return The accumulated time in seconds.
    double stop(){
      if(!running)
        throw std::runtime_error("Timer was already stopped!");
      running=false;

      const auto end_time = clock::now();

      accumulated_time += timediff(start_time,end_time);

      return accumulated_time;
    }

    ///Returns the timer's accumulated time. Throws an exception if the timer is
    ///running.
    ///
    ///@return The timer's accumulated time, in seconds.
    double accumulated(){
      if(running)
        throw std::runtime_error("Timer is still running!");
      return accumulated_time;
    }

    ///Returns the time between when the timer was started and the current
    ///moment. Throws an exception if the timer is not running.
    ///
    ///@return Time since the timer was started and current moment, in seconds.
    double lap(){
      if(!running)
        throw std::runtime_error("Timer was not started!");
      const auto lap_time = clock::now();
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
