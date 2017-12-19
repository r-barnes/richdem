/**
  @file
  @brief Defines functions for calculating memory usage.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _memory_hpp_
#define _memory_hpp_

#include <fstream>
#include <string>

namespace richdem {

/**
  @brief Return memory statistics of the process

  This code is drawn from "http://stackoverflow.com/a/671389/752843"

  @param[out]   vmpeak    Peak virtual memory size (kB)
  @param[out]   vmhwm     Peak resident set size (kB)
*/
void ProcessMemUsage(long &vmpeak, long &vmhwm){
  #if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    vmpeak = 0;
    vmhwm  = 0;

    std::ifstream fin("/proc/self/status");
    if(!fin.good())
      return;
    
    while (!vmpeak || !vmhwm){
      std::string line;

      if(!getline(fin,line)){ //Check if we could still read file
        return;
      }

      if(line.compare(0,7,"VmPeak:")==0){        //Peak virtual memory size
        vmpeak = std::stoi(line.substr(7));
      // } else if(line.compare(0,7,"VmSize:")==0){ //Virtual memory size
      //   std::cerr<<"T: "<<line.substr(7,10)<<std::endl;
      //   vmsize = std::stoi(line.substr(7,10));
      // } else if(line.compare(0,6,"VmRSS:")==0){  //Resident set size
      //   std::cerr<<"T: "<<line.substr(7,10)<<std::endl;
      //   vmrss = std::stoi(line.substr(7,10));
      } else if(line.compare(0,6,"VmHWM:")==0){  //Peak resident set size
        vmhwm = std::stoi(line.substr(6));
      }
    }
  #else
    #pragma message("Cannot check memory statistics for this OS.")
  #endif
}

}

#endif
