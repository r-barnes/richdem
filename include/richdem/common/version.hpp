/**
  @file
  @brief Defines RichDEM version, git hash, compilation time. Used for 
         program/app headers and for processing history entries.

  Richard Barnes (rbarnes@umn.edu), 2015
*/
#ifndef _richdem_version_hpp_
#define _richdem_version_hpp_

#include <string>
#include <iostream>

namespace richdem {

#ifndef RICHDEM_GIT_HASH
  #pragma message "Compiling without a git hash!"
  ///Git hash of program's source (used if RICHDEM_GIT_HASH is undefined)
  const std::string git_hash = "NO HASH SPECIFIED!";
#else
  ///Git hash of program's source (Used if RICHDEM_GIT_HASH is defined)
  const std::string git_hash = std::string(RICHDEM_GIT_HASH).substr(0,16);
#endif

#ifndef RICHDEM_COMPILE_TIME
  #pragma message "Compiling without UTC compile time falling back to local!"
  ///Date and time of when the program was compiled (used if RICHDEM_COMPILE_TIME is undefined)
  const std::string compilation_datetime = __DATE__ " " __TIME__;
#else
  ///Date and imte of when the program was compiled
  const std::string compilation_datetime = RICHDEM_COMPILE_TIME;
#endif

///Richdem vX.X.X
const std::string program_name = "RichDEM v2.2.0";

///Richard Barnes
const std::string author_name  = "Richard Barnes";

///Richard Barnes © 2018
const std::string copyright    = "Richard Barnes © 2018";

///Richdem vX.X.X (hash=GIT HASH, compiled=COMPILATION DATE TIME)
const std::string program_identifier = program_name + " (hash=" + git_hash + ", compiled="+compilation_datetime + ")";

std::string rdHash(){
  return git_hash;
}

std::string rdCompileTime() {
  return compilation_datetime;
}

///Takes the program's command line arguments and prints to stdout a header with
///a variety of useful information for identifying the particulars of what was
///run.
std::string PrintRichdemHeader(int argc, char **argv){
  std::string analysis;
  for(int i=0;i<argc;i++)
    analysis += std::string(argv[i])+" ";

  std::cout<<"c Program name       = " <<program_name        <<std::endl;
  std::cout<<"c Script compiled at = " <<rdCompileTime()     <<std::endl;
  std::cout<<"c Git hash           = " <<rdHash()            <<std::endl;
  std::cout<<"c Copyright          = " <<copyright           <<std::endl;
  std::cout<<"a Analysis command   = " <<analysis            <<std::endl;

  return analysis;
}

}

#endif
