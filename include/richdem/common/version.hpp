#ifndef _richdem_version_hpp_
#define _richdem_version_hpp_

#include <string>
#include <iostream>

#ifndef RICHDEM_GIT_HASH
  #message Compiling without a git hash!
  const std::string git_hash = "NO HASH SPECIFIED!";
#else
  const std::string git_hash = std::string(RICHDEM_GIT_HASH).substr(0,16);
#endif

#ifndef RICHDEM_COMPILE_TIME
  #message Compiling without UTC compile time falling back to local!
  const std::string compilation_datetime = __DATE__ " " __TIME__;
#else
  const std::string compilation_datetime = RICHDEM_COMPILE_TIME;
#endif

///Richdem vX.X.X
const std::string program_name = "RichDEM v0.0.0";

///Richard Barnes
const std::string author_name  = "Richard Barnes";

///Richard Barnes © 2016
const std::string copyright    = "Richard Barnes © 2016";

///Richdem vX.X.X (hash=GIT HASH, compiled=COMPILATION DATE TIME)
const std::string program_identifier = program_name + " (hash=" + git_hash + ", compiled="+compilation_datetime + ")";

std::string PrintRichdemHeader(int argc, char **argv){
  std::string analysis;
  for(int i=0;i<argc;i++)
    analysis += std::string(argv[i])+" ";

  std::cout<<"c Program name       = " <<program_name        <<std::endl;
  std::cout<<"c Script compiled at = " <<compilation_datetime<<std::endl;
  std::cout<<"c Git hash           = " <<git_hash            <<std::endl;
  std::cout<<"c Copyright          = " <<copyright           <<std::endl;
  std::cout<<"a Analysis command   = " <<analysis            <<std::endl;

  return analysis;
}

#endif