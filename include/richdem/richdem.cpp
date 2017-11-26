#include "richdem/common/logger.hpp"

#include <unordered_map>
#include <string>

namespace richdem {

std::unordered_map<LogFlag, std::string> log_flag_chars_begin = {
  {ALG_NAME, "\nA"},
  {CITATION, "C"},
  {DEBUG,    "\033[95md"},
  {ERROR,    "E"},
  {MEM_USE,  " "},
  {MISC,     "m"},
  {PROGRESS, "p"},
  {TIME_USE, "t"},
  {WARN,     "\033[91mW"}
};

std::unordered_map<LogFlag, std::string> log_flag_chars_end = {
  {ALG_NAME, ""},
  {CITATION, "\n"},
  {DEBUG,    ""},
  {ERROR,    ""},
  {MEM_USE,  ""},
  {MISC,     ""},
  {PROGRESS, ""},
  {TIME_USE, ""},
  {WARN,     ""}
};

void RDLOGfunc(LogFlag flag, const char* file, const char* func, unsigned line, std::string msg) {
  std::cerr<<log_flag_chars_begin.at(flag)<<" "<<msg
  #ifdef RICHDEM_DEBUG
  <<" ("<<file<<" : "<<func<<" : "<<line<<")"
  #endif
  <<"\033[39m"
  <<log_flag_chars_end.at(flag)
  <<std::endl;
}

}