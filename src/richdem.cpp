#include <richdem/common/logger.hpp>

#include <map>
#include <string_view>

namespace richdem {

// TODO(r-barnes): Could use vector for this since the enum just compiles to an
// int. But need to make sure things are in the right order.
std::string_view log_flag_chars_begin(LogFlag flag){
  switch(flag){
    case ALG_NAME:  return "\nA";
    case CITATION:  return "C";
    case CONFIG:    return "c";
    case DEBUG:     return "\033[95md";
    case ERROR_:    return "E";
    case MEM_USE:   return " ";
    case MISC:      return "m";
    case PROGRESS:  return "p";
    case TIME_USE:  return "t";
    case WARN:      return "\033[91mW";
    default:
      throw std::runtime_error("Unrecognized logging flag!");
  }
}

std::string_view log_flag_chars_end(LogFlag flag){
  switch(flag){
    case ALG_NAME:  return "";
    case CITATION:  return "\n";
    case CONFIG:    return "";
    case DEBUG:     return "";
    case ERROR_:    return "";
    case MEM_USE:   return "";
    case MISC:      return "";
    case PROGRESS:  return "";
    case TIME_USE:  return "";
    case WARN:      return "";
    default:
      throw std::runtime_error("Unrecognized logging flag!");
  }
}

void RDLOGfunc(LogFlag flag, const char* file, const char* func, unsigned line, const std::string &msg) {
  (void)file; // Suppress unused variable warning
  (void)func; // Suppress unused variable warning
  (void)line; // Suppress unused variable warning
  std::cerr<<log_flag_chars_begin(flag)<<" "<<msg
  #ifdef RICHDEM_DEBUG
  <<" ("<<file<<" : "<<func<<" : "<<line<<")"
  #endif
  <<"\033[39m"
  <<log_flag_chars_end(flag)
  <<std::endl;
}

}
