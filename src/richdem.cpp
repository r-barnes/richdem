#include <richdem/common/logger.hpp>

#include <map>
#include <string_view>

namespace richdem {

// TODO(r-barnes): Could use vector for this since the enum just compiles to an
// int. But need to make sure things are in the right order.
std::string_view log_flag_chars_begin(LogFlag flag){
  switch(flag){
    case LogFlag::ALG_NAME:  return "\nA";
    case LogFlag::CITATION:  return "C";
    case LogFlag::CONFIG:    return "c";
    case LogFlag::DEBUG:     return "\033[95md";
    case LogFlag::ERROR_:    return "E";
    case LogFlag::MEM_USE:   return " ";
    case LogFlag::MISC:      return "m";
    case LogFlag::PROGRESS:  return "p";
    case LogFlag::TIME_USE:  return "t";
    case LogFlag::WARN:      return "\033[91mW";
    default:
      throw std::runtime_error("Unrecognized logging flag!");
  }
}

std::string_view log_flag_chars_end(LogFlag flag){
  switch(flag){
    case LogFlag::ALG_NAME:  return "";
    case LogFlag::CITATION:  return "\n";
    case LogFlag::CONFIG:    return "";
    case LogFlag::DEBUG:     return "";
    case LogFlag::ERROR_:    return "";
    case LogFlag::MEM_USE:   return "";
    case LogFlag::MISC:      return "";
    case LogFlag::PROGRESS:  return "";
    case LogFlag::TIME_USE:  return "";
    case LogFlag::WARN:      return "";
    default:
      throw std::runtime_error("Unrecognized logging flag!");
  }
}

void RDLOGfunc(const LogFlag flag, const char* file, const char* func, unsigned line, const std::string &msg) {
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
