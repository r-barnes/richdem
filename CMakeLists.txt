cmake_minimum_required (VERSION 3.0)

project (librichdem 
  VERSION 2.2.9 
  DESCRIPTION "High-performance terrain analysis"
  LANGUAGES CXX
)

execute_process(
  COMMAND git log -1 --format=%H
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

STRING(TIMESTAMP TIMEZ UTC)

add_definitions(-DRICHDEM_GIT_HASH=${GIT_COMMIT_HASH})
add_definitions(-DRICHDEM_COMPILE_TIME=${GIT_COMMIT_HASH})

add_library(richdem
  include/richdem/richdem.cpp
)

# set_target_properties(richdem PROPERTIES VERSION ${PROJECT_VERSION})
# set_target_properties(richdem PROPERTIES SOVERSION 2)

target_include_directories(richdem 
  PUBLIC 
    $<INSTALL_INTERFACE:include>    
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(richdem 
  PUBLIC
    cxx_auto_type
    cxx_std_11
)