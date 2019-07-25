cmake_minimum_required (VERSION 3.0)

project (librichdem 
  VERSION 2.2.9 
  DESCRIPTION "High-performance terrain analysis"
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC RICHDEM_GIT_HASH)

STRING(TIMESTAMP RICHDEM_TIMESTAMP UTC)

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

target_compile_options(richdem
  PUBLIC
    -DRICHDEM_GIT_HASH="${RICHDEM_GIT_HASH}"
    -DRICHDEM_COMPILE_TIME="${RICHDEM_TIMESTAMP}"
)