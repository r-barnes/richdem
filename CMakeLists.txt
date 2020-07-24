cmake_minimum_required (VERSION 3.9)

project (librichdem
  VERSION 2.2.9
  DESCRIPTION "High-performance terrain analysis"
  LANGUAGES CXX
)

# add_subdirectory(apps)
# add_subdirectory(tests)

find_package(GDAL)
find_package(OpenMP)

option(USE_GDAL "Whether or not to compile with GDAL.")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC RICHDEM_GIT_HASH)

STRING(TIMESTAMP RICHDEM_TIMESTAMP UTC)

add_library(richdem
  src/richdem.cpp
  src/random.cpp
  src/gdal.cpp

  src/terrain_generation/terrain_generation.cpp
  src/terrain_generation/PerlinNoise.cpp
)

if(NOT OpenMP_CXX_FOUND)
  message(WARNING "Compiling without OpenMP!")
endif()

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

if(OpenMP_CXX_FOUND)
  target_link_libraries(richdem PUBLIC OpenMP::OpenMP_CXX)
endif()

target_compile_options(richdem
  PUBLIC
    -DRICHDEM_GIT_HASH="${RICHDEM_GIT_HASH}"
    -DRICHDEM_COMPILE_TIME="${RICHDEM_TIMESTAMP}"
)

if(USE_GDAL AND GDAL_FOUND)
  message(STATUS "Compiling RichDEM with GDAL.")
  target_link_libraries(richdem
    PUBLIC
      ${GDAL_LIBRARY}
  )
  target_include_directories(richdem
    PUBLIC
      ${GDAL_INCLUDE_DIR}
  )
  target_compile_options(richdem PUBLIC -DUSEGDAL)
else()
  message(WARNING "Compiling RichDEM without GDAL!")
endif()

add_executable(richdem_unittests
  tests/tests.cpp
)
target_link_libraries(richdem_unittests PRIVATE richdem)
target_compile_features(richdem_unittests PRIVATE cxx_std_17)