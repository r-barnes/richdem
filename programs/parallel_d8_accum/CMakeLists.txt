cmake_minimum_required (VERSION 3.9)

project (richdem_pfacc
  VERSION 2.2.9
  DESCRIPTION "RichDEM Parallel D8 Flow Accumulation"
  LANGUAGES CXX
)

find_package(MPI REQUIRED)
find_package(Boost COMPONENTS iostreams)

add_executable(parallel_d8_accum.exe main.cpp)
target_include_directories(parallel_d8_accum.exe PRIVATE .)
target_link_libraries(parallel_d8_accum.exe PRIVATE MPI::MPI_CXX richdem)
target_compile_features(parallel_d8_accum.exe
  PUBLIC
    cxx_auto_type
    cxx_std_11
)

if(Boost_FOUND)
  add_executable(parallel_d8_accum_with_compression.exe main.cpp)
  target_include_directories(parallel_d8_accum_with_compression.exe PRIVATE .)
  target_link_libraries(parallel_d8_accum_with_compression.exe PRIVATE MPI::MPI_CXX richdem)
  target_compile_definitions(parallel_d8_accum_with_compression.exe PRIVATE DWITH_COMPRESSION)
  target_compile_features(parallel_d8_accum_with_compression.exe
    PUBLIC
      cxx_auto_type
      cxx_std_11
  )
endif()

add_executable(pd8acc_test.exe test.cpp)
target_link_libraries(pd8acc_test.exe PRIVATE richdem)

#TODO
#timing:
#	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_d8_accum.exe main.cpp -lipm $(GDAL_LIBS)