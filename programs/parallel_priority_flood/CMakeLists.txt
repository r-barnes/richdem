cmake_minimum_required (VERSION 3.9)

project (richdem_ppf
  VERSION 2.2.9
  DESCRIPTION "RichDEM Parallel Priority-Flood"
  LANGUAGES CXX
)

find_package(MPI REQUIRED)
find_package(Boost COMPONENTS iostreams)

add_executable(parallel_pf.exe main.cpp)
target_include_directories(parallel_pf.exe PRIVATE .)
target_link_libraries(parallel_pf.exe PRIVATE MPI::MPI_CXX richdem)
target_compile_features(parallel_pf.exe
  PUBLIC
    cxx_auto_type
    cxx_std_11
)

if(Boost_FOUND)
  add_executable(parallel_pf_with_compression.exe main.cpp)
  target_include_directories(parallel_pf_with_compression.exe PRIVATE .)
  target_link_libraries(parallel_pf_with_compression.exe PRIVATE MPI::MPI_CXX richdem)
  target_compile_definitions(parallel_pf_with_compression.exe PRIVATE DWITH_COMPRESSION)
  target_compile_features(parallel_pf_with_compression.exe
    PUBLIC
      cxx_auto_type
      cxx_std_11
  )
endif()

add_executable(ppf_auth_gen.exe auth_gen.cpp)
target_link_libraries(ppf_auth_gen.exe PRIVATE richdem)

#TODO
#timing:
#	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_pf.exe main.cpp -lipm $(GDAL_LIBS)