export CXX=g++
export MPICXX=mpic++
export GDAL_LIBS=`gdal-config --libs`
export GDAL_CFLAGS=`gdal-config --cflags`
export CXXFLAGS=$(GDAL_CFLAGS) --std=c++11 -I../../include -I. -Wall -Wno-unknown-pragmas
export OPT_FLAGS=-O3
export DEBUG_FLAGS=-g
export COMPRESSION_LIBS=-lboost_iostreams -lz
export XSEDE_BOOST_INCLUDES=-I/opt/boost/intel/mvapich2_ib/include
export XSEDE_BOOST_LIBS=-L/opt/boost/intel/mvapich2_ib/lib

#-DNDEBUG -DSHOW_STAMPS

mpic++ -g --std=c++11  `gdal-config --cflags` -I/opt/boost/intel/mvapich2_ib/include/ -I../../include -I. -Wall -Wno-unknown-pragmas -DWITH_COMPRESSION main.cpp -L/opt/boost/intel/mvapich2_ib/lib/  -lboost_iostreams -lz `gdal-config --libs`  -o parallel_d8_accum.exe 


compile: main.cpp
	$(MPICXX) $(OPT_FLAGS) $(CXXFLAGS) -o parallel_d8_accum.exe main.cpp $(GDAL_LIBS)

compile_with_compression:
	$(MPICXX) $(OPT_FLAGS) $(CXXFLAGS) -o parallel_d8_accum.exe -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) $(COMPRESSION_LIBS)

timing:
	$(MPICXX) $(OPT_FLAGS) $(CXXFLAGS) -o parallel_d8_accum.exe main.cpp -lipm $(GDAL_LIBS)

xsede: main.cpp
	$(MPICXX) $(OPT_FLAGS) $(CXXFLAGS) -o parallel_d8_accum.exe main.cpp $(GDAL_LIBS) $(XSEDE_MPI_LIBS)

xsede_with_compression:
	$(MPICXX) $(OPT_FLAGS) $(CXXFLAGS) $(XSEDE_BOOST_INCLUDES) -o parallel_d8_accum.exe -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) $(XSEDE_BOOST_LIBS) $(COMPRESSION_LIBS)

debug: main.cpp
	$(MPICXX) $(DEBUG_FLAGS) $(CXXFLAGS) -o parallel_d8_accum.exe main.cpp $(GDAL_LIBS)
	#To run with debuggers: mpirun -n 2 xterm -hold -e gdb -ex run --args ./parallel_d8_accum.exe one @evict tests/dev/testdem10.d8 /z/out-%n.tif -w 5 -h 5

assemble_ascii:
	$(CXX) $(OPT_FLAGS) $(CXXFLAGS) -o assemble_ascii.exe assemble_ascii.cpp $(GDAL_LIBS)

merge_rasters_by_layout:
	$(CXX) $(OPT_FLAGS) $(CXXFLAGS) -o merge_rasters_by_layout.exe merge_rasters_by_layout.cpp $(GDAL_LIBS)

test:
	$(CXX) $(OPT_FLAGS) $(CXXFLAGS) -o test.exe test.cpp $(GDAL_LIBS)

clean:
	rm -f output* parallel_d8flow_accum.exe
