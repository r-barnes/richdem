export CXX=g++
export MPICXX=mpic++
export GDAL_LIBS=`gdal-config --libs`
export GDAL_CFLAGS=`gdal-config --cflags`
export CXXFLAGS=$(GDAL_CFLAGS) --std=c++11 -Irichdem/include -I. -Wall -Wno-unknown-pragmas
export OPT_FLAGS=-O3
export DEBUG_FLAGS=-g
export COMPRESSION_LIBS=-lboost_iostreams -lz
export XSEDE_MPI_LIBS=-L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/

#-DNDEBUG -DSHOW_STAMPS

compile: main.cpp
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_d8_accum.exe main.cpp $(GDAL_LIBS)

compile_with_compression:
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_d8_accum.exe -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) $(COMPRESSION_LIBS)

timing:
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_d8_accum.exe main.cpp -lipm $(GDAL_LIBS)

xsede: main.cpp
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_pit_fill.exe main.cpp $(GDAL_LIBS) $(XSEDE_MPI_LIBS)

xsede_with_compression:
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_pit_fill.exe -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) $(XSEDE_MPI_LIBS) $(COMPRESSION_LIBS)

debug: main.cpp
	$(MPICXX) $(CXXFLAGS) $(DEBUG_FLAGS) -o parallel_d8_accum.exe main.cpp $(GDAL_LIBS)
	#To run with debuggers: mpirun -n 2 xterm -hold -e gdb -ex run --args ./parallel_d8_accum.exe one @evict tests/dev/testdem10.d8 /z/out-%n.tif -w 5 -h 5

assemble_ascii:
	$(CXX) $(CXXFLAGS) $(OPT_FLAGS) -o assemble_ascii.exe assemble_ascii.cpp $(GDAL_LIBS)

test:
	$(CXX) $(CXXFLAGS) $(OPT_FLAGS) -o test.exe test.cpp $(GDAL_LIBS)

clean:
	rm -f output* parallel_d8flow_accum.exe
