#Note: This makefile is deprecated, please use cmake if possible
export OMPI_CXX=g++
export CXX=g++
export MPICXX=mpic++
export GDAL_LIBS=`gdal-config --libs` -DUSEGDAL
export GDAL_CFLAGS=`gdal-config --cflags`
RICHDEM_GIT_HASH=`git rev-parse HEAD`
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export CXXFLAGS=$(GDAL_CFLAGS) --std=c++17 -I../../include -I. -Wall -Wno-unknown-pragmas -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""
export OPT_FLAGS=-g -O3 -DNDEBUG
export DEBUG_FLAGS=-g
export COMPRESSION_LIBS=-lboost_iostreams -lz
export XSEDE_MPI_LIBS=-L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/

.PHONY: clean

compile: main.cpp
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_pf.exe main.cpp $(GDAL_LIBS)

compile_with_compression:
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_pf.exe -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) $(COMPRESSION_LIBS)

timing:
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_pf.exe main.cpp -lipm $(GDAL_LIBS)

xsede: main.cpp
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_pit_fill.exe main.cpp $(GDAL_LIBS) $(XSEDE_MPI_LIBS)

xsede_with_compression:
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_pit_fill.exe -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) $(XSEDE_MPI_LIBS) $(COMPRESSION_LIBS)

debug: main.cpp
	$(MPICXX) $(CXXFLAGS) $(DEBUG_FLAGS) -o parallel_pf.exe main.cpp $(GDAL_LIBS)
	#To run with debuggers: mpirun -n 2 xterm -hold  -e gdb -ex run --args ./parallel_pf.exe 0 beauford.tif

auth_gen:
	$(MPICXX) $(CXXFLAGS) $(OPT_FLAGS) -o auth_gen.exe auth_gen.cpp $(GDAL_LIBS)

clean:
	rm -f output* parallel_pf.exe
