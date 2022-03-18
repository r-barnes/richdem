export CXX=g++
export GDAL_LIBS=`gdal-config --libs`
export GDAL_CFLAGS=`gdal-config --cflags` -DUSEGDAL
export OPT_FLAGS=-g -O3 -DNDEBUG
export DEBUG_FLAGS=-g
RICHDEM_GIT_HASH=`git rev-parse HEAD`
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export CXXFLAGS=$(GDAL_CFLAGS) --std=c++17 -Wall -Wno-unknown-pragmas -I../../include -I. -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""

#-Wextra #-fsanitize=undefined #-Wextra -Wconversion

compile: main.cpp
	$(CXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_flats.exe main.cpp $(GDAL_LIBS)

compile_with_compression:
	$(CXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_flats.exe -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) -lboost_iostreams -lz

xsede_with_compression:
	$(CXX) $(CXXFLAGS) $(OPT_FLAGS) -o parallel_flats.exe -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) -lboost_iostreams -lz -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include

xsede_debug_with_compression:
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAGS) -o parallel_flats.exe -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) -lboost_iostreams -lz -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include

debug: main.cpp
	$(CXX) $(CXXFLAGS) $(DEBUG_FLAGS) -o parallel_flats.exe main.cpp $(GDAL_LIBS)

same_flowdirs: same_flowdirs.cpp
	$(CXX) $(CXXFLAGS) $(OPT_FLAGS) -o same_flowdirs.exe same_flowdirs.cpp $(GDAL_LIBS)

merge_rasters_by_layout: merge_rasters_by_layout.cpp
	$(CXX) $(CXXFLAGS) $(OPT_FLAGS) -o merge_rasters_by_layout.exe -g merge_rasters_by_layout.cpp $(GDAL_LIBS)

clean:
	rm -f output* parallel_flats.exe
