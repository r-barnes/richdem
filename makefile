export CXX=g++
export GDAL_LIBS=`gdal-config --libs`
export GDAL_CFLAGS=`gdal-config --cflags`
export CXXFLAGS=$(GDAL_CFLAGS) --std=c++11 -O3 -Wall -Wno-unknown-pragmas -Irichdem/include -I.

#-Wextra #-fsanitize=undefined #-Wextra -Wconversion #-DNDEBUG

compile: main.cpp
	$(CXX) $(CXXFLAGS) -o parallel_flats.exe -g -O3 -DSHOW_STAMPS main.cpp $(GDAL_LIBS)  

compile_with_compression:
	$(CXX) $(CXXFLAGS) -o parallel_flats.exe -g -O3 -DNDEBUG -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) -lboost_iostreams -lz

xsede_with_compression:
	$(CXX) $(CXXFLAGS) -o parallel_flats.exe -g -O3 -DNDEBUG -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) -lboost_iostreams -lz -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include

xsede_debug_with_compression:
	$(CXX) $(CXXFLAGS) -o parallel_flats.exe -g -DNDEBUG -DWITH_COMPRESSION main.cpp $(GDAL_LIBS) -lboost_iostreams -lz -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include

debug: main.cpp
	$(CXX) $(CXXFLAGS) -o parallel_flats.exe -g main.cpp $(GDAL_LIBS) 


clean:
	rm -f output* parallel_flats.exe
