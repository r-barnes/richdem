compile: main.cpp
	mpic++ -o parallel_d8flow_accum.exe -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization
	g++ -O3 -o generate_flowdirs.exe `gdal-config --cflags` `gdal-config --libs` generate_flowdirs.cpp -lgdal --std=c++11

debug: main.cpp
	mpic++ -o parallel_d8flow_accum.exe -g `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization
	#To run with debuggers: mpirun -n 3 xterm -hold  -e gdb --args ./parallel_d8flow_accum.exe tests/testdem1.d8

clean:
	rm -f output* parallel_d8flow_accum.exe generate_flowdirs.exe
