compile: main.cpp
	mpic++ -o parallel_d8flow_accum.exe -O3 -g  `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system -Wall
	g++ -O3 -o generate_flowdirs.exe `gdal-config --cflags` `gdal-config --libs` generate_flowdirs.cpp -lgdal --std=c++11

debug: main.cpp
	mpic++ -o parallel_d8flow_accum.exe -g `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
	#To run with debuggers: mpirun -n 3 xterm -hold  -e gdb -ex run --args ./parallel_d8flow_accum.exe one @offloadall ~/data/beauford.tif /z/out

clean:
	rm -f output* parallel_d8flow_accum.exe generate_flowdirs.exe
