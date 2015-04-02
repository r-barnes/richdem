compile: main.cpp
	mpic++ -o parallel_pit_fill.exe -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization

debug: main.cpp
	mpic++ -o parallel_pit_fill.exe -g `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization
	#To run with debuggers: mpirun -n 2 xterm -hold  -e gdb -ex run --args ./parallel_pit_fill.exe 0 beauford.tif

clean:
	rm -f output* parallel_pit_fill.exe
