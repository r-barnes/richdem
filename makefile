compile: main.cpp
	mpic++ -o parallel_pit_fill.exe -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization

debug: main.cpp
	mpic++ -o parallel_pit_fill.exe -g `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization
	#To run with debuggers: mpirun -n 3 xterm -hold  -e gdb --args ./parallel_pit_fill.exe tests/testdem1.d8

clean:
	rm -f output* parallel_pit_fill.exe
