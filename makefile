export OMPI_CXX=g++

compile: main.cpp
	mpic++ -o parallel_pf.exe -g -O3 -DNDEBUG `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

compile_with_compression:
	mpic++ -o parallel_pf.exe -g -O3 -DNDEBUG -DWITH_COMPRESSION `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall -lboost_iostreams -lz #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

timing:
	mpic++ -o parallel_pf.exe -g -O3 -DNDEBUG `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lipm -lgdal --std=c++11 -Wall #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

xsede: main.cpp
	mpic++ -o parallel_pit_fill.exe -g -O3 -DNDEBUG `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

xsede_with_compression:
	mpic++ -o parallel_pit_fill.exe -g -O3 -DNDEBUG -DWITH_COMPRESSION `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall -lboost_iostreams -lz -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

debug: main.cpp
	mpic++ -o parallel_pf.exe -g     `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall
	#To run with debuggers: mpirun -n 2 xterm -hold  -e gdb -ex run --args ./parallel_pf.exe 0 beauford.tif

test: test.cpp
	g++ -o test.exe -g -O3 -DNDEBUG `gdal-config --cflags` `gdal-config --libs` test.cpp -I. -lgdal --std=c++11 -Wall #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

auth_gen:
	mpic++ -o auth_gen.exe -g -O3 -DNDEBUG `gdal-config --cflags` `gdal-config --libs` auth_gen.cpp -I. -lgdal --std=c++11 -Wall

clean:
	rm -f output* parallel_pf.exe
