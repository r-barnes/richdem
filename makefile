compile: main.cpp
	mpic++ -o parallel_d8_accum.exe -g -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

compile_with_compression:
	mpic++ -o parallel_d8_accum.exe -g -O3 -DNDEBUG -DWITH_COMPRESSION `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall -lboost_iostreams -lz #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

timing:
	mpic++ -o parallel_d8_accum.exe -g -O3 -DNDEBUG `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lipm -lgdal --std=c++11 -Wall #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

xsede: main.cpp
	mpic++ -o parallel_pit_fill.exe -g -O3 -DNDEBUG `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

xsede_with_compression:
	mpic++ -o parallel_pit_fill.exe -g -O3 -DNDEBUG -DWITH_COMPRESSION `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall -lboost_iostreams -lz -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

debug: main.cpp
	mpic++ -o parallel_d8_accum.exe -g     `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall
	#To run with debuggers: mpirun -n 2 xterm -hold -e gdb -ex run --args ./parallel_d8_accum.exe one @evict tests/dev/testdem10.d8 /z/out-%n.tif -w 5 -h 5

auth_gen:
	mpic++ -o auth_gen.exe -g -O3 -DNDEBUG `gdal-config --cflags` `gdal-config --libs` auth_gen.cpp -I. -lgdal --std=c++11 -Wall

assemble_ascii:
	g++ -o assemble_ascii.exe -g `gdal-config --cflags` `gdal-config --libs` assemble_ascii.cpp -I. -lgdal --std=c++11 -Wall

clean:
	rm -f output* parallel_d8flow_accum.exe generate_flowdirs.exe
