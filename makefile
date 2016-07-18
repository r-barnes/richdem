export OMPI_CXX=g++

compile: main.cpp
	g++ -o parallel_flats.exe -g -O3 -DSHOW_STAMPS `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall #-Wextra #-fsanitize=undefined #-Wextra -Wconversion #-DNDEBUG

compile_with_compression:
	g++ -o parallel_flats.exe -g -O3 -DNDEBUG -DWITH_COMPRESSION `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall -lboost_iostreams -lz #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

xsede_with_compression:
	g++ -o parallel_flats.exe -g -O3 -DNDEBUG -DWITH_COMPRESSION `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall -lboost_iostreams -lz -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

xsede_debug_with_compression:
	g++ -o parallel_flats.exe -g -DNDEBUG -DWITH_COMPRESSION `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall -lboost_iostreams -lz -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include

debug: main.cpp
	g++ -o parallel_flats.exe -g `gdal-config --cflags` `gdal-config --libs` main.cpp -I. -lgdal --std=c++11 -Wall


clean:
	rm -f output* parallel_flats.exe
