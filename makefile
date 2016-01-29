export OMPI_CXX=g++

compile: main.cpp
	mpic++ -o parallel_pdf.exe -g -O3 -DNDEBUG `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

xsede: main.cpp
	mpic++ -o parallel_pdf.exe -g -O3 -DNDEBUG `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system -L/opt/boost/intel/mvapich2_ib/lib/ -I/opt/boost/intel/mvapich2_ib/include #-Wextra #-fsanitize=undefined #-Wextra -Wconversion

debug: main.cpp
	mpic++ -o parallel_pdf.exe -g     `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
	#To run with debuggers: mpirun -n 2 xterm -hold  -e gdb -ex run --args ./parallel_pdf.exe 0 beauford.tif

clean:
	rm -f output* parallel_pdf.exe
