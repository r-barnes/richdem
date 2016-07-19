CC = g++
CFLAGS = -O3

default:
	$(CC) $(CFLAGS) -o flat_resolution.exe -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall

clean:
	rm -rf *.exe
