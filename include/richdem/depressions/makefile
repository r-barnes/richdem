CC = g++
CFLAGS = -O3 --std=c++11

default:
	$(CC) $(CFLAGS) -o priority_flood.exe -O3 `gdal-config --cflags` `gdal-config --libs` main.cpp -lgdal --std=c++11 -Wall

clean:
	rm -f *exe
