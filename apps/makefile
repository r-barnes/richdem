CC = g++
CFLAGS = -O3

frances:
	$(CC) $(CFLAGS) -o frances.exe -O3 `gdal-config --cflags` `gdal-config --libs` frances.cpp -lgdal --std=c++11 -Wall

clean:
	rm -rf *.exe
