CC = g++
CFLAGS = -O3

test:
	$(CC) $(CFLAGS) -o test.exe -O3 `gdal-config --cflags` `gdal-config --libs` test.cpp -lgdal --std=c++11 -Wall

clean:
	rm -rf *.exe
