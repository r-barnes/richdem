all:
	$(CXX) --std=c++17 -O3 -g `gdal-config --cflags` -I../../include -DUSEGDAL -DRICHDEM_LOGGING main.cpp perlin.c ../../include/richdem/richdem.cpp -lgdal

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
