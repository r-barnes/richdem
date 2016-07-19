export CXX=g++
export GDAL_LIBS=`gdal-config --libs`
export GDAL_CFLAGS=`gdal-config --cflags`
export CXXFLAGS=$(GDAL_CFLAGS) --std=c++11 -O3 -Wall -Wno-unknown-pragmas -I../include

#-DNOPROGRESS -DNDEBUG

include makefile.vars