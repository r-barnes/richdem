export CXX=g++
export GDAL_LIBS=`gdal-config --libs`
export GDAL_CFLAGS=`gdal-config --cflags` -DUSEGDAL
RICHDEM_GIT_HASH=`git rev-parse HEAD`
RICHDEM_COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`
export LIBS=$(GDAL_LIBS) -lstdc++fs
export CXXFLAGS=$(GDAL_CFLAGS) --std=c++17 -O3 -Wall -Wno-unknown-pragmas -I../include -DRICHDEM_GIT_HASH="\"$(RICHDEM_GIT_HASH)\"" -DRICHDEM_COMPILE_TIME="\"$(RICHDEM_COMPILE_TIME)\""

#-DNOPROGRESS -DNDEBUG

tests:
	$(CXX) $(CXXFLAGS) tests.cpp -o tests.exe $(LIBS)