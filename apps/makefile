#frances_windows:
#	x86_64-w64-mingw32-g++ -o frances_win.exe -O3 `gdal-config --cflags` `gdal-config --dep-libs` frances.cpp -Wall -static-libgcc -static-libstdc++ --std=c++11 # -Wno-unknown-pragmas

frances:
	g++ -o frances.exe -O3 `gdal-config --cflags` `gdal-config --libs` frances.cpp -lgdal --std=c++11 -Wall # -Wno-unknown-pragmas

flood_for_flowdirs:
	g++ -o flood_for_flowdirs.exe -O3 `gdal-config --cflags` `gdal-config --libs` flood_for_flowdirs.cpp -lgdal --std=c++11 -Wall # -Wno-unknown-pragmas

flow_accumulation:
	g++ -o flow_accumulation.exe -g -O3 `gdal-config --cflags` `gdal-config --libs` flow_accumulation.cpp -lgdal --std=c++11 -Wall -DNOPROGRESS -DNDEBUG # -Wno-unknown-pragmas

fast_flood:
	g++ -o fast_flood.exe -g -O3 `gdal-config --cflags` `gdal-config --libs` fast_flood.cpp -lgdal --std=c++11 -Wall -DNOPROGRESS -DNDEBUG # -Wno-unknown-pragmas

d8_flowdirs:
	g++ -o d8_flowdirs.exe -g -O3 `gdal-config --cflags` `gdal-config --libs` d8_flowdirs.cpp -lgdal --std=c++11 -Wall # -Wno-unknown-pragmas

set_no_data:
	g++ -o set_no_data.exe -g -O3 `gdal-config --cflags` `gdal-config --libs` set_no_data.cpp -lgdal --std=c++11 -Wall # -Wno-unknown-pragmas

taudem_d8_to_richdem_d8:
	g++ -o taudem_d8_to_richdem_d8.exe -g -O3 `gdal-config --cflags` `gdal-config --libs` taudem_d8_to_richdem_d8.cpp -lgdal --std=c++11 -Wall # -Wno-unknown-pragmas

expand_dimensions:
	g++ -o expand_dimensions.exe -g -O3 `gdal-config --cflags` `gdal-config --libs` expand_dimensions.cpp -lgdal --std=c++11 -Wall # -Wno-unknown-pragmas

test:
	g++ -o test.exe -g -O3 `gdal-config --cflags` `gdal-config --libs` test.cpp -lgdal --std=c++11 -Wall # -Wno-unknown-pragmas

pitremove:
	g++ -o pitremove.exe -g -O3 `gdal-config --cflags` `gdal-config --libs` pitremove.cpp -lgdal --std=c++11 -Wall # -Wno-unknown-pragmas


clean:
	rm -rf *.exe
