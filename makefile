#CC=i586-mingw32msvc-g++
#CC=x86_64-w64-mingw32-g++
#CFLAGS=-Wall -Wno-unknown-pragmas -DARCGIS -static-libgcc -static-libstdc++ #-fopenmp

CC=g++
CFLAGS=-Wall -fopenmp #-DARCGIS #-lX11 -pthread #-ltbb
#-ffp-contract=fast enables forming of fused multiply-add instructions. Default is fast.

ODIR=obj
#PRE_FLAGS=-lgcov -g -fprofile-arcs -ftest-coverage
PRE_FLAGS=-O3
#PRE_FLAGS=-g

DEPS = d8_methods.h data_structures.h dinf_methods.h interface.h data_io.h pit_fill.h utility.h flat_resolution.h debug.h visualize.h unit_test.h interval_dinf.h

_OBJ = d8_methods.o dinf_methods.o data_io.o pit_fill.o utility.o debug.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp
	$(CC) $(PRE_FLAGS) -c -o $@ $< $(CFLAGS)

richdem: $(OBJ) obj/main.o
	$(CC) $(PRE_FLAGS) -o richdem.exe $^ $(CFLAGS)
	du -hs ./richdem.exe

unit: $(OBJ) obj/unit_test.o
	$(CC) $(PRE_FLAGS) -o richdem_unit.exe $^ $(CFLAGS)
	du -hs ./richdem_unit.exe

test: $(OBJ) obj/test.o
	$(CC) $(PRE_FLAGS) -o test.exe $^ $(CFLAGS)
	du -hs ./test.exe

interval: $(OBJ) obj/interval_dinf_main.o
	$(CC) $(PRE_FLAGS) -o richdem.exe $^ $(CFLAGS)
	du -hs ./richdem.exe

jake: $(OBJ) obj/jake.o
	$(CC) $(PRE_FLAGS) -o jake.exe $^ $(CFLAGS)
	du -hs ./jake.exe

jaked: $(OBJ) obj/jake_deliverables.o
	$(CC) $(PRE_FLAGS) -o jake_deliverables.exe $^ $(CFLAGS)
	du -hs ./jake_deliverables.exe

debug: $(OBJ) obj/main.o
	$(CC) $(PRE_FLAGS) -g -o richdem.exe $^ $(CFLAGS)
	du -hs ./richdem.exe

clean:
	rm -f $(ODIR)/*.o *~ core
