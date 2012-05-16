#CC=i586-mingw32msvc-g++
#CFLAGS=-Wall -Wno-unknown-pragmas -DARCGIS #-fopenmp

CC=g++
CFLAGS=-Wall -fopenmp #-DARCGIS #-lX11 -pthread #-ltbb

ODIR=obj
#PRE_FLAGS=-lgcov -g -fprofile-arcs -ftest-coverage
PRE_FLAGS=-O3
DEPS = d8_methods.h data_structures.h dinf_methods.h interface.h data_io.h pit_fill.h utility.h flat_resolution.h debug.h visualize.h watershed.h

_OBJ = d8_methods.o dinf_methods.o interface.o data_io.o main.o pit_fill.o utility.o debug.o pit_fill_test.o watershed.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) $(PRE_FLAGS) -c -o $@ $< $(CFLAGS)

richdem: $(OBJ)
	$(CC) $(PRE_FLAGS) -o richdem.exe $^ $(CFLAGS)
	du -hs ./richdem.exe

debug: 

clean:
	rm -f $(ODIR)/*.o *~ core
