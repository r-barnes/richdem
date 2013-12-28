#CC=i586-mingw32msvc-g++
#CFLAGS=-Wall -Wno-unknown-pragmas -DARCGIS -static-libgcc #-static-libstdc++ #-fopenmp
#In the 32-bit mingw32 compiler, the static-libstdc++ command is covered by the static-libgcc command

#CC=x86_64-w64-mingw32-g++
#CFLAGS=-Wall -Wno-unknown-pragmas -DARCGIS -static-libgcc -static-libstdc++ #-fopenmp

CC=g++
CFLAGS=-Wall -fopenmp -DNDEBUG #-DARCGIS #-lX11 -pthread #-ltbb
#-ffp-contract=fast enables forming of fused multiply-add instructions. Default is fast.

ODIR=obj
#PRE_FLAGS=-lgcov -g -fprofile-arcs -ftest-coverage
PRE_FLAGS=-O3
#PRE_FLAGS=-g

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

debug: $(OBJ) obj/test.o
	$(CC) $(PRE_FLAGS) -g -o richdem.exe $^ $(CFLAGS)
	du -hs ./richdem.exe

clean:
	rm -f $(ODIR)/*.o *~ core
