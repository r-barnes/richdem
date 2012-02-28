CC=g++
ODIR=obj
CFLAGS=-fopenmp -g
DEPS = d8_methods.h  data_structures.h  dinf_methods.h  interface.h  load_data.h pit_fill.h  utility.h

_OBJ = d8_methods.o data_structures.o dinf_methods.o interface.o load_data.o main.o pit_fill.o utility.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

richdem: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f $(ODIR)/*.o *~ core
