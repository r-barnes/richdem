IDIR =.
CC=g++

#Compilation flags. Feel free to add optimization options here
CFLAGS=-I$(IDIR) --std=c++11 -fpermissive -O3 `gdal-config --cflags` -I/home/rick/projects/watershed/richdem/include -DUSEGDAL

ODIR=.

LIBS=-lm -lgdal 

_DEPS = Node.h dem.h utils.h perlin.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ =dem.o fillDEM.o main.o utils.o

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

fill: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
