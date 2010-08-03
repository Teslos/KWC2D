# GNU make file
RUNFILE = sim
NAGDIR = /opt/NAG/cllux07db
CC = g++

SRC = extruder.cc \
	outcolor.cc \
	sim2d.cc    \
	solve2dnew.cc  \
	simmain.cc \
	nag.cc

OBJ = extruder.o \
	outcolor.o   \
	sim2d.o      \
	solve2dnew.o    \
	simmain.o  \
	nag.o

HEADERS = sim2d.h

# Options for debugger
CFLAGS = -g -c -Wno-unused-variable -I$(NAGDIR)/include 
LIBOPT = -L$(NAGDIR) 
# Options for Linux
#CFLAGS = -O2 -c 

LIB = -lc -lnagc-mkl -lpthread -lm -lblas 

$(RUNFILE) : $(HEADERS) $(OBJ)
	$(CC) $(LIBOPT) -g $(OBJ) -o $(RUNFILE) $(LIB)

extruder.o: extruder.cc
	$(CC) $(CFLAGS) $<
outcolor.o: outcolor.cc
	$(CC) $(CFLAGS) $<
sim2d.o: sim2d.cc
	$(CC) $(CFLAGS) $<
solve2dnew.o: solve2dnew.cc
	$(CC) $(CFLAGS) $<
simmain.o: simmain.cc
	$(CC) $(CFLAGS) $<
nag.o: nag.cc
	$(CC) $(CFLAGS) $<

clean:
	rm -f $(OBJ) sim
