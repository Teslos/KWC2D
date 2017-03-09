# GNU make file
RUNFILE = sim
NAGDIR = /opt/NAG/cllux07db
MKL = /opt/intel/mkl
CC = g++ 
FC = gfortran
LIBOPT = -L/usr/local/gfortran/lib

SRC = extruder.cc \
	outcolor.cc \
	sim2d.cc    \
	solve2dnew_iccg.cc  \
	simmain.cc \

OBJ = extruder.o \
	outcolor.o   \
	sim2d.o      \
	solve2dnew_iccg.o    \
	simmain.o  \


HEADERS = sim2d.h

# Options for debugger
#CFLAGS = -g -c -Wno-unused-variable -I$(NAGDIR)/include 
CFLAGS = -g -c -Wno-unused-variable 
FFLAGS = -g -c -cpp 
#LIBOPT = -L$(NAGDIR) 
# Options for Linux
#CFLAGS = -O2 -c 
LIBSOLVE = solvers/libsolve.a \
	   solvers/libslap.a
#LIB = -lc -lnagc-mkl -lpthread -lm -lblas 
LIB = -lc -lpthread -lm -lblas -lgfortran $(LIBSOLVE)

all:
	make -C solvers
	make $(RUNFILE)
	
$(RUNFILE) : $(HEADERS) $(OBJ)
	$(CC) $(LIBOPT) -g $(OBJ) -o $(RUNFILE) $(LIB)

extruder.o: extruder.cc
	$(CC) $(CFLAGS) $<
outcolor.o: outcolor.cc
	$(CC) $(CFLAGS) $<
sim2d.o: sim2d.cc
	$(CC) $(CFLAGS) $<
solve2dnew_iccg.o: solve2dnew_iccg.cc
	$(CC) $(CFLAGS) $<
simmain.o: simmain.cc
	$(CC) $(CFLAGS) $<

clean:
	rm -f $(OBJ) sim
