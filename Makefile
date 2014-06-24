# GNU make file
RUNFILE = sim
NAGDIR = /opt/NAG/cllux07db
MKL = /opt/intel/mkl
CC = g++ 
FC = gfortran

SRC = extruder.cc \
	outcolor.cc \
	sim2d.cc    \
	solve2dnew_iccg.cc  \
	simmain.cc \
	nag.cc \
	icslap.f \
	xersla.f \
	mach.f

OBJ = extruder.o \
	outcolor.o   \
	sim2d.o      \
	solve2dnew_iccg.o    \
	simmain.o  \
	nag.o \
	icslap.o \
	xersla.o \
	mach.o


HEADERS = sim2d.h

# Options for debugger
#CFLAGS = -g -c -Wno-unused-variable -I$(NAGDIR)/include 
CFLAGS = -g -c -Wno-unused-variable 
FFLAGS = -g -c -cpp 
#LIBOPT = -L$(NAGDIR) 
# Options for Linux
#CFLAGS = -O2 -c 

#LIB = -lc -lnagc-mkl -lpthread -lm -lblas 
LIB = -lc -lpthread -lm -lblas -lgfortran ./libslap.a

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
nag.o: nag.cc
	$(CC) $(CFLAGS) $<
icslap.o: icslap.f
	$(FC) $(FFLAGS) $<
xersla.o: xersla.f
	$(FC) $(FFLAGS) $<
mach.o: mach.f
	$(FC) $(FFLAGS) $<

clean:
	rm -f $(OBJ) sim
