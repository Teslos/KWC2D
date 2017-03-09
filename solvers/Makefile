LIBSOLVE = libsolve.a
CC = g++
FC = gfortran
AR = ar
CFLAGS = -g -c
FFLAGS = -g -c
LIBOPT = -L/usr/local/gfortran/lib

SRC = nag.cc  \
      iccg.cc \
      icslap.f \
      xersla.f \
      mach.f 

OBJ = nag.o \
      iccg.o \
      icslap.o \
      xersla.o \
      mach.o

LIB = ./libslap.a

all: 
	make libsolve.a
	
nag.o: nag.cc
	$(CC) $(CFLAGS) $<
icslap.o: icslap.f
	$(FC) $(FFLAGS) $<
xersla.o: xersla.f
	$(FC) $(FFLAGS) $<
mach.o: mach.f
	$(FC) $(FFLAGS) $<

$(LIBSOLVE) : $(OBJ)
	$(AR) $(ARFLAGS) $(LIBSOLVE) $(OBJ) $(LIB)

clean:
	rm -f *.o libsolve.a
