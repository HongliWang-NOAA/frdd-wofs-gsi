#! /bin/sh -v

# Define the fortran compiler and options

include ../../../dtc/configure.gsi

FC   = $(SFC)
FLAGS= $(FFLAGS_DEFAULT)   
INCLD=  -I../../../dtc/include
LIBS =  -L../../../src -lgsi

SRC = $(wildcard *.f90)
DEST = $(SRC:.f90=.exe)

all: $(DEST)

%.exe: %.f90 
	${FC} -o $(<:.f90=.exe) ${FLAGS} ${INCLD} $< ${LIBS}
	#rm $<.o

.SUFFIXES : .f90  .o

.f90.o : %.f90
	${FC} ${FLAGS} ${INCLD} -c $<

clean:
	/bin/rm -f *.o *.exe
