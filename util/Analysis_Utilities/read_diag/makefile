#! /bin/sh -v

# Define the fortran compiler and options

include ../../../dtc/configure.gsi

FC   = $(SFC)
FLAGS= $(FFLAGS_DEFAULT)   
INCLD=  -I../../../dtc/include
LIBS =  -L../../../src -lgsi

all: read_diag_conv read_diag_rad read_diag_conv_ens

read_diag_conv_ens: read_diag_conv_ens.o
	${FC} -o read_diag_conv_ens.exe ${FLAGS} read_diag_conv_ens.o ${LIBS} 
	rm read_diag_conv_ens.o

read_diag_conv: read_diag_conv.o
	${FC} -o read_diag_conv.exe ${FLAGS} read_diag_conv.o ${LIBS} 
	rm read_diag_conv.o

read_diag_rad: read_diag_rad.o
	${FC} -o read_diag_rad.exe ${FLAGS} read_diag_rad.o ${LIBS} 
	rm read_diag_rad.o

.SUFFIXES : .f90  .o

.f90.o :
	${FC} ${FLAGS} ${INCLD} -c $<

clean:
	/bin/rm -f *.o *.exe
