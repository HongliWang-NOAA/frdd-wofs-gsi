#! /bin/sh -v

# Define the fortran compiler and options
include ../../../dtc/configure.gsi

FC   = $(SF90)
FLAGS= $(FFLAGS_DEFAULT) -xAVX
INCLD=  -I../../../dtc/include
LIBS =  -L../../../src -lgsi

all: read_diag_conv_ens read_diag_conv read_diag_rad innov_mean_conv innov_mean_radiance

read_diag_conv_ens: read_diag_conv_ens.o
	${FC} -o read_diag_conv_ens.exe ${FLAGS} read_diag_conv_ens.o ${LIBS}
	#rm read_diag_conv_ens.o

read_diag_conv: read_diag_conv.o
	${FC} -o read_diag_conv.exe ${FLAGS} read_diag_conv.o ${LIBS}
	#rm read_diag_conv.o

read_diag_rad: read_diag_rad.o
	${FC} -o read_diag_rad.exe ${FLAGS} read_diag_rad.o ${LIBS}
	#rm read_diag_rad.o

innov_mean_conv: innov_mean_conv.o
	${FC} -o innov_mean_conv.exe ${FLAGS} innov_mean_conv.o ${LIBS}
	#rm innov_mean_conv.o

innov_mean_radiance: innov_mean_radiance.o
	${FC} -o innov_mean_radiance.exe ${FLAGS} innov_mean_radiance.o ${LIBS}
	#rm innov_mean_radiance.o

.SUFFIXES : .f90  .o

.f90.o :
	${FC} ${FLAGS} ${INCLD} -c $<

clean:
	/bin/rm -f *.o *.exe
