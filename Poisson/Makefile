# Makefile for poisson solver from chymera

include Makefile.DEFS

CC = gcc
FC = gfortran

FOPT = -O3
DATA_SIZE = -mcmodel=medium
OMPFLAGS = -fopenmp

CFLAGS =
TARCH  = -frecord-marker=4
CONSERVEOPT =

FFLAGS := $(FOPT) $(FINC)
FFLAGS += $(DATA_SIZE) 

OBJS = defines_mod.o fft.o pot3.o boundary.o blktri.o

#------------------------------------------------------------------

all: poisson_solver

poisson_solver: $(OBJS) poisson_solver.o
	$(FC) -o poisson_solver poisson_solver.o $(OBJS)

check: poisson_solver
	poisson_solver

clean:
	-rm -f *.o *.mod poisson_solver

########################################################################
#  OBJECTS

boundary.o: boundary.F
	$(FC) $(FINC) $(DATA_SIZE) $(CONSERVEOPT) -c $<
# DO NOT OPTIMIZE BOUNDARY UNLESS OTHERWISE NOTED ON WIKI. ONLY GCC PROBLEM.

%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

%.o: %.F
	$(FC) $(FFLAGS) -c $<

interp.o: interp.c
	$(CC) $(COPT) $(CLINK) $(CINC) -c $<

