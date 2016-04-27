# Makefile for poisson solver from chymera

#include Makefile.DEFS

# detect OS, 'linux' is output for Linux system, 'darwin' for OS X
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])


CC = mpicc
FC = mpif90

FOPT = -O3
OMPFLAGS = -fopenmp 
FINC = -fdefault-real-8

ifneq ($(OSLOWER),darwin)
   DATA_SIZE = -mcmodel=medium
endif

CFLAGS =
TARCH  = -frecord-marker=4
CONSERVEOPT =

FFLAGS := $(FOPT) $(FINC)
FFLAGS += $(DATA_SIZE) 

OBJS = defines_mod.o fft.o pot3.o boundary.o blktri.o
# Objects used for relaxation routine.
OBJS_RELAX = defines_mod.o io.o multigrid.o relaxMod.o 

#------------------------------------------------------------------

all: poisson_solver

poisson_solver: $(OBJS) poisson_solver.o
	$(FC) -o poisson_solver poisson_solver.o $(OBJS)

poisson_multigrid: io.o multigrid.o poisson_multigrid.o
	$(FC) -o poisson_multigrid poisson_multigrid.o multigrid.o io.o

poisson_relax: $(OBJS_RELAX) poissonRelax.o
	$(FC) -o poisson_relax poissonRelax.o $(OBJS_RELAX)

check: poisson_solver
	poisson_solver

check-multi: poisson_multigrid
	poisson_multigrid

check-relax: poisson_relax
	poisson_relax
clean:
	-rm -f *.o *.mod poisson_solver out_* poisson_multigrid poisson_relax

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

