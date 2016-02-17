# Makefile for poisson solver from chymera

include Makefile.DEFS

CC = gcc
FC = gfortran

OMPF95 = gfortran -fopenmp -mcmodel=medium
PREPROCESSOR = -x f95-cpp-input
CFLAGS =
TARCH  = -frecord-marker=4
FOPT =  -O3
CONSERVEOPT =

FFREE = -ffree-form
FFIXED = -ffixed-form
COPT = -m64
# CLINK = -pg 
CLINK =

FFLAGS := $(FOPT)
FFLAGS += $(FINC)

OBJS = defines_mod.o fft.o pot3.o  \
       boundary.o blktri.o         \

#------------------------------------------------------------------

all: poisson_solver

poisson_solver: $(OBJS) poisson_solver.o
	$(FC) -o poisson_solver poisson_solver.o $(OBJS)

chymera-omp: dummy
	make target "TARGET= $(EXEDIR)/chymera-omp"  \
 "FC=$(OMPF95)" "CFLAGS=$(CFLAGS) $(TARCH) $(DOPT)"  \
 "CC=$(CC)"  \
 "LFLAGS=$(LFLAGS)" \
 "OBJS=$(OBJS) $(FOBJ) $(WOBJ) $(POBJ) $(EOBJ)"

target: $(OBJS)
	$(FC) $(LFLAGS) -o $(TARGET) $(OBJS)

dummy:

check: poisson_solver
	poisson_solver

clean:
	-rm -f *.o *.mod poisson_solver
########################################################################
#  OBJECTS

boundary.o: boundary.F
	$(FC) $(FINC) $(CONSERVEOPT) -c $<
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

