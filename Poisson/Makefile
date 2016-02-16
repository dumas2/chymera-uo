# Makefile for poisson solver from chymera

include Makefile.DEFS

CC = gcc
FC = gfortran

OMPF95 = gfortran -fopenmp -mcmodel=medium
PREPROCESSOR = -x f95-cpp-input
CFLAGS =
TARCH  = -frecord-marker=4
FOPT =  -O3 $(PREPROCESSOR)
CONSERVEOPT = $(PREPROCESSOR)

FFREE = -ffree-form
FFIXED = -ffixed-form
COPT = -m64
# CLINK = -pg 
CLINK =

FFLAGS := -O3

OBJS = defines_mod.o fft.o pot3.o  \
       boundary.o blktri.o         \

#------------------------------------------------------------------

all: poisson_solver

poisson_solver: poisson_solver.o $(OBJS)
	$(FC) -o poisson_solver $< $(OBJS)

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
	-rm -f *.o *.mod
########################################################################
#  OBJECTS

boundary.o: boundary.f
	$(FC) $(CFLAGS) $(CONSERVEOPT) $(FFIXED) $(FINC) -c $<
# DO NOT OPTIMIZE BOUNDARY UNLESS OTHERWISE NOTED ON WIKI. ONLY GCC PROBLEM.

%.o: %.f
	$(FC) $(CFLAGS) $(FOPT) $(FFIXED) $(FINC) -c $<

%.o: %.f90
	$(FC) $(CFLAGS) $(FOPT) $(FFREE) $(FINC) -c $<

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

%.o: %.F
	$(FC) $(CFLAGS) $(FOPT) $(FFREE) $(FINC) -c $<

interp.o: interp.c
	$(CC) $(COPT) $(CLINK) $(CINC) -c $<

