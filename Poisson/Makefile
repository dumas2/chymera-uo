# Makefile for poisson solver from chymera

include Makefile.DEFS

 OMPF95 = gfortran -fopenmp -mcmodel=medium
 PREPROCESSOR = -x f95-cpp-input
 CFLAGS =
 TARCH  = -frecord-marker=4
 FOPT =  -O3 $(PREPROCESSOR)
 CONSERVEOPT = $(PREPROCESSOR)

 FFREE = -ffree-form
 FFIXED = -ffixed-form
 CC = gcc
 COPT = -m64
# CLINK = -pg 
 CLINK =

 OBJ = main.o fft.o pot3.o \
       boundary.o blktri.o\

#------------------------------------------------------------------
chymera-omp: dummy
	make target "TARGET= $(EXEDIR)/chymera-omp" \
 "FC=$(OMPF95)" "CFLAGS=$(CFLAGS) $(TARCH) $(DOPT)"  \
 "CC=$(CC)"  \
 "LFLAGS=$(LFLAGS)" \
 "OBJ=$(OBJ) $(FOBJ) $(WOBJ) $(POBJ) $(EOBJ)"

target: $(OBJ)
	$(FC) $(LFLAGS) -o $(TARGET) $(OBJ)

dummy:
	
clean:
	-rm -f *.o *.mod
########################################################################
#  OBJECTS

boundary.o: boundary.f
	$(FC) $(CFLAGS) $(CONSERVEOPT) $(FFIXED) $(FINC) -c $<
# DO NOT OPTIMIZE BOUNDARY UNLESS OTHERWISE NOTED ON WIKI. ONLY GCC PROBLEM.

%.o:%.f
	$(FC) $(CFLAGS) $(FOPT) $(FFIXED) $(FINC) -c $<

%.o:%.f90
	$(FC) $(CFLAGS) $(FOPT) $(FFREE) $(FINC) -c $<

%.o:%.F
	$(FC) $(CFLAGS) $(FOPT) $(FFREE) $(FINC) -c $<

interp.o: interp.c
	$(CC) $(COPT) $(CLINK) $(CINC) -c $<

