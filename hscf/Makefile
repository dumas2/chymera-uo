CC = cc
FC = gfortran

OPT_FLAG = -fdefault-real-8 


CFLAGS  = $(OPT_FLAG)
FFLAGS  = $(OPT_FLAG)
LDFLAGS =

OBJS = hscf.f
LIBS =

all : $(OBJS) hscf

hscf : $(OBJS) hscf.f 
	$(FC) -o $@ $(LDFLAGS) $(OBJS) $(OPT_FLAG)

clean : 
	rm *.o fort.41 fort.42 fort.47 fort.71 fort.93 anal.dat wdout fort.2 hscf.out polyout fact 
