# makefile driver for chymera
# SET SPECIFIC OPTIONS IN MAKEFILE.XX
# NOTE FOR GCC: Set GOMP_STACKSIZE to some very large number.

 include Makefile.DEFS
 include $(MAKEFILE)

 DOPT = -DPASSIVE=$(PASSIVE) -DVERBOSE=$(VERBOSE) -DPASSIVE_INIT=$(PASSIVE_INIT) -DFLUID=$(FLUID)\
        -DWIGGLE=$(WIGGLE) -DEXTERNAL_POT=$(EXTERNAL_POT) -DFLUID_RESTART=$(FLUID_RESTART)\
        -DWIGGLE_RESTART=$(WIGGLE_RESTART) -DSELFGRAVITY=$(SELFGRAVITY) -DINDIRECT=$(INDIRECT)\
        -DPARTICLE=$(PARTICLE) -DBINARY=$(BINARY)

 ifeq ($(NOLIMIT),1)
  DOPT+=-DNOLIMIT
 endif

 OBJ = particle_module.o eos.o 3dhyd-main.o fft.o flux.o\
       housekeeping.o io.o misc.o\
       pot3.o rad.o source.o\
       state.o boundary.o blktri.o\
       avisc.o hybrid.o \
       CoolingPrescriptions.o RadTran.o \
 
 ifeq ($(EXTERNAL_POT),1)
   OBJ+=ExternalPot.o
 endif

 ifeq ($(FLUID),1)
   FOBJ = fluid_mod.o fluid_start.o\
          fluid_setup.o fluid_interp.o \
          fluid_advance.o fluid_writeout.o\
          fluid_restart.o interp_project.o\
          interp.o
 else
   FOBJ =
 endif 

 ifeq ($(WIGGLE),1)
    WOBJ = wiggle.o
 else
    WOBJ = 
 endif

    POBJ = passive.o

 ifeq ($(EXTERNALS),1)
    EOBJ = ExternalRoutines.o
 else
    EOBJ = 
 endif

#*******************************************************************************
#                           THE EXECUTABLES

#--------------------  Serial versions  ----------------------------------------
chymera-omp: dummy
	@if [ $(WIGGLE) -eq "0" ] ; then\
         if [ $(WIGGLE_RESTART) -ne "0" ]; then\
           echo 'WIGGLE IS SET TO ZERO WHILE WIGGLE_RESTART=1';\
           exit 1;\
         fi\
        fi
	@if [ $(EXTERNAL_POT) -eq "0" ] ; then\
         if [ $(INDIRECT) -ne "0" ]; then\
           echo 'EXTERNAL_POT SET TO ZERO WHILE INDIRECT=1';\
           exit 1;\
          fi\
        fi
	make target "TARGET= $(EXEDIR)/chymera-omp" \
 "FC=$(OMPF95)" "CFLAGS=$(CFLAGS) $(TARCH) $(DOPT)"  \
 "CC=$(CC)"  \
 "LFLAGS=$(LFLAGS)" \
 "OBJ=$(OBJ) $(FOBJ) $(WOBJ) $(POBJ) $(EOBJ)"

chymera-ser: dummy
	@if [ $(WIGGLE) -eq "0" ] ; then\
         if [ $(WIGGLE_RESTART) -ne "0" ]; then\
           echo 'WIGGLE IS SET TO ZERO WHILE WIGGLE_RESTART=1';\
           exit 1;\
         fi\
        fi
	@if [ $(EXTERNAL_POT) -eq "0" ] ;then\
         if [ $(INDIRECT) -ne "0" ] ; then\
           echo 'EXTERNAL_POT SET TO ZERO WHILE INDRECT=1';\
           exit 1;\
         fi\
	fi
	make target "TARGET= $(EXEDIR)/chymera-ser" \
 "FC=$(F95)" "CFLAGS=$(CFLAGS) $(TARCH) $(DOPT)"  \
 "CC=$(CC)"  \
 "LFLAGS=$(LFLAGS)" \
 "OBJ=$(OBJ) $(FOBJ) $(WOBJ) $(POBJ) $(EOBJ)"

target: $(OBJ)
	$(FC) $(LFLAGS) -o $(TARGET) $(OBJ)

dummy:

clean:
	-rm -f *.o *.mod

realclean:
	-rm -f *.o *.mod *.lst chymera-omp chymera-ser

movie: 
	echo 'Making a movie from rho3d. files in $(MOVDIR)'
	cd $(MOVDIR) 
	cp ../nlss/sigplot $(MOVDIR)
	mkdir img
	./sigplot
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

