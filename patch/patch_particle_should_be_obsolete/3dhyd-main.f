***********************************************************************
C************************************************************************
C...THIS IS A NEW VERSION OF THE HYDRO CODE WHICH IS SECOND-ORDER IN BOTH
C...SPACE AND TIME.  S, T ARE FACE-CENTERED; A,RHO,& EPS ARE CELL-
C...CENTERED.   AN ENERGY EQUATION IS INCLUDED FOR ENERGY TRANSPORT.
!
C...Also, opacities added such that cooling can be calculated. Need to 
C...define real dimensions of the disk in the source routine.  acm 2002.
C...ARTIFICAL VISCOSITY IS USED TO TREAT SHOCKS AND HEAT THE MODEL.
!
! Substantial changes have been made to the hydrodynamics code.  These
! changes include fluxing eps directly, and not eps^(1/gamma).  The 
! internal energy of H2 is directly calculated at the beginning of the
! simulations, and the temperature is determined based on the eps, rho
! and the assumed ortho:para hydrogen statistics. The pressure is 
! calculated by p~rho T.  The mean molecular weight for the gas
! is calculated based on what is used for calculating the specific 
! internal energy.  This does create a slight difference when compared
! with D'Alessio's opacities, but it should be small.  In addition,
! a new radiative transfer routine is used, which uses vertical rays.
! Nora Bolig (ACB)
!
C...THE FOLLOWING SUBROUTINES ARE CALLED:
c
C     *SETUP   :  Read starting models, impose perts, etc. (io.f)
!     *initengtable : Creates specific energy table for calculating temp. (initengtable.f) ACB
!     *TempFind:  Finds temperature for each cell. (source.f) ACB
!     *TempFindSpec: Finds temperature based on input eps, rho. (source.f) ACB
C     *VLIMIT  :  Limit maximum velocity (usually < 2*SOUND). (housekeeping.f)
C     *DELTA   :  Calculate maximum safe delta time. (housekeeping.f)
C     *RITE    :  Write output information. (io.f)
C     *SLOPE   :  Calculate van Leer Slope. (flux.f)
C     *VLI     :  Calculate van Leer Monotonic Interpolation. (flux.f)
C     *FLUX    :  Advect S,T,A,RHO, and EPS. (flux.f)
C     *CLEANUP :  Fix velocities, densities and energy on boundaries.
c                 (housekeeping.f)
C     *SOURCE   :  Source S,T,A, and EPS. (source.f)
C     *VELOCITY :  From momentum densities, calculate velocities.
C     *CENTMASS :  Calculate Center of Mass. (housekeeping.f)
C     *STATE    :  Equation of State. (state.f)
C     *RAD     : Old routine for "radiative physics." (rad.f) 
C     *REALTR  :  -|Together These Perform a  (fft.f)
C     *FFT     :  -|Fast Fourier Transform.   (fft.f)
C     *POT3    :  Potential Solver. (pot3.f)
C     *ZAXPHI  :  -|                (pot3.f)
C     *BLKTRI  :  -|                (pot3.f)
C     *BLKTR1  :  -|Various Functions used in/with the (pot3.f)
C     *PRDCT   :  -|Potential Solvers. (pot3.f)
C     *COMBP   :  -|                   (pot3.f)
C     *TQLRAT  :  -|                   (pot3.f)
C     *SETBDY  :  Initialization before BDYGEN. (boundary.f)
C     *BDYGEN  :  Boundary Potential Solver.    (boundary.f)
C     *ExternalPot: External potential
C     *SORT    :  Sort. (housekeeping.f)
C     *CLEARUP :  Clearup after SETMODE or DAMP call. (housekeeping.f)
C     *TORQUE  :  Calculate instaneous torques.
c     *dumpphi :  dump the potential grid for comparison (io.f)
c     *cyl3d   :  dump density grid for use in gen. 3d images (io.f)
c     *lowresrandom: low resolution random pert. (io.f)
!     *hybrid  : Vertical ray radiative physics routine. (hybrid.f) ACB
!     *initengtable: initialize energy table. (initengtable.f) ACB
!     *tempfind: find temp based on energy table. (source.f) ACB
!     *tempfindspec: find the temp of just one cell based on the energy
!                    table. (source.f) ACB
C***********************************************************************
C     
      use particle
      IMPLICIT real*8 (a-h,o-z)
!#include "hydroparam.h"
!#include "globals.h"
!#include "units.h"

      COMMON /INSIDE/TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT
      COMMON /TIMEST/INDX,ISOADI,ALLOW,DMAX,CHGMAX
      COMMON /COEFS/COEF(POT3JMAX2,POT3KMAX2,LMAX2,2)
      COMMON /ITS/ITSTRT,ITSTOP,ITSTEP

      real*8 epsjr,rhojr,ommax,mirp,pdelt
      common /misc/ epsjr,rhojr,ommax

C The following arrays are local to the main program, and should not need 
C to be placed in common. However, some systems may try to put them on the
C stack, which may result in stack overflow. Putting them in a common block
C guarantees they will go on the heap. ! COMMON BLOCK REMOVED. ACB
      REAL*8 SS(JMAX2,KMAX2,LMAX),
     &       TT(JMAX2,KMAX2,LMAX),
     &       AA(JMAX2,KMAX2,LMAX),
     &       RRHO(JMAX2,KMAX2,LMAX),
     &       EEPS(JMAX2,KMAX2,LMAX)



      CHARACTER  tim*6,rhofile*80,tempfile*80,starfile*80,
     &           restart_star*80
      CHARACTER  coeffile*80,index*6,modefile*80,centfile*80
      CHARACTER  tcoeffile*80,tcoeftenfile*80,massfile*80,epsfile*80
      real*8     massc(jmax2), massf(jmax2),mdot,mout
      real*8     sumeven(8),sumodd(8),theta(lmax),tcool,theat,tflux,tirr
      real*8     volume(jmax2),totcoef(8)
      REAL*8     totcoefcym(10,8)
      real*8     start_time,finish_time,odelt
      DATA ISTRDN,CHANGD/25,2.00/,NCONS,DELCON/10,2.0/
!$      integer    OMP_GET_MAX_THREADS

      REAL*8 Oross(jmax2,kmax2,lmax),Oplck(jmax2,kmax2,lmax)
     &     ,Oabs(jmax2,kmax2,lmax),Otot(jmax2,kmax2,lmax)
      logical, save::kick_switch=.false.
      COMMON /COOLINGSHARED/Oross,Oplck,Oabs,Otot,sum

!
! FIRST TOUCH: INITIALIZE ALL ARRAYS.  THIS IS NOT ONLY GOOD PRACTICE, BUT
! IT CAN SIGNICANTLY SPEED UP THE OPENMP CODE.  THE ARRAYS ARE MORE 
! EFFICIENTLY DISTRIBUTED AMONG PROCESSORS IF TOUCHED ALL AT ONCE.

      call first_touch()      
      call CPU_TIME(start_time)
      print *,"---------------------------------------------------------&
     &-----"
      print *,"CCCCCCC  H     H  Y     Y  M     M  EEEEEEE  RRRRRR     A&
     &AA   "
      print *,"C        H     H  Y     Y  MM   MM  E        R     R   A &
     &  A  "
      print *,"C        H     H   Y   Y   M M M M  E        R    R   A  &
     &   A "
      print *,"C        HHHHHHH    YYY    M  M  M  EEEEEEE  RRRRR    AAA&
     &AAAA "
      print *,"C        H     H     Y     M     M  E        R    R   A  &
     &   A "
      print *,"C        H     H     Y     M     M  E        R     R  A  &
     &   A "
      print *,"CCCCCCC  H     H     Y     M     M  EEEEEEE  R     R  A  &
     &   A "
      print *,"---------------------------------------------------------&
     &-----"
      print *,"                        VERSION 2.1.0"
!$    write(6,*) 'MP enabled.  Nproc = ', OMP_GET_MAX_THREADS()

C....START EVOLUTION

      TEMPCD=ten*ten*CHANGD
      MSTORE=TEMPCD
      MAXTRM=10
      ISYM=2
      ENON=one
      CONSTP=0
      BDYTEM=1.d-3
      CORMAS=1.d-2
      A1NEWR=zero
      A1NEWZ=zero
      OMMAX=zero ! just define it.  Tracking down old uses still
      KWFW=int(log10(dble(KMAX))/log10(two))-1 ! this is used in pot3.f
      
cbkp..cs=0 for no shear viscosity (not implemented yet, so don't diddle!)
crpl..xxxtodo: make this parameter a read-in.      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ DATA AND SET UP GRID.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL SETUP(ITSTRT,ITSTOP,IDIAG,ISOADI,ISTOR,ITSTEP,
     &           ISYM,MAXTRM)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OPEN RECORDING FILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(index,'(i6.6)')itstop
      coeffile='coefs.'//index
      modefile='modes.'//index
      centfile='c_o_m.'//index
      tcoeffile='tcoef.'//index
      tcoeftenfile='tctot.'//index
      massfile='massflow.'//index
      epsfile='coolheat.'//index

      OPEN(unit=18,file=coeffile)
      OPEN(unit=19,file=centfile)
      open(unit=13,file=tcoeffile)
      open(unit=15,file=tcoeftenfile)
      open(unit=17,file=massfile,form='formatted')
      open(unit=21,file=epsfile,form='formatted')

#if WIGGLE>0
        starfile='starrysim.'//index//' '
        restart_star='starry_restart.'//index//' '
        open(unit=987,file=trim(starfile),form='FORMATTED')
        open(unit=988,file=trim(restart_star),form='FORMATTED')
#endif

      CALL RITE(3,-1,  1,1,1,  1,1,1,  1,1,1)
!$OMP PARALLEL DEFAULT(SHARED)
      CALL CLEARUP
!$OMP END PARALLEL
      CALL RAD(ISOADI)

      DMAX=nine*DEN/ten
      CHGMAX=0.0001
      MIRP=two*pi/ommax

#if FLUID>0
        call Fluid_Setup()
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c---------------------begin main loop-----------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(6,*) 'Beginning simulation'
      DO 90 ITSTEP=ITSTRT,ITSTOP
         INDX=ITSTEP-ITSTRT
         if(mod(ITSTEP,1).eq.0) write(6,10000) itstep,time,delt        !dkb
10000    format('step=',i8,'   time=',1pe26.13,'   delt=',1pe26.15)  !dkb
         print "(A,7(1X,1pe15.8))", 
     &       "TRACKPARTICLE1",time,x_p(1),y_p(1),z_p(1),
     &       vx_p(1),vy_p(1),vz_p(1)
         print "(A,7(1X,1pe15.8))", 
     &       "TRACKPARTICLE2",time,x_p(2),y_p(2),z_p(2),
     &       vx_p(2),vy_p(2),vz_p(2)


         IF(MOD(ITSTEP,IDIAG).EQ.0) then
            IHEAD=-1
            IPRINT=1
            if(MOD(ITSTEP,10000).EQ.0)
     &      cALL RITE(2,-1,ITSTOP,ITSTRT,ITSTOP,ISTOR,1,1,1,1,1)
         else 
            IPRINT=0
            IHEAD=0
         end if


C... Calculation of accreted mass and boundary outflow 
C... tmassini = initial disk mass 
C... tmass    = current disk mass
C... tmassadd = added mass by resetting bgrnd to rholmt
C... tmassacc = accreted mass
C... tmassout = mass flowing out of the top, bottom and outer boundary

         mdot=zero
         mout=zero
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(J,K,L) REDUCTION(+:mdot,mout)
!$OMP DO SCHEDULE(STATIC)
         do l=1,lmax
            do j=jmin,jmax
               mout=mout+T(j,kmax1,l)*dtheta*(R(j+1)**2-R(j)**2)
            enddo
         enddo
!$OMP END DO 
!$OMP DO SCHEDULE(STATIC) 
         do l=1,lmax
           do k=2,kmax
               mout=mout+S(jmax1,k,l)*two*dtheta*R(jmax1)*(zof3n)
            enddo
         enddo
!$OMP END DO 
      if (jmin.gt.2) then
!$OMP DO SCHEDULE(STATIC) 
         do l=1,lmax
            do k=2,kmax
               mdot=mdot+S(jmin,k,l)*two*dtheta*R(jmin)*(zof3n)
            enddo
         enddo
!$OMP ENDDO
      endif
!$OMP END PARALLEL
      tmassout=tmassout+max(mout*delt,zero)
      tmassacc=tmassacc+(abs(mdot)*delt)

         if(mod(indx,100)==0)
     &write(17,20)itstep,time,tmass,tmassadd,tmassout,tmassacc,mdot
 20      format(i8,1x,6(1pe24.16,1x))

         odelt=delt
         CALL DELTA(MOD(ITSTEP,IDIAG).EQ.0)
         call particle_timestep(pdelt)
         print *, delt,pdelt
         delt=min(delt,pdelt)

         TIME=TIME+DELT

C.............................................C
C.....START SECOND ORDER TIME INTEGRATION.....C
C.............................................C

C..(1) 1/2 SOURCE S, T, A, RHO, EPS.


         DELT=DELT*half
         ODELT=ODELT*half

#if FLUID>0
         call fluidtrace(itstep,itstop,.true.)
#endif

         CALL SOURCE
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,l)                           &
!$OMP&  SHARED(delt)

         CALL VELOCITY

         CALL VLIMIT

C..(2) STORE QUANTITIES.

!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
           DO K=1,KMAX2
             DO J=1,JMAX2
                SS(J,K,L)=S(J,K,L)
                TT(J,K,L)=T(J,K,L)
                AA(J,K,L)=A(J,K,L)
                RRHO(J,K,L)=RHO(J,K,L)
                EEPS(J,K,L)=EPS(J,K,L)
             ENDDO
           ENDDO
         ENDDO
C$OMP END DO
                      
C..(3) 1/2 FLUX S, T, A, RHO, EPS.

         CALL FLUX(S,T,A,RHO,EPS)

         CALL VELOCITY

         CALL VLIMIT

C..(4) 1 FLUX SS, TT, AA, RRHO, EEPS, 
C........BUT USE S, T, A, RHO, EEPS TO CALCULATE FLUXES.

!$OMP END PARALLEL

         DELT=two*DELT
         ODELT=two*ODELT

#if FLUID>0
         call fluidtrace(itstep,itstop,.false.)
#endif


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,l)                           &
!$OMP&  SHARED(delt)

         CALL FLUX(SS,TT,AA,RRHO,EEPS)

C..(5) UPDATE S, T, A, RHO, EPS.

!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
           DO K=1,KMAX2
             DO J=1,JMAX2
               S(J,K,L)=SS(J,K,L)
               T(J,K,L)=TT(J,K,L)
               A(J,K,L)=AA(J,K,L)
               RHO(J,K,L)=RRHO(J,K,L)
               EPS(J,K,L)=EEPS(J,K,L)
             ENDDO
           ENDDO
         ENDDO
C$OMP END DO NO WAIT

C..(6) 1/2 SOURCE S, T, A, RHO, EPS.

         CALL VELOCITY

         CALL VLIMIT

         CALL CLEARUP

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX
       do K = 1, KMAX2
        do J = 1, JMAX2
         phi(J,K,L) = zero
        enddo
       enddo
      enddo
!$OMP END DO 
!$OMP END PARALLEL

         call set_particle_density()

         CALL RAD(ISOADI)

         REDGE=R(JMAX1)

         CALL BDYGEN(MAXTRM,ISYM,REDGE)

         CALL POT3(8,IPRINT)

#if WIGGLE>0
         call wiggle(rhotot,rhf,zhf,vx,vy,fx,fy,delt,pi,rpstar,phi_star,
     &        rof3n,zof3n,dtheta,gsoft,JMAX,KMAX,LMAX,JMIN,.false.,
     &        restart_wiggle)
         if(mod(INDX,10)==0)then
           write(987,'(9(1X,1pe15.8))') time,vx,vy,fx,fy,rpstar,
     &     phi_star,rpstar*cos(phi_star),rpstar*sin(phi_star)
         endif
#endif

      tmass=zero
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,j,l)                           &
!$OMP&  SHARED(delt,tmassacc,zof3n,dtheta,tmass,mass_star)

!$OMP DO SCHEDULE(STATIC) REDUCTION(+:tmass)
            do l=1,lmax
               do k=2,kmax
                  do j=jmin,jmax
                     tmass=tmass+rho(j,k,l)*half*dtheta*zof3n
     &                    *(r(j+1)**2-r(j)**2)
                  enddo
               enddo
            enddo
!$OMP END DO
!$OMP MASTER
            tmass=tmass*two
!$OMP END MASTER 

#if EXTERNAL_POT > 0
         call ExternalPot()
#endif

!$OMP MASTER
         DELT=DELT*half
         ODELT=ODELT*half
!$OMP END MASTER

         CALL STATE
!$OMP END PARALLEL

         CALL SOURCE
         call particle_fullstep(odelt,delt)

C..(7) UPDATES DUE TO ENERGY EQUATION

!$OMP PARALLEL DEFAULT(SHARED)

         CALL STATE

         CALL VELOCITY

         CALL VLIMIT

         CALL CLEARUP

!$OMP END PARALLEL
         DELT=two*DELT
         ODELT=two*ODELT

C.............................................C
C.....COMPLETED ONE TIME STEP INTEGRATION.....C
C.............................................C




C...Written Output for result file...
c....Slices in j,k,l.................

          IF(MOD(ITSTEP,IDIAG).eq.0) then
         WRITE(3,10200) ITSTEP,TIME                                    !dkb
10200    FORMAT(//,' --------------------------------------  ',           !dkb
     &         'STEP=',i8,'   TIME= ',1pe15.8,                            !dkb
     &         '   -------------------------------------------------')    !dkb
!
! CALLING SUBROUTINE RITE(IWHAT,IHEAD,JST,JSP,JSK,KST,KSP,KSK,LST,LSP,LSK)
!
            LHALF=LMAX/2
            CALL RITE(1,-1,  2,JMAX1,1,  2,2,1,  1,1,1)
            CALL RITE(1, 1,  JMAX/4,JMAX/4,1,  2,KMAX1,1,  1,1,1)
            CALL RITE(1, 1,  JMAX/2,JMAX/2,1,  2,KMAX1,1,  1,1,1)
            CALL RITE(10,1,  1,1,1,  1,1,1,  1,1,1)
         END IF

         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENERATE OUTPUT FILES (RHO, TEMP), AND WRITE TO RECORD FILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (jmin.eq.2) then
            timc=int(10.d0*time/MIRP)
         else 
            timc=int(5.d0*time*omega(jmin+5,2,1)/pi)
         endif

         IF(MOD(itstep,IDIAG).EQ.0) THEN

            call dump_particles(itstep)

            write(tim,'(i6.6)')itstep

            rhofile='rho3d.'//tim
            OPEN(UNIT=14,FILE=rhofile,FORM='UNFORMATTED')
            write(14)rho
            write(14)time
            CLOSE(14)  


            tempfile='temperat3d.'//tim
            OPEN(UNIT=23,FILE=tempfile,FORM='UNFORMATTED')
            write(23)tempK
            write(23)time
            CLOSE(23)  

C	    call qcalc(tim)

            write(21,100)itstep,time,totcool,totdflux,totheat,totirr,
     &                   etotfl,eflufftot

            jstart=jmin
            do j=jstart,jmax
               if (lambda(j,2,1).ne.0.0) then
                  tcool=eps(j,2,1)/(abs(lambda(j,2,1))*torp)
               else
                  tcool=0.0
               endif
               if (divflux(j,2,1).ne.0.0) then
                  tflux=eps(j,2,1)/(abs(divflux(j,2,1))*torp)
               else
                  tflux=0.0
               endif
               if (hgamma(j,2,1).ne.0.0) then
                  theat=eps(j,2,1)/(hgamma(j,2,1)*torp)
               else
                  theat=0.0
               endif
               if (igamma(j,2,1).ne.0.0) then
                  tirr=eps(j,2,1)/(igamma(j,2,1)*torp)
               else
                  tirr=0.0
               endif

               write(21,101)j,eps(j,2,1),lambda(j,2,1),divflux(j,2,1)
     &              ,hgamma(j,2,1),igamma(j,2,1),tcool,tflux,theat,tirr
     &              ,tau(j,2,1,1),TempK(j,2,1),TeffK(j,1),TphK(j,1)
            enddo

            write(21,105)itstep,time,totcool,totdflux,totheat,totirr

            do i=1,3
               if (i.eq.1) j=30
               if (i.eq.2) j=100
               if (i.eq.3) j=200
               do k=2,kmax
                  if (lambda(j,k,1).ne.0.0) then
                     tcool=eps(j,k,1)/(abs(lambda(j,k,1))*torp)
                  else
                     tcool=0.0
                  endif
                  if (divflux(j,k,1).ne.0.0) then
                     tflux=eps(j,k,1)/(abs(divflux(j,k,1))*torp)
                  else
                     tflux=0.0
                  endif
                  if (hgamma(j,k,1).ne.0.0) then
                     theat=eps(j,k,1)/(hgamma(j,k,1)*torp)
                  else
                     theat=0.0
                  endif
                  if (igamma(j,k,1).ne.0.0) then
                     tirr=eps(j,k,1)/(igamma(j,k,1)*torp)
                  else
                     tirr=0.0
                  endif
                  write(21,106)j,k,eps(j,k,1),lambda(j,k,1),divflux(j
     &                 ,k,1),hgamma(j,k,1),igamma(j,k,1),tcool,tflux
     &                 ,theat,tirr,tau(j,k,1,1),TempK(j,k,1)
               enddo
            enddo

         END IF

 100     format(///,'STEP',i9,' TIME',1pe24.16,' TOTAL COOLING',1pe24.16
     &        ,' TOTAL DIVFLUX',1pe24.16,' TOTAL HEATING',1pe24.16,
     &        ' TOTAL IRRADIATION',1pe24.16,' TOTAL FLUXLOSS',1pe24.16
     &        ,' TOTAL EFLUFF LOSS',1pe24.16,
     &          //,3X,'J',6X,'EPS',9X,'LAMBDA',6X,'DIVFLUX',6X
     &        ,'HGAMMA',7X,'IGAMMA',5X,'Tcool(orp)',1X,
     &        'Tdflx(orp)',1X,'Theat(orp)',1X,'Tirr(orp)',5X,'TAU',6X,
     &        'Tmid(K)',4X,'Teff(K)',5X,'Tph(K)'
C    &        ,2X,'SDen(g/cm2)'
     &        )


 105     format(/,'STEP',i9,' TIME',1pe24.16,' TOTAL COOLING',1pe24.16
     &        ,' TOTAL DIVFLUX',1pe24.16,' TOTAL HEATING',1pe24.16
     &        ,' TOTAL IRRADIATION',1pe24.16,//,3X,'J',4X,'K',6X,'EPS'
     &        ,9X,'LAMBDA',6X,'DIVFLUX',6X,'HGAMMA',7X,'IGAMMA',5X
     &        ,'Tcool(orp)',1X,'Tdflx(orp)',1X,'Theat(orp)',2X
     &        ,'Tirr(orp)',4X,'TAU',6X,'Temp(K)')
                  
 101     format(i4,1x,5(1pe12.5,1x),8(1pe10.2,1x))
 106     format(i4,1x,i4,1x,5(1pe12.5,1x),6(1pe10.2,1x))

C...Store Coefficient Information....................
C...Use to extract modes and rho in Equatorial Plane.
         

         IF(MOD(ITSTEP,IDIAG).EQ.0) THEN

c...Equatorial fourier mode information (M=1-LMAX/2)

            CALL RITE(0,IHEAD,ITSTEP,INDX,1,1,1,1,ISTRDN,MSTORE,MS)
            LMWR9 = LMAX/2      
            WRITE(9,103)ITSTEP,TIME,(((COEF(J,2,L,I),L=1,LMWR9 )
     &          ,I=1,2),J=2,JMAX)
 103        FORMAT('COEF AT STEP ',I8,' TIME',1PE13.4,/,(1P10E13.5))
            WRITE(9,124)(RHF(I),I=2,JMAX)
 124        FORMAT(' HALF RADII:',/,(1P10E13.5))

c...Calculate total Fourier components m=1-8

            do j=2,jmax2
               massf(j)=0.d0
               massc(j)=0.d0
            end do
            
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,j,l,m)                         &
!$OMP&  SHARED(zof3n,dtheta)
!$OMP DO SCHEDULE(STATIC)
            do m=1,8
               sumeven(m)=0.d0
               sumodd(m)=0.d0
               totcoef(m)=0.d0
               do i=1,10
                  totcoefcym(i,m)=0.d0
               end do
            end do
!$OMP END DO NOWAIT
c
c....Calculate SQRT(a**2 +b**2)
c....don't forget normalizations! 2x for other half of
c....disk. Normalization by A(0) hidden because A(0)=
c....total mass = 1.
c
!$OMP DO SCHEDULE(STATIC)
             do l=1,lmax
                theta(l)=dtheta*(l-1)
             end do
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
             do j=2,jmax1
                volume(j)=0.5d0*dtheta*zof3n*(r(j+1)**2.-r(j)**2.)
             end do
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
             do m=1,8

               do j=2,jmax
               do k=2,kmax
               do l=1,lmax
                  sumeven(m)=sumeven(m)+rho(j,k,l)
     &                 *volume(j)*cos(m*theta(l))
                  sumodd(m)=sumodd(m)+rho(j,k,l)
     &                 *volume(j)*sin(m*theta(l))
               end do
               end do
               end do

               totcoef(m)=
     &              4.*SQRT(sumeven(m)**2. +sumodd(m)**2.)


            end do
!$OMP END DO nowait
!$OMP END PARALLEL

            write(13,104)time,(totcoef(m),m=1,8)

C...Diagnostic Information for COM and M=1 power...

            CALL CENTMASS(DISP)

            avgm = 0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j)                            &
!$OMP&  SHARED(jreq) REDUCTION(+:avgm)
            do j=2,JREQ
               avgm=avgm + coef(j,2,1,1)
            enddo
!$OMP END PARALLEL DO
            avgm = avgm/(JREQ-1)
            write(19,1111) itstep,time,disp,avgm
 1111              format(i8,1x,1pe13.4,1x,1pe13.4,1x,1pe13.4)



         END IF

 104     format(1pe13.4,8(1x,1pe13.4))

         
 90   CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c---------------end main loop-------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if FLUID>0
        call fluid_writerestart(itstop)
#endif

#if WIGGLE>0
           write(988,'(8(1X,1pe15.8))') time,vx,vy,rpstar,
     &     phi_star,rpstar*cos(phi_star),rpstar*sin(phi_star),delt
#endif

      call dump_particles(itstop)
      call clean_stop_particles
      ITSTEP = ITSTEP - 1

C.....Finished with Current Run Here.....C
C.........Run Diagnostics Follow.........C
     
      CALL CENTMASS(DISP)
      CALL RITE(2,-1,ITSTOP,ITSTRT,ITSTOP,ISTOR,1,1,1,1,1)


      close(17)
      close(21)
      close(987)
      close(988)

      call CPU_TIME(finish_time)
      write(6,"(A,1pe15.8)")" ELAPSED CPU TIME: ",finish_time-start_time
      write(6,"(A)")" ENDING SIMULATION"

      STOP
      END


      subroutine first_touch()
!
!     touch every field array  to organize data placement
!
#include "hydroparam.h"
#include "globals.h"
!
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        DO K=1,KMAX2
          DO J=1,JMAX2
            starphi(J,K,L)=zero
            tinphi(J,K,L)=zero
            rhotot(J,K,L)=zero       
            S(J,K,L)=zero       
            T(J,K,L)=zero       
            A(J,K,L)=zero       
            U(J,K,L)=zero       
            W(J,K,L)=zero       
            JN(J,K,L)=zero       
            OMEGA(J,K,L)=zero       
            P(J,K,L)=zero
            CV(J,K,L)=zero
            EPS(J,K,L)=zero
            poly_constant(J,K,L)=zero
            RHO(J,K,L)=zero
            PHI(J,K,L)=zero
            QRR(J,K,L)=zero
            QZZ(J,K,L)=zero
            QTT(J,K,L)=zero
            HGAMMA(J,K,L)=zero
            LAMBDA(J,K,L)=zero
            do i=1,4
              TAU(J,K,L,i)=zero
            enddo
            TEMPK(J,K,L)=zero
            DIVFLUX(J,K,L)=zero
            do i=1,3
              RADFLUX(J,K,L,i)=zero
            enddo
#if PASSIVE>0
            do i=1,PAS
              PASSFLUX(J,K,L,I)=zero
            enddo
#endif
            temporary(J,K,L)=zero
            dsdt(J,K,L)=zero
            sfunc(J,K,L)=zero
            l_tau_z(J,K,L)=zero
            dtau_z(J,K,L)=zero
            intensity_in_z(J,K,L)=zero
            intensity_z(J,K,L)=zero
            ddsdtt(J,K,L)=zero
            int_temp(J,K,L)=zero
          ENDDO
        ENDDO
      ENDDO
!$OMP ENDDO nowait

!$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        DO J=1,JMAX2
          TEFFK(J,L)=zero
          TPHK(J,L)=zero
          SURFCGS(J,L)=zero
          init_int_in(J,L)=zero
          KFITA(J,L)=zero
        ENDDO
      ENDDO
!$OMP ENDDO nowait
!$OMP END PARALLEL
      rpstar = zero
      phi_star=zero
      return
      END
