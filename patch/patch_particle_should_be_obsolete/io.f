C***********************************************************************

      SUBROUTINE SETUP(ITSTRT,ITSTOP,IDIAG,ISOADI,ISTOR,ITSTEP,
     &                 ISYM,MAXTRM)
      use particle
      implicit real*8(a-h,o-z)

!!#include "hydroparam.h"
!!#include "globals.h"
!!#include "units.h"

      PARAMETER (JMAX_s=JMAX/2,JMAX_s1=JMAX_s+1,JMAX_s2=JMAX_s+2)
      PARAMETER (KMAX_s=KMAX/2,KMAX_s1=KMAX_s+1,KMAX_s2=KMAX_s+2)

      real*8  ommax,omcen,them
      common  /coefs/coef(pot3jmax2,pot3kmax2,lmax2,2)
      COMMON  /MISC/EPSJR,RHOJR,OMMAX
      integer, PARAMETER::JKM1=2*POT3JMAX+POT3KMAX-1
      integer OMP_GET_MAX_THREADS

      character::date_date*8,date_time*10,date_zone*5

      real*8 limiter

c...arrays for various starts (see itype below).  Eventually,
c...only the arrays in POIS, STATE, and EOM are evolved.

c..axisymmetric start from the hachisu code (256x256) grid
      dimension denny(hj2,hk2),anggy(hj1,hk1)
c..axisymmetric start from formatted files (ie generated from the 2D hydrocode)
      real*8 s_axi(jmax2,kmax2),t_axi(jmax2,kmax2),a_axi(jmax2,kmax2),
     &     rho_axi(jmax2,kmax2),eps_axi(jmax2,kmax2)
c..start for a doubling of the radial grid
      real*8 s_r(jmax_s2,kmax2,lmax), t_r(jmax_s2,kmax2,lmax), 
     &     a_r(jmax_s2,kmax2,lmax), rho_r(jmax_s2,kmax2,lmax),
     &     eps_r(jmax_s2,kmax2,lmax)
c..start for a doubling of the vertical grid
      real*8 s_z(jmax2,kmax_s2,lmax), t_z(jmax2,kmax_s2,lmax), 
     &     a_z(jmax2,kmax_s2,lmax), rho_z(jmax2,kmax_s2,lmax),
     &     eps_z(jmax2,kmax_s2,lmax)
c..start for a doubling of the radial and vertical grids
      real*8 s_rz(jmax_s2,kmax_s2,lmax), t_rz(jmax_s2,kmax_s2,lmax), 
     &     a_rz(jmax_s2,kmax_s2,lmax), rho_rz(jmax_s2,kmax_s2,lmax),
     &     eps_rz(jmax_s2,kmax_s2,lmax)


      CHARACTER resfile*80,index*6
      DATA CURLYR,XMU/83.14,2.0/
 
C     ITYPE TELLS WHETHER INITIAL OR READ IN MODEL.
C
C     ITYPE=1-4 RESTARTED UNFORMATTED  MODELS
C
C        1 = READ FROM UNFORMATTED RESTART FILE
C        2 = READ FROM UNFORMATTED RESTART FILE AND INCREASE R GRID SIZE
C        3 = READ FROM UNFORMATTED RESTART FILE AND INCREASE Z GRID SIZE
C        4 = READ FROM UNFORMATTED RESTART FILE AND INCREASE RZ GRID SIZE
C
C     ITYPE=5-7 INITIAL HACHISU MODEL PLUS PERTURBATION
C
C        5 = READ MODEL CREATED BY HACHISU CODE, NO PERT.
C        6 = READ MODEL CREATED BY HACHISU CODE AND PERTURB IT WITH
C              MACLAURIN BAR MODE 
C        7 = READ MODEL CREATED BY HACHISU CODE AND PERTURB IT WITH
C              RANDOM PERTURBATION.
C
C     ITYPE GE 8 INITIAL MODEL FROM ANOTHER MACHINE PLUS PERTURBATION
C
C        8 = FROM FORMATTED DATA SET GENERATED ON ANOTHER MACHINE
C        9 = FROM AXISYMMETRIC  FORMATTED DATA SET GENERATED ON ANOTHER MACHINE.
C              APPLY NO PERTURBATION.
C       96 = FROM AXISYMMETRIC  FORMATTED DATA SET GENERATED ON ANOTHER MACHINE.
C              PERTURB IT WITH A MACLAURIN BAR MODE.
C       97 = FROM AXISYMMETRIC  FORMATTED DATA SET GENERATED ON ANOTHER MACHINE.
C              PERTURB IT WITH RANDOM PERTURBATION.


   98 FORMAT(1E25.15) 
  100 FORMAT(1P2E15.8)
  102 FORMAT(9(I6,1X))
  103 FORMAT('1THIS WILL BE AN ISOTHERMAL SIMULAT. STARTING AT TIMESTEP
     &NUMBER ',I8,' AND',/,' GOING THROUGH TIMESTEP NUMBER ',I8,'.  FULL
     &DIAGNOSTICS EVERY ',I5,' STEPS.',///)
  104 FORMAT('1THIS WILL BE AN ISENTROPIC SIMULAT. STARTING AT TIMESTEP
     &NUMBER ',I8,' AND',/,' GOING THROUGH TIMESTEP NUMBER ',I8,'.  FULL
     &DIAGNOSTICS EVERY ',I5,' STEPS.',///)
 1044 FORMAT('1THIS WILL BE AN ENERGY EQUATION EVOLUTION, FROM TIMESTEP
     &NUMBER ',I8,' AND',/,' GOING THROUGH TIMESTEP NUMBER ',I8,'.  FULL
     &DIAGNOSTICS EVERY ',I5,' STEPS.',///)


c-------------------------------------------------------------------------------
c  Read run parameter file.
c

      WRITE(6,10000) 'fort.5'
10000 FORMAT('Reading run parameter file ',a)
      OPEN(UNIT=5,FILE='fort.5',STATUS='OLD')

      READ(5,100) KONST,XN
      WRITE(6,10010) KONST,XN
10010 FORMAT('  KONST  ',1PE20.12,/,'  XN     ',1PE20.12)
      READ(5,98)GAMMA
      WRITE(6,10020) GAMMA
10020 FORMAT('  GAMMA  ',1PE20.12)
                  
      READ(5,102)ITSTRT,ITSTOP,IDIAG,ISOADI,ITYPE,NMODL,ISTOR,IGRID
      WRITE(6,10030)ITSTRT,ITSTOP,IDIAG,ISOADI,ITYPE,NMODL,ISTOR,IGRID
     &     ,JMIN
10030 FORMAT('  ITSTRT ',I8,/,'  ITSTOP ',I8,/,'  IDIAG  ',I8,
     &     /,'  ISOADI ',I8,/,'  ITYPE  ',I8,/,'  NMODL  ',I8,
     &     /,'  ISTOR  ',I8,/,'  IGRID  ',I8,/,'  JMIN   ',I8)

      READ(5,100) NPRIME
      WRITE(6,10040) NPRIME
10040 FORMAT('  NPRIME ',1PE20.12)
      NTAPE=7
      CLOSE(5)

 136  FORMAT(///,5X,'PINDEX =',1PE11.3,/,5X,
     &  'NPRIME=',1pe11.3,/,5x,'CON2 =',1PE11.3,/,5X,
     &  'R2 =',1PE11.3,/,5X,'OMCEN =',1PE11.3,/,5X,'DENCEN =',
     &  1PE11.3,/,5X,'TOVERW =',1PE11.3,/,5X,'DELR =',1PE11.3,
     &  5X,'DELZ =',1PE11.3,/,5X,'A1NEWZ =',1PE11.3,/,5X,'JREQ =',I5,/,
     &  5X,'KZPOL =',I5,/,5X,'LMAX=',I5,/,5X,'SOUND SPEED =',1PE11.3,/,
     &  5x,'CS=',1pe11.3,/,5X,'CQ=',1pe11.3,//)

      TMASS=zero
      ENEW=zero
      ELOST=zero
      EDIF=zero
      PHICHK=zero
      KLOCAT=0
      PI=3.1415926535897932d0
      GRAV=one
       
      CVHEAT=CURLYR/(XMU*(GAMMA-1.0))
      DTHETA=two*PI/dble(LMAX)

      
c
c*******************************************************************************
c*******************************************************************************
c
c     READ IN THE MODEL:  DIFFERENT OPTIONS
c
c-------------------------------------------------------------------------------
c  Initial axisymmetric Hachisu model
c
      IF ((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.7)) THEN
         IF(ITYPE.EQ.5) THEN
            WRITE(6,10100) ITYPE, ' no perturbations'
         ELSE IF(ITYPE.EQ.6) THEN
            WRITE(6,10100) ITYPE, ' McClaurin bar mode perturbation'
         ELSE
            WRITE(6,10100) ITYPE, ' random perturbation'
         ENDIF
10100    FORMAT(/,'Reading model type',I4,':  Axisymmetric Hachisu,',a)

         OPEN(UNIT=2,FILE='fort.2',STATUS='OLD')    

         READ(2,1685)PINDEX,CON2,RRR2,OMCEN,DENCEN,TOVERW,
     &        ROF3N,ZOF3N,A1NEWZ,JREQ,KZPOL

 1685    FORMAT(3X,1PE22.15,2X,8(1PE22.15,2X),2I4)

         write(*,1685)PINDEX,CON2,RRR2,OMCEN,DENCEN,TOVERW,
     &        ROF3N,ZOF3N,A1NEWZ,JREQ,KZPOL
         
         READ(2,1617) DENNY
         READ(2,1617) ANGGY
 1617    FORMAT(8(1PE22.15,2X))
         CLOSE(2)

         den=dencen
         rholmt=dencen*gridlim

         rholmt_p=rholmt*dust_to_gas
         epslmt=(1.d0/(gamma-1.0))*rholmt*tbgrnd/tconv*bkmpcode
         dumlmt=epslmt**(1.d0/gamma)
         OMMAX=OMCEN
         sound=SQRT(gamma*(DENNY(2,2)**(gamma-1.d0)))

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(LP,j,k,l)                        &
!$OMP&  SHARED(gamma,konst,den)
!$OMP DO SCHEDULE(STATIC)
         do l=1,lmax
            do k=2,kmax1
               do j=2,jmax1
                  rho(j,k,l)=denny(j,k)
                  jn(j,k,L)=anggy(j,k)
                  if(rho(j,k,l).lt.1.e-10*den) then
                    rho(j,k,l)=1.e-10*den
                    jn(j,k,L)=0.d0
                  endif
                  s(j,k,l)=0.d0
                  t(j,k,l)=0.d0
                  a(J,K,L)=rho(j,k,l)*jn(j,k,l)
                  u(j,k,l)=0.d0
                  w(j,k,l)=0.d0
                  p(j,k,l)=konst*rho(j,k,l)**gamma
                  eps(j,k,l)=p(j,k,l)/(gamma-1.d0)
               end do
            end do
         end do
!$OMP END DO NOWAIT

C....Set quantities on outside boundary.
!$OMP DO SCHEDULE(STATIC)
         do l=1,lmax
           do k=2,kmax1
             jn(jmax2,k,l)  = jn(jmax1,k,l)
             s(jmax2,k,l)   = s(jmax1,k,l)
             t(jmax2,k,l)   = t(jmax1,k,l)
             a(jmax2,k,l)   = a(jmax1,k,l)
             u(jmax2,k,l)   = u(jmax1,k,l)
c            omega(jmax2,k,l)= omega(jmax1,k,l)
             w(jmax2,k,l)   = w(jmax1,k,l)
             rho(jmax2,k,l) = rho(jmax1,k,l)
             p(jmax2,k,l)   = p(jmax1,k,l)
             eps(jmax2,k,l) = eps(jmax1,k,l)
           enddo
         enddo
!$OMP ENDDO NOWAIT
!$OMP DO SCHEDULE(STATIC)
         do l=1,lmax
           do j=2,jmax2
             jn(j,kmax2,l)  = jn(j,kmax1,l)
             s(j,kmax2,l)   = p(j,kmax1,l)
             t(j,kmax2,l)   = t(j,kmax1,l)
             a(j,kmax2,l)   = a(j,kmax1,l)
             u(j,kmax2,l)   = u(j,kmax1,l)
c            omega(j,kmax2,l)= omega(j,kmax1,l)
             w(j,kmax2,l)   = w(j,kmax1,l)
             rho(j,kmax2,l) = rho(j,kmax1,l)
             p(j,kmax2,l)   = p(j,kmax1,l)
             eps(j,kmax2,l) = eps(j,kmax1,l)
           enddo
         enddo
!$OMP ENDDO

C....Set quantities around z-axis.
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            LP=1+MOD(L-1+LMAX/2,LMAX)
            DO K=2,KMAX2
               JN(1,K,L)  = JN(2,K,LP)
               S(1,K,L)   = -S(3,K,LP)
               T(1,K,L)   = T(2,K,LP)
               A(1,K,L)   = A(2,K,LP)
               U(1,K,L)   = -U(3,K,LP)
c              OMEGA(1,K,L) = OMEGA(2,K,LP)
               W(1,K,L)   = W(2,K,LP)
               RHO(1,K,L) = RHO(2,K,LP)
               P(1,K,L)   = P(2,K,LP)
               EPS(1,K,L) = EPS(2,K,LP)
            ENDDO
         ENDDO
!$OMP ENDDO NOWAIT

C....Set quantities below the equatorial plane.
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            DO J=1,JMAX2
               JN(J,1,L)  = JN(J,2,L)
               S(J,1,L)   = S(J,2,L)
               T(J,1,L)   = -T(J,3,L)
               A(J,1,L)   = A(J,2,L)
               U(J,1,L)   = U(J,2,L)
               W(J,1,L)   = -W(J,3,L)
c              OMEGA(J,1,L) = OMEGA(J,2,L)
               RHO(J,1,L) = RHO(J,2,L)
               P(J,1,L)   = P(J,2,L)
               EPS(J,1,L) = EPS(J,2,L)
            ENDDO
         ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      END IF


C
C-------------------------------------------------------------------------------
C  Restarted model
C
      IF((ITYPE.LE.4).OR.(ITYPE.GE.8)) THEN


c...  Standard read in .........................................................
10200    FORMAT(/,'Reading model type',I4,
     &          ':  Unformatted Restart File',a)


         IF (ITYPE.EQ.99) THEN

            OPEN(UNIT=7,FILE='fort.7',FORM='UNFORMATTED',STATUS='OLD')
            WRITE(6,10200) ITYPE,'.'
            read(7) S(1:JMAX2,1:KMAX2,1:8)
            read(7) T(1:JMAX2,1:KMAX2,1:8)
            read(7) A(1:JMAX2,1:KMAX2,1:8)
            read(7) RHO(1:JMAX2,1:KMAX2,1:8)
            read(7) EPS(1:JMAX2,1:KMAX2,1:8)
            read(7)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,
     &        JREQ,OMMAX
            read(7,IOSTAT=ios) tmassini,tmass,tmassadd,
     &         tmassout,tmassacc,totcool,totdflux,totheat,totirr,etotfl,
     &         eflufftot  !ACB

            if (ios /= 1) then 
               print *, "Last set of data missing. Check input."
            endif

            tmassini = tmass
            tmassadd = zero
            tmassout = zero
            tmassacc = zero
            totcool  = zero
            totdflux = zero
            totheat  = zero
            totirr   = zero
            etotfl   = zero
            eflufftot= zero
            time     = zero

            sound=half

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(J,K,L)
      do L = 1, LMAX
       do K = 1, KMAX2
        do J = 1, JMAX2
          s(J,K,L) = s(J,K,1)            
          a(J,K,L) = a(J,K,1)            
          t(J,K,L) = t(J,K,1)            
          rho(J,K,L) = rho(J,K,1)            
          eps(J,K,L) = eps(J,K,1)            
        enddo
       enddo
      enddo
!$OMP END PARALLEL DO

            ommax=omega(jmin,2,1)
            ommax=max(1.d-14,ommax)     !dkb - temp fix for underflow problem
            dencen=den
            rholmt=dencen*gridlim
            epslmt=(1.d0/(gamma-1.0))*rholmt**gamma*gridlim
            rholmt_p=rholmt*dust_to_gas
            CLOSE(7)

         END IF

         IF (ITYPE.EQ.1.OR.ITYPE.EQ.98) THEN
            
            OPEN(UNIT=7,FILE='fort.7',FORM='UNFORMATTED',STATUS='OLD')
            WRITE(6,10200) ITYPE,'.'
            read(7) S
            read(7) T
            read(7) A
            read(7) RHO
            read(7) EPS
            read(7)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,
     &        JREQ,OMMAX
            read(7,IOSTAT=ios) tmassini,tmass,tmassadd,
     &         tmassout,tmassacc,totcool,totdflux,totheat,totirr,etotfl,
     &         eflufftot  !ACB

            if (ios /= 1) then 
               print *, "Last set of data missing. Check input."
            endif

            sound=sqrt(1.5d0/ten)

            ommax=omega(jmin,2,1)
            ommax=max(1.d-14,ommax)     !dkb - temp fix for underflow problem
            dencen=den
            rholmt=dencen*gridlim
            epslmt=(1.d0/(gamma-1.0))*rholmt**gamma*gridlim
            CLOSE(7)

         END IF


c...  Read in data, then increase the radial grids by 2x .......................

         IF (ITYPE.EQ.2) THEN 

            OPEN(UNIT=7,FILE='fort.7',FORM='UNFORMATTED',STATUS='OLD')
            WRITE(6,10200) ITYPE,'. Increasing R-grid size.'
            read(7) S_R
            read(7) T_R
            read(7) A_R
            read(7) RHO_R
            read(7) EPS_R
            read(7)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,
     &        JREQ,OMMAX
            read(7,IOSTAT=ios) tmassini,tmass,tmassadd, 
     &         tmassout,tmassacc,totcool,totdflux,totheat,totirr,etotfl, 
     &         eflufftot  !ACB 
 
            if (ios /= 1) then  
               print *, "Last set of data missing. Check input." 
            endif 

            dencen=den
            rholmt=dencen*gridlim
            epslmt=(1.d0/(gamma-1.0))*rholmt**gamma*gridlim
            ommax=7.483E-02            

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k,l)                        &
!$OMP&  SHARED(zof3n,dtheta,epslmt,rholmt) REDUCTION(+:tmassadd)
            do l=1,lmax
               do k=1,kmax2
                  do j=1,jmax_s
                     rho(j,k,l)=rho_r(j,k,l)
                     s(j,k,l)=s_r(j,k,l)
                     t(j,k,l)=t_r(j,k,l)
                     a(j,k,l)=a_r(j,k,l)
                     eps(j,k,l)=eps_r(j,k,l)
                  end do
                  do j=jmax_s1,jmax2
                     rho(j,k,l)=rholmt
                     eps(j,k,l)=epslmt
                     s(j,k,l)=0.d0
                     t(j,k,l)=0.d0
                     a(j,k,l)=0.d0
                     if ((j.le.jmax).and.(k.le.kmax).and.(k.ge.2)) then
                        tmassadd=tmassadd+rholmt*dtheta*zof3n*(r(j+1)**2
     &                       -r(j)**2)
                     endif
                  end do
               end do
            end do
!$OMP END PARALLEL DO

            CLOSE(7)

         END IF

c...  Read in data, then increase vertical grids by x2 .........................

         IF (ITYPE.EQ.3) THEN 

            OPEN(UNIT=7,FILE='fort.7',FORM='UNFORMATTED',STATUS='OLD')
            WRITE(6,10200) ITYPE,'. Increasing Z-grid size.'
            read(7) S_Z
            read(7) T_Z
            read(7) A_Z
            read(7) RHO_Z
            read(7) EPS_Z
            read(7)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,
     &        JREQ,OMMAX
            read(7,IOSTAT=ios) tmassini,tmass,tmassadd, 
     &         tmassout,tmassacc,totcool,totdflux,totheat,totirr,etotfl, 
     &         eflufftot  !ACB 
 
            if (ios /= 1) then  
               print *, "Last set of data missing. Check input." 
            endif 
       
            dencen=den
            rholmt=dencen*gridlim
            epslmt=(1.d0/(gamma-1.0))*rholmt**gamma*gridlim

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,l)                        &
!$OMP&  SHARED(zof3n,dtheta,epslmt,rholmt) REDUCTION(+:tmassadd)
            do l=1,lmax
               do k=1,kmax_s                  
                 do j=1,jmax2
                     rho(j,k,l)=rho_z(j,k,l)
                     s(j,k,l)=s_z(j,k,l)
                     t(j,k,l)=t_z(j,k,l)
                     a(j,k,l)=a_z(j,k,l)
                     eps(j,k,l)=eps_z(j,k,l)
                  end do
                enddo
                do k=kmax_s1,kmax2
                  do j=1,JMAX2
                     rho(j,k,l)=rholmt
                     eps(j,k,l)=epslmt
                     s(j,k,l)=zero
                     t(j,k,l)=zero
                     a(j,k,l)=zero
                     if ((j.le.jmax).and.(j.ge.jmin).and.(k.le.kmax))
     &                    then
                        tmassadd=tmassadd+rholmt*dtheta*zof3n*(r(j+1)**2
     &                       -r(j)**2)
                     endif
                  end do
               end do
            end do
!$OMP END PARALLEL DO

            CLOSE(7)

         END IF

c...  Read in data, then increase both radial and vertical grids by 2x .........

         IF (ITYPE.EQ.4) THEN 

            WRITE(6,10200) ITYPE,'. Increasing RZ-grid size.'
            read(7) S_RZ
            read(7) T_RZ
            read(7) A_RZ
            read(7) RHO_RZ
            read(7) EPS_RZ
            read(7)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,
     &        JREQ,OMMAX

            read(7,IOSTAT=ios) tmassini,tmass,tmassadd, 
     &         tmassout,tmassacc,totcool,totdflux,totheat,totirr,etotfl, 
     &         eflufftot  !ACB 
 
            if (ios /= 1) then  
               print *, "Last set of data missing. Check input." 
            endif 
       
            ommax=omega(jmin,2,1)
            ommax=max(1.d-14,ommax)
            dencen=den
            rholmt=dencen*gridlim
            epslmt=(1.d0/(gamma-1.0))*rholmt**gamma*gridlim

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,l)                        &
!$OMP&  SHARED(zof3n,dtheta,epslmt,rholmt) REDUCTION(+:tmassadd)
            do l=1,lmax
               do k=1,kmax_s                  
                 do j=1,jmax_s
                     rho(j,k,l)=rho_rz(j,k,l)
                     s(j,k,l)=s_rz(j,k,l)
                     t(j,k,l)=t_rz(j,k,l)
                     a(j,k,l)=a_rz(j,k,l)
                     eps(j,k,l)=eps_rz(j,k,l)
                  enddo
                 enddo
                 do k=kmax_s1,kmax2
                   do J=1,jmax_s 
                     rho(j,k,l)=rholmt
                     eps(j,k,l)=epslmt
                     s(j,k,l)=zero
                     t(j,k,l)=zero
                     a(j,k,l)=zero
                     if ((j.ge.jmin).and.(k.le.kmax)) then
                        tmassadd=tmassadd+rholmt*dtheta*zof3n*(r(j+1)
     &                       **2-r(j)**2)
                     endif   
                  enddo
               enddo
               do k=1,kmax2
                 do j=jmax_s1,jmax2
                     rho(j,k,l)=rholmt
                     eps(j,k,l)=epslmt
                     s(j,k,l)=zero
                     t(j,k,l)=zero
                     a(j,k,l)=zero
                     if ((j.le.jmax).and.(k.le.kmax).and.(k.ge.2)) then
                        tmassadd=tmassadd+rholmt*dtheta*zof3n*(r(j+1)**2
     &                       -r(j)**2)
                     endif
                  end do
               end do
            end do
!$OMP END PARALLEL DO

            CLOSE(7)

         END IF

C...............................................................................
C...  Continue a job from cthulhu (or somewhere else) on E10K

         IF (ITYPE.EQ.8) THEN 
            
            open(unit=99,file='cthustart.dat')

            READ(99,666)S
            READ(99,666)T
            READ(99,666)A
            READ(99,666)RHO
            READ(99,666)EPS
            READ(99,667)ROF3N
            read(99,667)ZOF3N
            read(99,667)DELT
            read(99,667)TIME
            read(99,667)ELOST
            read(99,667)den
            read(99,667)sound
            read(99,667)ommax
            read(99,668)jreq
            read(99,668)kzpol

            dencen=den
            rholmt=dencen*gridlim
            epslmt=(1.d0/(gamma-1.0))*rholmt**gamma*gridlim


 666        format(1pe22.15)
 667        format(1pe22.15)
 668        format(i4)
            close(99)

         END IF

C...............................................................................
c...  Continue an axisymmetric job on the E10K in 3D.  This is useful
c...  for evolving the output from the 2D hydro code
         
         IF (ITYPE.GE.9.AND.ITYPE<98) THEN 
            
            open(unit=99,file='cthustart_axi.dat')
            
            READ(99,666)S_AXI
            READ(99,666)T_AXI
            READ(99,666)A_AXI
            READ(99,666)RHO_AXI
            READ(99,666)EPS_AXI
            READ(99,667)ROF3N
            read(99,667)ZOF3N
            read(99,667)DELT
            read(99,667)TIME
            read(99,667)ELOST
            read(99,667)den
            read(99,667)sound
            read(99,667)ommax
            read(99,668)jreq
            read(99,668)kzpol
            if (jmin.gt.2) then
            read(99,667)tmassini
            read(99,667)tmass
            read(99,667)tmassadd
            read(99,667)tmassout
            read(99,667)tmassacc
            endif
            
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,l)                        &
            do j=1,jmax2
               do k=1,kmax2
                  do L=1,LMAX
                     s(j,k,l)=s_AXI(j,k)
                     t(j,k,l)=t_AXI(j,k)
                     a(j,k,l)=a_AXI(j,k)
                     rho(j,k,l)=rho_AXI(j,k)
                     eps(j,k,l)=eps_AXI(j,k)
                  end do
               end do
            end do
!$OMP END PARALLEL DO

            dencen=den
            rholmt=dencen*gridlim
            epslmt=(1.d0/(gamma-1.0))*rholmt**gamma*gridlim
            totcool=0.0
            totdflux=0.0
            totheat=0.0
            totirr=0.0
            etotfl=0.d0
            eflufftot=0.d0
            close(99)
      
         END IF


c...  Verfication we read the restarted model ok ...............................

      call DATE_AND_TIME(date_date,date_time,date_zone)
      write(6,"(A,3(1X,A))")"DATE CCYYMMDD TIME HHMMSS.SSS ZONE +-HHMM U&
     &TC: ",date_date,date_time,date_zone
      write(6,*)
         
      write(6,114) ' ROF3N   ',rof3n
      write(6,114) ' ZOF3N   ',zof3n
      write(6,114) ' DELT    ',delt
      write(6,114) ' TIME    ',time
      write(6,114) ' DEN     ',den
      write(6,114) ' SOUND   ',sound
      write(6,114) ' RHOLMT  ',rholmt
      write(6,114) ' EPSLMT  ',epslmt
      write(6,114) ' TBGRND  ',tbgrnd
      write(6,114) ' TENVK   ',tenvk
      write(6,114) ' XABUN   ',xabun
      write(6,114) ' YABUN   ',yabun
      write(6,114) ' ZABUN   ',zabun
      write(6,114) ' OMMAX   ',ommax
      write(6,117) ' JREQ    ',jreq
      if (external_pot) then
        write(6,114) ' TMASSINI',tmassini
        write(6,114) ' TMASS   ',tmass
        write(6,114) ' TMASSADD',tmassadd
        write(6,114) ' TMASSOUT',tmassout
        write(6,114) ' TMASSACC',tmassacc
        write(6,114) ' TOTCOOL ',totcool
        write(6,114) ' TOTDFLUX',totdflux
        write(6,114) ' TOTHEAT ',totheat
        write(6,114) ' TOTIRR  ',totirr
      endif
      write(6,*)
      select case(H2STAT)
      case(1)
       write(6,"( ' H2 TREATMENT: PURE ORTHO')")
       write(6,"( ' Do you know what you are doing?')")
      case(2)
       write(6,"( ' H2 TREATMENT: SET MIXTURE')")
       write(6,114) ' ORTHO->   ',bc
       write(6,114) ' PARA ->   ',ac
      case(3)
       write(6,"( ' H2 TREATMENT: EQUILIBRIUM')")
      case(4)
       write(6,"( ' H2 TREATMENT: SINGLE GAMMA')")
       write(6,114) " GAMMA ->", gamma
      end select
      write(6,*)
      if(cq.eq.0.)write(6,21)
      if(cq.ne.0.)write(6,210)CQ
      if(cs.eq.0.)write(6,22)
      if(cs.ne.0.)write(6,220)CS
      if(ictype.eq.0) write(6,"(' No Cooling')")
      if(ictype.eq.1) write(6,23) cct
      if(ictype.eq.2) write(6,230)
      if(ictype.eq.3) write(6,231)
      if(ictype.eq.4) write(6,232)
      if(ictype.eq.5) write(6,"(' Tcool Omega = ',1pe15.8)")cct
      if(ictype.eq.69) write(6,233)
      if(fgsoft.gt.0d0)write(6,"(' Softened Potential. FGSOFT= ',1pe10.3&
     &                         )")fgsoft
      if(fgsoft.eq.0d0)write(6,"(' NO Potential Softening.')")
      write(6,'(" EXTERNAL POTENTIAL  ->",L)')external_pot
#if WIGGLE>0
      write(6,'(" INTEGRATE STAR      ->",L)').true.
      write(6,'(" RESTART STAR        ->",L)')restart_wiggle
#else 
      write(6,'(" INTEGRATE STAR      ->",L)').false.
      write(6,'(" RESTART STAR        ->",L)').false.
#endif
#if FLUID>0
      write(6,'(" TRACE FLUID ELEMENTS->",L)').true.
      write(6,'(" RESTART FLUID ELE.  ->",L)')restart_fluid
#else
      write(6,'(" TRACE FLUID ELEMENTS->",L)').false.
      write(6,'(" RESTART FLUID ELE.  ->",L)').false.
#endif
      write(6,*)
 114     format(a,1pe15.8)
 115     format(a,i7)
         
      END IF



c
c*******************************************************************************
c*******************************************************************************
c
c     MODEL READ.  Now setup grid, choose perturbations, find phi, read
c     opacities, etc
c
      
      write(index,'(i6.6)')itstop
      resfile='rslts.'//index
      OPEN(unit=3,file=resfile)

      write(3,"(A,3(1X,A))")"DATE CCYYMMDD TIME HHMMSS.SSS ZONE +-HHMM U&
     &TC: ",date_date,date_time,date_zone
!$    write(3,"(A22,1X,I4)")'OMP enabled. Nproc=',OMP_GET_MAX_THREADS()

      IF(ISOADI.EQ.1) WRITE(3,103) ITSTRT,ITSTOP,IDIAG
      IF(ISOADI.EQ.2) WRITE(3,104) ITSTRT,ITSTOP,IDIAG
      IF(ISOADI.EQ.3) WRITE(3,1044) ITSTRT,ITSTOP,IDIAG
c         WRITE(3,116) XN,GAMMA,KONST,NPRIME
c 116     FORMAT(//,10X,'ROTATING POLYTROPE...  N =',0PF5.2,/,27X,
c     &        '1 + 1/N =',1PE10.3,/,33X,'K =',1PE10.3,/,28X,
c     &        'NPRIME =',1pe10.3,/)

 116  format(a,1pe10.3)
 117  format(a,i4)

c
c...verification on the rslts file
c
      write(3,*)
      write(3,116) ' ROF3N   ',rof3n
      write(3,116) ' ZOF3N   ',zof3n
      write(3,116) ' DELT    ',delt
      write(3,116) ' TIME    ',time
      write(3,116) ' DEN     ',den
      write(3,116) ' SOUND   ',sound
      write(3,116) ' RHOLMT  ',rholmt
      write(3,116) ' EPSLMT  ',epslmt
      write(3,116) ' TBGRND  ',tbgrnd
      write(3,116) ' TENVK   ',tenvk  
      write(3,116) ' XABUN   ',xabun  
      write(3,116) ' YABUN   ',yabun  
      write(3,116) ' ZABUN   ',zabun  
      write(3,116) ' OMMAX   ',ommax
      write(3,117) ' JREQ    ',jreq
      if (external_pot) then
        write(3,116) ' TMASSINI',tmassini
        write(3,116) ' TMASS   ',tmass
        write(3,116) ' TMASSADD',tmassadd
        write(3,116) ' TMASSOUT',tmassout
        write(3,116) ' TMASSACC',tmassacc
        write(3,116) ' TOTCOOL ',totcool
        write(3,116) ' TOTDFLUX',totdflux
        write(3,116) ' TOTHEAT ',totheat
        write(3,116) ' TOTIRR  ',totirr
      endif
      write(3,*)
      select case(H2STAT)
      case(1)
       write(3,"( ' H2 TREATMENT: PURE ORTHO')")
       write(3,"( ' Do you know what you are doing?')")
      case(2)
       write(3,"( ' H2 TREATMENT: SET MIXTURE')")
       write(3,116) ' ORTHO->   ',bc
       write(3,116) ' PARA ->   ',ac
      case(3)
       write(3,"( ' H2 TREATMENT: EQUILIBRIUM')")
      case(4)
       write(3,"( ' H2 TREATMENT: SINGLE GAMMA')")
       write(3,116) " GAMMA ->", gamma
      end select
      write(3,*)

      if(cq.eq.0.)write(3,21)
      if(cq.ne.0.)write(3,210)CQ
      if(cs.eq.0.)write(3,22)
      if(cs.ne.0.)write(3,220)CS
      if(ictype.eq.0) write(3,"(' No Cooling')")
      if(ictype.eq.1) write(3,23) cct
      if(ictype.eq.2) write(3,230)
      if(ictype.eq.3) write(3,231)
      if(ictype.eq.4) write(3,232)
      if(ictype.eq.5) write(3,"(' Tcool Omega = ',1pe15.8)")cct
      if(ictype.eq.69) write(3,233)
      if(fgsoft.gt.0d0)write(3,"(' Softened Potential. FGSOFT= ',1pe10.3&
     &                         )")fgsoft
      if(fgsoft.eq.0d0)write(3,"(' NO Potential Softening.')")
      write(3,'(" EXTERNAL POTENTIAL  ->",L)')external_pot
#if WIGGLE>0
      write(3,'(" INTEGRATE STAR      ->",L)').true.
      write(3,'(" RESTART STAR        ->",L)')restart_wiggle
#else 
      write(3,'(" INTEGRATE STAR      ->",L)').false.
      write(3,'(" RESTART STAR        ->",L)').false.
#endif
#if FLUID>0
      write(3,'(" TRACE FLUID ELEMENTS->",L)').true.
      write(3,'(" RESTART FLUID ELE.  ->",L)')restart_fluid
#else
      write(3,'(" TRACE FLUID ELEMENTS->",L)').false.
      write(3,'(" RESTART FLUID ELE.  ->",L)').false.
#endif
 21   format(' Artificial Viscosity is OFF.')
 210  format(' Artificial Viscosity is ON, CQ=',1pe10.3)
 22   format(' Shear Viscosity is OFF.')
 220  format(' Shear Viscosity is ON, CS=',1pe10.3)
 23   format(' Constant Cooling Time =', 1pe10.3)
 230  format(' Eddington Atmosphere')
 231  format(' Flux Limited Diffusion')
 232  format(' Flux Limited Diffusion and Atmosphere')  
 233  format(' Vertical Rays and Diffusion')
      write(3,*)

c
c...grid setup
c
      DELR=ROF3N
      DELZ=ZOF3N

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,j)                             &
!$OMP&  SHARED(DELZ,DELR)
!$OMP DO SCHEDULE(STATIC)
      DO 448 J=1,JMAX2
         R(J)=(J-2)*DELR
         RHF(J)=R(J)+DELR/2.0
  448 CONTINUE
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
      DO 450 K=1,KMAX2
         Z(K)=(K-2)*DELZ
         ZHF(K)=Z(K)+DELZ/2.0
  450 CONTINUE
!$OMP END DO nowait
!$OMP END PARALLEL

C...............................................................................
c...  Get a few conversions done. These will be used in the source
c routine.

      select case(SETUNITS)

      case (0)

      Msyscgs=Mstar*Msuncgs*(1.0+(tmassini/(1.0-tmassini)))
      mass_star = 1.d0-tmassini
      
      PKcgs = ( (Rdiskau*AUcgs/r(jreq))**(3.0-xn) * 
     &     Gcgs**xn * Msyscgs**(xn-1.0) )**(1.0/xn)
      
      sigma = (sigmacgs/(bkmpcgs**4)) * 
     &     ( Gcgs**(15.0/2.0 - 3.0*xn) * 
     &     Msyscgs**(5.0 - 2.0*xn) * PKcgs**(xn/2.0)
     &     )**(1.0/(3.0-xn)) 
      
      Tconv = (1.0/bkmpcgs) *( Gcgs**3 * Msyscgs**2 * PKcgs**(-xn)
     &     )**(1.0/(3.0-xn))
      
      Sconv = (Gcgs**(2.0*xn) * Msyscgs**(1.0+xn) * PKcgs**(-2
     &     .0*xn))**(1.0/(3.0-xn))
      
      Pconv = (Gcgs**(3.0+3.0*xn) * Msyscgs**(2.0+2.0*xn) * PKcgs**(-4
     &     .0*xn))**(1.0/(3.0-xn))

      Dconv = (Gcgs**(-xn) * Msyscgs**(1.0-xn) * PKcgs**(xn))**(1.0/(3.0
     &     -xn))

      rhoconv = (Gcgs**3*Msyscgs**2*PKcgs**(-3))
     &        **(1.d0/(3.d0*gamma-4.d0))

      engconv = pconv/rhoconv

      bkmpcode = 1.d0

      case(1)

      mass_star = mstar
      msyscgs   = (mstar + tmassini)*msuncgs
      bkmpcode = bkmpcgs*aucgs/(gcgs*msuncgs)
      rhoconv = msuncgs/aucgs**3
      sconv = rhoconv*aucgs
      sigma = sigmacgs*aucgs**4.5d0/(msuncgs**2.5d0*gcgs**1.5d0)
      tconv = 1.d0
      pconv = bkmpcgs*rhoconv/bkmpcode
      engconv = bkmpcgs/bkmpcode

      end select

      epslmt=(1.d0/(gamma-1.0))*rholmt*tbgrnd/tconv*bkmpcode
      dumlmt=epslmt**(1.d0/gamma)
      rholmt_p=rholmt*dust_to_gas
      rhoa=rhoacgs/rhoconv

      print "(A,1X,1pe15.8)", ' Msyscgs ', msyscgs
      print "(A,1X,1pe15.8)", ' PKcgs ', pkcgs  
      print "(A,1X,1pe15.8)", ' sigma ', sigma  
      print "(A,1X,1pe15.8)", ' Tconv ', tconv  
      print "(A,1X,1pe15.8)", ' Sconv ', sconv  
      print "(A,1X,1pe15.8)", ' Pconv ', pconv  
      print "(A,1X,1pe15.8)", ' rhoconv ', rhoconv
      print "(A,1X,1pe15.8)", ' engconv ', engconv
      print "(A,1X,1pe15.8)", ' Time convert (yr) ',
     &                        1.d0/(sqrt(gcgs*rhoconv)*3.155760d7)

      write(3,114)' Star Mass Code ', mass_star
      write(3,114)' Msyscgs ', msyscgs
      write(3,114) ' bkmpcode ', bkmpcode
      write(3,114) ' PKcgs ', pkcgs
      write(3,114) ' sigma ', sigma
      write(3,114) ' Tconv ', tconv
      write(3,114) ' Sconv ', sconv
      write(3,114) ' Pconv ', pconv
      write(3,114) ' rhoconv ', rhoconv
      write(3,114) ' engconv ', engconv
      write(3,114) ' Time convert (yr) ',
     &                        1.d0/(sqrt(gcgs*rhoconv)*3.155760d7)
      write(3,*)

      call Initengtable()
      epslmt=engtable(1)/engconv*rholmt
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC)
      do L=1,LMAX
       do K=1,KMAX2
        do J=1,JMAX2
          poly_constant(J,K,L)=(poly_factor*rof3n)**2/pi
!          poly_constant(J,K,L)=(poly_factor*
!     &      max(rhf(J)*dtheta,rof3n,zof3n))**2/pi
!!          eps(J,K,L)=eps(J,K,L)+poly_constant(J,K,L)*rho(J,K,L)**2
        enddo
       enddo
      enddo
!$OMP ENDDO
      call State()
!$OMP END PARALLEL

#if ROTATING_SHIFT>0
      sound=max(sound,omega_frame*rhf(JMAX2)*1.5)
!$OMP SINGLE 
      print *, "*********************************************"
      print *, " POSSIBLE sound ADJUST DUE TO ROTATING FRAME."
      print *, " COMPARE WITH PREVIOUS OUTPUT. sound = ", sound
      print *, "*********************************************"
      write(3,*)"*********************************************"
      write(3,*)" POSSIBLE sound ADJUST DUE TO ROTATING FRAME."
      write(3,*)" COMPARE WITH PREVIOUS OUTPUT. sound = ", sound
      write(3,*)"*********************************************"
!$OMP END SINGLE NO WAIT
#endif


C...............................................................................
c...  The next few lines open the mean molecular weight and opacities files
c...  for mean rosseland, mean planck, planck weighted by T* absoption
c...  only and absorption plus scattering.


      print*, ' READING MOLECULAR WEIGHTS AND OPACITIES'
      write(3,"(a)")" Reading opacities and MMW from the following:"
      write(3,"(a,a)")" MMW   -> ",mmwfile
      write(3,"(a,a)")" ROSS  -> ",rosstablefile
      write(3,"(a,a)")" PLANCK-> ",plancktablefile
      write(3,"(a,a)")" IRR   -> ", irrtablefile
      write(3,"(a,a)")" SCA   -> ",scatablefile
      write(3,"(' ')")
      open(unit=70,file=trim(mmwfile),status='old',form
     &     ='formatted')
      do i=1,itable
         do n=1,itable
            read(70,*) XMMWtable(i,n)
         enddo
      enddo
      close(70)
      
      open(unit=71,file=trim(rosstablefile),status='old'
     &     ,form='formatted')
      do i=1,itable
         do n=1,itable
            read(71,*) ROStable(i,n)
         enddo
      enddo
      close(71)
      
      open(unit=72,file=trim(plancktablefile), status='old',
     &     form='formatted')
      do i=1,itable
         do n=1,itable
            read(72,*) PLKtable(i,n)
         enddo
      enddo
      close(72)
      
      open(unit=73,file=trim(irrtablefile), status
     &     ='old',form='formatted')
      do i=1,itable
         do n=1,itable
            read(73,*) IRRtable(i,n)
         enddo
      enddo
      close(73)

      open(unit=74,file=trim(scatablefile),
     &     status='old',form='formatted')
      do i=1,itable
         do n=1,itable
            read(74,*) SCAtable(i,n)
         enddo
      enddo
      close(74)

      Pmin=-12.0
      Pmax=9.0
      Tmin=0.5
      Tmax=7.0
      dP=(Pmax-Pmin)/(itable-1)
      dT=(Tmax-Tmin)/(itable-1)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ip,it)                           &
!$OMP&  SHARED(dP,Pmin,dT,Tmin)
!$OMP DO SCHEDULE(STATIC)
      do it=1,itable
         Ttab(it)=Tmin+(it-1)*dT
      enddo
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
      do ip=1,itable
         Ptab(ip)=Pmin+(ip-1)*dP
      enddo
!$OMP END DO nowait
!$OMP END PARALLEL

C 
C....FOR RESTARTED MODELS; Not for hachisu models!  Note the EOS is different
C
      IF ((ITYPE.LE.4).OR.(ITYPE.GE.8)) THEN

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k,l)                           &
!$OMP&  SHARED(zof3n,dtheta,rholmt,epslmt,CVHEAT) REDUCTION(+:tmassadd)
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            DO K=1,KMAX2
               DO J = 1, JMAX2
                  CV(J,K,L)=CVHEAT
                  EPS(j,k,l)=MAX(EPSLMT,EPS(j,k,l))
                  if ((rho(j,k,l).lt.rholmt).and.((j.le.jmax).and.(j.ge
     &                 .jmin)).and.((k.le.kmax).and.(k.ge.2))) then 
                     tmassadd=tmassadd+(rholmt-rho(j,k,l))
     &                    *dtheta*zof3n*(r(j+1)**2-r(j)**2)
                  endif   
                  rho(j,k,l)=MAX(RHOLMT,RHO(j,k,l))
               enddo
            enddo
         enddo
!$OMP END DO
         call State()
c
c...find u,w,omega, and limit them
c

         CALL VELOCITY
         CALL VLIMIT
!$OMP END PARALLEL
      END IF


C...FIND POTENTIALS FOR RESTARTED EVOLUTION or axisymmetric

c
c...Initial grav potential. Only restarted or unperturbed models here.
c...Perturbed models wait until after perturbation.
c
      IF((ITYPE.LE.5).OR.(ITYPE.EQ.8).OR.(ITYPE.EQ.9)) THEN
         IPRINT=0
         REDGE=0.d0
         gsoft = fgsoft*rof3n
         if (gsoft>0d0) tmassacc=0d0

        if(read_particle_file)then
         call read_particles(P_FILEID)
        else
         call initialize_particles()
        endif
        call set_particle_density()


#if SELFGRAVITY>0
!$OMP PARALLEL DEFAULT(SHARED)                                          &
!$OMP&  SHARED(JKMAX,ISYM)
         CALL SETBDY(0,ISYM)
!$OMP END PARALLEL
         CALL BDYGEN(MAXTRM,ISYM,REDGE)
         CALL POT3(8,IPRINT)
#endif

#if WIGGLE>0
         call delta(.false.)
         call wiggle(rhotot,rhf,zhf,vx,vy,fx,fy,delt,pi,rpstar,phi_star,
     &        rof3n,zof3n,dtheta,gsoft,JMAX,KMAX,LMAX,JMIN,.true.,
     &        restart_wiggle)
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
#if EXTERNAL_POT>0
      call ExternalPotInit()
#endif

!$OMP END PARALLEL

      if(.not.read_particle_file)call set_particle_vel()
      END IF
 
C...Smooth the edge of the initial model.


      limiter = den*phylim 
      rholmt_p=rholmt*dust_to_gas

   
      DO 313 J=2,JMAX2
         DO 314 K=2,KMAX2
            IF(rho(J,K,1).LT.limiter) THEN
               WRITE(3,315) J,K
 315           FORMAT(' JTOP=',I3,' KTOP=',I3)
               GO TO 313
            END IF
 314     CONTINUE
 313  CONTINUE

c
c...Axisymmetric models done.  Restarted models done.  Return and hydro.
c

      IF((ITYPE.LE.5).OR.(ITYPE.EQ.8).OR.(ITYPE.EQ.9))
     &     RETURN

c-------------------------------------------------------------------------------
c     PERTURBATIONS
c
c...Maclaurin Bar Mode
      
      IF((ITYPE.EQ.6).OR.(ITYPE.EQ.96)) THEN

c...  AMP3 is amplitude of hit         
C...  SIGMA IS THE EIGENFREQUENCY

         AMP3=0.001
         SIGMA=0.037
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(AMP2,AMP1,AMP,j,k,them,l)     &
!$OMP&  SHARED(AMP3,sigma,DELR,dtheta)
         DO 712 L=1,LMAX
         THEM=2.0*(L-1)*DTHETA
         DO 712 K=2,KMAX1
         DO 712 J=3,JMAX
            AMP=-RHF(J)*COS(THEM)*(RHO(J+1,K,L)-RHO(J-1,K,L))/(2*DELR)
            AMP1=(SIGMA-2*OMEGA(J,K,L)-
     &           RHF(J)*(OMEGA(J+1,K,L)-OMEGA(J-1,K,L))
     &           /(2*DELR))*COS(THEM)*RHF(J)**2
            AMP2=R(J)*SIN(THEM)*(SIGMA-OMEGA(J,K,L)-OMEGA(J-1,K,L))
            AMP=AMP*AMP3
            AMP1=AMP1*AMP3
            AMP2=AMP2*AMP3
            IF(ABS(AMP).LE.RHO(J,K,L)) THEN
               RHO(J,K,L)=RHO(J,K,L)+AMP
            ELSE
               IF(AMP.GE.0.)RHO(J,K,L)=2.*RHO(J,K,L)
               IF(AMP.LT.0.)RHO(J,K,L)=0.d0
            END IF
            JN(J,K,L)=JN(J,K,L)+AMP1
            IF(RHO(J,K,L).EQ.0.)JN(J,K,L)=0.d0
            U(J,K,L)=AMP2
            IF(RHO(J,K,L).EQ.0.)U(J,K,L)=0.d0
 712     CONTINUE
!$OMP END PARALLEL DO
      END IF

c
c...Random initial hit
c
      IF((ITYPE.EQ.7).OR.(ITYPE.GE.97)) THEN
         print*,' RANDOM PERTURBATION'
c
c...amplitude of initial hit
c
         WRITE(3,1001)AMP0
 1001    FORMAT(5X,' RANDOM PERT, AMP0 =',1PE15.4)
         DO 462 L=1,LMAX
            !PERT=amp0*cos(16d0*dtheta*dble(L))
         DO 462 K=2,KMAX1
         DO 462 J=jmin,jmax1
            !PERT2=PERT+AMP0*(2.0*ran4(j+k+l)-1.)*1e-2
            PERT=ran4(j+k+l)
            PERT=AMP0*(2.0*PERT-1.0)
            RHO(J,K,L)=(1.0+PERT)*RHO(J,K,L)
 462     CONTINUE
      END IF


C....FIND POTENTIALS FOR INITIAL MODEL.
      IPRINT=1
      REDGE=0.d0
      gsoft = fgsoft*rof3n
      if(gsoft>zero)tmassacc=zero

        if(read_particle_file)then
         call read_particles(P_FILEID)
        else
         call initialize_particles()
        endif
        call set_particle_density()

#if SELFGRAVITY>0
!$OMP PARALLEL DEFAULT(SHARED)                                          &
!$OMP&  SHARED(JKMAX,ISYM)
      CALL SETBDY(0,ISYM)
!$OMP END PARALLEL

      CALL BDYGEN(MAXTRM,ISYM,REDGE)
      CALL POT3(8,IPRINT)
#endif
!!
#if WIGGLE>0
         call delta(.false.)
         call wiggle(rhotot,rhf,zhf,vx,vy,fx,fy,delt,pi,rpstar,phi_star,
     &        rof3n,zof3n,dtheta,gsoft,JMAX,KMAX,LMAX,JMIN,.true.,
     &        restart_wiggle)
#endif
      print *, "Set STAR"

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

#if EXTERNAL_POT>0
      call ExternalPotInit()
#endif

!$OMP END PARALLEL
      if(.not.read_particle_file)call set_particle_vel()
         
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(LP,LHF,j,k,l,LHAF)                &
!$OMP&  SHARED(rholmt,itype)
      if(itype.eq.7)then 
C....FROM U, W, AND JN, FIND S, T, AND A.
C....S ARE CENTERED AT THE CENTER OF R-DIRECTIONED CELL SURFACES.  T ARE
C....CENTERED AT THE CENTER OF Z-DIRECTIONED CELL SURFACES.  A ARE
C....CENTERED AT CELL CENTERS.
!$OMP DO SCHEDULE(STATIC)
      DO 630 K=2,KMAX1
      DO 630 J=2,JMAX1
      DO 630 L=1,LMAX
         S(J,K,L)=U(J,K,L)*(RHO(J,K,L)+RHO(J-1,K,L))/2.0
         T(J,K,L)=W(J,K,L)*(RHO(J,K,L)+RHO(J,K-1,L))/2.0
         A(J,K,L)=JN(J,K,L)*RHO(J,K,L)
630   CONTINUE
!$OMP END DO NOWAIT

      end if

!$OMP DO SCHEDULE(STATIC)
      DO 635 L=1,LMAX
      DO 635 K=2,KMAX1
      DO 635 J=2,JMAX1
         OMEGA(J,K,L) = JN(J,K,L)/(RHF(J)**2)
635   CONTINUE
!$OMP END DO NOWAIT

C...............CLEAN UP BY SETTING QUANTITIES AROUND Z-AXIS.

      IF (JMIN.GT.2) THEN 
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            DO K=1,KMAX
               DO J=1,JMIN2-1
                  S(J,K,L) = 0.0
                  T(J,K,L) = 0.0
                  A(J,K,L) = 0.0
                  U(J,K,L) = 0.0
                  OMEGA(J,K,L) = 0.0
                  JN(J,K,L)= 0.0
                  W(J,K,L) = 0.0
                  RHO(J,K,L) = RHOLMT
               ENDDO
            ENDDO
         ENDDO
!$OMP END DO         
      ELSE   
         LHAF=LMAX/2
!$OMP DO SCHEDULE(STATIC)         
         DO 680 K=2,KMAX1
            DO 680 L=1,LHAF
               LP=L+LHAF
               S(1,K,LP) = -S(3,K,L)
               T(1,K,LP) = T(2,K,L)
               A(1,K,LP) = A(2,K,L)
               U(1,K,LP) = -U(3,K,L)
               OMEGA(1,K,LP) = OMEGA(2,K,L)
               JN(1,K,LP)= JN(2,K,L)
               W(1,K,LP) = W(2,K,L)
               
               S(1,K,L) = -S(3,K,LP)
               T(1,K,L) = T(2,K,LP)
               A(1,K,L) = A(2,K,LP)
               U(1,K,L) = -U(3,K,LP)
               OMEGA(1,K,L) = OMEGA(2,K,LP)
               JN(1,K,L)= JN(2,K,LP)
               W(1,K,L) = W(2,K,LP)
 680     CONTINUE
!$OMP ENDDO
            
      ENDIF
      
C..............SET QUANTITIES BELOW THE EQUATORIAL PLANE.
!$OMP DO SCHEDULE(STATIC)
      DO 690 J=1,JMAX1
      DO 690 L=1,LMAX
         S(J,1,L) = S(J,2,L)
         T(J,1,L) = -T(J,3,L)
         A(J,1,L) = A(J,2,L)
         U(J,1,L) = U(J,2,L)
         OMEGA(J,1,L) = OMEGA(J,2,L)
         JN(J,1,L)= JN(J,2,L)
         W(J,1,L) = -W(J,3,L)
690   CONTINUE
!$OMP END DO NOWAIT

C....SET RHO AROUND Z-AXIS AND BELOW THE EQUATORIAL PLANE
      LHF=LMAX/2
!$OMP DO SCHEDULE(STATIC)
      DO 52 L=1,LHF
      LP=L+LHF
      DO 52 K=2,KMAX2
         RHO(1,K,LP)=RHO(2,K,L)
         RHO(1,K,L)=RHO(2,K,LP)
   52 CONTINUE
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
      DO 55 J=1,JMAX2
      DO 55 L=1,LMAX
   54   RHO(J,1,L)=RHO(J,2,L)
   55 CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
      
      RETURN
      END


C***********************************************************************
 

      SUBROUTINE RITE(IWHAT,IHEAD,  JST,JSP,JSK,  KST,KSP,KSK,
     &                LST,LSP,LSK)
      IMPLICIT real*8 (a-h,o-z)      

#include "hydroparam.h"
#include "globals.h"

      COMMON /INSIDE/TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT
      COMMON /TIMEST/INDX,ISOADI,ALLOW,DMAX,CHGMAX
      COMMON /ITS/ITSTRT,ITSTOP,ITSTEP


      REAL*8 FAREA(JMAX,KMAX),
     &       HAREA(JMAX,KMAX),
     &       RD(JMAX2),
     &       ZD(KMAX2),
     &       RINV(JMAX2),
     &       Tauross(jmax2,kmax2,lmax)
      CHARACTER np*2,tw*3
      CHARACTER index*6,savedfile*80,epsfull*80
      save ekold,egold,pdvold
      integer jstart

      real*8 OMMAX

      npr=int(10.0*nprime)
      tovw=int(1000.d0*toverw)
      write(index,'(i6.6)') itstep
      write(np,'(i2.2)') npr
      if(nprime.gt.3.0) np="in"
      write(tw,'(i2.2)') int(tovw)
      OMMAX=zero
 
C     IWHAT = 1  PRINTS ALL VARIABLES OUT IN 1PE12.4 FORMAT.
C           = 0  BRIEF DIAGNOSTICS (2 LINES) ONLY.
C           = 2  STORES THIS MODEL.
C           = 10 PRINTS INTEGRATED PARAMETERS.
C     IHEAD = NEGATIVE SKIPS TO A NEW PAGE?  IHEAD.GE.0 DOES NOT.
C     ABS(IHEAD) = 1  PRINTS HEADING ONCE PER CALL TO SUBROUTINE.
      
c     Initialization of previously uninitialized variables.  I make no
c     guarantee as to the correctness of these values. -rpl

      xn1 = 0.0

c     end of bogus initializations.
      
 
 100  FORMAT('1')
 101  FORMAT(//,' J  K  L ',6X,'S',11X,'T',11X,'A',11X,'U',11X,'W',10X,
     &     'JN',9X,'OMEGA',9X,'EPS',9X,'P',10X,'RHO',9X,'PHI',/)
 102  FORMAT(3I3,1P11E12.4)
 103  FORMAT(//)
 104  FORMAT('   TSTEP',4X,'TIME',8X,'DELT',8X,'ETOT/JT',4X,'EGRAV/ROT',
     &     2X,'EKIN/RZKIN',4X,'ENEW/CD',4X,'EDIF/DMAX',3X,'ELOST/JKL',
     &     4X,'TMASS',7X,'PHICHK',4X,'K',/)
 105  FORMAT(I8,1P10E12.4,I4)
 106  FORMAT(4X,1PE12.4,' CODETIME',6X,1P5E12.4,I4,2I3,2X,1P2E12.4,/)
 118  FORMAT(///)
 119  FORMAT(2I4,1P5E11.3,2X,1P5E11.3)
 121  FORMAT('INITIAL, TOTAL, ADDED, OUTFLOW, ACCRETED MASSES',1P5E16.6)

C-----------------------------------------------------------------------

      jstart=jmin
 
      IF(IWHAT.NE.1) GO TO 50
      IF(IHEAD.LT.0) WRITE(3,100)
      IF(IABS(IHEAD).EQ.1) WRITE(3,101)
      IF(IABS(IHEAD).NE.1) WRITE(3,103)
      DO L=LST,LSP,LSK
        DO K=KST,KSP,KSK
          DO J=JST,JSP,JSK
            WRITE(3,102) J,K,L,S(J,K,L),T(J,K,L),A(J,K,L),
     &           U(J,K,L),W(J,K,L),JN(J,K,L),OMEGA(J,K,L),
     &           EPS(J,K,L),P(J,K,L),RHO(J,K,L),PHI(J,K,L)
          ENDDO
        ENDDO
      ENDDO
      GO TO 999



C-----------------------------------------------------------------------
 50   IF(IWHAT.NE.0) GO TO 95
      IF(IHEAD.LT.0) WRITE(3,100)
      IF(IHEAD.GT.0) WRITE(3,103)
      IF(IABS(IHEAD).EQ.1) WRITE(3,104)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rhomax_cap1,LD_CAP1,KD_CAP1,     &
!$OMP& JD_CAP1,AREA1,VOL,AREA2,l,AREA,k,j)                              &
!$OMP&  SHARED(EKIN,RZKIN,EROT,TOTJN,PI,RZK,EGRAV,jstart,ER,TOTJ,EG,    &
!$OMP& DTHETA,LD,KD,JD,rhomax)
!$OMP DO SCHEDULE(STATIC)
      DO J=2,JMAX1
        RD(J)=R(J+1)**2-R(J)**2
      ENDDO
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
      DO K=2,KMAX1
        ZD(K)=Z(K+1)-Z(K)
      ENDDO
!$OMP ENDDO

C  FIND EGRAV, TOTJ, EROT.
      AREA=0.5*DTHETA
!$OMP MASTER
      EG=0.d0
      TOTJ = 0.d0
      ER = 0.d0
!$OMP END MASTER
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ER,TOTJ,EG)
      DO L=1,LMAX
        DO K=2,KMAX
          AREA2=AREA*ZD(K)
          DO J=jstart,JMAX
            VOL = RD(J)*AREA2
            EG = EG + VOL*(PHI(J,K,L)+starphi(J,K,L)+tinphi(J,K,L))     &
     &          * RHO(J,K,L)
            TOTJ = TOTJ + A(J,K,L)*VOL
            ER = ER + 0.5*VOL*A(J,K,L)*OMEGA(J,K,L)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP MASTER
      EGRAV=EG
!$OMP END MASTER
!$OMP DO SCHEDULE(STATIC)
      DO J=3,JMAX1
        RD(J)=RHF(J)**2-RHF(J-1)**2
        RINV(J)=1.0/R(J)**2
      ENDDO
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
      DO K=2,KMAX1
        ZD(K)=ZHF(K)-ZHF(K-1)
      ENDDO
!$OMP END DO
!$OMP MASTER
      ZD(2)=0.5*ZD(2)

C  FIND EROT, RZKIN, TOTJN.
      RZK=0.d0
!$OMP END MASTER
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RZK)
      DO L=1,LMAX
        DO K=2,KMAX1
          AREA1=AREA*ZD(K)
          DO J=3,JMAX1
            VOL=RD(J)*AREA1
            RZK=RZK + 0.5*VOL*(U(J,K,L)*S(J,K,L) + W(J,K,L)*T(J,K,L))
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
      AREA1=PI*RHF(2)**2
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RZK)
      DO K=2,KMAX1
        VOL=AREA1*ZD(K)
        RZK=RZK + 0.5*VOL*W(2,K,1)*T(2,K,1)
      ENDDO
!$OMP END DO
!$OMP MASTER
      TOTJN=2.0*TOTJ
      EROT=2.0*ER
      RZKIN=2.0*RZK
      EKIN=EROT+RZKIN
      JD=2
      KD=2
      LD=1
      RHOMAX=RHO(jmin,2,1)
!$OMP END MASTER

!$OMP BARRIER
      rhomax_cap1=rhomax
      JD_CAP1=JD
      KD_CAP1=KD
      LD_CAP1=LD
!$OMP DO SCHEDULE(STATIC)
        do l=1,lmax,1
          do k=2,kmax,1
            do j=2,jmax,1
            if (rho(j,k,l).gt.rhomax_cap1) then
              rhomax_cap1=rho(j,k,l)
              JD_CAP1=j
              KD_CAP1=k
              LD_CAP1=l
            endif
            enddo
          enddo
        enddo
!$OMP END DO 
!$OMP CRITICAL
      if (rhomax_cap1.gt.rhomax) then
        rhomax=rhomax_cap1
        JD=JD_CAP1
        KD=KD_CAP1
        LD=LD_CAP1
      endif
!$OMP END CRITICAL
!$OMP END PARALLEL
      DMAX=rhomax
      CD=rho(2,2,1)
      PDV=0.d0
      ETOT=EGRAV+EKIN+ENEW+ELOST-PDV
      ECHECK=0.d0
      if (JSP.gt.0) then
        ECHECK=(EKIN-ekold+EDIF)/(egold-EGRAV+PDV-pdvold)-1.0
      endif
      egold=EGRAV
      pdvold=PDV
      ekold=EKIN
      write(3,105)JST,time,delt,ETOT,EGRAV,EKIN,ENEW,EDIF,ELOST,tmass,  &
     &PHICHK,KLOCAT
      write(3,106)TIME,TOTJN,EROT,RZKIN,CD,DMAX,JD,KD,LD,ECHECK
      write(3,121)tmassini,tmass,tmassadd,tmassout,tmassacc
      GOTO 999

C-----------------------------------------------------------------------
C  Store the model.
C
   95 CONTINUE
      IF(IWHAT.NE.2) GO TO 300
      NTAPES=8
      if (ITSTEP .eq. ITSTOP) ITSTEP = JST
      NUM=1
      
      savedfile='saved.'//index
      open(unit=8,file=savedfile,form='unformatted')
      WRITE(8) S
      WRITE(8) T
      WRITE(8) A
      WRITE(8) RHO
      WRITE(8) EPS
      WRITE(8) ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,JREQ,OMMAX
      if (jmin.gt.2) write(8) tmassini,tmass,tmassadd,tmassout
     &     ,tmassacc,totcool,totdflux,totheat,totirr,etotfl,eflufftot
      close(8)
      WRITE(3,110) JST
      WRITE(3,111)
  110 FORMAT('1',///,
     &  ' SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS',/,
     &  ' S                                            S',/,
     &  ' S  MODEL FROM TIME STEP NUMBER ',I8,' HAS  S'  ,/,
     &  ' S  BEEN STORED IN FILE FORT.8                S',/,
     &  ' S                                            S',/,
     &  ' SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS',/////)
  111 FORMAT(6X,'DATE STORED:',/,6X,'DATE PURGED:',/,10X,'PF NAME:',/,
     &5X,'CYCLE NUMBER:',/,12X,'OTHER:',////)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,j,l)
      do l=1,lmax
         do j=1,jmax2
            do k=1,kmax2
               tauross(j,k,l)=tau(j,k,l,1)
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO

      epsfull='coolheat_full.'//index
      open(unit=9,file=epsfull,form='unformatted')
      write(9) divflux
      write(9) lambda
      write(9) hgamma
      write(9) igamma
      write(9) tauross
      write(9) TempK
      write(9) TeffK
      write(9) TphK
      write(9) time
      close(9)

      open(unit=9,file="gamma1."//index,form='unformatted')
      write(9) gamma1
      write(9) time
      close(9)

C-----------------------------------------------------------------------
C  Print some grid spacing information.
C
  300 CONTINUE
      IF(IWHAT.NE.3) GO TO 400
      WRITE(3,115) ROF3N,A1NEWR,CORMAS
      WRITE(3,116) R,RHF
      WRITE(3,117) ZOF3N,A1NEWZ
      WRITE(3,116) Z,ZHF
  115 FORMAT(//,'1',5X,'R-DIRECTION',/,' ROF3N   ',1PE16.6,
     &  /,' A1NEWR  ',1PE16.6,/,' CORMAS  ',1PE16.6)
  116 FORMAT(1P8E15.6)
  117 FORMAT(////,6X,'Z-DIRECTION',/,' ZOF3N   ',1PE16.6,
     &  /,' A1NEWZ  ',1PE16.6)



C-----------------------------------------------------------------------
C  Not sure what the following does.
C
  400 CONTINUE
      IF(IWHAT.NE.4) GO TO 450
      WRITE(3,118)
      LL=IHEAD
      RHOMAX=0.d0
      DO L=1,LMAX
        DO K=2,KMAX
          DO J=2,JMAX
            IF(RHO(J,K,L).GT.RHOMAX) THEN
              RHOMAX=RHO(J,K,L)
              JRHO=J
              KRHO=K
              LRHO=L
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      IF(LL.EQ.0) LL=LRHO
      IF(LL.LT.1.OR.LL.GT.LMAX)LL=1
      WRITE(3,125)
  125 FORMAT('  J  KRHO',2X,'J...R',7X,'RHO',8X,'PHI',8X,'OMG',9X,'U',
     &  10X,'K...Z',7X,'RHO',8X,'OMG',9X,'W')
      DO KK=2,JMAX1
        WRITE(3,119)KK,JRHO,R(KK),RHO(KK,2,LL),PHI(KK,2,LL),
     &        OMEGA(KK,2,LL),U(KK,2,LL),Z(KK),RHO(JRHO,KK,LL),
     &        OMEGA(JRHO,KK,LL),W(JRHO,KK,LL)
      ENDDO
      KSTRTJ=KMAX1+1
      IF(JMAX1.GT.KMAX1) THEN
        DO KK=KSTRTJ,JMAX1
          WRITE(3,119)KK,KRHO,R(KK),RHO(KK,2,LL),PHI(KK,2,LL),
     &            OMEGA(KK,2,LL),U(KK,2,LL)
        ENDDO
      ENDIF
      IF(LL.GE.0) GO TO 450
      DO KK=2,KMAX1
        JJ=KK/2+1
        D1=SQRT(RHF(KK)**2+ZHF(KK)**2)
        D2=SQRT(R(KK)**2+Z(KK)**2)
        D3=SQRT(RHF(JJ)**2+ZHF(KK)**2)
        D4=SQRT(R(JJ)**2+Z(KK)**2)
        V1=SQRT(U(KK,KK,LL)**2+W(KK,KK,LL)**2)
        V2=SQRT(U(JJ,KK,LL)**2+W(JJ,KK,LL)**2)
        WRITE(3,119)JJ,KK,D1,RHO(KK,KK,LL),P(KK,KK,LL),D2,V1,D3,
     &          RHO(JJ,KK,LL),P(JJ,KK,LL),D4,V2
      ENDDO
      WRITE(3,118)


C-----------------------------------------------------------------------
  450 CONTINUE
      IF(IWHAT.LT.10) GOTO 999

C  The following section prints out various parameters of a model
C  integrated over cylinders about the rotation axis. All quantities
C  except loss rates are calculated on half-integer grid points
C  including the translational kinetic energy. Loss rates of mass and
C  angular momentum are at the integer grid points.
 
C  Set up needed areas and volumes.
C  In the following 4.*DTHETA has been changed to 2.*DTHETA. 4.*DTHETA was
C  used for zero to pi calculations.
        open(unit=33, file='mfrac.dat')
        R1=2.*DTHETA
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,R4,R2,j)                    &
!$OMP&  SHARED(R1)
        do j=2,jmax,1
        R2=R1*r(j)
        R4=R1*(r(j+1)-r(j))*rhf(j)
          do k=2,kmax,1
          FAREA(j,k)=R2*(z(k+1)-z(k))
          HAREA(j,k)=R4*(z(k+1)-z(k))
          enddo
        enddo
!$OMP END PARALLEL DO
C  Print output header.
      write(3,502)
  502   FORMAT(/,3X,'MODEL PROPERTIES INTEGRATED OVER CYLINDERS'/)
      write(3,503)
  503   FORMAT(2X,'J',4X,'R  ',4X,'MFRAC ',4X,'MDIFF ',5X,              &
     &  'MDOT ',3X,'ANGSUM ',2X,'ANGDIFF ',3X,'ANGDOT ',4X,             &
     &  'EGRAV ',5X,'EROT ',3X,'ETRANS ',3X,'TOTKIN ',5X,'EINT '        &
     &  ,5X,'ETOT '/)
!  Calculate and print mass and angular momentum distributions and flow
!  and various energies.
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ETOTT,TKINT,ETOT,TKIN,W1,JW,U1,  &
!$OMP& L1,l,EI1,ET1,ER1,EG1,ADOT1,dot1,ADIFF1,DIFF1,K1,k,J1,j,EIT,ETT,  &
!$OMP& ERT,EGT,R2,R1,AFRAC,FRAC)                                        &
!$OMP&  SHARED(EI,ET,ER,EG,ADOT,DOT,ADIFF,DIFF,jstart,gamma)
      FRAC=0.d0
      AFRAC=0.d0
      R1=1./128.
      R2=1./(gamma-1.)
      EGT=0.d0
      ERT=0.d0
      ETT=0.d0
      EIT=0.d0
        do 504  j=jstart,jmax,1
        J1=j+1
!$OMP BARRIER
!$OMP MASTER
        DIFF=0.d0
        ADIFF=0.d0
        DOT=0.d0
        ADOT=0.d0
        EG=0.d0
        ER=0.d0
        ET=0.d0
        EI=0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:ADOT,DOT,EI,ET,ER,EG,ADIFF,DIFF)
          do 505  k=2,kmax,1
          K1=k+1
          DIFF1=0.d0
          ADIFF1=0.d0
          dot1=0.d0
          ADOT1=0.d0
          EG1=0.d0
          ER1=0.d0
          ET1=0.d0
          EI1=0.d0
            do 506  l=1,lmax,1
            L1=l+1
            if (l.eq.lmax) then
              L1=1
            endif
            DIFF1=DIFF1+rho(j,k,l)
            ADIFF1=ADIFF1+a(j,k,l)
            if (jmin.gt.2) then
              EG1=EG1+rho(j,k,l)*(phi(J,K,L)+starphi(J,K,L)             &
     &           + tinphi(J,K,L))
            else
            endif
            ER1=ER1+a(j,k,l)*omega(j,k,l)*0.5
            EI1=EI1+p(j,k,l)*R2
            U1=u(J1,k,l)+u(J1,K1,l)+u(J1,k,L1)+u(J1,K1,L1)
            if (U1.gt.0.) then
              JW=j
            endif
            if (U1.LE.0.) then
              JW=J1
            endif
Cacm          DOT1=DOT1+U1*RHO(JW,K,L)*0.25
            dot1=dot1+s(j,k,l)
            ADOT1=ADOT1+U1*a(JW,k,l)*0.25
            U1=U1+u(j,k,l)+u(j,K1,l)+u(j,k,L1)+u(j,K1,L1)
            U1=U1*U1
            W1=w(j,k,l)+w(j,K1,l)+w(j,k,L1)+w(J1,K1,L1)+w(J1,k,l)+w(J1, &
     &      K1,l)+w(J1,k,L1)+w(J1,K1,L1)
            W1=W1*W1
            ET1=ET1+rho(j,k,l)*(U1+W1)*R1
  506       CONTINUE
          DIFF=DIFF+DIFF1*HAREA(j,k)
          ADIFF=ADIFF+ADIFF1*HAREA(j,k)
          EG=EG+EG1*HAREA(j,k)
          ER=ER+ER1*HAREA(j,k)
          ET=ET+ET1*HAREA(j,k)
          EI=EI+EI1*HAREA(j,k)
          DOT=DOT+dot1*FAREA(j,k)
          ADOT=ADOT+ADOT1*FAREA(j,k)
  505     CONTINUE
!$OMP END DO
        TKIN=ER+ET
        ETOT=EG+TKIN+EI
        FRAC=FRAC+DIFF
        AFRAC=AFRAC+ADIFF
        EGT=EGT+EG*0.5
        ERT=ERT+ER
        ETT=ETT+ET
        TKINT=ERT+ETT
        EIT=EIT+EI
        ETOTT=EGT+TKINT+EIT
!$OMP MASTER
        write(3,507)j,r(j),FRAC,DIFF,DOT,AFRAC,ADIFF,ADOT,EG,ER,ET,TKIN,&
     &  EI,ETOT
        write(33,667)j,r(j),FRAC,DIFF,DOT,AFRAC,ADIFF,ADOT
        write(3,508)EGT,ERT,ETT,TKINT,EIT,ETOTT
!$OMP END MASTER
  504   CONTINUE
!$OMP END PARALLEL
  507   FORMAT(1X,I4,F7.3,1P12E10.2/)
  667   FORMAT(1X,I4,F7.3,1P7E10.2/)
  508   FORMAT(45X,'TOTAL ENERGIES =    ',1P6E10.2/)
C-----------------------------------------------------------------------
  999 RETURN
      END
      
C***********************************************************************
      subroutine dumpphi(a,i1max,i2max,i3max,iunit)
      implicit real*8(a-h,o-z)
#include "hydroparam.h"
      real*8 a(pot3jmax2,pot3kmax2,lmax)
      integer i1,i2,i3
      
      do i3=1,i3max
         write(iunit,666) i3
         do i2=1,i2max
            write(iunit,*) (a(i1,i2,i3),i1=1,i1max)
            write(iunit,*) ' '
         end do
      end do
 666  format('l= ',I5)
      close(iunit)
      return
      end

      

C***********************************************************************
      subroutine cyl3d
      IMPLICIT real*8 (a-h,o-z)

#include "hydroparam.h"
#include "globals.h"

      integer jm,km,km2
      parameter (jm=jmax-1,km=kmax-1,km2=2*(kmax-1))
      real*4 rhoout(jm,km2,lmax)

c     (don't use j=1 or k=1, flip k=2--kmax)
      integer j,k,l,jp,kp
      
      do l=1,lmax
        do j=1,jm
          jp=j+1
          do k=1,km
            kp=k+1
            rhoout(j,kmax-k,l) = log10(rho(jp,kp,l))
            rhoout(j,km+k,l) = log10(rho(jp,kp,l))
          end do
        end do
      end do
      
      write (14) rhoout
      return
      end
      


C***********************************************************************
      subroutine lowresrandom(rho,jn,denny,anggy)
c     perturb the model with the random perturbation from a low
c     resolution model.  -rpl feb. 1996
      IMPLICIT real*8 (a-h,o-z)      
#include "hydroparam.h"
      parameter (amp0=0.005)
      parameter (lrkmax1=lrkmax+1)
      parameter (jratio=jmax/lrjmax)
      parameter (kratio=kmax/lrkmax)
      parameter (lratio=lmax/lrlmax)
      real*8 rho(jmax2,kmax2,lmax)
      real*8 jn(jmax2,kmax2,lmax)
      real*8 denny(hj2,hk2)
      real*8 anggy(hj1,hk1)
      real*8 pert
      integer j,k,l,jj,kk,ll
      integer jmap,kmap,lmap
      
      write (3,901) lrjmax,lrkmax,lrlmax,amp0
 901  format(5X,'Random perturbation, lowresmdel: ',I5,'x',I5,'x',
     $     I5,'  amp0 = ',1PE15.4) 
      do l=1,lmax
        do k=2,kmax1
          do j=2,jmax1
            jn(j,k,l)=anggy(j,k)
          end do
        end do
      end do
      
      do l=1,lrlmax
        do k=1,lrkmax
          do j=1,lrjmax
c           pert=drand48()
            pert=rand()
            pert=amp0*(2.0*pert-1.0)
            lmap = lratio*(l-1) + 1
            kmap = 1+ kratio*(k-1) + 1 ! extra +1 because k goes from 2 to kmax+1 
            jmap = 1+ jratio*(j-1)+1 ! ditto
            do ll=lmap,lmap+lratio
              do kk=kmap,kmap+kratio
                do jj=jmap,jmap+jratio
                  if ((ll .gt. lmax) .or. (kk.gt.kmax1) .or. 
     $                 (jj .gt. jmax1)) then
                    continue
                  else
                    rho(jj,kk,ll) = (1.0+pert)*denny(jj,kk)
                  end if
                end do
              end do
            end do
          end do
        end do
      end do
      return
      end
      


      Real*8 FUNCTION RAN4(IDUM)
C* math can be done in integer if two comments Cs are moved
C* see numerical recipes
      Implicit None
      Integer iff,idum,i,inext,inextp,k,ii
      Real*8 mbig, mseed, mz, fac, ma(55), mj, mk
      save inext,INEXTP, ma
      
      PARAMETER (MBIG=4000000.D0,MSEED=1618033.D0,MZ=0.D0,FAC=2.5D-7)
CC      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DATA IFF /0/

      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
         IFF=1
         MJ=MSEED-dble(IABS(IDUM))
         MJ=MOD(MJ,MBIG)
         MA(55)=MJ
         MK=1
         DO I=1,54,1
            II=MOD(21*I,55)
            MA(II)=MK
            MK=MJ-MK
            IF(MK.LT.MZ)MK=MK+MBIG
            MJ=MA(II)
         END DO
         DO K=1,4,1
            DO I=1,55,1
               MA(I)=MA(I)-MA(1+MOD(I+30,55))
               IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
            END DO
         END DO
         INEXT=0
         INEXTP=31
         IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.ge.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.ge.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN4=MJ*FAC

      RETURN
      END
