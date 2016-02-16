      program poissonSolver
      implicit real*8(a-h,o-z)
#include "hydroparam.h"
#include "globals.h"
#include "units.h"


      CHARACTER potfile*80,x1*6
      common /coefs/coef(pot3jmax2,pot3kmax2,lmax2,2)

      integer, PARAMETER::JKM1=2*POT3JMAX+POT3KMAX-1
!$      integer OMP_GET_MAX_THREADS
      dimension denny(hj2,hk2)

   98 FORMAT(1E25.15)
  100 FORMAT(1P2E15.8)
  102 FORMAT(9(I6,1X))

!...  Read in some run parameters.
      OPEN(UNIT=5,FILE='fort.5',STATUS='OLD')

      READ(5,100) KONST,XN
      WRITE(6,10010) KONST,XN
10010 FORMAT('  KONST  ',1PE20.12,/,'  XN     ',1PE20.12)
      READ(5,98)GAMMA
      WRITE(6,10020) GAMMA
10020 FORMAT('  GAMMA  ',1PE20.12)

      READ(5,102)ITSTRT,ITSTOP,IDIAG,ISOADI,ITYPE,NMODL,ISTOR,IGRID
      WRITE(6,10030)ITSTRT,ITSTOP,IDIAG,ISOADI,ITYPE,NMODL,ISTOR,IGRID  &
     &     ,JMIN
10030 FORMAT('  ITSTRT ',I8,/,'  ITSTOP ',I8,/,'  IDIAG  ',I8,          &
     &     /,'  ISOADI ',I8,/,'  ITYPE  ',I8,/,'  NMODL  ',I8,          &
     &     /,'  ISTOR  ',I8,/,'  IGRID  ',I8,/,'  JMIN   ',I8)          

      READ(5,100) NPRIME
      WRITE(6,10040) NPRIME
10040 FORMAT('  NPRIME ',1PE20.12)
      NTAPE=7
      CLOSE(5)





!...Define some things...
      MAXTRM=10
      ISYM  =2
      KWFW=int(log10(dble(KMAX))/log10(two))-1 

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

      if(ITYPE.EQ.7) then
! Read in the density array, denny
      OPEN(UNIT=2,FILE='fort.2',STATUS='OLD')
      READ(2,*)PINDEX,CON2,RRR2,OMCEN,DENCEN,TOVERW,ROF3N,ZOF3N,        &
     &        A1NEWZ,JREQ,KZPOL
 1685    FORMAT(3X,1PE22.15,2X,8(1PE22.15,2X),2I4)

      write(*,1685)PINDEX,CON2,RRR2,OMCEN,DENCEN,TOVERW,                &
     &        ROF3N,ZOF3N,A1NEWZ,JREQ,KZPOL
      READ(2,*) DENNY
 1617    FORMAT(8(1PE22.15,2X))
      CLOSE(2)
      DEN=DENCEN

!! Define rho array from read in data.
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(LP,j,k,l)
!&
!$OMP&  SHARED(gamma,konst,den)
!$OMP DO SCHEDULE(STATIC)
         do l=1,lmax
            do k=2,kmax1
               do j=2,jmax1
                  rho(j,k,l)=denny(j,k)
                  if(rho(j,k,l).lt.gridlim*den) then
                    rho(j,k,l)=gridlim*den
                  endif
               end do
            end do
         end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      elseif(ITYPE.EQ.1) then

            OPEN(UNIT=7,FILE='fort.7',FORM='UNFORMATTED',STATUS='OLD')
            read(7) S
            read(7) T
            read(7) A
            read(7) RHO
            read(7) EPS
            read(7)ROF3N,ZOF3N,DELT,TIME,ELOST,DEN,SOUND,               &
     &        JREQ,OMMAX
            if (jmin.gt.2) then
            read(7) tmassini,tmass,tmassadd,                            &
     &         tmassout,tmassacc,totcool,totdflux,totheat,totirr,etotfl,&
     &         eflufftot  !ACB
            endif
            dencen=den
            rholmt=dencen*gridlim
            epslmt=(1.d0/(gamma-1.0))*rholmt**gamma*gridlim
            CLOSE(7)


      endif
!...Set up the grid. (From radhydro/io.f)
!...grid setup
      DELR=ROF3N
      DELZ=ZOF3N
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,j)
!&
!$OMP&  SHARED(DELZ,DELR)
!$OMP DO SCHEDULE(STATIC)
      DO J=1,JMAX2
         R(J)=(J-2)*DELR
         RHF(J)=R(J)+DELR/2.0
      END DO
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
      DO K=1,KMAX2
         Z(K)=(K-2)*DELZ
         ZHF(K)=Z(K)+DELZ/2.0
      END DO
!$OMP END DO nowait
!$OMP END PARALLEL

!...Calling the potential solver now.
      IPRINT = 0
      REDGE  = 0.d0
!$OMP PARALLEL DEFAULT(SHARED)
!&
!$OMP&  SHARED(JKMAX,ISYM)
         CALL SETBDY(0,ISYM)
!$OMP END PARALLEL
         CALL BDYGEN(MAXTRM,ISYM,REDGE)
         CALL POT3(8,IPRINT)

!...Write gravitational potential to output file.
         write(x1,'(i6.6)')ITSTOP
            potfile='phi3d.'//trim(x1)
            OPEN(UNIT=23,FILE=potfile,FORM='UNFORMATTED')
            write(23)phi
            write(23)time
            CLOSE(23)


      end program 

