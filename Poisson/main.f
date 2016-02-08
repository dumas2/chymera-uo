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

! Read in the density array, denny
      OPEN(UNIT=2,FILE='fort.2',STATUS='OLD')
      READ(2,*)PINDEX,CON2,RRR2,OMCEN,DENCEN,TOVERW,ROF3N,ZOF3N,
     &        A1NEWZ,JREQ,KZPOL
 1685    FORMAT(3X,1PE22.15,2X,8(1PE22.15,2X),2I4)

      write(*,1685)PINDEX,CON2,RRR2,OMCEN,DENCEN,TOVERW,
     &        ROF3N,ZOF3N,A1NEWZ,JREQ,KZPOL
      READ(2,*) DENNY
 1617    FORMAT(8(1PE22.15,2X))
      CLOSE(2)
!! Define rho array from read in data.
      DEN=DENCEN
      print *, 'gridlim,den = ', gridlim, den

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

        print *, 'rho = ',rho(2,65,75)


! Setup the grid. (From radhydro/io.f)
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
        print *, 'here'
!$OMP PARALLEL DEFAULT(SHARED)
!&
!$OMP&  SHARED(JKMAX,ISYM)
         CALL SETBDY(0,ISYM)
!$OMP END PARALLEL
        print *, 'there'
        print *, 'here'
         CALL BDYGEN(MAXTRM,ISYM,REDGE)
        print *, 'there'
        print *, 'here'
         CALL POT3(8,IPRINT)
        print *, 'there'

!...Write gravitational potential to output file.
         write(x1,'(i6.6)')ITSTOP
            potfile='phi3d.'//trim(x1)
            OPEN(UNIT=23,FILE=potfile,FORM='UNFORMATTED')
            write(23)phi
            write(23)time
            CLOSE(23)


      end program 

