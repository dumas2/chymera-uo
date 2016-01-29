C***********************************************************************
      SUBROUTINE RAD(I)
      IMPLICIT real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"
      common /INSIDE/TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT
C Arrays rd and zd are used only in this subroutine, and thus should not
C need to be put in a common block. But the OpenMP version of the code
C might fail if they are not.
      real*8 RD(jmax2),ZD(kmax2)
      common /RADPRIVATE/RD,ZD
C     EVENTUALLY THIS ROUTINE WILL DO THE RADIATIVE TRANSFER.
C     IF I=1, ISOTHERMAL CLOUD.
C     IF I=2, ADIABATIC COLLAPSE.
C     IF I=3, INTERNAL ENERGY TRANSPORT ON
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(K,J)
!$OMP DO SCHEDULE(STATIC)
        DO J=2,jmax1,1
        RD(J)=r(J+1)**2-r(J)**2
        ENDDO
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
        DO K=2,kmax1,1
        ZD(K)=z(K+1)-z(K)
        ENDDO
!$OMP END DO nowait
!$OMP END PARALLEL
      AREA2=0.5*DTHETA
C-----------------------------------------------------------------------
      IF (I.NE.1) THEN
        GOTO 200
      ENDIF
C  Isothermal case.
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(E2,E1,VOL,J,K,L,EN)              &
!$OMP&  SHARED(AREA2,bdytem,EPSNEW,EPSSUM)
!$OMP MASTER
      EPSSUM=0.0
      EPSNEW=0.0
!$OMP END MASTER
!xOMP END SINGLE
      EN=cv(2,2,1)*bdytem
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:EPSNEW,EPSSUM)
        DO L=1,lmax,1
          DO K=2,kmax1,1
            DO J=2,jmax1,1
            VOL=AREA2*ZD(K)*RD(J)*rho(J,K,L)
            E1=eps(J,K,L)*VOL
            E2=EN*VOL
            EPSSUM=EPSSUM+(E1-E2)
            EPSNEW=EPSNEW+E2
            eps(J,K,L)=EN
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
        DO L=1,lmax,1
          DO K=2,kmax1,1
          eps(1,K,L)=eps(2,K,L)
          ENDDO
          DO J=1,jmax1,1
          eps(J,1,L)=eps(J,2,L)
          ENDDO
        ENDDO
!$OMP END DO nowait
!$OMP END PARALLEL
C  EDIF is energy lost via radiation during this timestep alone.
      EDIF=2.0*EPSSUM
C  Total energy radiated away since initial model is in ELOST.
      ELOST=ELOST+EDIF
C  ENEW should remain constant in time. It is total internal energy.
      ENEW=2.0*EPSNEW
      GOTO 999
C-----------------------------------------------------------------------
  200 IF (I.NE.2) THEN
        GOTO 300
      ENDIF
C.... Adiabatic collapse
C.... ISOAD=2 means isentropic fluid.
      EPSSUM=0.0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(VOL,J,K,L)                    &
!$OMP&  SHARED(AREA2) REDUCTION(+:EPSSUM)
        DO L=1,lmax,1
          DO K=2,kmax1,1
            DO J=2,jmax1,1
            VOL=AREA2*ZD(K)*RD(J)
c           EPS2=EPS2+RHO(J,K,L)**(1.0/XN)*VOL
            EPSSUM=EPSSUM+eps(J,K,L)*VOL
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
CxOMP END DO
CxOMP SINGLE
      EPSSUM=2.0*EPSSUM
C The following calculates the internal energy.
C     ENEW=EPS2*2./(GAMMA-1.)
C EDIF is increase in internal energy this time step.
      EDIF=EPSSUM-ENEW
C Total internal energy is in ENEW.
      ENEW=EPSSUM
C-----------------------------------------------------------------------
  300 IF(I.NE.3) GO TO 999
      EPSSUM=0.0
      EOSSUM=0.0

      DO L=1,LMAX
        DO K=2,KMAX1
          DO J=2,JMAX1
            VOL=AREA2*ZD(K)*RD(J)
C           EPSSUM=EPSSUM+EPS(J,K,L)*VOL
            EOSSUM=EOSSUM+P(j,k,l)*VOL/(gamma-1.0)
          ENDDO
        ENDDO
      ENDDO

      EOSSUM=2.0*EOSSUM
C EDIF is increase in internal energy this time step.
      EDIF=EOSSUM-ENEW
C Total internal energy is in ENEW.
      ENEW=EOSSUM

C-----------------------------------------------------------------------
  999 CONTINUE

      RETURN
      END


