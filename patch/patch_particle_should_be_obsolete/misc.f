C***********************************************************************
      SUBROUTINE VELOCITY
      IMPLICIT real*8 (a-h,o-z)

#include "units.h"
#include "hydroparam.h"
#include "globals.h"
      real*8 rhox(kmax2)
      save rhox
      integer jstart,LP,K

C...FROM MOMENTA FIND VELOCITIES

      if (jmin.gt.2) then
         jstart=jmin
      else
         jstart=2
      endif

!$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            U(J,K,L)=two*S(J,K,L)/(RHO(J,K,L)+RHO(J-1,K,L))
            W(J,K,L)=two*T(J,K,L)/(RHO(J,K,L)+RHO(J,K-1,L))
            JN(J,K,L)=A(J,K,L)/RHO(J,K,L)
#if ROTATING>0
            OMEGA(J,K,L) = JN(J,K,L)/(RHF(J)**2)-omega_frame
#else
            OMEGA(J,K,L) = JN(J,K,L)/(RHF(J)**2)
#endif
          ENDDO
        ENDDO
      ENDDO
!$OMP ENDDO nowait

!$OMP DO SCHEDULE(STATIC)
      DO K=2,KMAX2
        RHOX(K)=zero
      ENDDO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:RHOX)
        do L=1,lmax,1
          do K=2,kmax2,1
          RHOX(K)=RHOX(K)+rho(2,K,L)
          ENDDO
        ENDDO
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
        do K=2,kmax2,1    
        RHOX(K)=RHOX(K)/lmax
        ENDDO
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
        do L=1,lmax,1
          do K=2,kmax2,1    
          u(2,K,L)=s(2,K,L)/RHOX(K)
          ENDDO
        ENDDO
!$OMP ENDDO nowait

C...SET VELOCITIES AROUND Z-AXIS.
      if (jmin.gt.2) then
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            DO K=2,KMAX1
               do j=1,jmin1
                  U(1,K,L)    = zero
                  W(1,K,L)    = zero
                  OMEGA(1,K,L)= zero
                  JN(1,K,L)   = zero
               ENDDO
            ENDDO
         ENDDO
!$OMP END DO NOWAIT        
      else  
!$OMP DO SCHEDULE(STATIC) 
         DO L=1,LMAX
            LP=1+MOD(L-1+LMAX/2,LMAX)
            DO K=2,KMAX1
               IF(L.LE.LMAX/2) U(2,K,L)=zero !dkb- what's this for?
               U(1,K,L)    = -U(3,K,LP)
               W(1,K,L)    = W(2,K,LP)
               OMEGA(1,K,L)= OMEGA(2,K,LP)
               JN(1,K,L)   = JN(2,K,LP)
            ENDDO
         ENDDO
!$OMP END DO
      endif

C...SET VELOCITES BELOW THE EQUATORIAL PLANE.

      if (jmin.gt.2) then
         jstart=jmin
      else
         jstart=1
      endif

!$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        DO J=jstart,JMAX1
          U(J,1,L)    = U(J,2,L)
          W(J,2,L)    = zero
          W(J,1,L)    = -W(J,3,L)
          OMEGA(J,1,L)= OMEGA(J,2,L)
          JN(J,1,L)   = JN(J,2,L)
        ENDDO
      ENDDO
!$OMP END DO


      RETURN
      END
 
