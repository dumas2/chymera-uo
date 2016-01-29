C*******************************************************************************
C...ARTIFICIAL VISCOSITY FROM NORMAN AND WINKLER 1986 ARH.
C
C...VISCOSITY COEFFICIENTS
C
      SUBROUTINE AVISC
      implicit real*8(a-h,o-z)
#include "hydroparam.h"
#include "globals.h"
#include "units.h"
      integer jstart

      DR=ROF3N
      DZ=ZOF3N
      RHOLMTAV=DEN*phylim
      
      if (jmin.gt.2) then
         jstart=jmin
      else
         jstart=2
      endif   
      

!$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LM=L-1
        LP=L+1
        IF(L.EQ.1) LM=LMAX
        IF(L.EQ.LMAX) LP=1
        DO K=2,KMAX
          DO J=jstart,JMAX
            IF ((RHO(J,K,L).LE.RHOLMTAV)) THEN
              QRR(J,K,L)=0.d0
              QZZ(J,K,L)=0.D0
              QTT(J,K,L)=0.D0
            ELSE
              DELVR=(U(J+1,K,L)-U(J,K,L))
              DELVZ=(W(J,K+1,L)-W(J,K,L))
              DELVT=0.5d0*(OMEGA(J,K,LP)-OMEGA(J,K,LM))*RHF(J)

              QRR(J,K,L)=(CQ*DELVR*MIN(DELVR,0.d0))  
              QZZ(J,K,L)=(CQ*DELVZ*MIN(DELVZ,0.d0))          
              QTT(J,K,L)=(CQ*DELVT*MIN(DELVT,0.d0))
            END IF
               
          END DO
          QRR(2,K,L)=0.d0
        END DO
      END DO
!$OMP END DO

      RETURN
      END
