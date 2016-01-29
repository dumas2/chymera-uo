      subroutine State() ! this routine calculates the pressure
        implicit none
     
#include "hydroparam.h"
#include "globals.h"
#include "units.h"

        integer JSTART,J,K,L

        if (jmin.gt.2) then
           jstart=jmin2
        else   
           jstart=1
        endif

      call TempFind()
   
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
         DO K=1,KMAX2
            DO J=jstart,jmax2

              p(J,K,L) = max(bkmpcgs*rho(J,K,L)*tempk(J,K,L)*rhoconv
     &            / (muc*pconv),poly_constant(J,K,L)*rho(J,K,L)**2)

            end do
         ENDDO
      ENDDO
C$OMP END DO

      RETURN
      END
 

