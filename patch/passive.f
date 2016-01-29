! THIS IS AN EXAMPLE OF HOW TO INTIALZE THE PASSIVE ARRAY.
! THIS EXAMPLE JUST SETS IT TO THE DENSITY ARRAY.
      subroutine set_passive()
      implicit none
#include "hydroparam.h"
#include "globals.h"
#include "units.h"
      integer::J,K,L,I,JMAXC,KMAXC,LMAXC,PASC
      real*8::pf(JMAX2,KMAX2,LMAX,4),timec

#if VERBOSE>0
      print *, " INTIALIZING PASSIVE ARRAY "
#endif

#if PASSIVE_INIT>0
!$OMP PARALLEL DO SCHEDULE(STATIC)
       do L=1,LMAX
        do K=1,KMAX2
         do J=1,JMAX2
         enddo
        enddo
       enddo
!$OMP END PARALLEL DO
#else
      open(unit=472,file="passive.restart",form="UNFORMATTED")
      read(472)JMAXC,KMAXC,LMAXC,PASC
      if(JMAX/=JMAXC.or.KMAX/=KMAXC.or.LMAX/=LMAXC.or.PASC/=4)then
        print *, "HEADER OF passive.restart DOES NOT MATCH WHAT IS IN"
        print *, "hydroparam.h.  STOPPING"
        stop
      endif
      read(472)pf
      read(472)timec
      close(472)
#if VERBOSE>0
      print *,"HEADER OF passive.restart:",JMAXC,KMAXC,LMAXC,PASC
      print *,"TIME FOR passive.restart:",timec
#endif
#endif
      return
      end
