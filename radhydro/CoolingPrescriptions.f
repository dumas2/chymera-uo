      subroutine ConstantTcool()
      implicit none
#include "hydroparam.h"
#include "units.h"
#include "globals.h"
      integer::j,k,l
      real*8::tbgrnd4,limiter

      tbgrnd4=(tbgrnd/tconv)**four
      limiter = den*phylim

!$OMP DO SCHEDULE(STATIC)  
      do L=1,LMAX
        do J=JCOOL,JMAX1
          teffk(J,L)=zero
        enddo
      enddo
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)  
      do l=1,lmax
        do k=2,kmax
         do j=jcool,jmax

               if (rho(j,k,l).ge.limiter) then
                  lambda(j,k,l)=eps(j,k,l)/(cct*torp)
               else
                  lambda(J,K,L)=zero
               endif

! An effective cooling temperature.

               Teffk(J,L) = Teffk(J,L) + lambda(j,k,l)*zof3n/sigma

C     Weighted average for photospheric temperature (at tau=2/3)
               
               if ((tau(j,k,l,1).ge.(twothree)).and.(tau(j,k+1,l,1).lt
     &              .(twothree))) TphK(j,l)=ten**((log10(TempK(j,k,l))+
     &              log10(TempK(j,k+1,l)))/two)
               
            enddo
         enddo  
      enddo
!$OMP END DO 
!$OMP DO SCHEDULE(STATIC)  
      do L = 1, LMAX
       do J=2,JMAX1
         if (teffk(J,L)<tbgrnd4)teffk(J,L)=tbgrnd4 ! tbgrnd4 is in code units
         teffk(J,L)=sqrt(sqrt(teffk(J,L)))*tconv
       enddo
      enddo
!$OMP ENDDO NOWAIT
      return
      end subroutine ConstantTcool

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine TcoolOmega()
      implicit none
#include "hydroparam.h"
#include "units.h"
#include "globals.h"
      integer::j,k,l
      real*8::tbgrnd4

      tbgrnd4=(tbgrnd/tconv)**four

!$OMP DO SCHEDULE(STATIC)
      do L=1,LMAX
       do J=2,JMAX1
        teffk(J,L)=zero
       enddo
      enddo
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
      do L=1,LMAX
       do K=2,KMAX1
        do J=2,JMAX1
          if (rho(J,K,L)>=rholmt) then
            lambda(J,K,L)=eps(J,K,L)*omega(J,K,L)/cct
          else
            lambda(J,K,L)=zero
          endif
          teffk(J,L)=teffk(J,L)+lambda(J,K,L)*zof3n/sigma
        enddo
       enddo
      enddo
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
      do L=1,LMAX
       do J=2,JMAX1
        if (teffk(J,L)<tbgrnd4)teffk(J,L)=tbgrnd4
        teffk(J,L)=sqrt(sqrt(teffk(J,L)))*tconv
       enddo
      enddo
!$OMP ENDDO NO WAIT

      return
      end subroutine TcoolOmega
