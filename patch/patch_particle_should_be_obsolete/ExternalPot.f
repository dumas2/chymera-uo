      subroutine ExternalPot()
      implicit none
#include "hydroparam.h"
#include "globals.h"
#include "units.h"

      integer::J,K,L

#if WIGGLE>0

!$OMP DO SCHEDULE(STATIC) 
            do l=1,lmax
                do k=1,kmax2
                  do j=1,jmax2
                     starphi(j,k,l) = -(mass_star+tmassacc)
     &                       /(sqrt(rhf(J)*rhf(J)+rpstar*rpstar-two
     & *rpstar*rhf(J)*cos((dble(L)-half)*dtheta-phi_star)
     & +zhf(K)*zhf(K)+gsoft**2))
                     phi(j,k,l)= phi(j,k,l)+starphi(j,k,l) ! ELIMINATED diskphi. ACB.
                  enddo
               enddo
            enddo
!$OMP END DO

#else

!$OMP DO SCHEDULE(STATIC) 
            do l=1,lmax
               do k=1,kmax2
                  do j=1,jmax2
                     if(tmassacc.gt.zero) then
C                    add potential of accreted mass 
                        tinphi(j,k,l) = -(tmassacc)
     &                       /(sqrt(RHF(j)*RHF(j)+zHF(k)*zHF(k)         &
     &                       +gsoft**2))
                     else
                        tinphi(j,k,l)=zero
                     endif
                     phi(j,k,l)= phi(j,k,l)+starphi(j,k,l)!+ ! ELIMINATED diskphi. ACB.
!     &                    tinphi(j,k,l)
                  enddo
               enddo
            enddo
!$OMP END DO
#endif

      return
      end subroutine ExternalPot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ExternalPotInit()
      implicit none
#include "hydroparam.h"
#include "globals.h"
#include "units.h"

      integer::J,K,L

#if WIGGLE>0

!$OMP DO SCHEDULE(STATIC) 
            do l=1,lmax
                do k=1,kmax2
                  do j=1,jmax2
                     starphi(j,k,l) = -(mass_star+tmassacc)
     &                       /(sqrt(rhf(J)*rhf(J)+rpstar*rpstar-two
     & *rpstar*rhf(J)*cos((dble(L)-half)*dtheta-phi_star)
     & +zhf(K)*zhf(K)+gsoft**2))
                     phi(j,k,l)= phi(j,k,l)+starphi(j,k,l) ! ELIMINATED diskphi. ACB.
                  enddo
               enddo
            enddo
!$OMP END DO

#else

!$OMP DO SCHEDULE(STATIC) 
            do l=1,lmax
               do k=1,kmax2
                  do j=1,jmax2
                     starphi(j,k,l)=-mass_star
     &                             / sqrt(rhf(J)**2+zhf(K)**2+gsoft**2)
                     if(tmassacc.gt.zero) then
C                    add potential of accreted mass 
                        tinphi(j,k,l) = -(tmassacc)
     &                       /(sqrt(RHF(j)*RHF(j)+zHF(k)*zHF(k)         &
     &                       +gsoft**2))
                     else
                        tinphi(j,k,l)=zero
                     endif
                     phi(j,k,l)= phi(j,k,l)+starphi(j,k,l)!+ ! ELIMINATED diskphi. ACB.
!     &                    tinphi(j,k,l)
                  enddo
               enddo
            enddo
!$OMP END DO
#endif

      return
      end subroutine ExternalPotInit
