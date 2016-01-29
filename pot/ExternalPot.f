      subroutine ExternalPot(ITSTEP)
      implicit none
#include "hydroparam.h"
#include "globals.h"
#include "units.h"

      integer::J,K,L,ITSTEP
      real*8::thisMass,omegaBin2,rcom,theta,x,y
      real*8::rBin,thetaBin,thisTime,omegaBin,massDisk,muMass,rdist
      thisTime=time

#if BINARY>0
      call getBinaryPosition(thisTime,rBin,thetaBin,omegaBin,
     &   massDisk,muMass)
      omegaBin2=omegabin**2
      rcom=rbin-rbin*mass_star/muMass
#endif
             indirectx = rpstar*cos(phi_star)
             indirecty = rpstar*sin(phi_star)

#if INDIRECT>0
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:indirectx,indirecty)
c Then calculate the position of the star is such that the center of
c mass of the system remains zero.
      do L=1,LMAX
        do K=2,KMAX1
          do J=JMIN,JMAX1
             thisMass=rho(J,K,L)*rhf(J)*rof3n*zof3n*two*dtheta ! take into account both sides

             indirectx=indirectx+thisMass*cos((dble(L-1)+half)*dtheta)
     &                *rhf(J)
             indirecty=indirecty+thisMass*sin((dble(L-1)+half)*dtheta)
     &                *rhf(J)


          enddo
        enddo
      enddo
!$OMP ENDDO
             indirectx = -indirectx/mass_star
             indirecty = -indirecty/mass_star
#endif

!$OMP DO SCHEDULE(STATIC) PRIVATE(theta,rdist,x,y)
            do l=1,lmax
       theta=(dble(L)-half)*dtheta
                do k=1,kmax2
                  do j=1,jmax2
                     starphi(j,k,l) = -(mass_star+tmassacc)
     & /sqrt(rhf(J)*rhf(J)+zhf(K)*zhf(K) + indirectx*indirectx
     &  + indirecty*indirecty - two*rhf(J)
     & *(cos((dble(L-1)+half)*dtheta)*indirectx 
     & + sin((dble(L-1)+half)*dtheta)*indirecty))    
                     phi(j,k,l)= phi(j,k,l)+starphi(j,k,l) ! ELIMINATED diskphi. ACB.
                    
#if BINARY>0
       !!!!!!!!!!!!!!
       !! Add Binary!
       !!!!!!!!!!!!!!

       rdist=sqrt(rBin*rBin+rhf(J)*rhf(J)-two*rBin*rhf(J)*
     &       cos( theta- pi ) +
     &       zhf(K)*zhf(K) )

       x=rhf(J)*cos( theta )
       y=rhf(J)*sin( theta )

       phi(J,K,L)=phi(J,K,L)-massBin/rdist-half*omegaBin2*
     &            ( (x+rcom)**2+y**2 )
#endif

c#if INDIRECT>0
c       !!!!!!!!!!!!!!!!
c       !! Add Indirect!
c       !!!!!!!!!!!!!!!!
c                     phi(J,K,L)=phi(J,K,L)
c     &+indirectx*rhf(J)*cos((dble(L-1)+half)*dtheta)
c     &+indirecty*rhf(J)*sin((dble(L-1)+half)*dtheta)


c#endif 

                  enddo
               enddo
            enddo
!$OMP END DO
      return
      end subroutine ExternalPot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ExternalPotInit() ! this should be obsolete now
      implicit none
#include "hydroparam.h"
#include "globals.h"
#include "units.h"

      real*8::rBin,thetaBin,thisTime,omegaBin,massDisk,muMass,theta
      integer::J,K,L
      real*8::thisMass,omegaBin2,rcom,x,y,rdist

#if BINARY>0
      call getBinaryPosition(thisTime,rBin,thetaBin,omegaBin,
     &   massDisk,muMass)
      omegaBin2=omegaBin**2
      rcom=rBin-rBin*(mass_star)/muMass
#endif

#if INDIRECT>0
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:indirectx,indirecty) 
!$OMP&PRIVATE(thisMass)
      do L=1,LMAX
        do K=2,KMAX1
          do J=JMIN,JMAX1
             thisMass=rho(J,K,L)*rhf(J)*rof3n*zof3n*two*dtheta ! take into account both sides
             indirectx=indirectx+thisMass*cos((dble(L-1)+half)*dtheta)
     &                *rhf(J)/(sqrt(rhf(J)*rhf(J)+zhf(K)*zhf(K)))**3
             indirecty=indirecty+thisMass*sin((dble(L-1)+half)*dtheta)
     &                *rhf(J)/(sqrt(rhf(J)*rhf(J)+zhf(K)*zhf(K)))**3
          enddo
        enddo
      enddo
!$OMP ENDDO
#endif

 
!$OMP DO SCHEDULE(STATIC) private(rdist,x,y)
            do l=1,lmax
       theta=(dble(L)-half)*dtheta
                do k=1,kmax2
                  do j=1,jmax2
                     starphi(j,k,l) = -(mass_star+tmassacc)
     &                       /(sqrt(rhf(J)*rhf(J)+zhf(K)*zhf(K)))
                     phi(j,k,l)= phi(j,k,l)+starphi(j,k,l) ! ELIMINATED diskphi. ACB.
#if BINARY>0
       !!!!!!!!!!!!!!
       !! Add Binary!
       !!!!!!!!!!!!!!

       rdist=sqrt(rBin*rBin+rhf(J)*rhf(J)-two*rBin*rhf(J)*
     &       cos( theta- pi ) +
     &       zhf(K)*zhf(K) )

       x=rhf(J)*cos( theta )
       y=rhf(J)*sin( theta )

       phi(J,K,L)=phi(J,K,L)-massBin/rdist-half*omegaBin2*
     &            ( (x+rcom)**2+y**2 )
#endif

#if INDIRECT>0
       !!!!!!!!!!!!!!!!
       !! Add Indirect!
       !!!!!!!!!!!!!!!!
                     phi(J,K,L)=phi(J,K,L)
     &+indirectx*rhf(J)*cos((dble(L-1)+half)*dtheta)
     &+indirecty*rhf(J)*sin((dble(L-1)+half)*dtheta)
#endif
                  enddo
               enddo
            enddo
!$OMP END DO

      return
      end subroutine ExternalPotInit
