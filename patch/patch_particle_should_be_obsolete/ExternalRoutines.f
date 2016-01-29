!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THIS IS AN ACCRETION ROUTINE
      subroutine TsigConstant(tmass)
      implicit none
#include "hydroparam.h"
#include "units.h"
#include "globals.h"
      integer::J,K,L
      integer::JRIN=32,JROUT=152
      real*8::limiter,grow,ran4,nrho,madd,avgNrho,madd1,tmass
      limiter=den*phylim

! this is stupid, but I am keeping it in for stupid reasons.  
! you should change it anyway for your own use.
      grow=1d-1*.1592d0/ten/pi/
     &             (r(JROUT)**1.5-r(int(JRIN*1.5))**1.5)*delt
      grow=grow/zof3n/1.7d0*two

      madd=zero
      madd1=zero
      avgNrho=zero
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nrho) 
!$OMP&REDUCTION(+:madd,avgNrho)
      K=KMAX*3/4
!$OMP DO SCHEDULE(STATIC)
      do L=1,LMAX,2
        do J=JRIN,JROUT
           nrho=grow/rhf(J)**1.5d0
           avgNrho = avgNrho + nrho
           madd=madd+nrho*rhf(J)*zof3n*rof3n*dtheta*two
           eps(J,K,L)=eps(J,K,L)+nrho*bkmpcode*1.5d0/muc*
     &                TenvK/Tconv
!           t(J,K,L) = zero
           t(J,K,L)  = t(J,K,L)-
     &              nrho*sqrt(two*mass_star/rhf(J))
           rho(J,K,L)=rho(J,K,L)+nrho
           a(J,K,L)  =rho(J,K,L)*sqrt(mass_star/(rhf(J)**2+gsoft**2)
     &               *rhf(J))*rhf(J)
           eps(J,K,L)=max(eps(J,K,L),epslmt)
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      print *, "Actual Mdot is", madd/delt/.1592d0

      return
      end subroutine TsigConstant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine simpleCool
      implicit none
#include "hydroparam.h"
#include "units.h"
#include "globals.h"
      integer::J,K,L
      real*8::limiter,dtau,area,volume,tau_fac,coolTime,etenv
      real*8::ttenv

      REAL*8 Oross(jmax2,kmax2,lmax),Oplck(jmax2,kmax2,lmax)
     &     ,Oabs(jmax2,kmax2,lmax),Otot(jmax2,kmax2,lmax),
     &      oplck_env(JMAX2,KMAX2,LMAX),sum

      COMMON /COOLINGSHARED/Oross,Oplck,Oabs,Otot,sum,oplck_env

      limiter=den*phylim   

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(dtau,area,volume,tau_fac,
!$OMP&coolTime,J,K,L,ttenv,etenv)
!$OMP DO SCHEDULE(STATIC)
      do L=1,LMAX
       do K=2,KMAX1
         do J=JCOOL,JMAX1
           ttenv=tenvk/tconv/sqrt(rhf(J))+10d0/tconv
           etenv = ttenv*bkmpcode/(muc*tconv)
           if (rho(J,K,L)>limiter) then

           area=two*(rhf(J)*rof3n*dtheta+rof3n*zof3n)+(r(J)+r(J+1)
     &         )*dtheta*zof3n
           volume=rhf(J)*rof3n*zof3n*dtheta
           dtau=sconv*rho(J,K,L)*oross(J,K,L)*volume**(one/three)
!           dtau=tau(J,K,L,1)
!           dtau=sconv*rho(J,K,L)*( (one-exp(-two*dtau))*oross(J,K,L)
!     &         +oplck(J,K,L)*exp(-two*dtau) )*volume**(one/three)
           dtau=max(dtau,1d-15)
           
           tau_fac=dtau+one/dtau
           divflux(J,K,L)=area*sigma*((TempK(J,K,L)/Tconv)**4
     &                   -(ttenv/Tconv)**4)/(volume*tau_fac)
           coolTime = eps(J,K,L)/abs(divflux(J,K,L))
           divflux(J,K,L)=( eps(J,K,L)-etenv*rho(J,K,L)
     &                   /(gamma1(J,K,L)-one) )/delt
     &                   *(one-exp(-(delt/coolTime)**2)) 
     &                   + divflux(J,K,L)*exp(-(delt/coolTime)**2)
          
           else

           divflux(J,K,L)=zero

           endif 
         enddo
        enddo
       enddo
!$OMP ENDDO
!$OMP END PARALLEL

      return
      end subroutine simpleCool




