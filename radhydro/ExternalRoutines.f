!!!BASIC UTILITIES
      module utilities
      implicit none
      real*8,parameter::tol_kep=1d-3
       
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is essentially atan2, but with some additional control for arg=0,0
      real*8 function getAngle(x,y) 
       real*8,parameter::pi=3.141592653589793d0
       real*8,intent(in):: x,y
       if(x==0d0)then
        if(y==0d0)then
          getAngle=0d0
        elseif(y>0d0)then
          getAngle=pi*0.5d0
        else
          getAngle=pi*1.5d0
        endif 
       else
        getAngle=atan(y/x)
        if (y>0d0)then
         if(x<0d0)getAngle=getAngle+pi 
        else
         if(x<0d0)then
          getAngle=getAngle+pi
         else
          getAngle=getAngle+2.d0*pi
         endif
        endif
       endif  
       return 
      end function getAngle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!
      ! radial part of Kepler
      !!!!!!!!!!!!!!!!!!!!!!!

      real*8 function radialCoorKep(a,ee,eAnom)
       real*8,intent(in)::a,ee,eAnom
       radialCoorKep=a*(1d0-ee**2)/(1d0+ee*cos(eAnom))
       return
      end function radialCoorKep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!
      ! Solve Kepler's problem
      !!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine trueAnomaly(time,mu,a,ee,theta,eAnom)
       real*8,parameter::pi=3.141592653589793d0
       real*8,intent(in)::time,mu,a,ee
       real*8,intent(out)::theta,eAnom
       real*8::r,calc,errMag,err,orbits,val

       val=time*sqrt(mu/a**3)
       orbits=dble(int(val/(2.d0*pi)))
       !print *,"Anomaly r,val,orbits,mu", r,val,orbits,mu
       val=val-orbits*2d0*pi
       if(val<1d-32)then ! don't do work unless we must
        theta=0d0;eAnom=0d0;return
       endif

       eAnom=1d0
       errMag=1d0
       do while(errMag>tol_kep)
        calc=eAnom-ee*sin(eAnom)
        err=(val-calc)/val
        errMag=abs(err)
        eAnom=eAnom+0.6d0*min(errMag,1d0)*eAnom*err/errMag
        !print *, eAnom,val,err,time
       end do
       theta=2d0*getAngle(sqrt(1d0-ee)*cos(eAnom/2d0),
     &  sqrt(1d0+ee)*sin(eAnom/2d0)) ! arg(x,y)
       return
      end subroutine trueAnomaly

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This routine calculates the position of a wide-orbit binary.
!!!
!!! For the position of the binary, only consider r and phi because
!!! we are still in mirror symmetry land for the disk. 
      subroutine getBinaryPosition(tBin,rBin,thetaBin,omegaBin,
     &  massDisk,muBin)
      use utilities
      implicit none
#include "hydroparam.h"
#include "units.h"
#include "globals.h"
      real*8,intent(in)::tBin,massDisk
      real*8,intent(out)::rBin,thetaBin,omegaBin,muBin
      real*8::eAnom

! Currently set up for Roche potential in a rotating frame

      muBin=massBin+mass_star
      omegaBin=sqrt((muBin)/(aBin**3))
      rbin=abin
 
      return

! if you want to do something fancier such as follow the binary in the 
! frame of the planet, then the above return can be removed.

      call trueAnomaly(tBin,muBin,rbin,eBin,thetaBin,eAnom)
      rBin=radialCoorKep(abin,eBin,eAnom)
      !!print *, "GetBINARY",tBin,muBin,omegaBin,rBin,thetaBin,massDisk
      return
      end subroutine

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
! this should only be thought of as a template.  Perhaps a bad one.

      print *, "You are attempting to run an accretion routine."
      print *, "You must be aware of the consequences of doing so."
      print *, "The routine is designed to dump mass on the grid from"
      print  *,"the upper boundary. You will need to remove the STOP"
      print *, "command in the source of "
      print *, "ExternalRoutines.f: TsigConstant(). You will likely"
      print *, "Need to edit the code for your purpose."

      stop
      
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
! radiative cooling routine.  Very simple and easy.  Don't fool yourself
! into thinking this is terribly accurate. It gets the job done, and it is stable.
! Note that ttenv is user specified.
      use eos, only : get_gamma_from_tk
      implicit none
#include "hydroparam.h"
#include "units.h"
#include "globals.h"
      integer::J,K,L
      real*8::limiter,dtau,area,volume,tau_fac,coolTime,etenv
      real*8::ttenv,opacfac,fluxfac,ds,tacc,mmw,gam

      REAL*8 Oross(jmax2,kmax2,lmax),Oplck(jmax2,kmax2,lmax)
     &     ,Oabs(jmax2,kmax2,lmax),Otot(jmax2,kmax2,lmax),
     &      oplck_env(JMAX2,KMAX2,LMAX),sum

      COMMON /COOLINGSHARED/Oross,Oplck,Oabs,Otot,sum,oplck_env

      logical :: include_accretion=.false.

      limiter=den*phylim   

      ds=rof3n
      opacfac=sqrt(3.)/2.
      fluxfac=8./sqrt(3.)

      if(include_accretion)then
      if (mdot <= zero ) then
              
      tacc=SQRT(SQRT(((Mstar+tmassacc)*(-mdot))   !referenced at 1AU
     &   /((ten+six)*pi*sigma)))
      else 
      print*, 'WARNING: positive mdot expecting neg value'
       
      endif 

!$OMP MASTER
      print*, 'Lacc (Lsun)', (mstar+tmassacc)*(-mdot)/(rstar*two)
     &*engconv*msuncgs/5d6/3.8d33,
     & tmassacc, -mdot, tacc

!$OMP END MASTER
      else
        tacc=zero
      endif 
 
!$OMP DO SCHEDULE(STATIC) private(mmw,gam,ttenv,etenv)
      do L=1,LMAX
       do K=2,KMAX1
         do J=JCOOL,JMAX1
!!! ttenv is user specified.  
           ttenv=((tenvk/tconv/(rhf(J))**(.75)+10.d0/tconv)**4
     &          +(tacc/rhf(J)**.75)**4)**.25
           call get_gamma_from_tk(etenv,rho(J,K,L),ttenv,mmw,gam)
           if (rho(J,K,L)>limiter) then

           dtau=max(tau(J,K,L,1),1d-15) !vertically integrated Rossland mean opacity
           
           tau_fac=dtau+one/dtau
           divflux(J,K,L)=fluxfac*sigma*((TempK(J,K,L)/Tconv)**4
     &                   -(ttenv/Tconv)**4)/(ds*tau_fac)
           coolTime = eps(J,K,L)/abs(divflux(J,K,L))
           divflux(J,K,L)=( eps(J,K,L)-etenv )/delt
     &                   *(one-exp(-(delt/coolTime)**2)) 
     &                   + divflux(J,K,L)*exp(-(delt/coolTime)**2)
          
           else

           tempk(J,K,L)=ttenv
           eps(J,K,L)=etenv
           divflux(J,K,L)=zero

           endif 
         enddo
        enddo
       enddo
!$OMP ENDDO

      return
      end subroutine simpleCool
