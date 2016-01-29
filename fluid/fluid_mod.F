# 1 "fluid_mod.F"
      module Fluid_Mod
      implicit  none
#include "units.h" 
#include "hydroparam.h"
! define grid coordinates and grid motion
      integer,parameter :: NPART=1000
      integer,parameter :: R_END=200,R_START=20,Z_TOP=3,TSTOP=10
      integer :: TSTEP
      real*8 :: dphi,dr,dz,delta_t
! define fluid element coordinates and element motion
      real*8,dimension(-1:TSTOP,0:NPART-1) :: r_e,z_e,angle_e
      real*8,dimension(         0:NPART-1) :: ro,zo,phio
      real*8,dimension(-1:TSTOP,0:NPART-1) :: vr_e,vz_e,omega_e,  &
     &eps_e,rho_e,temp_e,gamma1_e,tau_e
      real*8,dimension(1:lmax) :: angle
      logical,dimension(0:NPART-1) :: partskip
!***********************************************************************
!
      real*8 dtheta,pi,grav,bgden,KWFW
      COMMON /BLOK6/dtheta,pi,grav,bgden,KWFW
      real*8 etotfl,efl,eflufftot,gamma1
      COMMON /etally/etotfl,efl,eflufftot,gamma1(jmax2,kmax2,lmax)
      real*8 rholmt,epslmt,dumlmt,sound
      COMMON /RELIMITS/rholmt,epslmt,dumlmt,sound
      real*8 rcloud,constp,delt,bdytem,den,time,cormas,epscen
      COMMON /BLOK7/rcloud,constp,delt,bdytem,den,time,cormas,epscen
      real*8 JN,s,t,a,u,w,omega
      COMMON /EOM/s(jmax2,kmax2,lmax),t(jmax2,kmax2,lmax),a(jmax2,kmax2,&
     &lmax),u(jmax2,kmax2,lmax),w(jmax2,kmax2,lmax),JN(jmax2,kmax2,lmax)&
     &,omega(jmax2,kmax2,lmax)
      real*8 p,cv,eps
      real*8 phi,rho
      COMMON /STATES/enon,p(jmax2,kmax2,lmax),cv(jmax2,kmax2,lmax),     &
     &eps(jmax2,kmax2,lmax)
      COMMON /POIS/phi(pot3jmax2,pot3kmax2,lmax),rho(pot3jmax2,pot3kmax2&
     &,lmax)
      real*8 Msyscgs,PKcgs,Tconv,Sconv,Dconv,Pconv,sigma,rhoconv,       &
     &       engconv,bkmpcode
      COMMON /CONVERT/ Msyscgs,PKcgs,Tconv,Sconv,                       &
     &     dconv,Pconv,sigma,rhoconv,engconv,bkmpcode
      real*8 lambda,tau,TempK,TeffK,TphK,Surfcgs,divflux,radflux
      real*8 totcool,totdflux
      COMMON /COOLING/lambda(jmax2,kmax2,lmax),tau(jmax2,kmax2,lmax,4), &
     &TempK(jmax2,kmax2,lmax),TeffK(jmax2,lmax),TphK(jmax2,lmax),Surfcgs&
     &(jmax2,lmax),divflux(jmax2,kmax2,lmax),radflux(jmax2,kmax2,lmax,3)&
     &,totcool,totdflux
      real*8 r,z,rhf,zhf,rof3n,zof3n,enon
      integer jreq,kzpol
      COMMON /GRID/jreq,KZPOL,r(pot3jmax2),z(pot3kmax2),                &
     &rhf(pot3jmax2),zhf(pot3kmax2),rof3n,zof3n
      endmodule
