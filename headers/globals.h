!***********************************************************************
!     
      logical radtrigger

      real*8 :: vx,vy,fx,fy,rpstar,phi_star
      COMMON /movingstar/ vx,vy,fx,fy,rpstar,phi_star

      integer :: KWFW
      real*8 dtheta, pi, grav, bgden, gsoft
      common /blok6/dtheta, pi, grav, bgden, gsoft, KWFW

      real*8 etotfl, efl, eflufftot, gamma1
      COMMON /etally/ etotfl, efl, eflufftot,gamma1(JMAX2,KMAX2,LMAX)

      integer HHIT, HCOUNT, LO, CCOUNTER, HCOUNTER
      common /kant/radtrigger,HHIT,HCOUNT,CCOUNTER,HCOUNTER,LO

      REAL*8 KONST,NPRIME,xn,gamma,toverw
      COMMON /PTROPE/XN,GAMMA,KONST,NPRIME,TOVERW

      real*8 rcloud,constp,delt,bdytem,den,time,cormas,epscen
      COMMON /BLOK7/RCLOUD,CONSTP,DELT,BDYTEM,DEN,TIME,CORMAS,epscen 

      real*8 rholmt,epslmt,dumlmt,sound,rholmt_p,rhoa
      COMMON /RELIMITS/RHOLMT,EPSLMT,DUMLMT,sound,rholmt_p,rhoa

      integer jin,itype
      real*8 tmassini,tmassadd,tmassout,tmassacc,mass_star, mdot
      real*8 starphi,tinphi
      COMMON /GAP/starphi(jmax2,kmax2,lmax),                            &
     &     tinphi(jmax2,kmax2,lmax),                                    &
     &     tmassini,tmassadd,tmassout,tmassacc, mdot, mass_star,        &
     &     jin,itype


      REAL*8 JN,s,t,a,u,w,omega
      COMMON /EOM/                                                      &
     &     S(JMAX2,KMAX2,LMAX),                                         &
     &     T(JMAX2,KMAX2,LMAX),                                         &
     &     A(JMAX2,KMAX2,LMAX),                                         &
     &     U(JMAX2,KMAX2,LMAX),                                         &
     &     W(JMAX2,KMAX2,LMAX),                                         &
     &     JN(JMAX2,KMAX2,LMAX),                                        &
     &     OMEGA(JMAX2,KMAX2,LMAX)

      real*8 p,cv,eps,poly_constant
      real*8 phi,rho,rhotot,indirectx,indirecty
      COMMON /STATES/ENON,                                              &
     &     P(JMAX2,KMAX2,LMAX),                                         &
     &     CV(JMAX2,KMAX2,LMAX),                                        &
     &     EPS(JMAX2,KMAX2,LMAX),                                       &
     &     poly_constant(JMAX2,KMAX2,LMAX)


      COMMON /POIS/ indirectx,indirecty,                                &
     &   PHI(POT3JMAX2,POT3KMAX2,LMAX),                                 &
     &   RHO(POT3JMAX2,POT3KMAX2,LMAX),                                 &
     &   RHOTOT(POT3JMAX2,POT3KMAX2,LMAX)

      real*8 Msyscgs,PKcgs,Tconv,Sconv,Dconv,Pconv,sigma,rhoconv,       &
     &       engconv,bkmpcode
      COMMON /CONVERT/                                                  &
     &     Msyscgs,PKcgs,Tconv,Sconv,                                   &
     &     dconv,Pconv,sigma,rhoconv,engconv,bkmpcode

      real*8 qrr,qzz,qtt,hgamma,cs,totheat
      COMMON /AVIS/                                                     &
     &     QRR(JMAX2,KMAX2,LMAX),                                       &
     &     QZZ(JMAX2,KMAX2,LMAX),                                       &
     &     QTT(JMAX2,KMAX2,LMAX),                                       &
     &     HGAMMA(JMAX2,KMAX2,LMAX),                                    &
     &     CS,totheat

      real*8 lambda,tau,TempK,TeffK,TphK,Surfcgs,divflux,radflux
      real*8 totcool,totdflux
      COMMON /COOLING/                                                  &
     &     LAMBDA(JMAX2,KMAX2,LMAX),                                    &
     &     TAU(JMAX2,KMAX2,LMAX,4),                                     &
     &     TEMPK(JMAX2,KMAX2,LMAX),                                     &
     &     TEFFK(JMAX2,LMAX),                                           &
     &     TPHK(JMAX2,LMAX),                                            &
     &     SURFCGS(JMAX2,LMAX),                                         & 
     &     DIVFLUX(JMAX2,KMAX2,LMAX),                                   &
     &     RADFLUX(JMAX2,KMAX2,LMAX,3),                                 &
     &     totcool,totdflux

      real*8 Igamma,totirr
      COMMON /IRRAD/                                                    &
     &     IGAMMA(JMAX2,KMAX2,LMAX),                                    &
     &     totirr
   
      real*8 xmmwtable,rostable,plktable,irrtable,scatable,ptab,ttab    
      COMMON /OPACITY/                                                  & 
     &     XMMWtable(ITABLE,ITABLE),                                    &
     &     ROStable(ITABLE,ITABLE),                                     &
     &     PLKtable(ITABLE,ITABLE),                                     &
     &     IRRtable(ITABLE,ITABLE),                                     &
     &     SCAtable(ITABLE,ITABLE),                                     &
     &     Ptab(ITABLE),                                                &
     &     Ttab(ITABLE)

      real*8 r,z,rhf,zhf,rof3n,zof3n,enon
      integer jreq,kzpol
      COMMON /GRID/JREQ,KZPOL,                                          &
     &     R(pot3jmax2),                                                &
     &     Z(POT3KMAX2),                                                &
     &     RHF(POT3JMAX2),                                              &
     &     ZHF(POT3KMAX2),                                              &
     &     ROF3N,ZOF3N

     
      integer :: KFITA
      real*8 :: temporary,dsdt,sfunc,l_tau_z,dtau_z,                    &
     &          intensity_in_z,intensity_z,ddsdtt,                      &
     &          int_temp,init_int_in 
      common /intensity/ KFITA(JMAX+2,LMAX),                            &
     &          temporary(JMAX+2,KMAX+2,LMAX),                          &
     &          dsdt(JMAX+2,KMAX+2,LMAX),                               &
     &          sfunc(JMAX+2,KMAX+2,LMAX),                              &
     &          l_tau_z(JMAX+2,KMAX+2,LMAX),                            &
     &          dtau_z(JMAX+2,KMAX+2,LMAX),                             &
     &          intensity_in_z(JMAX+2,KMAX+2,LMAX),                     &
     &          intensity_z(JMAX+2,KMAX+2,LMAX),                        &
     &          ddsdtt(JMAX+2,KMAX+2,LMAX),                             &
     &          int_temp(JMAX+2,KMAX+2,LMAX),                           &
     &          init_int_in(JMAX+2,LMAX)



      real*8 :: temptable,engtable,gammatable,muc
      common /engtables/temptable(TTABLE),engtable(TTABLE),             &
     &        gammatable(TTABLE),muc

      real*8 :: a1newr,a1newz
      common/obsolete/a1newr,a1newz

#if PASSIVE>0
      real*8 :: passflux
      common/passivearray/passflux(JMAX+2,KMAX+2,LMAX,PAS)
#endif

