!    The following units: Mstar, Rdiskau, Tstar are model
!    dependent.  Must be set by the user.   
!    The rest are astronomical units.

!    Tstar is used in RadTran.f to calculate Igamma.

!    Mstar is used in io.f for normilization purposes to calculate Msyscgs,
!    as well as in hybrid to calculate the accretion time,
!    as well as in ExternalRoutines.f to calculate tacc.

!    Rstar is used to calculate accretion luminosity in ExternalRoutines.f
!    as well as in hybrid.f for the same purpose,
!    as well as in RadTran.f to calculate Igamma.

!    Rdiskau is used in io.f for normalization purposes.

      REAL*8 Mstar,Rstar,Rdiskau,Tstar,gridlim,phylim
      REAL*8 Msuncgs,Rsuncgs,AUcgs,sigmacgs,Gcgs,bkmpcgs
      PARAMETER (Mstar=1.d0,Rstar=1.d0,Rdiskau=1.d2,Tstar=2.d3)
      PARAMETER (Msuncgs=1.989d33, Rsuncgs=6.96d10, AUcgs=1.496d13)
      PARAMETER (sigmacgs=5.670d-5, Gcgs=6.672d-8, phylim=1.d-4) 
      PARAMETER (bkmpcgs=8.254d7,gridlim=1.d-10)
      real*8,parameter::psize=1d1/AUcgs,dust_to_gas=0.01d0,rhoacgs=3d0
      real*8,parameter::poly_factor=4d0
      real*8,parameter::massBin=1.,ebin=0.,abin=100.

!    Most of the following variables are used in the source routine only. 
!    torp is one ORP in code time units, defined as 2*pi/omega(200,2,1) 
!    of the starting model.
!    cq = 3 for AV on, 0 for AV off. 
!    irtype = 0 for no irr, 1 for irradiation.
!    ictype is the cooling type: 1 for constant Tcool, 2 for 
!    Edd grey atmosphere only, 3 for diffusion approx, 4 for diffusion
!    approx with shining atmosphere.
!    cct is only used when using type 1.  It is the cooling
!    time in ORPS.      
!    tcoollmt is the lower limit of cooling time in ORP.
!    tirrlmt  is the lower limit of irradiation time in ORP.
!    Tbgrnd is the lower limit of the temperature
!    irop is used to spread the irradiation to a few more cells by dividing 
!    the opacities by some number. 
!    jirr is the first first j zone to be irradiated.
!    jcool is the j zone at which cooling starts.  Disabled in most ictypes.
!    fk = metallicity factor: fraction of standard metallicity to keep.
!    tenvk = constant envelop temperature
!    SETUNITS = 0 for polytropic units (G=K=M=1) and 1 for Aaron's units
!       (G=1, 1 R = 1 AU, 1 M = Msun)
!    Use H2STAT to select what type of mixture you want.
!    0 = pure para
!    1 = pure ortho -- this is not really physical for astro applications
!    2 = set mixture
!    3 = equilibrium
!    -1 = single gamma law
!    Be sure to set ac and bc for mixture (ac: para component, bc: ortho component)
!    also pick a metallicity by adjusting xabun, yabun, and zabun.


      REAL*8 torp,cct,Tbgrnd,tcoollmt,theatlmt,cq,irop,fk,tenvk,amp0
      real*8, parameter :: fgsoft=0d0 ! set softening for star
      real*8, parameter :: ac = 1.d0, bc = 3.d0
      real*8, parameter :: xabun = .73d0,yabun=.25d0,zabun=.02d0
      real*8, parameter :: tk_eos_cutoff=1d3
      real*8, parameter :: tk_bgrnd=3d0
      real*8, parameter :: dtk_eos=3d1
      real*8, parameter :: mp=1.67d-24
      integer, parameter :: NEOS_RHO=500
      real*8,parameter::rho_eos_high=1d-4
      real*8,parameter::rho_eos_low=1d-15

      
      INTEGER ictype,irtype,jcool,jirr,SETUNITS,H2STAT
      integer,parameter::JCONST=80
      PARAMETER(torp=1.446d2,cct=5d0,cq=3d0,tcoollmt=5d1,theatlmt=5d1)   
      PARAMETER(ictype=1,irtype=0,jcool=2)
      PARAMETER(jirr=11,irop=1.0,Tbgrnd=3d0)
      PARAMETER(tenvk=150d0,fk=1d0)
      PARAMETER(SETUNITS=0,H2STAT=-1) 
      logical,parameter::find_eos_rho=.true.
      PARAMETER(amp0=0.005) ! amplitude of initial, random perturbation

#if EXTERNAL_POT>0
      logical, parameter :: external_pot = .true.
#else
      logical, parameter :: external_pot = .false.
#endif
#if FLUID_RESTART>0
      logical, parameter :: restart_fluid = .true.
#else
      logical, parameter :: restart_fluid = .false.
#endif
#if WIGGLE_RESTART>0
      logical, parameter :: restart_wiggle = .true.
#else
      logical, parameter :: restart_wiggle = .false.
#endif
#if ROTATING>0
      real*8, parameter :: omega_frame=5d-2  
#endif

!     leave a space after the file name, or set the exact character length.
!     This will help to ensure that all compilers do the right thing.
      character*11, parameter :: mmwfile='mmw_new.dat'
      character*33, parameter :: rosstablefile                          &
     &                          ='rosseland.abpoll.p3p5.amax1.dat'
      character*30, parameter :: plancktablefile                        &
     &                          ='planck.abpoll.p3p5.amax1.dat'
      character*37, parameter :: irrtablefile                           &
     &                          ='planck.abpoll.p3p5.amax1.ts4000.dat'
      character*37, parameter :: scatablefile                           &
     &                          ='planck.abpoll.p3p5.amax1.ts4000.dat'

