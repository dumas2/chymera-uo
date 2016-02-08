!     This file includes the parameter statements for the 3dhyd.f code.
!     These can be modified as necessary to run different sized grids;
!     however, parameters starting with 'pot3' must be powers of 2.
!     
!     -rpl, Sep, 1996

! Do not currently support different resolutions in the l-direction.
      integer  jmax,jmax1,jmax2
      integer  kmax,kmax1,kmax2
      integer  lmax,lmax2,lrjmax,lrkmax,lrlmax,lrjmax1
      integer  pot3jmax,pot3jmax1,pot3jmax2
      integer  pot3kmax,pot3kmax1,pot3kmax2
      integer  jmin,jmin1,jmin2

!  High resolution problem.
      parameter (jmax=256, jmax1=jmax+1, jmax2=jmax+2)
      parameter (kmax=128, kmax1=kmax+1, kmax2=kmax+2)
      parameter (lmax=256,lmax2=lmax/2)
      parameter(pot3jmax=256,pot3jmax1=pot3jmax+1,pot3jmax2=pot3jmax+2)
      parameter(pot3kmax=128,pot3kmax1=pot3kmax+1,pot3kmax2=pot3kmax+2)

!  Minimum radial grid point, for cutting out central star.
      parameter (jmin=2,jmin1=jmin-1,jmin2=jmin-2)

!  Other parameters.
      integer hj,hk,hj1,hj2,hk1,hk2,itable
      parameter (hj=256,hk=256,hj1=hj+1,hk1=hk+1,hj2=hj+2,hk2=hk+2)
      parameter (lrlmax=128,lrjmax=64,lrkmax=8,lrjmax1=lrjmax+1)
      parameter (itable=100)

      integer, parameter :: TTABLE=400
      integer MAXTHREADS
      parameter (MAXTHREADS=200)

! change the following for passive array size modification.
#if PASSIVE>0
      integer, parameter :: PAS=PASSIVE
#endif

      real*8,parameter::one=1d0,two=2d0,three=3d0,four=4d0,five=5d0,    &
     & six=6d0,seven=7d0,eight=8d0,nine=9d0,ten=10d0,zero=0d0,          &
     & twothree=0.666666666666667d0,half=0.5d0,quarter=0.25d0
