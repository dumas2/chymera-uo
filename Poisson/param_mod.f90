module param_mod

 use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD

!Contains

!subroutine hydro_mod

!     This file includes the parameter statements for the 3dhyd.f code.
!     These can be modified as necessary to run different sized grids;
!     however, parameters starting with 'pot3' must be powers of 2.
!     
!     -rpl, Sep, 1996


! Do not currently support different resolutions in the l-direction.
      integer  jmax,jmax1,jmax2
      integer  kmax,kmax1,kmax2,ksplit
      integer  lmax,lmax2,lrjmax,lrkmax,lrlmax,lrjmax1
      integer  pot3jmax,pot3jmax1,pot3jmax2
      integer  pot3kmax,pot3kmax1,pot3kmax2
      integer  jmin,jmin1,jmin2

!  High resolution problem.
      parameter (jmax=32, jmax1=jmax+1, jmax2=jmax+2)
      parameter (kmax=16, kmax1=kmax+1, kmax2=kmax+2)
      parameter (lmax=32,lmax2=lmax/2)
      parameter(pot3jmax=32,pot3jmax1=pot3jmax+1,pot3jmax2=pot3jmax+2)
      parameter(pot3kmax=16,pot3kmax1=pot3kmax+1,pot3kmax2=pot3kmax+2)

!  Minimum radial grid point, for cutting out central star.
      parameter (jmin=4,jmin1=jmin-1,jmin2=jmin-2)

!  Other parameters.
      integer hj,hk,hj1,hj2,hk1,hk2,itable
      parameter (hj=32,hk=16,hj1=hj+1,hk1=hk+1,hj2=hj+2,hk2=hk+2)
      parameter (lrlmax=128,lrjmax=64,lrkmax=8,lrjmax1=lrjmax+1)
      parameter (itable=100)

      integer, parameter :: TTABLE=400
      integer MAXTHREADS
      parameter (MAXTHREADS=200)

      real*8,parameter::one=1d0,two=2d0,three=3d0,four=4d0,five=5d0,    &
     & six=6d0,seven=7d0,eight=8d0,nine=9d0,ten=10d0,zero=0d0,          &
     & twothree=0.666666666666667d0,half=0.5d0,quarter=0.25d0

!end subroutine hydro_mod

!subroutine read_mod

real    :: konst,xn,gamma  ! from fort.5
real    :: pindex,con2,rrr2,omcen,dencen,toverw,rof3n,zof3n,a1newz

integer :: jreq,kzpol
integer :: itstrt,itstop,idiag,isoadi,itype,nmodl,istor
real*8, parameter:: dr = 5.103179588499317d-1
!end subroutine read_mod

real    :: edif, elost, enew, grav , nprime
integer :: iprint, klocat

! Some common blocks from the boundary solver
      PARAMETER(JKM1=2*POT3JMAX+POT3KMAX-1)
 
      COMMON /BDY/ COSM(LMAX,10),SINM(LMAX,10),RBDY(JKM1),BDYCHK
      COMMON /BDY1/ JPOS(JKM1),KPOS(JKM1),JKMAX 
contains
subroutine dime

integer :: rank, numRanks
call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)

ksplit=kmax/numRanks
end subroutine dime
end module
