#undef USE_NETCDF

Program PoissonRelax
!=============================================================================
! Name        : poissonRelax.F90
! Version     : 1.0
! Description : Solves the Poisson equation in two dimensions.
!
! Method      :
!
!   Iterative (over)relaxation. 
!==============================================================================
use MultiGrid, only : Relax, Residual, Prolongate, Restrict, RestrictRho
use io       , only : readDensity, writeData
#ifdef USE_NETCDF
use io       , only : initFileNetCDF, writeDataNetCDF
#endif
use MPI_F08  , only : MPI_Init, MPI_Finalize, MPI_Barrier
use MPI_F08  , only : MPI_Send, MPI_Recv, MPI_DOUBLE_PRECISION, MPI_Status
use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
implicit none

real   , parameter :: tol =  1e-5

integer, parameter :: Nj  = 256    ! number of r interior elements
integer, parameter :: NTk =  64    ! number of z interior elements total
integer, parameter :: Nl  =  1     ! number of phi interior elements
integer            :: Nk           ! number of z interior elements per rank in z
integer, parameter :: m   =  1     ! consider m=1 Fourier mode for testing


integer   :: i,k,mn,ir,iz,p,je,ke,j
integer   :: rank, numRanks, nr, loop,debug
integer   :: nsteps = 5

integer   :: msteps = 1000
integer   :: diag = 10099

real, allocatable :: V1h(:,:), Tmp(:,:),V2h(:,:),error(:,:),Tmp2(:,:)
real, allocatable :: rho(:,:),Resid(:,:),rho2h(:,:),resid2(:,:)
 real, parameter   :: dr = 0.09316442489377
 real, parameter   :: dz = 0.09316442489377

real              :: w, errmax, buf, residTot(32)

character(len=6) :: iter
character(len=6) :: numrlx

type(MPI_Status) :: status

! Initialize MPI library
call MPI_Init()

call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)

! Calculate the number of z elements Nk for each rank (should be evenly distributed)
Nk = NTk/numRanks
if (numRanks*Nk /= NTk) then
   if (rank == 0) then
      print *, "ERROR: z dimension can't be evenly distributed"
      print *, "numRanks is", numRanks, "Nk total is", NTk
   end if
   call MPI_Finalize()
   stop 1
end if

! write nsteps to string for textual output
write(numrlx,"(i6.6)") nsteps

! Allocate arrays
allocate(V1h(-1:Nj+1,-1:Nk+1),   Tmp(-1:Nj+1,-1:Nk+1), V2h(-1:Nj/2+1,-1:Nk/2+1))
allocate(rho(-1:Nj+1,-1:Nk+1), Resid(-1:Nj+1,-1:Nk+1),rho2h(-1:Nj/2+1,-1:Nk/2+1))
allocate(error(-1:Nj+1,-1:Nk+1),Resid2(-1:Nj/2+1,-1:Nk/2+1),Tmp2(-1:Nj/2+1,-1:Nk/2+1))
!! Initialize
rho    = 0.0
rho2h  = 0.0
resid  = 0.0
resid2 = 0.0
V1h    = 0.0
V2h    = 0.0

!! Read in source term and boundary data from file
call readDensity(rho,Nj,Nk)
!Exchange halo once to receive upper boundaries of source term
if (numRanks>1) then
call ExchangeHalo(rank,Nj,Nk,rho)
end if
call RestrictRho(Nj  ,  Nk  ,  rho  ,  rho2h)

call writeData(-1,Nj,Nk,rho, "rhoInit")
call writeData(-1,Nj/2,Nk/2,rho2h, "rhoInit2h")

#ifdef USE_NETCDF
call writeDataNetCDF(rho, 1, Nj, Nk, Ntk, Nl)
#endif

errmax = 1e6
mn = 0
do while(mn.lt.msteps)
!do while(errmax.gt.tol)

!    -------------------------
mn = mn + 1
!Update boundary before relaxing.
!! Dirichlet boundary conditions in the midplane.
if (rank == 0)then
  do ir = 0,Nj+1
    V1h(ir, 0) = V1h(ir,1) 
    V1h(ir,-1) = V1h(ir,2) 
  end do
end if

!Top boundary, calculated directly before call to potential solve
if (rank==numRanks-1)then
  do ir = 0,Nj+1
    V1h(ir,Nk) = rho(ir,Nk)
  end do
end if
!! Left and right boundary.
! Left boundary uses symmetry, right boundary is calculated directly.
 do iz = 0,Nk
   V1h(Nj,iz) = rho(Nj,iz)
!  V1h(0 ,iz) = -V1h(1,iz)  ! Left boundary actually not required for relaxation.
 end do
do j = 1,nsteps
!Relax on V1h using rho as source
if (numRanks>1) then
 call ExchangeHalo(rank,Nj,Nk,V1h)
end if
call Relax(Nj, Nk, m, V1h, Tmp, rho, dr, dz)
end do

call Residual(Nj, Nk, m, V1h, Resid, rho, dr, dz)
call RestrictRho(Nj,Nk,Resid,V2h)

!call writeData(-1,Nj,Nk,V1h, "sol")
!call writeData(-1,Nj/2,Nk/2,V2h, "sol.2")

!call writeData(-1,Nj/2,Nk/2,V2h, "out_sol.2h.")
do j=1,nsteps
call ExchangeHalo(rank,Nj/2,Nk/2,resid2)
call Relax(Nj/2, Nk/2, m, resid2, Tmp2, V2h, 2.0*dr, 2.0*dz)
end do
!write data here

error = 0.0
call prolongate(Nj, Nk, error, resid2)

do j = -1,Nj+1
 do k = -1,Nk+1
  V1h(j,k) = V1h(j,k) - error(j,k)
 end do
end do
do j = 1,nsteps
call ExchangeHalo(rank,Nj,Nk,V1h)
call Relax(Nj, Nk, m, V1h, Tmp, rho, dr, dz)
end do

! Calculate residual. Max value is current error.

call Residual(Nj, Nk, m, V1h, Resid, rho, dr, dz)
 if(mod(mn,diag)==0)then
   write(iter,"(i6.6)") mn
   call writeData(0,Nj,Nk,V1h, "sol." // iter )
call writeData(-1,Nj,Nk,error, "error." // iter )
!  call writeData(0,Nj/2,Nk/2,V2h, "sol.2." // iter )
   print *, "Rank",rank,"Wrote to file at iteration " // iter // ",", maxval(abs(Resid))

  call writeData(0,Nj, Nk, Resid, "resid." // iter ) 
 end if

if (rank ==0 ) then
 residTot(1) = maxval(abs(error))
 do nr = 1,numRanks-1
  call MPI_Recv(buf,1,MPI_DOUBLE_PRECISION, nr, 17, MPI_COMM_WORLD, status)
  if (buf>residTot(1)) then
   residTot(1) = buf
  end if
 end do
else 
 residTot(rank+1) = maxval(abs(error))
 call MPI_Send(residTot(rank+1), 1, MPI_DOUBLE_PRECISION, 0, 17, MPI_COMM_WORLD)
end if

 errmax=maxval(residTot)
 write(75,"(i6.6,1x,f30.24)") mn,errmax
end do ! While loop

print *, "rank",rank," believes we are done!!!"
! Write results
call writeData(1,Nj,Nk,V1h, "final" )
deallocate(V1h,Tmp,rho,Resid,rho2h,resid2,tmp2,v2h)

! Shutdown MPI
call MPI_Finalize()

CONTAINS

Subroutine ExchangeHalo(rank, Nj, Nk, A)
!
! Exchange halo information between neighboring processes
!
use MPI_F08, only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
use MPI_F08, only : MPI_Send, MPI_Recv, MPI_DOUBLE_PRECISION, MPI_Status
   Implicit None
   Integer, intent(in   ) :: Nj,Nk, rank
   Real,    intent(inout) :: A(-1:Nj+1,-1:Nk+1)

   Integer :: above, below, tag,bufSize
   Real    :: topBC, bottomBC   
   type(MPI_Status) :: status
   real  , allocatable :: bufFromAbove(:),bufFromBelow(:)

Call MPI_Comm_size( MPI_COMM_WORLD, numRanks)

   allocate(bufFromAbove(Nj+3),bufFromBelow(Nj+3))
   bufSize=NJ+3
   tag    = 156
   above  = rank + 1
   below  = rank - 1
    
   !! MPI halo exchange for parallel version
!Rank=0, bottom of grid, no below
if (rank == 0) then
  Call MPI_Send(A(:,Nk-1), bufSize , MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD)
  Call MPI_Recv(bufFromAbove, bufSize, MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD, status)
  A(:,Nk+1) = bufFromAbove

!Rank=numRanks-1, top of grid, no above
else if (rank == numRanks-1) then

  Call MPI_Send(A(:,1), bufSize , MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD)
  Call MPI_Recv(bufFromBelow, bufSize, MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD, status)
  A(:,-1) = bufFromBelow
!Middle of grid
else
  Call MPI_Send(A(:,Nk-1), bufSize , MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD)

  Call MPI_Send(A(:,1), bufSize , MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD)
  Call MPI_Recv(bufFromBelow, bufSize, MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD, status)
  Call MPI_Recv(bufFromAbove, bufSize, MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD, status)
  A(:,Nk+1) = bufFromAbove
  A(:, -1) = bufFromBelow

end if   

deallocate(bufFromAbove,bufFromBelow)
End Subroutine ExchangeHalo

end program PoissonRelax

