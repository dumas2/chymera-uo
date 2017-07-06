#undef USE_NETCDF

subroutine PoissMultigrid(rho,phi,Nj,Nk,Nl,l,dr,dz)
!=============================================================================
! Name        : poissonRelax.F90
! Version     : 1.0
! Description : Solves the Poisson equation in two dimensions.
!
! Method      :
!
!   Multigrid method.   
!==============================================================================
use MultiGrid, only : Relax, Residual, Prolongate, Restrict, RestrictRho
use io       , only : readDensity, writeData
#ifdef USE_NETCDF
use io       , only : initFileNetCDF, writeDataNetCDF
#endif
use MPI_F08  , only : MPI_Init, MPI_Finalize, MPI_Barrier,MPI_Reduce
use MPI_F08  , only : MPI_Send, MPI_Recv, MPI_Status,mpi_bcast
use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
use MPI_F08  , only : MPI_Integer, MPI_DOUBLE_PRECISION, mpi_max
implicit none

real   , parameter :: tol =  1e-7

real, intent(inout) :: rho(-1:Nj+1, -1:Nk+1)
real, intent(inout) :: phi(-1:Nj+1, -1:Nk+1)
real, intent(in) :: dz, dr
integer, intent(in) :: Nj, Nk, Nl, l


integer   :: i,k,mn,ir,iz,p,j,m
integer   :: rank, numRanks, nr, loop
integer   :: nsteps = 2, debug = 1

integer   :: msteps = 200
integer   :: diag = 100
real, allocatable :: V1h(:,:)  , Tmp(:,:)  , V2h(:,:)   , error(:,:)
real, allocatable :: Resid(:,:), rho2h(:,:), resid2(:,:), Tmp2(:,:)

real              :: w, errmax, buf, residTmp

character(len=6) :: iter
character(len=6) :: numrlx
character(len=3) :: fm

type(MPI_Status) :: status

! Initialize MPI library
call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)

! write nsteps to string for textual output
write(numrlx,"(i6.6)") nsteps

! Allocate arrays
allocate(V1h(-1:Nj+1,-1:Nk+1),   Tmp(-1:Nj+1,-1:Nk+1), V2h(-1:Nj/2+1,-1:Nk/2+1))
allocate(Resid(-1:Nj+1,-1:Nk+1),rho2h(-1:Nj/2+1,-1:Nk/2+1))
allocate(error(-1:Nj+1,-1:Nk+1),Resid2(-1:Nj/2+1,-1:Nk/2+1),Tmp2(-1:Nj/2+1,-1:Nk/2+1))

!! Initialize
phi    = 0.0
rho2h  = 0.0
resid  = 0.0
resid2 = 0.0
V1h    = 0.0
V2h    = 0.0

!Exchange halo once to receive upper boundaries of source term
call RestrictRho(Nj  ,  Nk  ,  rho  ,  rho2h)
!call writeData(-1,Nj,Nk,rho, "rhoInit")
!call writeData(-1,Nj/2,Nk/2,rho2h, "rhoInit2h")
#ifdef USE_NETCDF
call writeDataNetCDF(rho, 1, Nj, Nk, Ntk, Nl)
#endif

! this is the fourier mode. 
if(l.lt.Nl/2+1) then
m  = l-1
else
m = l - Nl/2-1
end if

!! Main while loop
errmax = 1e6
mn = 0
!do while(mn.lt.msteps) ! used for testing
do while(errmax.gt.tol) ! Main while loop

!    -------------------------
mn = mn + 1
!Update boundary before relaxing.
resid = 0.0
!! Dirichlet boundary conditions in the midplane.
if (rank == 0)then
  do ir = 0,Nj+1
    V1h(ir, 0) = V1h(ir,1) 
!   V1h(ir,-1) = V1h(ir,2) 
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
 do iz = -1,Nk+1
  V1h(Nj,iz) = rho(Nj,iz)
 end do

!Relax on V1h using rho as source
do j = 1,nsteps
  if(numRanks>1) then
   call ExchangeHalo(rank,Nj,Nk,V1h)
  end if
 call Relax(Nj, Nk, m, V1h, Tmp, rho, dr, dz)
end do

 do iz = -1,Nk+1
  if(l==1) then
    V1h(0 ,iz) = V1h(1,iz)  ! Left boundary actually not required for relaxation
  else if(l.le.Nl/2) then
    V1h(0 ,iz) = -V1h(1,iz) ! Left boundary actually not required for relaxation
  else
    V1h(0 ,iz) = V1h(1,iz)  ! Left boundary actually not required for relaxation
  end if
 end do

!Calculate the residual, and restrict to coarse mesh.
call Residual(Nj, Nk, m, V1h, Resid, rho, dr, dz)
call RestrictRho(Nj,Nk,Resid,V2h)

! Relax on coarse mesh
do j=1,nsteps
call ExchangeHalo(rank,Nj/2,Nk/2,resid2)
call Relax(Nj/2, Nk/2, m, resid2, Tmp2, V2h, 2.0*dr, 2.0*dz)
end do

! Interpolate error to fine mesh.
error = 0.0
call ExchangeHalo(rank,Nj/2,Nk/2,resid2)
call prolongate(Nj, Nk, error, resid2)
call ExchangeHalo(rank,Nj,Nk,error)

! Correct
if(rank==numRanks-1) then
 do j = 0,Nj-1
  do k = 0,Nk-1
   V1h(j,k) = V1h(j,k) - error(j,k)
  end do
 end do
else
 do j = 0,Nj-1
  do k = 0,Nk
   V1h(j,k) = V1h(j,k) - error(j,k)
  end do
 end do
end if

!update boundary again
if (rank == 0)then
  do ir = 0,Nj+1
    V1h(ir, 0) = V1h(ir,1) 
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
 do iz = -1,Nk+1
  V1h(Nj,iz) = rho(Nj,iz)
 end do

! Relax again on fine mesh.
do j = 1,nsteps
call ExchangeHalo(rank,Nj,Nk,V1h)
call Relax(Nj, Nk, m, V1h, Tmp, rho, dr, dz)
end do
call ExchangeHalo(rank,Nj,Nk,V1h)



! Calculate residual. Max value is current error.
call Residual(Nj, Nk, m, V1h, Resid, rho, dr, dz)
call Relax(Nj, Nk, m, resid, Tmp, V1h, dr, dz)

! Write data for debug
 if((mod(mn,diag)==0).and.(debug==1))then
   write(iter,"(i6.6)") mn
   write(fm,"(i3.3)") l
     call writeData(0,Nj,Nk,V1h, "sol." // iter // fm)
     call writeData(-1,Nj,Nk,error, "error." // iter // fm)
     call writeData(-1,Nj/2,Nk/2,resid2, "error.2h." // iter // fm)
   print *, "Rank",rank,"Wrote to file at iteration " // iter // ",", maxval(abs(error))

 end if
 ! Calculate max error on any process. 
 residTmp = maxval(abs(error))
call MPI_Reduce(residTmp,errmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD)
call MPI_Bcast(errmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD)

!write max error to file
if(rank==0) then
write(75,"(i2.2,1x,i6.6,1x,f30.24)") rank,mn,errmax
end if

end do ! While loop

! set output
phi = V1h

deallocate(V1h,Tmp,Resid,rho2h,resid2,tmp2,v2h)


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

end subroutine PoissMultigrid

