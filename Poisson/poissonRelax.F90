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
use MultiGrid, only : Relax, Residual
use io       , only : readDensity, writeData
#ifdef USE_NETCDF
use io       , only : initFileNetCDF, writeDataNetCDF
#endif
use MPI_F08  , only : MPI_Init, MPI_Finalize, MPI_Barrier
use MPI_F08  , only : MPI_Send, MPI_Recv, MPI_DOUBLE_PRECISION, MPI_Status
use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
implicit none

real   , parameter :: tol =  1e-5

integer, parameter :: Nj  =  256     ! number of r interior elements
integer, parameter :: NTk =  128     ! number of z interior elements total
integer, parameter :: Nl  =    1     ! number of phi interior elements
integer            :: Nk             ! number of z interior elements per rank in z
integer, parameter :: m   =  1       ! consider m=1 Fourier mode for testing


integer   :: i,k,mn,ir,iz,p,je,ke,j
integer   :: rank, numRanks, nr
integer   :: nsteps = 1000
integer   :: msteps = 400
integer   :: diag = 100

real, allocatable :: V1h(:,:), Tmp(:,:)
real, allocatable :: rho(:,:),Resid(:,:)
real, parameter   :: dr = 0.09316442463639991
real, parameter   :: dz = 0.09316442463639991

real              :: w, errmax, buf, residTot(32)

character(len=6) :: iter
character(len=6) :: numrlx

type(MPI_Status) :: status
! Initialize MPI library
call MPI_Init()

call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)
print *, "numRanks = ",numRanks
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
allocate(V1h(-1:Nj+1,-1:Nk+1),   Tmp(-1:Nj+1,-1:Nk+1))
allocate(rho(-1:Nj+1,-1:Nk+1), Resid(-1:Nj+1,-1:Nk+1))

!! Initialize
V1h =  0.0
rho =  0.0
resid= 0.0

!! Read in source term and boundary data from file
call readDensity(rho,Nj,Nk)
!Exchange halow once to receive upper boundaries
call ExchangeHalo(rank,Nj,Nk,rho)

call writeData(-1,Nj,Nk,rho, "rhoInit")

#ifdef USE_NETCDF
call writeDataNetCDF(rho, 1, Nj, Nk, Ntk, Nl)
#endif

errmax = 1e6
mn = 0
do while(mn.lt.msteps)
!do while(errmax.gt.tol)
!... Relax solution on 1h mesh

!    -------------------------
mn = mn + 1
!!! Update boundary before relaxing.
!! Top boundary, calculated directly before call to potential solve
if (rank==numRanks-1)then
  do ir = 0,Nj
    V1h(ir,Nk) = rho(ir,Nk)
  end do
end if
!! Left and right boundary.
! Left boundary uses symmetry, right boundary is calculated directly.
 do iz = 0,Nk
   V1h(Nj,iz) = rho(Nj,iz)
!   V1h(0 ,iz) = -V1h(1,iz)  ! Left boundary actually not required for relaxation.
 end do
! Dirichlet boundary conditions in the midplane.
if (rank == 0)then
  do ir = 0,Nj
    V1h(ir,0) = V1h(ir,1) 
  end do
end if

do j = 1,5
! Relax on V1h using rho as source
 call Relax(Nj, Nk, m, V1h, Tmp, rho, dr, dz)
 call ExchangeHalo(rank,Nj,Nk,V1h)
! Calculate residual. Max value is current error. 
end do

call Residual(Nj, Nk, m, V1h, Resid, rho, dr, dz)

 if(mod(mn,diag)==0)then
   write(iter,"(i6.6)") mn
   call writeData(0,Nj,Nk,V1h, "sol." // iter )
   print *, "Rank",rank,"Wrote to file at iteration " // iter // ",", maxval(abs(Resid))

  call writeData(-1,Nj, Nk, Resid, "resid." // iter ) 
 end if

if (rank ==0 ) then
 residTot(1) = maxval(abs(Resid))
 do nr = 1,numRanks-1
  call MPI_Recv(buf,1,MPI_DOUBLE_PRECISION, nr, 17, MPI_COMM_WORLD, status)
  if (buf>residTot(1)) then
   residTot(1) = buf
  end if
 end do
 errmax=maxval(residTot)
 write(75,"(i6.6,1x,f30.24)") mn,errmax
else
 residTot(rank+1) = maxval(abs(Resid))
 call MPI_Send(residTot(rank+1), 1, MPI_DOUBLE_PRECISION, 0, 17, MPI_COMM_WORLD)
end if

end do ! While loop
print *, "rank",rank," believes we are done!!!"
! Write results
call writeData(1,Nj,Nk,V1h, "final" )
    do iz=-1,Nk+1
      do ir = -1,Nj+1
        write(30+rank,*) ir, iz, V1h(ir,iz)
      end do
    end do
deallocate(V1h,Tmp,rho,Resid)

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
!Call MPI_Comm_rank( MPI_COMM_WORLD, rank)

   allocate(bufFromAbove(Nj+3),bufFromBelow(Nj+3))
   bufSize=NJ+3
   tag    = 156
   above  = rank + 1
   below  = rank - 1
    
   !! MPI halo exchange for parallel version
   !
if (rank == 0) then
  Call MPI_Send(A(:,Nk-1), bufSize , MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD)
! print *, "rank",rank," sends to ",above
  Call MPI_Recv(bufFromAbove, bufSize, MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD, status)
! print *, "rank", rank,"receives from",above
  A(:,Nk+1) = bufFromAbove
    ! write(30+rank,*) rank, above, below
    ! do iz=-1,Nk+1
    !   do ir = -1,Nj+1
    !     write(30+rank,*) ir, iz, A(ir,iz)
    !   end do
    ! end do
else if (rank == numRanks-1) then

  Call MPI_Send(A(:,1), bufSize , MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD)
! print *, "rank",rank," sends to ",below
  Call MPI_Recv(bufFromBelow, bufSize, MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD, status)
! print *, "rank", rank,"receives from",below
  A(:,-1) = bufFromBelow
    ! write(30+rank,*) rank, above, below
    ! do iz=-1,Nk+1
    !   do ir = -1,Nj+1
    !     write(30+rank,*) ir, iz, A(ir,iz)
    !   end do
    ! end do
else

  Call MPI_Send(A(:,Nk-1), bufSize , MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD)
! print *, "rank",rank," sends to ",above
  Call MPI_Recv(bufFromAbove, bufSize, MPI_DOUBLE_PRECISION, above, tag, MPI_COMM_WORLD, status)
  A(:,Nk+1) = bufFromAbove

  Call MPI_Send(A(:,1), bufSize , MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD)
! print *, "rank",rank," sends to ",below
  Call MPI_Recv(bufFromBelow, bufSize, MPI_DOUBLE_PRECISION, below, tag, MPI_COMM_WORLD, status)
! print *, "rank", rank,"receives from",below

! print *, "rank", rank,"receives from",above
  A(:, -1) = bufFromBelow

    ! write(30+rank,*) rank, above, below
    ! do iz=-1,Nk+1
    !   do ir = -1,Nj+1
    !     write(30+rank,*) ir, iz, A(ir,iz)
    !   end do
    ! end do
end if   

deallocate(bufFromAbove,bufFromBelow)
End Subroutine ExchangeHalo

end program PoissonRelax

