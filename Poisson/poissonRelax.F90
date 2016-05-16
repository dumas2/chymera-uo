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
use MultiGrid, only : Relax, Residual,RelaxB
use io       , only : readDensity, writeData
#ifdef USE_NETCDF
use io       , only : initFileNetCDF, writeDataNetCDF
#endif
use MPI_F08  , only : MPI_Init, MPI_Finalize
use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
implicit none

real   , parameter :: tol  =  1e-4

integer, parameter :: Nj   =  256     ! number of r interior elements
integer, parameter :: NTk  =   64     ! number of z interior elements total
integer, parameter :: Nl   =    1     ! number of phi interior elements
integer            :: Nk              ! number of z interior elements per rank in z

integer   :: i,m,k,mn,ir,iz,p,je,ke
integer   :: rank, numRanks
integer   :: nsteps = 1000
integer   :: msteps = 300
integer   :: diag = 500  

real, allocatable :: V1h(:,:), Tmp(:,:)
real, allocatable :: rho(:,:),Resid(:,:)
real      :: dr, dz, w, errmax

character(len=6) :: iter
character(len=6) :: numrlx

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

! Define dr and dz, these should be passed in eventually
dr = 0.09316442463639991
dz = 0.09316442463639991

! We are considering the m=1 Fourier mode for testing.
m = 1

! Allocate arrays
allocate(V1h(-1:Nj+1,-1:Nk+1),   Tmp(-1:Nj+1,-1:Nk+1))
allocate(rho(-1:Nj+1,-1:Nk+1), Resid(-1:Nj+1,-1:Nk+1))

!! Initialize
V1h =  0.0
rho =  0.0
resid= 0.0

!! Read in source term and boundary data from file
call readDensity(rho,Nj,Nk)

!! CHANGE THIS so that we start on the coursest mesh, and go up!!!
!call writeData(1,Nj,Nk,V1h, numrlx // ".000")

#ifdef USE_NETCDF
call writeDataNetCDF(rho, 1, Nj, Nk, Ntk, Nl)
#endif

errmax = 1e6
mn = 0
!do while(mn.lt.msteps)
do while(errmax.gt.tol)
!... Relax solution on 1h mesh

!print *, "mn= ",mn
!    -------------------------
! do i = 1, nsteps
mn = mn + 1
!!! Update boundary before relaxing.
!! Top boundary, calculated directly before call to potential solve
  do ir = 1,Nj-1
    V1h(ir,Nk) = rho(ir,Nk)
 !   V1h(ir,1) = rho(ir,0)
  end do

!! Left and right boundary.
! Left boundary uses symmetry, right boundary is calculated directly.
  do iz = 1,Nk-1
    V1h(Nj,iz) = rho(Nj,iz)
 !   V1h(0 ,iz) = -V1h(1,iz)  ! Left boundary actually not required for relaxation.
  end do
! Dirichlet boundary conditions in the midplane.
do ir = 1,Nj-1
   V1h(ir,0) = V1h(ir,1) 
end do
 

! Relax on V1h using rho as source
  call RelaxB(Nj, Nk, m, V1h, Tmp, rho, dr, dz)

! end do

!! Call residual
! do i = 1, nsteps
!mn=mn+1
!!! Update boundary before relaxing.
!! Top boundary, calculated directly before call to potential solve
!  do ir = 1,Nj-1
!    V1h(ir,Nk) = rho(ir,Nk)
 !   V1h(ir,1) = rho(ir,0)
!  end do

!! Left and right boundary.
! Left boundary uses symmetry, right boundary is calculated directly.
!  do iz = 1,Nk-1
!    V1h(Nj,iz) = rho(Nj,iz)
 !   V1h(0 ,iz) = -V1h(1,iz)  ! Left boundary actually not required for relaxation.
!  end do
! Dirichlet boundary conditions in the midplane.
!do ir = 1,Nj-1
!   V1h(ir,0) = V1h(ir,1) 
!end do
 

! Relax on V1h using rho as source
! call RelaxB(Nj, Nk, m, V1h, Tmp, rho, dr, dz)

!end do


 call Residual(Nj, Nk, m, V1h, Resid, rho, dr, dz)

 if(mod(mn,diag)==0)then
   write(iter,"(i6.6)") mn
   call writeData(0,Nj,Nk,V1h, "sol." // iter )
   print *, "Wrote to file at iteration " // iter // ",", maxval(abs(Resid))

  call writeData(-1,Nj, Nk, Resid, "resid." // iter ) 
 end if


  errmax = maxval(abs(Resid))
  write(75,"(i6.6,1x,f30.24)") mn,errmax



end do ! While loop

! Write results
call writeData(1,Nj,Nk,V1h, "final" )

deallocate(V1h,Tmp,rho,Resid)

! Shutdown MPI
call MPI_Finalize()

end program PoissonRelax

