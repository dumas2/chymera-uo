#undef USE_NETCDF

module io

  ! dimension ids
  integer :: rDimId, zDimId, phiDimId, modeDimId

  ! variable ids
  integer :: jmaxVarId, kmaxVarId, lmaxVarId, ktmaxVarId, iterVarId  ! scalars
  integer :: rhoModeVarId, residVarId                                ! 2D arrays
  integer :: rhoVarId                                                ! 3D arrays

CONTAINS

subroutine readBoundary(A,J,K)
! Reads in the fourier transformed boundary data, calculated by boundary.F
implicit none
integer, intent(in)          :: J,K
real   , intent(inout)       :: A(1:J,1:K)
real                         ::  bndryValues(J+K-2)
integer, parameter           :: fd=12
character(len=*), parameter  :: fmt="(i3,1x,i3,1x,i3,1x,1e16.5)"
character(len=*), parameter  :: fmt2="(i3,1x,i3,1x,1e16.5)"
integer                      :: jj,kk,ll,i,imax
imax=J+K-2

!open file to read
open(unit=fd, file="boundaryConditions.dat")

! read in the data and assign it to array
do i = 1,imax
   read(fd,fmt)  jj,kk,ll,bndryValues(i)
   A(jj-1,kk-1)   =  bndryValues(i)
!   print *, jj, kk, jj-1, kk-1, A(jj-1,kk-1)
end do
!do i = 1,255
!   read(57,fmt2)  jj,kk,bndryValues(i)
!   A(jj,kk)   =  bndryValues(i)
!   print *, jj, kk, jj-1, kk-1, A(jj-1,kk-1)
!end do
!do i = 1,63
!   read(57,fmt2)  jj,kk,bndryValues(i)
!   A(jj,kk)   =  bndryValues(i)
!   print *, jj, kk, jj-1, kk-1, A(jj-1,kk-1)
!end do
close(fd)

end subroutine readBoundary

subroutine readDensity(A,jmax,kmax)
! Reads in the fourier transformed density for a given mode,
! as calculated by fft.F by a call in pot3.F
!
! Note: jmax and kmax must be powers of 2 for multigrid to work
!       interior for multigrid is (1:jmax-1,1:kmax-1)
!
use MPI_F08, only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
use MPI_F08, only : MPI_Send, MPI_Recv, MPI_DOUBLE_PRECISION, MPI_Status
implicit none
integer, intent(in)  :: jmax,kmax
real   , intent(out) :: A(-1:jmax+1,-1:kmax+1) ! dimensions for multigrid, odd interior pts, 2 halo cells
real   , allocatable :: buf(:,:)
integer, parameter   :: fd=13
integer              :: ir,iz,jj,kk,nr
integer              :: rank, numRanks, bufSize, tag=11
type(MPI_Status)     :: status
character(len=*), parameter :: fmt="(i3,1x,i3,1x,1e16.5)"

call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)

bufSize = (jmax+3)*(kmax+3)
allocate(buf(-1:jmax+1,-1:kmax+1))

! initialize so that boundaries will be zero
buf = 0.0

!! Rank 0 should read in data and then send it to cohorts
!    - z dimension is broken into numRanks partitions
!
if (rank == 0) then
  open(unit=fd, file="density.dat", FORM="FORMATTED")

  buf(:,-1) = 0.0    ! extra boundaries only needed in neighbor below

  ! read in data for rank 0
  do iz = 0, kmax
    do ir = 0, jmax
      read(fd,fmt) jj, kk, buf(ir,iz)
    end do
  end do
  if (numRanks > 1) then   ! extra boundary cells (kmax+1) extends into neighbor above
    do ir = 0, jmax
      read(fd,fmt) jj, kk, buf(ir,kmax+1)
    end do
  else
      buf(:,kmax+1) = 0.0     ! extra boundary cells (kmax+1) not needed
  end if

  A(:,:) = buf(:,:)      ! copy entire array

  buf(:,-1) = buf(:,kmax-1)  ! lower bc for neigbor copied from last interior cell of previous
  buf(:, 0) = buf(:,kmax  )  ! lower shared boundary for neigbor from shared upper cell of previous
  buf(:, 1) = buf(:,kmax+1)  ! first interior cell for neigbor upper boundary cell of previous

  ! read in data for other ranks and send the buffer
  do nr = 1, numRanks-1
    do iz = 2, kmax          ! halo and shared cell (iz=-1:0) already copied
      do ir = 0, jmax
        read(fd,fmt) jj, kk, buf(ir,iz)
      end do
    end do
    if (nr /= numRanks-1) then
      do ir = 0, jmax
        read(fd,fmt) jj, kk, buf(ir,kmax+1)
      end do
    else
      buf(kmax+1,:) = 0.0     ! extra boundaries only needed in neighbor above
    end if
    call MPI_Send(buf, bufSize, MPI_DOUBLE_PRECISION, nr, tag, MPI_COMM_WORLD)
    buf(:,-1) = buf(:,kmax-1)  ! lower bc for neigbor copied from last interior cell of previous
    buf(:, 0) = buf(:,kmax  )  ! lower shared boundary for neigbor from shared upper cell of previous
    buf(:, 1) = buf(:,kmax+1)  ! first interior cell for neigbor upper boundary cell of previous
  end do

  close(fd)

else  ! receiving ranks (other than rank 0)

  call MPI_Recv(buf, bufSize, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, status)
  A(:,:) = buf(:,:)      ! copy entire array

end if

deallocate(buf)

end subroutine readDensity


subroutine writeData(b,J,K,A,id)
! Writes array as text to file for pretty plotting.
implicit none

integer, intent(in)               :: J,K,b
real                              :: A(:,:)
integer                           :: i,l,fd
character(len=*), parameter       :: fmt = "(i3,1x,i3,1x,1e22.10)"
character(len=*), intent(in)      :: id

fd = 14
!open file for writing
open(fd,file="out_" // id // ".txt")

! write data to file
do l = b,K
 do i = b,J
    write(fd,fmt)  i, l, A(i,l)
 end do
end do

! close file
close(fd)

end subroutine writeData

#ifdef USE_NETCDF
subroutine writeDataNetCDF(rhoMode, id, jmax, kmax, ktmax, lmax)
  ! Writes array as a netCDF file
  use netcdf
  use MPI_F08, only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
  use MPI_F08, only : MPI_Send, MPI_Recv, MPI_DOUBLE_PRECISION, MPI_Status
  implicit none
  integer, intent(in) :: id, jmax, kmax, ktmax, lmax
  real, intent(in)    :: rhoMode(-1:,-1:)

  ! local variables
  integer             :: rank, numRanks, bufSize, nr, tag
  real, allocatable   :: buf(:,:)
  character(len=64)   :: filename
  type(MPI_Status)    :: status

  ! netCDF variables
  integer, parameter :: NDIMS = 4
  integer :: ncid, dimids(NDIMS)

  call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank)

  tag = 13
  bufSize = (jmax+3)*kmax             ! everyone sends k=1:kmax

  !! Rank 0 should recv data from cohorts and then write it
  !    - z dimension is broken into numRanks partitions
  !
  if (rank == 0) then
     ! allocate memory for the interior (plus upper shared boundary) in z
     allocate(buf(-1:jmax+1,1:kmax))

    ! initialize netCDF metadata
    call initFileNetCDF(id, jmax, kmax, ktmax, lmax)

    write(filename,'("rho3d_", i7.7, ".nc")') id
    call nc_check( nf90_open(trim(filename), NF90_WRITE, ncid) )

    ! write lower (in z) boundary values plus interior of rank 0
    call nc_check( nf90_put_var(ncid, rhoModeVarId, rhoMode(:,-1:kmax), start=[1,1]) )

    ! recv data from other ranks and write it to the netCDF file
    do nr = 1, numRanks-1
      call MPI_Recv(buf, bufSize, MPI_DOUBLE_PRECISION, nr, tag, MPI_COMM_WORLD, status)
      call nc_check( nf90_put_var(ncid, rhoModeVarId, buf, start=[1,1+nr*(kmax+2)]) )
    end do

    ! top boundary layer (ktmax+1) is zero (unused)

    buf(:,1) = 0.0d0
    call nc_check( nf90_put_var(ncid, rhoModeVarId, buf(:,1), start=[1,ktmax+3]) )

    ! close the file
    call nc_check( nf90_close(ncid) )
    deallocate(buf)

  else

    ! send data from interior plus upper shared boundary (in z) to rank 0 (-1:jmax+1,1:kmax)
    call MPI_Send(rhoMode(-1,1), bufSize, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD)

  end if

end subroutine writeDataNetCDF


subroutine initFileNetCDF(id, jmax, kmax, ktmax, lmax)
  use netcdf
  implicit none
  integer, intent(in) :: id, jmax, kmax, ktmax, lmax

  integer, parameter :: NDIMS = 2
  integer :: ncid, j_dimid, k_dimid, dimids(NDIMS), varid
  character(len=64) :: filename

  !! define meta data for the file
  !
  ! create the netcdf file, overwriting the file if it already exists
  write(filename,'("rho3d_", i7.7, ".nc")') id
  call nc_check( nf90_create(trim(filename), NF90_CLOBBER, ncid) )

  ! define the dimensions
  call nc_check( nf90_def_dim(ncid, "r", jmax+3, rDimId) )          ! range (-1:jmax+1)
  call nc_check( nf90_def_dim(ncid, "z", ktmax+3, zDimId) )         ! range (-1:ktmax+1)
  call nc_check( nf90_def_dim(ncid, "phi", lmax+3, phiDimId) )      ! range (-1:lmax+1)
  call nc_check( nf90_def_dim(ncid, "mode", kmax+3, modeDimId) )    ! range (-1:kmax+1)

  ! scalar variables
  call nc_check( nf90_def_var(ncid, "jmax", NF90_INT, jmaxVarId) )
  call nc_check( nf90_put_att(ncid, jmaxVarId, "long_name", "r dim size is jmax+3, indices (-1:jmax+1)") )
  call nc_check( nf90_def_var(ncid, "kmax", NF90_INT, kmaxVarId) )
  call nc_check( nf90_put_att(ncid, kmaxVarId, "long_name", "z dim size per rank is kmax+3, indices (-1:kmax+1)") )
  call nc_check( nf90_def_var(ncid, "ktmax", NF90_INT, ktmaxVarId) )
  call nc_check( nf90_put_att(ncid, ktmaxVarId, "long_name", "z dim size total is ktmax+3, indices (-1:ktmax+1)") )
  call nc_check( nf90_def_var(ncid, "lmax", NF90_INT, lmaxVarId) )
  call nc_check( nf90_put_att(ncid, lmaxVarId, "long_name", "phi dim size is lmax+3, indices (-1:lmax+1)") )
  call nc_check( nf90_def_var(ncid, "iter", NF90_INT, iterVarId) )
  call nc_check( nf90_put_att(ncid, iterVarId, "long_name", "iteration number") )

  !! note that fortran arrays are stored in column-major order
  !
  ! 3D arrays
  call nc_check( nf90_def_var(ncid, "rho", NF90_FLOAT, [rDimId, zDimId, phiDimId], rhoVarId) )
  call nc_check( nf90_put_att(ncid, rhoVarId, "long_name", "density distribution") )
  call nc_check( nf90_put_att(ncid, rhoVarId, "units", "g/cm^3?") )
  ! 2D arrays
  call nc_check( nf90_def_var(ncid, "rho_mode", NF90_FLOAT, [rDimId, zDimId], rhoModeVarId) )
  call nc_check( nf90_put_att(ncid, rhoModeVarId, "long_name", "transformed density distribution") )
  call nc_check( nf90_def_var(ncid, "resid", NF90_FLOAT, [rDimId, zDimId], residVarId) )
  call nc_check( nf90_put_att(ncid, residVarId, "long_name", "relaxation residual") )

  ! finished creating metadata
  call nc_check( nf90_enddef(ncid) )

  ! write data to file
  call nc_check( nf90_put_var(ncid, jmaxVarId,  jmax) )
  call nc_check( nf90_put_var(ncid, kmaxVarId,  kmax) )
  call nc_check( nf90_put_var(ncid, lmaxVarId,  lmax) )
  call nc_check( nf90_put_var(ncid, ktmaxVarId, ktmax) )
  call nc_check( nf90_put_var(ncid, iterVarId,  id) )

  ! close the file
  call nc_check( nf90_close(ncid) )

end subroutine initFileNetCDF

subroutine nc_check(status)
  use netcdf
  implicit none
  integer, intent(in) :: status

  if (status /= nf90_noerr) then
     print *, "ERROR in netCDF call, status is:", trim(nf90_strerror(status))

     call MPI_Finalize
     STOP 1
  end if
end subroutine nc_check

! end USE_NETCDF
#endif

end module io
