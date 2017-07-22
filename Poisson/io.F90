#undef USE_NETCDF

module io

  ! dimension ids
  integer :: rDimId, zDimId, phiDimId, modeDimId

  ! variable ids
  integer :: jmaxVarId, kmaxVarId, lmaxVarId, ktmaxVarId, iterVarId  ! scalars
  integer :: rhoModeVarId, residVarId                                ! 2D arrays
  integer :: rhoVarId                                                ! 3D arrays

CONTAINS

subroutine readData(A,jmax,kmax,lmax)
! Reads in density and calculated boundary conditions.      
!
! Note: jmax and kmax must be powers of 2 for multigrid to work
!       interior for multigrid is (1:jmax-1,1:kmax-1)
!
use MPI_F08, only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
use MPI_F08, only : MPI_Send, MPI_Recv, MPI_DOUBLE_PRECISION, MPI_Status
implicit none
integer, intent(in)  :: jmax,kmax,lmax
real   , intent(inout) :: A(-1:jmax+1,-1:kmax+1,1:lmax) ! dimensions for multigrid, odd interior pts, 2 halo cells
real   , allocatable :: buf(:,:,:),buf2(:,:,:)
integer, parameter   :: fd=14
integer              :: ir,iz,jj,kk,ll,nr
integer              :: rank, numRanks, bufSize, tag=11
type(MPI_Status)     :: status
character(len=*), parameter :: fmt="(i3,1x,i3,1x,1pe22.5)"

call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)

bufSize = (jmax+2)*(kmax+3)*lmax
allocate(buf(1:jmax+2,1:numRanks*kmax+2,1:lmax),buf2(1:jmax+2,1:kmax+3,1:lmax))

! initialize so that boundaries will be zero
buf = 0.0
buf2 = 0.0
!! Rank 0 should read in data and then send it to cohorts
!    - z dimension is broken into numRanks partitions
!


if (rank == 0) then
  open(unit=fd, FORM="UNFORMATTED")

  ! read in data for rank 0
  read(14) buf

 ! write buf for master to array  
do ll = 1,lmax
 do kk = 1,kmax+2
  do jj = 1,jmax+2
  A(jj-1,kk-1,ll) = buf(jj,kk,ll)      ! copy entire array
  end do
 end do
end do

do nr = 1,numRanks-1
      call MPI_Send(buf(:,nr*kmax:kmax*(nr+1)+2,:),bufSize,MPI_DOUBLE_PRECISION,nr,tag,MPI_COMM_WORLD)
     print *, "sent k= ",nr*kmax,"-",kmax*(nr+1)+2,"data to rank ",nr
end do

  ! read in data for other ranks and send the buffer

close(fd)

else  ! receiving ranks (other than rank 0)

  call MPI_Recv(buf2, bufSize, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, status)

do ll = 1,lmax
 do kk = 1,kmax+3     ! contains cell above and below shared cell.
  do jj = 1,jmax+2
  A(jj-1,kk-2,ll) = buf2(jj,kk,ll)  
  end do
 end do
end do

end if

deallocate(buf,buf2)

end subroutine readData


subroutine writeData(b,J,K,A,id)
! Writes array as text to file for easy plotting using gnuplot.
use MPI_F08, only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD
use MPI_F08, only : MPI_Send, MPI_Recv, MPI_DOUBLE_PRECISION, MPI_Status
implicit none
integer, intent(in)               :: J,K,b
real, intent(in)                  :: A(-1:J+1,-1:K+1)
character(len=*), intent(in)      :: id

!local variables
integer                           :: i,l,fd,tag
integer                           :: bufSize,nr,rank,numRanks
real, allocatable                 :: buf(:,:)
character(len=*), parameter       :: fmt = "(i3,1x,i3,1x,1e22.10)"
!character(len=*), parameter       :: fmt = "(7g8.3)"
type(MPI_Status)                  :: status

call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)

tag = 13
bufSize = (J+3)*(K+3)

if (rank == 0) then
 ! allocate memory for the interior (plus upper shared boundary) in z
  allocate(buf(-1:J+1,-1:K+1))

  fd = 14
 !open file for writing
  open(fd,file="out_" // id // ".txt")
 ! write lower (in z) boundary values plus interior of rank 0

 do l = b,K
   do i = b,J
     write(fd,fmt)  i, l, A(i,l)
   end do
 end do
 ! recv data from other ranks and write it to the file
  do nr = 1, numRanks-1
    call MPI_Recv(buf, bufSize, MPI_DOUBLE_PRECISION, nr, tag, MPI_COMM_WORLD, status)
    do l = 1,K
     do i = b,J
       write(fd,fmt)  i, K*nr+l, buf(i,l)
     end do
    end do
 ! close file
  end do
  close(fd)
  deallocate(buf)

else 

 ! send data from interior plus upper shared boundary (in z) to rank 0 (-1:jmax+1,1:kmax)
  call MPI_Send(A(-1:J+1,-1:K+1), bufSize, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD)
end if

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
