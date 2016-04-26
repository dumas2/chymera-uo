module io

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
real   , intent(out) :: A(0:jmax,0:kmax) ! dimensions set up for multigrid, odd interior pts
real   , allocatable :: buf(:,:)
integer, parameter   :: fd=13
integer              :: ir,iz,jj,kk,nr
integer              :: rank, numRanks, bufSize, tag=11
type(MPI_Status)     :: status
character(len=*), parameter :: fmt="(i3,1x,i3,1x,1e16.5)"

call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)

bufSize = (jmax+1)*(kmax-1)
allocate(buf(jmax+1,kmax-1))  ! temporary buffer (interior points in z only)

!! Rank 0 should read in data and then send it to cohorts
!    - z dimension is broken into numRanks partitions
!
if (rank == 0) then
  open(unit=fd, file="density.dat", FORM="FORMATTED")

  ! read lower bc in z (global)
  do ir = 0, jmax
    read(fd,fmt) jj, kk, A(ir,0)
  end do

  ! read in data for rank 0 (interior only for z)
  do iz = 1, kmax-1
    do ir = 0, jmax
      read(fd,fmt) jj, kk, A(ir,iz)
    end do
  end do

  ! read in data for other ranks and send the buffer (interior only in z)
  do nr = 1, numRanks-1
    do iz = 1, kmax-1      ! interior points only
      do ir = 0, jmax
        read(fd,fmt) jj, kk, buf(ir,iz)
      end do
    end do
    call MPI_Send(buf, bufSize, MPI_DOUBLE_PRECISION, nr, tag, MPI_COMM_WORLD)
  end do

  ! read upper bc in z (global)
  do ir = 0, jmax
    read(fd,fmt) jj, kk, A(ir,kmax)
  end do

  A(:,0   ) = 0.0d0   ! initialize lower boundary in z (0 for now)
  A(:,kmax) = 0.0d0   ! initialize upper boundary in z (0 for now)

  close(fd)

else  ! receiving ranks (other than rank 0)

  call MPI_Recv(buf, bufSize, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, status)

  A( :    ,0       ) = 0.0d0   ! initialize lower boundary in z (0 for now)
  A(0:jmax,1:kmax-1) = buf
  A( :    ,  kmax  ) = 0.0d0   ! initialize upper boundary in z (0 for now)

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

end module io
