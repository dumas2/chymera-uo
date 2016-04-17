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

subroutine readDensity(A,J,K)
! Reads in the fourier transformed density for a given mode,
! as calculated by fft.F by a call in pot3.F
implicit none
integer, intent(in)  :: J,K
real   , intent(out) :: A(1:J,1:K)
integer, parameter   :: fd=13
character(len=*), parameter :: fmt="(i3,1x,i3,1x,1e16.5)"
integer              :: ir,iz,jj,kk
!imax=J+K-1

!open file to read
open(unit=fd, file="density.dat", FORM="FORMATTED")
! read in the data and assign it to array
 do iz = 1,K-1
  do ir = 1,J-1
      read(fd,fmt) jj,kk,A(ir,iz)
!      print *, jj,kk,ir,iz,A(ir,iz)
  end do
 end do
close(fd)
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
