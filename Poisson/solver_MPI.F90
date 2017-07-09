program poissonSolverMPI
 use defines_mod
 use param_mod
 use io, only : readData, writeData
 use MPI_F08  , only : MPI_Init, MPI_Finalize, MPI_Barrier
 use MPI_F08  , only : MPI_Send, MPI_Recv, MPI_DOUBLE_PRECISION, MPI_Status
 use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD,MPI_INTEGER

implicit none 

real, allocatable :: buf2(:,:,:),buf(:,:,:),array(:,:,:),phi(:,:,:)
! common /coefs/coef(pot3jmax2,pot3kmax2,lmax2,2),coef2(pot3jmax2,pot3kmax2,lmax)
!  integer, PARAMETER::JKM1=2*POT3JMAX+POT3KMAX-1
  real    :: denny(hj2,hk2)
  integer :: bufSize
  real    :: cvheat, curlyr, delr, delz, den
  real    :: kwfw,phichk,redge,tmass,xmu

integer :: rank, numRanks, nr,tag=17,jread
integer :: maxtrm, isym ,igrid,j,k

!! start MPI here
!Initialize MPI
type(MPI_Status) :: status

! Initialize MPI library
call MPI_Init()

call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)

print *, "numRanks = ",numRanks
ksplit = kmax/numRanks
bufSize = (jmax+2)*(ksplit)*lmax

if (numRanks*ksplit /= kmax) then
   if (rank == 0) then
      print *, "ERROR: z dimension can't be evenly distributed"
      print *, "numRanks is", numRanks, "Nk total is", kmax
   end if
   call MPI_Finalize()
   stop 1
end if

! Allocate arrays
allocate(array(-1:jmax+1,-1:ksplit+1,1:lmax),phi(-1:jmax+1,-1:ksplit+1,1:lmax))
! rank 0 reads in run parameters
if(rank == 0) then
!...  Read in some run parameters.
      OPEN(UNIT=5,FILE='fort.5',STATUS='OLD')

      READ(5,100)    KONST,XN
      WRITE(6,10010) KONST,XN
100   FORMAT(1P2E15.8)
10010 FORMAT('  KONST  ',1PE20.12,/,'  XN     ',1PE20.12)

      READ(5,98)     GAMMA
      WRITE(6,10020) GAMMA
 98   FORMAT(1E25.15)
10020 FORMAT('  GAMMA  ',1PE20.12)

      READ(5,102)    ITSTRT,ITSTOP,IDIAG,ISOADI,ITYPE,NMODL,ISTOR,IGRID
      WRITE(6,10030) ITSTRT,ITSTOP,IDIAG,ISOADI,ITYPE,NMODL,ISTOR,IGRID,JMIN
102   FORMAT(9(I6,1X))
10030 FORMAT('  ITSTRT ',I8,/,'  ITSTOP ',I8,/,'  IDIAG  ',I8,      &
     &     /,'  ISOADI ',I8,/,'  ITYPE  ',I8,/,'  NMODL  ',I8,      &
     &     /,'  ISTOR  ',I8,/,'  IGRID  ',I8,/,'  JMIN   ',I8)

      READ(5,100) NPRIME
      WRITE(6,10040) NPRIME
10040 FORMAT('  NPRIME ',1PE20.12)

      CLOSE(5)
!send itype to other ranks
  do nr = 1,numranks-1
   call MPI_Send(itype,1,MPI_INTEGER,nr,tag,MPI_COMM_WORLD)
  end do
else
  call MPI_Recv(itype,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,status)
end if
!...Define some things..

  MAXTRM=10
  ISYM  =2
  KWFW=int(log10(dble(KMAX))/log10(two))-1 
  CVHEAT=CURLYR/(XMU*(GAMMA-1.0))

  !TMASS=zero
  !ENEW=zero
  !ELOST=zero
  !EDIF=zero
  !PHICHK=zero
  !KLOCAT=0
  !GRAV=one


! rank 0 reads in the data, then distributes as necessary to other ranks.
! use call to seperate subroutine to read in data and pass to other ranks
! call readDensity

  if (ITYPE.EQ.7) then
! Read in the density array, denny
     OPEN(UNIT=2,FILE='fort.2',STATUS='OLD')
     READ(2,*)PINDEX,CON2,RRR2,OMCEN,DENCEN,TOVERW,ROF3N,ZOF3N,A1NEWZ,JREQ,KZPOL

1685 FORMAT(3X,1PE22.15,2X,8(1PE22.15,2X),2I4)
     write(*,1685)PINDEX,CON2,RRR2,OMCEN,DENCEN,TOVERW,ROF3N,ZOF3N,A1NEWZ,JREQ,KZPOL
     READ(2,*) DENNY
1617 FORMAT(8(1PE22.15,2X))
     CLOSE(2)
     DEN=DENCEN

  else if(ITYPE.EQ.1234) then
    call readData(array,jmax,ksplit,lmax)
    call writeData(-1,jmax,ksplit,array(:,:,1),'test')
end if

!!...Calling the potential solver now.
  IPRINT = 1
  REDGE  = 0.d0

! call the potential solver
CALL SETBDY(0,ISYM)
CALL BDYGEN(MAXTRM,ISYM,REDGE,phi,array)

!do k = -1,ksplit+1
!  do j=jmax-2,jmax+1
!    write(100+rank,'(i3,1x,i3,1x,e12.3)') k,j,phi(j,k,1)
!  end do
!end do
CALL POT3MPI(8,IPRINT,jmax,ksplit,lmax,array,dr,phi)


  call writeData(1,jmax,ksplit,phi(:,:,1),'fin')


! Shutdown MPI
call MPI_Finalize()

end program poissonSolverMPI

