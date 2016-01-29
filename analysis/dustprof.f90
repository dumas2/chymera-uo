! use simple bilinear interpolation to make surface density maps.
! takes density file and creates output for the
! the shell script makesig. Search for SWITCH to
! look for any hardcoded switch in this program.
program sigplotter
 implicit none
!
! creating the map.  rconv converts from cells to AU.
! dx is zof3n
!
 integer :: JMAX, KMAX, LMAX, YMAX, XMAX
 integer :: ISTART, IEND, ISKIP,IDUM,JMAX1,KMAX1,LMAX1
 integer :: II, JJ, IRR, IAN, IAN2, I, L, J, K,ERR

 real(KIND=8), allocatable, dimension(:,:,:) :: array,dust,tempk
 real(KIND=8), allocatable, dimension(:) :: array1d,dust1d,press

 real(KIND=8) :: dx,rconv,torp,sconv,mass,dmass
 real(KIND=8) :: y1, y2, y3, y4, t, u, angle, dphi, pi, ir, jr, rr,time, bg
 character :: filein*72, fileout*72,filenum*6,dum(15)*72,dustin*72,tempin*72
 pi = acos(-1d0)

  CALL GetArg(1,dum(1))
  CALL GetArg(2,dum(2))
  CALL GetArg(3,dum(3))
  CALL GetArg(4,dum(4))
  CALL GetArg(5,dum(5))
  CALL GetArg(6,dum(6))
  CALL GetArg(7,dum(7))
  CALL GetArg(8,dum(8))
  CALL GetArg(9,dum(9))
  CALL GetArg(10,dum(10))
  CALL GetArg(11,dum(11))
  CALL GetArg(12,dum(12))
  CALL GetArg(13,dum(13))
  CALL GetArg(14,dum(14))
  CALL GetArg(15,dum(15))

  if ( len(trim(dum(1))) == 0 .or. &
       len(trim(dum(2))) == 0 .or. &
       len(trim(dum(3))) == 0 .or. &
       len(trim(dum(4))) == 0 .or. &
       len(trim(dum(5))) == 0 .or. &
       len(trim(dum(6))) == 0 .or. &
       len(trim(dum(7))) == 0 .or. &
       len(trim(dum(8))) == 0 .or. &
       len(trim(dum(9))) == 0 .or. &
       len(trim(dum(10)))== 0 .or. &
       len(trim(dum(11)))== 0 .or. &
       len(trim(dum(12)))== 0 .or. &
       len(trim(dum(13)))== 0 .or. &
       len(trim(dum(14)))== 0 .or. &
       len(trim(dum(15)))== 0        ) then
    print *, "A CHYMTOOL SPONSORED EVENT"
    print *, " "
    print *, " "
    print *, "sigplot v 1. 11.10.2007. Nora Bolig"
    print *, " "
    print *, " "
    print *, "SIGPLOT USAGE: sigplot {stepbegin} {stepend} {stepskip}"
    print *, "                       {jmax} {kmax} {lmax} {torp} {dz}"
    print *, "                       {cell-to-size unit conversion  }"
    print *, "                       {min surface density background}"
    print *, "                       {code to cgs surfden conversion}"
    print *, "                       {filein prefix} {dustin prefix} {fileout prefix}" 
    print *, "                       {tempin } "
    print *, " "
    print *, " "
    print *, "Calling Cthulu...(run)"
    stop
  endif

  print *, "A CHYMTOOL SPONSORED EVENT"
  print *, " "
  print *, " "
  print *, "sigplot v 1. 11.10.2007. Nora Bolig"
  print *, " "
  print *, " "
 
  dum(1) = trim(dum(1)) 
  dum(2) = trim(dum(2)) 
  dum(3) = trim(dum(3)) 
  dum(4) = trim(dum(4)) 
  dum(5) = trim(dum(5)) 
  dum(6) = trim(dum(6)) 
  dum(7) = trim(dum(7)) 
  dum(8) = trim(dum(8)) 
  dum(9) = trim(dum(9)) 
  dum(10) = trim(dum(10)) 
  dum(11) = trim(dum(11)) 

  read(dum(1),"(i8)")ISTART
  read(dum(2),"(i8)")IEND
  read(dum(3),"(i8)")ISKIP
  read(dum(4),"(i8)")JMAX
  read(dum(5),"(i8)")KMAX
  read(dum(6),"(i8)")LMAX
  read(dum(7),"(f15.8)")torp
  read(dum(8),"(f15.8)")dx
  read(dum(9),"(f15.8)")rconv
  read(dum(10),"(f15.8)")bg
  read(dum(11),"(f15.8)")sconv

  print *, " SIGPLOT OUT: ISTART -> ", ISTART
  print *, " SIGPLOT OUT: IEND -> ", IEND
  print *, " SIGPLOT OUT: ISKIP -> ", ISKIP
  print *, " SIGPLOT OUT: JMAX -> ", JMAX
  print *, " SIGPLOT OUT: KMAX -> ", KMAX
  print *, " SIGPLOT OUT: LMAX -> ", LMAX
  print *, " SIGPLOT OUT: torp -> ", torp
  print *, " SIGPLOT OUT: dz -> ", dx
  print *, " SIGPLOT OUT: r conversion -> ", rconv
  print *, " SIGPLOT OUT: background -> ", bg
  print *, " SIGPLOT OUT: surfden conversion -> ", sconv
  print *, " SIGPLOT OUT: filein prefix -> ", trim(dum(12))
  print *, " SIGPLOT OUT: dustin prefix -> ", trim(dum(13))
  print *, " SIGPLOT OUT: fileout prefix-> ", trim(dum(14))
  print *, " SIGPLOT OUT: tempin  prefix-> ", trim(dum(15))

  allocate(array(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(dust(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(tempk(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(array1d(-1:JMAX)      )
  allocate(dust1d(-1:JMAX)      )
  allocate(press(-1:JMAX)      )

 I = ISTART
 do while ( I <= IEND )
    write (filenum,'(I6.6)')I
    filein=trim(dum(12))//filenum//' ' 
    dustin=trim(dum(13))//filenum//' ' 
    fileout=trim(dum(14))//filenum//' ' 
    tempin=trim(dum(15))//filenum//' ' 
    print "(4(1x,a))", trim(filein), trim(dustin), trim(fileout) , trim(tempin)
    OPEN(UNIT=12, FILE=trim(filein),FORM='UNFORMATTED', STATUS="OLD")

 read(12,iostat=ERR)array
 read(12)time
 close(12)
 print "(a,1pe9.2,a,I6.6,a)", " THE TIME IS ", (time/torp), " FOR FILE ",I,"."
if (ERR/=0)then
   print *, "found error ",ERR," during read of ",filein
   exit
 endif

    OPEN(UNIT=15, FILE=trim(tempin),FORM='UNFORMATTED', STATUS="OLD")
read(15,iostat=ERR)tempk
 read(15)time
 close(15)
 print "(a,1pe9.2,a,I6.6,a)", " THE TIME IS ", (time/torp), " FOR FILE ",I,"."
if (ERR/=0)then
   print *, "found error ",ERR," during read of ",tempin
   exit
 endif

 OPEN(UNIT=13, FILE=trim(dustin),FORM='UNFORMATTED', STATUS="OLD")

 read(13)jmax1,kmax1,lmax1,idum
 print *,jmax1,kmax1,lmax1,idum
 call assert(jmax1==jmax)
 call assert(kmax1==kmax)
 call assert(lmax1==lmax)
 read(13,iostat=ERR)dust,dust,dust,dust
 close(13)

 if (ERR/=0)then
   print *, "found error ",ERR," during read of ",dustin 
   exit
 endif

 open(unit=14,file=trim(fileout))
 dphi=2d0*acos(-1d0)/dble(LMAX)

 array1d = 0.d0
 dust1d =0d0 
  press=0d0
 do L = 0, LMAX-1
  do K = 0,  KMAX-1
    do J = 0, JMAX-1
      array1d(J) = array(J,K,L)*dx*2.d0/dble(LMAX)*sconv + array1d(J) ! take into account both sides of disk
      dust1d(J) = dust(J,K,L)*dx*2.d0/dble(LMAX)*sconv + dust1d(J) ! take into account both sides of disk
      press(J) = (sconv*array(J,K,L))**2*dx*2d0/dble(LMAX)*tempk(J,K,L)+press(J)
    enddo
  enddo
 enddo

 dmass=0d0
 mass=0d0
 do J=0,JMAX-1
   press(J)=press(J)/array1d(J)
   write(14,"(4(1X,1pe15.8))") (dble(J)+.5)*dx,dust1d(J),array1d(J),press(J)
   dmass=dmass+dust1d(J)*2d0*pi*(dble(J)+.5d0)*dx*dx
   mass=mass+array1d(J)*2d0*pi*(dble(J)+.5d0)*dx*dx
 enddo
 I=I+ISKIP
 print *, "Total mass is in gas is ",dmass," total mass in gas is ", mass, " and ratio ",dmass/mass

 enddo ! end while loop

 deallocate(array,array1d,dust,dust1d,press)

 stop
end
  
subroutine assert(logic)
 implicit none
 logical::logic

 if(.not.logic)then
  print *, "assert not met. stopping."
  stop
 endif

 return
end subroutine

    
