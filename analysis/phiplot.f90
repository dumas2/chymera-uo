! use simple bilinear interpolation to make surface density maps.
! takes density file and creates output for the
! the shell script makesig. Search for SWITCH to
! look for any hardcoded switch in this program.
program phiplotter
 implicit none
!
! creating the map.  rconv converts from cells to AU.
! dx is zof3n
!
 integer :: JMAX, KMAX, LMAX, YMAX, XMAX
 integer :: ISTART, IEND, ISKIP
 integer :: II, JJ, IRR, IAN, IAN2, I, L, J, K

 real(KIND=8), allocatable, dimension(:,:,:) :: array
 real(KIND=8), allocatable, dimension(:,:) :: array2d

 real(KIND=8) :: dx,rconv,torp,sconv
 real(KIND=8) :: y1, y2, y3, y4, t, u, angle, dphi, pi, ir, jr, rr,time, bg
 character :: filein*72, fileout*72,filenum*6,dum(13)*72

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
       len(trim(dum(13)))== 0        ) then
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
    print *, "                       {filein prefix} {fileout prefix}" 
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
  print *, " SIGPLOT OUT: fileout prefix-> ", trim(dum(13))

  allocate(array(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(array2d(-1:JMAX,0:LMAX-1)      )

  YMAX = 2*JMAX; XMAX = 2*JMAX
 pi = acos(-1.d0)
 I = ISTART
 do while ( I <= IEND )
    write (filenum,'(I6.6)')I
    filein=trim(dum(12))//filenum//' ' 
    fileout=trim(dum(13))//filenum//' ' 
    print "(a,1x,a)", trim(filein), trim(fileout) 
    OPEN(UNIT=12, FILE=trim(filein),FORM='UNFORMATTED', STATUS="OLD")
    open(unit=13,file=trim(fileout))

 read(12)array
 read(12)time
 close(12)

 print "(a,1pe9.2,a,I6.6,a)", " THE TIME IS ", (time*torp), " FOR FILE ",I,"."
 k=2
 do L = 0, LMAX-1
    do J = 0, JMAX-1
     write(13,"(1pe11.4,1x,1pe11.4,1x,1pe15.8)") j*dx, l*2.d0*pi/(LMAX-1), array(j,k,l)
     enddo
     write(13,"(a)") " "
 enddo
 close(13)

 I = I+ISKIP
 enddo ! end while loop
!
 deallocate(array,array2d)

 stop
end
  

    
