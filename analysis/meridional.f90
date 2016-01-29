! takes density file and creates output for the
! gnuplot utility. Meridional slice.
program meridional
 implicit none
!
 integer :: JMAX, KMAX, LMAX, YMAX, XMAX
 integer :: ISTART, IEND, ISKIP, USE_LOG
 integer :: II, JJ, IRR, IAN, IAN2, I, L, J, K,LEVAL

 real(KIND=8), allocatable, dimension(:,:,:) :: array

 real(KIND=8) :: rconv,torp,val
 real(KIND=8) :: y1, y2, y3, y4, t, u, angle, dphi, pi, ir, jr, rr,time,iir,jjr
 character :: filein*72, fileout*72,filenum*6,dum(12)*72

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
       len(trim(dum(12)))== 0        ) then
    print *, "MERIDIONAL USAGE: meridional {stepbegin} {stepend} {stepskip}"
    print *, "                             {jmax}{kmax}{lmax}{leval}{torp} "
    print *, "                             {cell-to-size unit conversion  }"
    print *, "                             {use logaritmic scaling [0|1]  }"
    print *, "                             {filein prefix} {fileout prefix}" 
    stop
  endif

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

  read(dum(1),"(i8)")ISTART
  read(dum(2),"(i8)")IEND
  read(dum(3),"(i8)")ISKIP
  read(dum(4),"(i8)")JMAX
  read(dum(5),"(i8)")KMAX
  read(dum(6),"(i8)")LMAX
  read(dum(7),"(i8)")LEVAL
  read(dum(8),"(f15.8)")torp
  read(dum(9),"(f15.8)")rconv
  read(dum(10),"(i8)")USE_LOG

  print *, " MERIDIONAL OUT: ISTART -> ", ISTART
  print *, "  OUT: IEND -> ", IEND
  print *, " MERIDINAL OUT: ISKIP -> ", ISKIP
  print *, " MERIDINAL OUT: JMAX -> ", JMAX
  print *, " MERIDINAL OUT: KMAX -> ", KMAX
  print *, " MERIDINAL OUT: LMAX -> ", LMAX
  print *, " MERIDINAL OUT: LEVAL -> ", LEVAL
  print *, " MERIDINAL OUT: torp -> ", torp
  print *, " MERIDINAL OUT: conversion -> ", rconv
  print *, " MERIDINAL OUT: USE_LOG -> ", USE_LOG
  print *, " MERIDINAL OUT: filein prefix -> ", trim(dum(11))
  print *, " MERIDINAL OUT: fileout prefix-> ", trim(dum(12))

  allocate(array(-1:JMAX,-1:KMAX,0:LMAX-1))

  YMAX = KMAX; XMAX = JMAX

 I = ISTART

 array = 0.d0
 
 do while ( I <= IEND )
    write (filenum,'(I6.6)')I
    filein=trim(dum(11))//filenum//' ' 
    fileout=trim(dum(12))//filenum//' ' 
    print "(a,1x,a)", trim(filein), trim(fileout) 
    OPEN(UNIT=12, FILE=trim(filein),FORM='UNFORMATTED', STATUS="OLD")
    open(unit=13,file=trim(fileout))

 read(12)array
 read(12)time
 close(12)

 print "(a,I9,a,I6.6,a)", " THE TIME IS ", int(time/torp), " FOR FILE ",I,"."

 pi = acos(-1.d0)
 dphi = 2.d0*pi/dble(LMAX)

 ir = -.5
 do II = 0, XMAX-1
  ir = ir + 1.
   jr = -.5
  do JJ = 0, YMAX-1
   jr = jr+1.

   iir = ir
   jjr = jr
   if (jr < 0 ) jjr = jr +1.
   if (ir < 0 ) iir = ir +1.
   IRR = int(ir)
   IAN = int(jr)

   select case (USE_LOG)
   case(0)
   y1 = (array(IRR,IAN,LEVAL))
   y2 = (array(IRR,IAN+1,LEVAL))
   y3 = (array(IRR+1,IAN+1,LEVAL))
   y4 = (array(IRR+1,IAN,LEVAL))
   case(1)
   y1 = LOG10(array(IRR,IAN,LEVAL))
   y2 = LOG10(array(IRR,IAN+1,LEVAL))
   y3 = LOG10(array(IRR+1,IAN+1,LEVAL))
   y4 = LOG10(array(IRR+1,IAN,LEVAL))
   end select

   t = (jjr - IAN*1.d0)
   u = (iir - IRR*1.d0)
   val =  (1.-t)*(1.-u)*y1+t*(1.-u)*y2+t*u*y3+(1.-t)*u*y4
   write(13,"(3(1pe15.8,1X))") (ir)*rconv, (jr)*rconv, val
  enddo
  write(13,"(a)")
 enddo
 close(13)

 I = I+ISKIP
 enddo ! end while loop

 deallocate(array)

 stop
end
  

    
