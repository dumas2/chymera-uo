! takes density file and creates output for the
! gnuplot utility. Meridional slice.
program meridional
 implicit none
!
 integer :: JMAX=256, KMAX=64, LMAX=8, YMAX, XMAX
 integer :: ISTART=010000, IEND, ISKIP, USE_LOG=1
 integer :: II, JJ, IRR, IAN, IAN2, I, L, J, K,LEVAL=1

 real(KIND=8), allocatable, dimension(:,:,:) :: array,array3
 real(KIND=8), allocatable, dimension(:,:,:,:) :: array2

 real(KIND=8) :: rconv=2.,torp,val
 real(KIND=8) :: y1, y2, y3, y4, t, u, angle, dphi, pi, ir, jr, rr,time,iir,jjr
 character :: filein*72, filein2*72, fileout*72,filenum*6,dum(12)*72


 allocate(array(-1:JMAX,-1:KMAX,0:LMAX-1 ))
 allocate(array3(-1:JMAX,-1:KMAX,0:LMAX-1))

  write (filenum,'(I6.6)')ISTART
  filein='../run/rho3d.'//filenum//' ' 
  filein2='../run/passive.'//filenum//' ' 
  fileout='g2d.'//filenum//' ' 
  print "(a,1x,a)", trim(filein), trim(fileout) 
  OPEN(UNIT=12, FILE=trim(filein),FORM='UNFORMATTED', STATUS="OLD")
  OPEN(UNIT=13, FILE=trim(filein2),FORM='UNFORMATTED', STATUS="OLD")

 read(12)array
 read(12)time
 close(12)

 read(13)J,K,L,II
 print *, J,K,L,II
 allocate(array2(-1:J,-1:K,0:L-1,0:II-1))
 read(13)array2
 close(13)

 do L=0,LMAX-1
 do K=-1,KMAX
 do J=-1,JMAX 
 array3(J,K,L)=array2(J,K,L,3)!/array(J,K,L)
 enddo
 enddo
 enddo  

 open(unit=14,file=trim(fileout))

 print "(a,I9,a,I6.6,a)", " THE TIME IS ", int(time/torp), " FOR FILE ",I,"."

 pi = acos(-1.d0)
 dphi = 2.d0*pi/dble(LMAX)

 XMAX=JMAX
 YMAX=KMAX

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
   y1 = (array3(IRR,IAN,LEVAL))
   y2 = (array3(IRR,IAN+1,LEVAL))
   y3 = (array3(IRR+1,IAN+1,LEVAL))
   y4 = (array3(IRR+1,IAN,LEVAL))
   case(1)
   y1 = LOG10(array3(IRR,IAN,LEVAL))
   y2 = LOG10(array3(IRR,IAN+1,LEVAL))
   y3 = LOG10(array3(IRR+1,IAN+1,LEVAL))
   y4 = LOG10(array3(IRR+1,IAN,LEVAL))
   end select

   t = (jjr - IAN*1.d0)
   u = (iir - IRR*1.d0)
   val =  (1.-t)*(1.-u)*y1+t*(1.-u)*y2+t*u*y3+(1.-t)*u*y4
   write(14,"(3(1pe15.8,1X))") (ir)*rconv, (jr)*rconv, val
  enddo
  write(14,"(a)")
 enddo
 close(14)

 deallocate(array,array2,array3)

 stop
end
  

    
