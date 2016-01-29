! takes density file and creates output for the
! gnuplot utility. Meridional slice.
program meridional
 implicit none
!
 integer :: JMAX=256, KMAX=64, LMAX=512, YMAX, XMAX
 integer :: ISTART=056000, IEND, ISKIP, USE_LOG=1
 integer :: II, JJ, IRR, IAN, IAN2, I, L, J, K,LEVAL=1,count

 real(KIND=8), allocatable, dimension(:,:,:) :: rho,temp
 real(KIND=8), allocatable, dimension(:) :: longrho,longtemp,size

 real(KIND=8) :: rconv=2.,torp,val,limit=1d-9,bkmpcode,rof3n=1.,dtheta,factor=4d0
 real(KIND=8) :: y1, y2, y3, y4, t, u, angle, dphi, pi, ir, jr, rr,time,iir,jjr
 character :: filein*72, filein2*72, fileout*72,filenum*6,dum(12)*72

 bkmpcode=8.254e7*1.5e13/(6.67e-8*2e33)
 bkmpcode=9.30476132E-06

 allocate(rho (-1:JMAX,-1:KMAX,0:LMAX-1 ))
 allocate(temp(-1:JMAX,-1:KMAX,0:LMAX-1 ))

  write (filenum,'(I6.6)')ISTART
  filein='../run/rho3d.'//filenum//' ' 
  filein2='../run//temperat3d.'//filenum//' ' 
  fileout='rhoT.'//filenum//' ' 
  print "(a,1x,a)", trim(filein), trim(fileout) 
  OPEN(UNIT=12, FILE=trim(filein),FORM='UNFORMATTED', STATUS="OLD")
  OPEN(UNIT=13, FILE=trim(filein2),FORM='UNFORMATTED', STATUS="OLD")

 read(12)rho
 read(12)time
 close(12)

 read(13)temp
 close(13)

 open(unit=14,file=trim(fileout))

 print "(a,I9,a,I6.6,a)", " THE TIME IS ", int(time/torp), " FOR FILE ",I,"."

 pi = acos(-1.d0)
 dphi = 2.d0*pi/dble(LMAX)

 count=0
 do L=0,LMAX-1
  do K=0,KMAX
   do J=0,JMAX
      if(rho(J,K,L)>limit)count=count+1
   enddo
  enddo
 enddo

 allocate(longrho(count))
 allocate(longtemp(count))
 allocate(size(count))

 count=0
 do L=0,LMAX-1
  do K=0,KMAX
   do J=0,JMAX
      if(rho(J,K,L)>limit)then 
         count=count+1
         longrho(count)=rho(J,K,L)
         longtemp(count)=temp(J,K,L)
         !size(count)=max(rof3n,dphi*(dble(J)+.5d0)*rof3n)*factor
         size(count)=rof3n*factor
      endif
   enddo
  enddo
 enddo

 do L=1,count
   write(14,'(5(1X,1pe15.8))')longrho(L),longtemp(L), &
      size(L)**2/pi*longrho(L)*2.33/bkmpcode, &
      size(L)**2/pi*longrho(L)**2,longrho(L)*longtemp(L)&
      *bkmpcode/2.33
 enddo

 close(14)

 deallocate(rho,temp,longrho,longtemp,size)

 stop
end
  

    
