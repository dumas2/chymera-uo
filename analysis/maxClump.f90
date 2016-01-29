program maxClump
 implicit none
!
! creating the map.  rconv converts from cells to AU.
! dx is zof3n
!
 integer :: JMAX=256, KMAX=64, LMAX=512, YMAX, XMAX
 integer :: II, JJ, IRR, IAN, IAN2, I, L, J, K
 integer :: block=10,IJ,IK,IL,LL

 real(KIND=8), allocatable, dimension(:,:,:) :: array
 real(KIND=8), allocatable, dimension(:) :: den,rad,massr

 real(KIND=8) :: dx,rconv,torp,sconv,maxden,mass
 real(KIND=8) :: y1, y2, y3, y4, t, u, angle, dphi, pi, ir, jr, rr,time, bg
 character :: filein*72, fileout*72,filenum*6,dum(13)*72

  CALL GetArg(1,dum(1))
  CALL GetArg(2,dum(2))

  filein = trim(dum(1)) 
  dum(2) = trim(dum(2)) 

  read(dum(2),"(f15.8)")dx

  print *, " SIGPLOT OUT: dx -> ", dx
  print *, " SIGPLOT OUT: filein prefix -> ", trim(dum(1))

  allocate(array(-1:JMAX,-1:KMAX,0:LMAX-1))
  allocate(rad(0:(2*block+1)**3))
  allocate(den(0:(2*block+1)**3))
  allocate(massr(0:(2*block+1)**3))

 OPEN(UNIT=12, FILE=trim(filein),FORM='UNFORMATTED', STATUS="OLD")

 read(12)array
 read(12)time
 close(12)

 print "(a,1pe9.2,a,I6.6,a)", " THE TIME IS ", (time), " FOR FILE ",I,"."

 pi = acos(-1d0)
 dphi=2d0*pi/dble(LMAX)

 maxden = 0d0
 mass=0d0
 do L = 0, LMAX-1
  do K = 0, KMAX-1
    do J = 0, JMAX-1
      mass=mass+array(J,K,L)*2d0*(J*1.+.5)*dx**3*dphi
      if (array(J,K,L)>maxden)then
        maxden = array(J,K,L)
        IJ=J; IK=K; IL=L
      endif
    enddo
  enddo
 enddo

 print *, "Total mass", mass

 mass = 0d0
 I=0
 do L=IL-block,IL+block
  LL=L
  if (L<0)LL=L+LMAX
  if (L>LMAX-1)LL=L-LMAX
  do K= 0,block
    do J=IJ-block,IJ+block
      if (array(J,K,LL)>0.0*maxden) mass=mass+array(J,K,LL)*2d0*(J*1.+.5)*dx**3*dphi
      den(I)=array(J,K,LL)
      rad(I)=sqrt(J*J*1.+K*K*1.+IJ*IJ*1.-2.*J*IJ*cos((IL-LL)*dphi))
      massr(I)=mass
      I=I+1
    enddo
  enddo
 enddo     

 print *, "#",mass*1d3,IJ,IK,IL
 print *, "#",IJ*cos(IL*dphi)+256,IJ*sin(IL*dphi)+256
 do J=0,I-1
   print *, rad(J),den(J),massr(J)
 enddo


 deallocate(array,rad,den,massr)

 stop
end
  

    
