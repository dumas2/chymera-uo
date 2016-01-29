program chym2part
  implicit none
  
  integer, parameter :: jmax=256,kmax=64,lmax=512,nn=jmax*lmax
  integer :: j,k,l,I
  
  real*8, parameter :: msun=1.989d33,au=1.496d13,rgas=8.254d7,grav=6.67e-8
  real*8, dimension(-1:jmax,-1:kmax,0:lmax-1) :: s,t,a,rho,eps
  real*8, dimension(-1:jmax,0:lmax-1) :: mass,xpos,ypos,zpos,vx,vy,vz,cangle
  real*8 :: rhf(-1:jmax),angle(0:lmax-1)

  real*8 :: den,dummy,dx,rholmt,tmass,rhoconv,vconv,tconv,dtheta

  rhoconv = msun/au**3
  tconv = 1.d0/sqrt(rhoconv*grav)
  vconv = au/tconv 

  open(unit=11,file="../run/ic.w4",form='unformatted')

  read(11)s
  read(11)t
  read(11)a
  read(11)rho
  read(11)eps
  read(11)dx,dummy,dummy,dummy,dummy,den
  read(11)dummy,tmass

  rholmt=den*1e-7
  print *, den

  vx=0.;vy=0.;vz=0.;xpos=0.;ypos=0.;zpos=0.


  do j = -1,jmax
   rhf(j) = (dble(j)+0.5d0)*dx 
  enddo

  dtheta = 2.d0*acos(-1.)/dble(lmax)
  print *, "# ",dx,rholmt,dtheta,tmass
  do l = 0,lmax-1
   angle(l) = (dble(l)+.5d0)*dtheta
  enddo

  tmass = 0.d0
  i = 0;k=0
  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
    mass(j,l)=0.
    if(rho(J,K,L)<rholmt)cycle
    mass(j,l) = rho(j,k,l)*rhf(j)*dx**2*dtheta
    tmass = tmass + mass(j,l)
    i = i+1
  enddo
  enddo
  enddo

  i = 0
  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
     xpos(j,l) = rhf(j)*cos(angle(l))
     i = i + 1
  enddo
  enddo
  enddo

  i = 0
  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
     ypos(j,l) = rhf(j)*sin(angle(l))
     i = i + 1
  enddo
  enddo
  enddo

  i = 0
  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
     zpos(j,l) = (dble(k)+.5)*dx
     i = i + 1
  enddo
  enddo
  enddo

  i = 0
  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
       dummy = (s(j,k,l)+s(j+1,k,l))/(2.d0*rho(j,k,l))*cos(angle(l))
       vx(j,l) = dummy - a(j,k,l)/(rhf(j)*rho(j,k,l))*sin(angle(l))
     i = i + 1
  enddo
  enddo
  enddo

  i = 0
  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
       dummy = (s(j,k,l)+s(j+1,k,l))/(2.d0*rho(j,k,l))*sin(angle(l))
       vy(j,l) = dummy + a(j,k,l)/(rhf(j)*rho(j,k,l))*cos(angle(l))
     i = i + 1
  enddo
  enddo
  enddo

  i = 0
  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
       vz(j,l) = (t(j,k,l)+t(j,k+1,l))/(rho(j,k,l)*2.d0)
     i = i + 1
  enddo
  enddo
  enddo

  i = 0
  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
       dummy = eps(j,k,l)/rho(j,k,l)
     i = i + 1
  enddo
  enddo
  enddo

  do l = 0,lmax-1
  do j = 0,jmax-1
    if (rho(J,K,L)>rholmt) then
    print "(9(1X,1pe15.8))" , mass(j,l),rhf(J)*cos(angle(l)),rhf(J)*sin(angle(l)),zpos(j,l),vx(j,l),vy(j,l),vz(j,l), &
         -sqrt(.7d0)/sqrt(rhf(J))*sin(angle(l)), &
          sqrt(.7d0)/sqrt(rhf(J))*cos(angle(l))
    endif
  enddo
  enddo

  print *, "# half and total disk mass: ",tmass,2d0*tmass
  print *, "# number of data points: ", i
!  print *, "# Data order and UNITS: "
!  print *, "#  mass (Msun)"
!  print *, "#  x (AU)"
!  print *, "#  y (AU)"
!  print *, "#  z (AU)"
!  print *, "#  vx (cm/s)"
!  print *, "#  vy (cm/s)"
!  print *, "#  vz (cm/s)"
!  print *, "#  e (cm^2/s^2)"
   

  stop
end

