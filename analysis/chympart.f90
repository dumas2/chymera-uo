program chym2part
  implicit none
  
  integer, parameter :: jmax=256,kmax=64,lmax=512
  integer :: j,k,l,I
  
  real*8, parameter :: msun=1.989d33,au=1.496d13,rgas=8.254d7,grav=6.67e-8
  real*8, dimension(-1:jmax,-1:kmax,0:lmax-1) :: s,t,a,rho,eps
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

  do j = -1,jmax
   rhf(j) = (dble(j)+0.5d0)*dx 
  enddo

  dtheta = acos(-1.)/dble(lmax)*2.d0
  rholmt = den*1.d-7

  do l = 0,lmax-1
   angle(l) = (dble(l)+.5d0)*dtheta
  enddo

  tmass = 0.d0
  i = 0
  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
    if (rho(j,k,l) > rholmt) then
      dummy = rho(j,k,l)*rhf(j)*dx**2*dtheta
      print "(1X,1pe15.8)",dummy
      tmass = tmass + dummy
      i = i+1
    endif
  enddo
  enddo
  enddo

  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
     if (rho(j,k,l) > rholmt) then
       dummy = rhf(j)*cos(angle(l))
       print "(1X,1pe15.8)",dummy
     endif
  enddo
  enddo
  enddo

  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
     if (rho(j,k,l) > rholmt) then
       dummy = rhf(j)*sin(angle(l))
       print "(1X,1pe15.8)",dummy
     endif
  enddo
  enddo
  enddo

  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
     if (rho(j,k,l) > rholmt) then
       dummy = (dble(k)+.5)*dx
       print "(1X,1pe15.8)",dummy
     endif
  enddo
  enddo
  enddo

  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
     if (rho(j,k,l) > rholmt) then
       dummy = (s(j,k,l)+s(j+1,k,l))/(2.d0*rho(j,k,l))*cos(angle(l))
       dummy = dummy - a(j,k,l)/(rhf(j)*rho(j,k,l))*sin(angle(l))
       print "(1X,1pe15.8)",dummy*vconv
     endif
  enddo
  enddo
  enddo

  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
     if (rho(j,k,l) > rholmt) then
       dummy = (s(j,k,l)+s(j+1,k,l))/(2.d0*rho(j,k,l))*sin(angle(l))
       dummy = dummy + a(j,k,l)/(rhf(j)*rho(j,k,l))*cos(angle(l))
       print "(1X,1pe15.8)",dummy*vconv
     endif
  enddo
  enddo
  enddo

  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
     if (rho(j,k,l) > rholmt) then
       dummy = (t(j,k,l)+t(j,k+1,l))/(rho(j,k,l)*2.d0)
       print "(1X,1pe15.8)",dummy*vconv
     endif
  enddo
  enddo
  enddo

  do l = 0,lmax-1
  do k = 0,kmax-1
  do j = 0,jmax-1
     if (rho(j,k,l) > rholmt) then
       dummy = eps(j,k,l)/rho(j,k,l)
       print "(1X,1pe15.8)",dummy*vconv**2
     endif
  enddo
  enddo
  enddo

  print *, "# Total mass (half disk) and full disk mass: ",tmass,tmass*2.d0
  print *, "#",I," Data ponts. order and UNITS as follows: "
  print *, "#  mass (Msun)"
  print *, "#  x (AU)"
  print *, "#  y (AU)"
  print *, "#  z (AU)"
  print *, "#  vx (cm/s)"
  print *, "#  vy (cm/s)"
  print *, "#  vz (cm/s)"
  print *, "#  e (cm^2/s^2)"
   

  stop
end

