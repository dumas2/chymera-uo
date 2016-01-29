! quick hack for freeing the star.  It'll wiggle -- hence the name of the routine.
! to avoid using include files, I just pass what is required.  
! This is supposed to be a second-order leap frog method.  The idea behind 
! the leap frog is that it is time-reversible and conserves energy.
! Hopefully I wrote it correctly.
! rho = density array
! rhf = half cell radius
! vx, vy, fx, and fy are the velocities and forces.  Orginal implemenation required
! them to be arrays, but this method does not.  The arrays could be removed.
! r = star's radial position
! phi = star's phi position
! simply account for the star's potential by adding
! phi_star = -mass_star/sqrt(r^2 + r'^2 -2 r r' cos (phi-phi') + z^2), where the ' denotes
! the star's position and the other coordinates are for the point of interest.
subroutine wiggle(rho,rhf,zhf,vx,vy,fx,fy,delt,pi,r,phi,dr,dz,dphi,gsoft,JMAX,KMAX,LMAX,JMIN,start_star,restart_wiggle)
  implicit none
  integer, intent(in) :: JMAX, KMAX, LMAX, JMIN
  real(KIND=8), intent(in) :: rho(JMAX+2,KMAX+2,LMAX),rhf(JMAX+2),zhf(KMAX+2)
  real(KIND=8), intent(in) :: delt,pi,dr,dz,dphi,gsoft
  real(KIND=8) :: vx,vy,fx, fy, r,phi,x,y,dummy
  real(KIND=8) :: corx,cory
  real(KIND=8), save :: odelt
  logical :: start_star,restart_wiggle

  if (START_STAR) then
    if (restart_wiggle) then
       open(unit=989,file='starry_restart.dat')
          read(989,'(8(1X,1pe15.8))') dummy,vx,vy,r, &
          phi,x,y,odelt
       close(989)
    else
     odelt = delt
     x = 0.d0
     y = 0.d0
     r = x*x + y*y
     r = sqrt(r)
     if (x > 0.d0)then
       phi = atan( (y)/(x))
     else if (x < 0.d0 ) then
       phi = atan( (y)/(x))+pi
     else if (x == 0.d0 .and. y == 0.d0) then
       phi = 0.d0
     else if (x == 0.d0) then
         if (y > 0.d0 ) phi = .5d0*pi
         if (y < 0.d0 ) phi = 1.5d0*pi
     else
         if (x > 0.d0 ) phi = 0.d0
         if (x < 0.d0 ) phi = pi
     endif
     if (phi < 0.d0 ) phi = phi + 2.d0*pi
     call find_force(rho,rhf,zhf,x,y,dr,dz,dphi,fx,fy,gsoft,JMAX,KMAX,LMAX,JMIN)
     vy =  0.
     vx =  0.
!     PRINT *, "STOP. THINK. DO YOU KNOW WHAT YOU ARE DOING? SET UP THE STAR'S IC."
!     STOP
    endif
    return
  endif

  ! find Cart coordinates
  x = r*cos(phi)
  y = r*sin(phi)
  ! find current fx and fy

  call find_force(rho,rhf,zhf,x,y,dr,dz,dphi,fx,fy,gsoft,JMAX,KMAX,LMAX,JMIN)
  vx = vx + fx*odelt*.5d0
  vy = vy + fy*odelt*.5d0
  vx = vx + fx*delt*.5d0
  vy = vy + fy*delt*.5d0
  corx  = x + vx*delt
  cory  = y + vy*delt
  odelt = delt
  
  r = sqrt(corx*corx+cory*cory) ! convert back to polar coordinates

       if (corx > 0.d0)then
         phi = atan( (cory)/(corx))
       else if (corx < 0.d0 ) then
         phi = atan( (cory)/(corx))+pi
       else if (corx == 0.d0 .and. cory == 0.d0) then
         phi = 0.d0
       else if (corx == 0.d0) then
           if (cory > 0.d0 ) phi = .5d0*pi
           if (cory < 0.d0 ) phi = 1.5d0*pi
       else
           if (corx > 0.d0 ) phi = 0.d0
           if (corx < 0.d0 ) phi = pi
       endif
       if (phi < 0.d0 ) phi = phi + 2.d0*pi

  return
  
end subroutine wiggle

! just assume that the force is well-described by approximating
! all of the mass in a cell as a point mass at the cell center.
! in retrospect, maybe it would be better to consider the mass
! at the half-mass radial center, but that should be a minor difference.
subroutine find_force(rho,rhf,zhf,x,y,dr,dz,dphi,fx,fy,gsoft,JMAX,KMAX,LMAX,JMIN)
  implicit none
  integer :: JMAX, KMAX, LMAX, JMIN
  integer :: J,K,L
  real(KIND=8), intent(in) :: rho(JMAX+2,KMAX+2,LMAX),rhf(JMAX+2),zhf(KMAX+2)
  real(KIND=8), intent(in) :: dr,dz,dphi,x,y,gsoft
  real(KIND=8) :: vol,lphi,tmass,xx,yy,fx,fy,dd,zz

  vol = dr*dz*dphi*2.d0 ! The extra 2 accounts for both sides of the disk.
  fx = 0.d0
  fy = 0.d0
  tmass = 0.d0

!$OMP PARALLEL
!$OMP DO SCHEDULE(STATIC) PRIVATE(tmass,dd,lphi,xx,yy,zz,J,K,L) REDUCTION(+:fy,fx)
  do L = 1, LMAX
   lphi = L*dphi-.5d0*dphi
   do K = 2, KMAX+1
     zz = zhf(K)
     do J = JMIN-3, JMAX+1
       xx = rhf(J)*cos(lphi)-x
       yy = rhf(J)*sin(lphi)-y
                  
       tmass = rho(J,K,L)*rhf(J)*vol
       dd = 1.d0/sqrt( xx**2 + yy**2 +zz**2+gsoft**2)**3
       fx = fx + tmass*dd*xx
       fy = fy + tmass*dd*yy

     enddo
   enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  return

end subroutine find_force 
 
 
