      module particle
      implicit none
#include  "units.h"
#include "hydroparam.h"
#include "globals.h"  

      logical,save::read_particle_file=.true.
      integer,save::P_FILEID=40000
      
      integer,save::NPARTICLE=100000
      logical,save,allocatable,dimension(:)::particle_skip

      real*8,save,allocatable,dimension(:)::x_p,y_p,z_p ,mass_p
      real*8,save,allocatable,dimension(:)::vx_p,vz_p,vy_p
!!      real*8,save,allocatable,dimension(:)::Ugas,Wgas,Ogas,r_ph
!!      real*8,save,allocatable,dimension(:)::new_dv_a,new_dv_r,new_dv_z


      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_cylinder(x,y,r,an)

      real*8, intent(in):: x,y
      real*8, intent(out)::r,an

      r=sqrt(x*x+y*y)
      if(x==zero)then
       an=pi*half 
       if(y<zero)an=1.5d0*pi
      elseif(y==zero)then
       an=zero
       if(x<zero)an=pi
      else
       an=atan(y/x)
       if(x<zero)then
         an=an+pi
       elseif(y<zero)then
         an=an+two*pi
       endif
      endif

      return

      end subroutine get_cylinder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_cyl_vel(vx,vy,r,an,vr,o)
 
      real*8,intent(in):: vx,vy,an,r
      real*8,intent(out)::vr,o

      vr=vx*cos(an)+vy*sin(an)
      o =(vy*cos(an)-vx*sin(an))/r

      return
      end subroutine get_cyl_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_cart_vel(vr,o,vx,vy,r,an)

      real*8,intent(in) ::vr,o,r,an
      real*8,intent(out)::vx,vy

      vx=vr*cos(an)-o*r*sin(an)
      vy=vr*sin(an)+o*r*cos(an)

      return
      end subroutine get_cart_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_cartesian(r,an,x,y)

      real*8,intent(in )::r,an
      real*8,intent(out)::x,y

      x=r*cos(an)
      y=r*sin(an)

      return
      end subroutine get_cartesian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine set_particle_vel()

      integer::thisj,thisk,thisllow,thislhi,thisl,I
      real*8::tmp,ran4,mag,r_p,angle_p,vr_p,omega_p
  
      do i=1,nparticle

         if(particle_skip(I))cycle

         call get_cylinder(x_p(I),y_p(I),r_p,angle_p)
         if(i==1)print *, x_p(I),y_p(I),r_p,angle_p


         thisJ=int(r_p/rof3n)+2
         thisK=int(z_p(I)/zof3n)+2
         thisL=int(angle_p/dtheta)+1

         omega_p= ( phi(thisJ+1,thisK,thisL)-phi(thisJ-1,thisK,thisL) )
     &          / ( two*rof3n )*r_p

         omega_p=sqrt(omega_p)/r_p

         tmp=(two*ran4(i) -one)
         mag=0.1*ran4(i)
         if(tmp/=zero)then
           vr_p=abs(tmp)/tmp*mag*omega_p*r_p
           vz_p(I)=abs(tmp)/tmp*mag*omega_p*r_p
         else
           vz_p(i)=zero
           vr_p=zero
         endif 
         call get_cart_vel(vr_p,omega_p,vx_p(I),vy_p(I),r_p,angle_p)

      enddo 
      return
      end subroutine set_particle_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine initialize_particles() 

      integer::SEED
      integer*8::N,NF,J,K,L,I,ITER,NSTART
      real*8::dev_r,tmass,ran4,limiter,max_cell_mass
      real*8::r_p,angle_p


      logical,allocatable,dimension(:)::flag
      integer,allocatable,dimension(:)::JCELL,KCELL,LCELL
      real*8, allocatable,dimension(:)::mass_cell

      limiter=den*phylim*1d4

      allocate(x_p    (NPARTICLE))
      allocate(y_p    (NPARTICLE))
      allocate(z_p    (NPARTICLE))
      allocate(mass_p (NPARTICLE))
      allocate(vx_p   (NPARTICLE))
      allocate(vy_p   (NPARTICLE))
      allocate(vz_p   (NPARTICLE))
!      allocate(Ugas   (NPARTICLE))
!      allocate(Wgas   (NPARTICLE))
!      allocate(Ogas   (NPARTICLE))
!      allocate(new_dv_a(NPARTICLE))
!      allocate(new_dv_r(NPARTICLE))
!      allocate(new_dv_z(NPARTICLE))

      allocate(particle_skip(NPARTICLE))
      
      !-----------------------
      ! calculate mass in cell
      !-----------------------
      N=0
      tmass=zero
      do L = 1, LMAX
        do K = 2, KMAX1
          do J = 2, JMAX1
            if(rho(J,K,L)>limiter)then
                N=N+1
                tmass=tmass+rho(J,K,L)*rhf(J)*rof3n*zof3n*dtheta!*LMAX ! only half mass
             endif
          enddo
        enddo
      enddo
      print *," Half mass for placing particles is ",tmass

      allocate(jcell(N)    )
      allocate(kcell(N)    )
      allocate(lcell(N)    )
      allocate(flag(N)     )
      allocate(mass_cell(N))

      I=1
      max_cell_mass=zero
      do L = 1,LMAX
        do K = 2,KMAX1
          do J = JMIN,JMAX1
            if(rho(J,K,L)>limiter)then
              jcell(I)=J
              kcell(I)=K
              lcell(I)=L 
              flag(I)=.true.
              mass_cell(I)=rho(J,K,L)*rhf(J)*rof3n*zof3n*dtheta!*LMAX
              max_cell_mass=max(max_cell_mass,mass_cell(I))
              I=I+1
            endif
          enddo
        enddo
      enddo
         
      SEED=0
      NF=NPARTICLE
      do while(NF>0)
        ITER=ran4(seed)*N
        I=ITER;if(ITER>N)I=I-N
        dev_r= (ran4(SEED)); SEED=SEED+1
        if(dev_r<=mass_cell(I)/tmass*1d2)then
            mass_p(NF)=dust_to_gas*tmass/dble(NPARTICLE)
            dev_r= (ran4(SEED)); SEED=SEED+1
            !print *, dev_r*rof3n,r(jcell(I))+dev_r*rof3n
            r_p=r(jcell(I))+dev_r*rof3n
            dev_r= (ran4(SEED)); SEED=SEED+1
            z_p(NF)=z(kcell(I))+dev_r*zof3n
            dev_r= (ran4(SEED)); SEED=SEED+1
            !angle_p=dev_r*dtheta*LMAX!dtheta*(dev_r+dble(lcell(I)-1))
            angle_p=dtheta*(dev_r+dble(lcell(I)-1))
            particle_skip(NF)=.false.
            call get_cartesian(r_p,angle_p,x_p(NF),y_p(NF))
            print *, NF, " Particles to go."
            NF=NF-1
         endif 
      enddo
!      r_p=r(jmin+10)+.1*rof3n
!      z_p(1)=z(4)+.1*zof3n
!      angle_p=pi
!      call get_cartesian(r_p,angle_p,x_p(1),y_p(1))
!      r_p=r(75)+.1*rof3n
!      z_p(2)=z(4)+.1*zof3n
!      angle_p=pi
!      call get_cartesian(r_p,angle_p,x_p(2),y_p(2))
      deallocate(jcell,kcell,lcell,flag,mass_cell)
      call set_particle_vel() 

      return
      end subroutine initialize_particles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine update_particle_vel(dt)

      integer::thisJ,thisK,thisLlow,thisLhi,thisL,I,J,K,L
      real*8::dt,frbottom,frtop,fzbottom,fztop,fabottom,fatop
      real*8::accr,accz,acca,accx,accy
      real*8::r_p,angle_p,zalpha,xstar,ystar

#if WIGGLE>0
      xstar=rpstar*cos(phi_star)
      ystar=rpstar*sin(phi_star)
#else
      xstar=zero
      ystar=zero
#endif


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,thisj,thisk,thisl,thisLlow,
!$OMP& thislhi,frbottom,frtop,fzbottom,fztop,fabottom,fatop,accr,accz,
!$OMP& acca,r_p,angle_p,accx,accy,zalpha)
!$OMP DO SCHEDULE(STATIC)
      do L=1,LMAX
      do K=1,KMAX2
      do J=1,JMAX2
        phi(J,K,L)=phi(J,K,L)-starphi(J,K,L)
      enddo
      enddo
      enddo
!$OMP ENDDO
!$OMP DO SCHEDULE(DYNAMIC)
      do I=1,NPARTICLE

         if(particle_skip(I))cycle

         call get_cylinder(x_p(I),y_p(I),r_p,angle_p)

         thisJ=int(r_p/rof3n)+2
         thisK=int(z_p(I)/zof3n)+2
         thisL=int(angle_p/dtheta)+1 
         thisLlow=thisL-1
         thisLhi =thisL+1
         if(thisLlow<1  )thisLlow=LMAX
         if(thisLhi>LMAX)thisLhi =1
      
         frbottom=-(phi(thisJ,thisK,thisL)-phi(thisJ-1,thisK,thisL))
     &           /rof3n
         frtop   =-(phi(thisJ+1,thisK,thisL)-phi(thisJ,thisK,thisL))
     &           /rof3n
         fzbottom=-(phi(thisJ,thisK,thisL)-phi(thisJ,thisK-1,thisL))
     &           /zof3n
         fztop   =-(phi(thisJ,thisK+1,thisL)-phi(thisJ,thisK,thisL))
     &           /zof3n
         fabottom=-(phi(thisJ,thisK,thisL)-phi(thisJ,thisK,thisLlow))
     &           /(dtheta*rhf(thisJ))
         fatop   =-(phi(thisJ,thisK,thisLhi)-phi(thisJ,thisK,thisL))
     &           /(dtheta*rhf(thisJ))
 
         accr=frbottom+(frtop-frbottom)*(r_p-r(thisJ))/rof3n

         accz=fzbottom+(fztop-fzbottom)*(z_p(I)-z(thisK))/zof3n

         acca=fabottom+(fatop-fabottom)
     &      *(angle_p-dble(thisL-1)*dtheta)
     &      /(dtheta)

         accx=accr*cos(angle_p)-acca*sin(angle_p)
         accy=accr*sin(angle_p)+acca*cos(angle_p)

         accr=-mass_star/sqrt((x_p(I)-xstar)**2+(y_p(I)-ystar)**2
     &        + z_p(I)**2+gsoft**2)**3
         accx=accr*(x_p(I)-xstar)+accx
         accy=accr*(y_p(I)-ystar)+accy
         accz=accr*z_p(I)+accz

         vx_p(I)=vx_p(I)+accx*dt
         vy_p(I)=vy_p(I)+accy*dt
         vz_p(I)=vz_p(I)+accz*dt

      enddo 
!$OMP ENDDO
!$OMP DO SCHEDULE(STATIC)
      do L=1,LMAX
      do K=1,KMAX2
      do J=1,JMAX2
        phi(J,K,L)=phi(J,K,L)+starphi(J,K,L)
      enddo
      enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL
      return
      end subroutine update_particle_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine update_particle_pos(dt)

      integer :: I
      real*8::dt

!$OMP PARALLEL DEFAULT(SHARED) 
!$OMP DO SCHEDULE(DYNAMIC)
      do I=1,NPARTICLE

       if(particle_skip(I))cycle
 
       x_p(I)=x_p(I)+vx_p(I)*dt
       y_p(I)=y_p(I)+vy_p(I)*dt
       z_p(I)=z_p(I)+vz_p(I)*dt

       !----------------------
       ! boundary condtion
       !----------------------  
       if(z_p(I)<zero)then
          z_p(I)=-z_p(I)
          vz_p(I)=-vz_p(I)
       endif 

      enddo
!$OMP ENDDO
!$OMP END PARALLEL
      return
      end subroutine update_particle_pos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine clean_particles

      integer::I
      real*8::r_p,angle_p

      do I=1,NPARTICLE

        call get_cylinder(x_p(I),y_p(I),r_p,angle_p)

        if(z_p(I)<zero)then
          z_p(I)=-z_p(I)
          vz_p(I)=-vz_p(I)
        endif 
        if(r_p>r(JMAX1))particle_skip(I)=.true.
        if(z_p(I)>z(KMAX1))particle_skip(I)=.true.

      enddo 
      return
      end subroutine clean_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine set_particle_density()

      integer::thisj,thisk,thisL,I
      real*8::r_p,angle_p
      rhotot=rho
  
      do i=1,NPARTICLE

         if(particle_skip(I))cycle

         call get_cylinder(x_p(I),y_p(I),r_p,angle_p)

         thisj=int(r_p/rof3n)+2
         thisk=int(z_p(i)/zof3n)+2
         thisl=int(angle_p/dtheta)+1 

!         print *, thisj,thisk,thisl,I
!         print *,r_p(I),z_p(I),angle_p(I),mass_p(I) 

         rhotot(thisj,thisk,thisl)=rhotot(thisj,thisk,thisl)
     &     + mass_p(I)/(rof3n*zof3n*rhf(thisj)*dtheta)

      enddo 
      return
      end subroutine set_particle_density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine particle_timestep(tmin)

      integer::I
      real*8::tmin,r_p,angle_p,omega_p,vr_p

      tmin=1d6
      do I=1,NPARTICLE

        if(particle_skip(I))cycle
         
        call get_cylinder(x_p(I),y_p(I),r_p,angle_p)
        call get_cyl_vel(vx_p(I),vy_p(I),r_p,angle_p,vr_p,omega_p)

        if(omega_p/=zero)then
           tmin=min(tmin,one/abs(omega_p)*.001)
           tmin=min(tmin,dtheta/abs(omega_p)*.1)
        endif
        if(vz_p(I)/=zero)then
           tmin=min(tmin,zof3n/abs(vz_p(I))*.1)
        endif
        if(vr_p/=zero)then
           tmin=min(tmin,rof3n/abs(vr_p)*.1)
        endif

      enddo
      return
      end subroutine particle_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine dump_particles(INDX)

      integer::INDX,I
      character::cindex*6,outfile*15
      real*8:: r_p,an,vr_p,o_p

      write(cindex,'(i6.6)')INDX
      outfile="particle."//cindex
      open(unit=176,FILE=outfile)
 
      do I=1,NPARTICLE
        call get_cylinder(x_p(I),y_p(I),r_p,an)
        call get_cyl_vel(vx_p(I),vy_p(I),r_p,an,vr_p,o_p)
        write(176,'(I6,7(1X,1pe15.8),1X,L)')I,x_p(I),y_p(I),z_p(I),
     &     vx_p(I),vy_p(I),vz_p(I),mass_p(I),particle_skip(I)
      enddo
      close(176)
      return
      end subroutine dump_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine read_particles(INDX)

      integer::INDX,I,II,file_eof
      logical::check_particle
      character::cindex*6,infile*15,dum*119

      write(cindex,'(i6.6)')INDX
      infile="particle."//cindex
      open(unit=176,FILE=infile)
      II=0
      do while (.true.)
        read(176,'(A119,L)',iostat=file_eof)dum,check_particle
        if(file_eof<0)exit
        if(.not.check_particle)II=II+1
      enddo
      rewind(176)
      print *, "Found ",II," particles in file ",infile
      NPARTICLE=II
      allocate(x_p    (NPARTICLE))
      allocate(y_p    (NPARTICLE))
      allocate(z_p    (NPARTICLE))
      allocate(mass_p (NPARTICLE))
      allocate(vx_p   (NPARTICLE))
      allocate(vy_p   (NPARTICLE))
      allocate(vz_p   (NPARTICLE))
!      allocate(Ugas   (NPARTICLE))
!      allocate(Wgas   (NPARTICLE))
!      allocate(Ogas   (NPARTICLE))
!      allocate(new_dv_a(NPARTICLE))
!      allocate(new_dv_r(NPARTICLE))
!      allocate(new_dv_z(NPARTICLE))

      allocate(particle_skip(NPARTICLE))
       
      do I=1,NPARTICLE
       read(176,'(I6,7(1X,1pe15.8),1X,L)')II,x_p(I),y_p(I),z_p(I),
     &    vx_p(I),vy_p(I),vz_p(I),mass_p(I),particle_skip(I)
      enddo
      close(176)
      return
      end subroutine read_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine clean_stop_particles()

      deallocate(x_p,y_p,z_p,mass_p,vx_p,vy_p,vz_p)
!      deallocate(Ugas,Ogas,Wgas,new_dv_a,new_dv_r,new_dv_z)
      deallocate(particle_skip)

      end subroutine clean_stop_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8 function updateDvWithDrag
     &       (dv,kn,csound,rhogas,rhoa,asize,dt)
      IMPLICIT none
      real*8::dv,kn,csound,rhogas,rhoa,asize,mach,kd,Re
      real*8::beta2,beta,alpha,dt,BB,dv_eps,dv_sto,magDv,PQ,PQINV


            if(dv==zero)then
              updateDvWithDrag=zero
              return
            endif
            magDv=abs(dv)
            mach=magDv/csound
            Re=three*sqrt(pi/eight)*mach/kn
            beta2=128d0*csound**2/(nine*pi)
            beta=sqrt(beta2)
            alpha=rhogas*csound/(rhoa*asize)*sqrt(eight/pi)
            BB=beta+sqrt(dv**2+beta2)
            dv_eps=two*BB*dv*beta*exp(-alpha*dt)
     &            /(BB*BB-dv**2*exp(-two*alpha*dt))

!!!!!!!!! NOW WORK ON THE STOKES LIMIT

            if(Re<=500d0)then
              PQ=0.687d0
              PQINV=one/PQ
              beta=0.15d0*(three*sqrt(pi/eight)/(csound*Kn))**PQ
              dv_sto=exp(-three*alpha*kn*dt)
              if (dv_sto>zero)then
                dv_sto=dv_sto/( (beta*magDv**PQ+one)*magDv**(-PQ)
     &                -beta*exp(-three*PQ*alpha*kn*dt))**PQINV
              endif
              dv_sto=dv_sto*magDv/dv

            elseif(Re<=1500d0)then
              PQ=2.4d0
              PQINV=one/PQ
              beta=3.96d-6*(three*sqrt(pi/eight)/(csound*Kn))**PQ
              dv_sto=(magDv**(-PQ)+7.2d0*Kn*alpha*beta*dt)
     &              **(-PQINV)*dv/magDv
            else
               dv_sto=dv/(one+magDv*0.99d0*sqrt(pi/eight)
     &               *alpha/csound*dt)
            endif

            updateDvWithDrag=dv_eps*(three*kn)**2/((three*kn)**2+one)
     &            +one/((three*kn)**2+one)*dv_sto

      return

      end function updateDvWithDrag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine account_for_drag(dt)
     
      integer:: I,thisj,thisk,thisl,thisLlow,thisLhi

      real*8::kn,oldTotmom,oldGasmom,Gasmom,ftop,fbottom
      real*8::new_vel,limiter
      real*8:: dv_r,dv_a,dv_z,dustRho,csound1,dt
      real*8::jn_p,omega_p,vr_p,angle_p,r_p
      real*8::ogas,wgas,ugas,new_dv_a,new_dv_r,new_dv_z

      logical::trigger=.false.
  
      limiter=den*phylim
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(kn,oldTotmom,oldGasmom,gasmom,
!$OMP& ftop,fbottom,new_vel,dv_r,dv_a,dv_z,dustRho,csound1,jn_p,
!$OMP& omega_p,vr_p,angle_p,r_p,thisj,thisk,thisl,thisLlow,thislHi,
!$OMP& ugas,wgas,ogas,new_dv_r,new_dv_a,new_dv_z)
!$OMP DO SCHEDULE(DYNAMIC)
      do I=1,NPARTICLE

       if(particle_skip(I))cycle

       call get_cylinder(x_p(I),y_p(I),r_p,angle_p)
       call get_cyl_vel(vx_p(I),vy_p(I),r_p,angle_p,vr_p,omega_p)

       thisj=int(r_p/rof3n)+2
       thisk=int(z_p(i)/zof3n)+2
       thisl=int(angle_p/dtheta)+1 
       thisLlow=thisL-1
       thisLhi =thisL+1
       if(thisLlow<1  )thisLlow=LMAX
       if(thisLhi>LMAX)thisLhi =1

       if(rho(thisj,thisk,thisl)<limiter)cycle

      !-------------------------
      ! get velocity differences
      !-------------------------

       Ugas= (s(thisj,thisk,thisl)*(rof3n-r_p+r(thisj))/rof3n +
     &        s(thisj+1,thisk,thisl)*(r_p-r(thisj))/rof3n) /
     &        rho(thisj,thisk,thisl)

       Wgas= (t(thisj,thisk,thisl)*(zof3n-z_p(I)+z(thisk))/zof3n +
     &        t(thisj,thisk+1,thisl)*(z_p(I)-z(thisk))/zof3n) /
     &        rho(thisj,thisk,thisl)

        Ogas=omega(thisj,thisk,thisl)+half*( 
     &    omega(thisj+1,thisk,thisl)-omega(thisj-1,thisk,thisl) ) /
     &    rof3n*(r_p-rhf(thisj))

       dv_r = Ugas - vr_p
       dv_z = Wgas - vz_p(I)
       if(r_p/=zero)then
         dv_a = (Ogas-omega_p)*r_p
       else
         dv_a = zero
         omega_p=zero
       endif

      !-----------------------------
      ! find new velocity difference
      !-----------------------------

       kn=half*muc*1.67d-24 / 
     &    (rho(thisj,thisk,thisl)*rhoconv*pi*1d-16*psize*AUcgs)
       csound1=sqrt(gamma1(thisj,thisk,thisl)
     &        *bkmpcode/muc*tempk(thisj,thisk,thisl))

       new_dv_r=updateDvWithDrag(dv_r,kn,csound1,
     &             rho(thisj,thisk,thisl),rhoacgs/rhoconv,
     &             psize,dt)

       new_dv_z=updateDvWithDrag(dv_z,kn,csound1,
     &             rho(thisj,thisk,thisl),rhoacgs/rhoconv,
     &             psize,dt)

       new_dv_a=updateDvWithDrag(dv_a,kn,csound1,
     &             rho(thisj,thisk,thisl),rhoacgs/rhoconv,
     &             psize,dt)


      !----------------------------------------
      ! update velocities and conserve momentum
      !----------------------------------------

       dustRho=mass_p(I)/(rof3n*zof3n*dtheta*rhf(thisj))
       
      !-----------------
      ! radial direction
      !-----------------     

       fbottom=(rof3n-r_p+r(thisj))/rof3n
       ftop   =(r_p-r(thisj))/rof3n
  
       oldGasmom=Ugas*rho(thisj,thisk,thisl)

       oldTotmom=oldGasmom+dustRho*vr_p

       new_vel = (oldTotmom-rho(thisj,thisk,thisl)*new_dv_r)/ 
     &           (rho(thisj,thisk,thisl) + dustRho)

       vr_p=new_vel

       Gasmom = oldTotmom-dustRho*vr_p

       s(thisj,thisk,thisl)=(Gasmom-oldgasmom)*fbottom
     &                     +s(thisj,thisk,thisl)

       s(thisj+1,thisk,thisl)=(Gasmom-oldGasmom)*ftop
     &                       +s(thisj+1,thisk,thisl)

      !------------------
      ! vertical direction
      !------------------     

       fbottom=(zof3n-z_p(I)+z(thisk))/zof3n
       ftop   =(z_p(I)-z(thisk))/zof3n

       oldGasmom=Wgas*rho(thisj,thisk,thisl)

       oldTotmom=oldGasmom+dustRho*vz_p(I)

       new_vel = (oldTotmom-rho(thisj,thisk,thisl)*new_dv_z)/ 
     &           (rho(thisj,thisk,thisl) + dustRho)
 
       !print *, dv_z,new_dv_z,new_vel,vz_p(I)
       vz_p(I)=new_vel

       Gasmom = oldTotmom-dustRho*vz_p(I)

       if(thisk==2)then
         t(thisj,thisk+1,thisl)=t(thisj,thisk+1,thisl)+
     &        (Gasmom-oldGasmom)*ftop
       else
         t(thisj,thisk,thisl)=t(thisj,thisk,thisl)+
     &        (Gasmom-oldGasmom)*fbottom
         t(thisj,thisk+1,thisl)=t(thisj,thisk+1,thisl)+
     &        (Gasmom-oldGasmom)*ftop
       endif

      !--------------------
      ! azimuthal direction
      !--------------------     

       oldGasmom=Ogas*rho(thisj,thisk,thisl)*r_p

       oldTotmom=oldGasmom+dustRho*omega_p*r_p

       new_vel = (oldTotmom-rho(thisj,thisk,thisl)*new_dv_a)/ 
     &           (rho(thisj,thisk,thisl) + dustRho)

       jn_p=new_vel*r_p
       omega_p=new_vel/r_p

       Gasmom = oldtotmom*r_p-dustRho*jn_p

       a(thisj,thisk,thisl)=a(thisj,thisk,thisl) + 
     &                      (Gasmom-oldGasmom*r_p)
       jn(thisj,thisk,thisl)=a(thisj,thisk,thisl)/rho(thisj,thisk,thisl)
       omega(thisj,thisk,thisl)=jn(thisj,thisk,thisl)/rhf(thisj)**2

       call get_cart_vel(vr_p,omega_p,vx_p(I),vy_p(I),r_p,angle_p)

      enddo
!$OMP ENDDO
!$OMP END PARALLEL

      return
      end subroutine account_for_drag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine particle_fullstep(odelt,delt)

      real*8 odelt,delt

      call update_particle_vel(odelt)
      call account_for_drag(odelt)
      call update_particle_vel(delt)
      call account_for_drag(delt)
      call update_particle_pos(two*delt)
      call clean_particles

      return
      end subroutine particle_fullstep



      end module

