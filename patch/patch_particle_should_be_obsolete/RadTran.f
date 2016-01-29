      SUBROUTINE RADTRANS
      IMPLICIT real*8 (a-h,o-z)      
#include "hydroparam.h"
#include "globals.h"
#include "units.h"      
      
      integer jstart,kfit, kfitmax
      REAL*8 sum,Tfit,Teff4,Tbd4, Tenv4, limiter
      REAL*8 Oross(jmax2,kmax2,lmax),Oplck(jmax2,kmax2,lmax)
     &     ,Oabs(jmax2,kmax2,lmax),Otot(jmax2,kmax2,lmax),absfrac
      REAL*8 Oplck_env(jmax2,kmax2,lmax), Tfit_old
     &     ,Oplck_ph(jmax2,kmax2,lmax), Oirr, tau_face(kmax2,4)
      REAL*8 Tbdry4, tau_f(jmax2,3),tbgnrd4,fakk(jmax2,kmax2,lmax)
      real*8 dtau,dtau1,fkk
      COMMON /COOLINGSHARED/Oross,Oplck,Oabs,Otot,sum, Oplck_env, 
     &      Oplck_ph
      
      DR=ROF3N
      DZ=ZOF3N
      
      limiter = den*phylim 

      if(jmin.gt.2) then
         jstart=jmin
      else   
         jstart=2
      endif 

C... New cooling (acm 2002)   

C    To convert from polytropic units of temperature to Kelvin
C    Use T(K)=T(pu)*Tconv
C    Opacities in cgs, so always use Oross or Oabs or Oplck times Sconv.

C    In this next loop, read opacities from table using the pollack 
C    function or the dalessio subroutine.
      tbgrnd4=(tbgrnd/tconv)**4
      Tenv4 = (TenvK/Tconv)**4
c... Added by Kai: Envelope heating term

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Teff4,absfrac,dtau1,dtau,k,j,l,  &
!$OMP&Oirr,fkk)                                                         &
!$OMP&  SHARED(Dconv,Tconv,sigma,Sconv,dz,Pconv,                        &
!$OMP& limiter,jstart)
!$OMP DO SCHEDULE(STATIC)
      do l=1,lmax
         do j=jstart,jmax
            TeffK(j,l)=zero
            TphK(j,l)=zero
            surfcgs(j,l)=zero
         enddo

         do k=1,kmax1
           do J = JSTART, JMAX 
               tau(j,k,l,:)=zero
               l_tau_z(J,K,L) = zero
               Lambda(j,k,l)=zero
               Hgamma(j,k,l)=zero
               Oross(j,k,l)=zero
               Oplck(j,k,l)=zero
               Oplck_env(j,k,l)=zero
               Oplck_ph(j,k,l)=zero
               Oabs(j,k,l)=zero
               Otot(j,k,l)=zero
               Radflux(j,k,l,1)=zero
               Radflux(j,k,l,2)=zero
               Radflux(j,k,l,3)=zero
               Divflux(j,k,l)=zero

               if (rho(j,k,l).ge.limiter) then
                  
!                  call dalessio(TenvK, (P(j,k,l)*Pconv), Oross(j,k,l),
!     &                 Oplck_env(j,k,l), Oirr, Otot(j,k,l))
        call dalessio_simp(TempK(j,k,l),(P(j,k,l)*Pconv),Oross(j,k,l))
        !ffk=fk*rho_p(J,K,L)/(rho(J,K,L)*dust_to_gas)
        ffk=fk*(dust_to_gas/0.02d0*0.2d0
     &     +(rhotot(J,K,L)-rho(j,k,l))/rho(j,k,l)/0.02d0*0.8d0)
!          if(K==2.and.L==1) print *, J,K,L,fakk(j,k,l),
!     &       (rhotot(J,K,L)-rho(J,K,L))/rho(J,K,L)
!                 Oross(j,k,l) = ffk*pollack(tempk(J,K,L))
                 Oross(j,k,l) = ffk*oross(J,K,L)
               endif
            enddo
          enddo

C     The following loop calculates vertical optical depths and the surface
c     density.  All taus are cell centered.
! the commented lines could be used to control maximum dtau change if such
! an occasion becomes necessary

          dtau=zero;dtau1=zero
          do k=kmax,2,-1
            do J = JSTART,JMAX
               if (tau(j,k+1,l,1).eq.zero) then
                  if (rho(j,k,l).ge.limiter) then
                 dtau = max(oross(J,K,L)*rho(J,K,L)*sconv*dz*half,1d-15)
                     tau(j,k,l,1)=dtau

                     surfcgs(j,l)=rho(j,k,l)*Sconv*dz*half
                  endif
               else      

         dtau = oross(J,K+1,L)*rho(J,K+1,L)*sconv*dz*half
         dtau1= oross(J,K,L)*rho(J,K,L)*sconv*dz*half
                  tau(j,k,l,1)=tau(J,K+1,L,1) + dtau + dtau1

                  surfcgs(j,l)=surfcgs(j,l) + ((rho(j,k+1,l)+rho(j,k,l))
     &                 *Sconv*dz*half)
               endif 
            enddo
         enddo

! I have taken out the following code because it is obsolete.  I am
! leaving it commented because it might need to be resurrected in
! a much more efficient format.

      enddo    
!$OMP ENDDO ! END BIG L LOOP

      if (ICTYPE.EQ.1) call ConstantTcool
      if (ICTYPE.EQ.5) call TcoolOmega

!$OMP END PARALLEL
      
!         CALL CPU_TIME(time_in)
      if ( ictype == 69 ) call Hybrid()
!         CALL CPU_TIME(time_out)
!         print *, "HYBRID ",time_out-time_in
      if ( ictype == 100 ) call simpleCool

      RETURN
      END

      FUNCTION POLLACK(TK)
      IMPLICIT real*8 (a-h,o-z)
      real*8 TK

C     Pollack et al. (1994) rosseland opacities

      if(TK.lt.80.0) then
         pollack=(TK**2)/3200.0
      else if(TK.lt.170.0) then
         pollack=-2.0 + 0.050*TK
      else if(TK.lt.180.0) then
         pollack=62.60 - 0.330*TK
      else if(TK.lt.270.0) then
         pollack=-1.0 + 0.0233*TK
      else if(TK.lt.300.0) then
         pollack=8.0 - 0.010*TK
      else if(TK.lt.425.0) then
         pollack=1.88 + 0.0104*TK
      else if(TK.lt.440.0) then
         pollack=128.13 - 0.2867*TK
      else if(TK.lt.670.0) then
         pollack=0.57 + 0.0033*TK
      else if(TK.lt.700.0) then
         pollack=19.50 - 0.0250*TK
      else if(TK.lt.1300.0) then
         pollack=-0.33 + 0.0033*TK
      else if(TK.lt.1350.0) then
         pollack=24.80 - 0.0160*TK
      else if(TK.lt.1450.0) then
         pollack=46.40 - 0.0320*TK
      else
         pollack=0.01
      endif

      return

      END


      FUNCTION XMMW(TK,Pcgs)

      implicit real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"

C     By courtesy of Paola D'Alessio (2002)  
c     limits in T&P: 10**0.5 to 1.e7 K, 1.e-12 to 1.e9 din/cm^2
c     tables in log T, log P, 
      real*8  o1(2),o2(2)
      real*8 TK, Pcgs, xmmw
      
      xlt=dlog10(TK)
      xlp=dlog10(Pcgs)
     
c     Out of bounds
      if (xlp.lt.Ptab(1)) xlp=Ptab(1)
      if (xlp.gt.Ptab(itable)) xlp=Ptab(itable)
      if (xlt.lt.Ttab(1)) xlt=Ttab(1)
      if (xlt.gt.Ttab(itable)) xlt=Ttab(itable)
 
      it1=dint((xlt-Ttab(1))/(Ttab(itable)-Ttab(1))*
     &     dfloat(itable-1)+1.)
      it2=it1+1
      if (it2.gt.itable) it2=it1
      
      ip1=dint((xlp-Ptab(1))/(Ptab(itable)-Ptab(1))*
     &     dfloat(itable-1)+1.)
      ip2=ip1+1
      if (ip2.gt.itable) ip2=ip1
      
      l=1

      do it=it1,it2
         i=1
         do ip=ip1,ip2
            o1(i)=XMMWtable(ip,it)
            i=i+1
C           P
         enddo 
         
         if (ip2.ne.ip1) then
            o2(l)=o1(1)+(o1(2)-o1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
         else
            o2(l)=o1(1)
         endif
         
         l=l+1
c	 T
      enddo 

      if (it2.ne.it1) then
         xmmw=o2(1)+(o2(2)-o2(1))*(xlt-Ttab(it1))/(Ttab(it2)-Ttab(it1))
      else
         xmmw=o2(1)
      endif

      RETURN
      END
      
      SUBROUTINE DALESSIO_simp(TK,Pcgs,Ors)

      implicit real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"

C     By courtesy of Paola D'Alessio (2002)  
c     limits in T&P: 10**0.5 to 1.e7 K, 1.e-12 to 1.e9 din/cm^2
c     tables in log T, log P, 
      real*8 o1(2),o2(2),p1(2),p2(2),g1(2),g2(2),s1(2),s2(2)
      real*8 TK,Pcgs,Ors,Opl,Oab,Oas

      Ors=0.0
            
      xlt=log10(TK)
      xlp=log10(Pcgs)

c     Out of bounds
      
      if (xlt.gt.Ttab(itable)) xlt=Ttab(itable)
      if (xlp.lt.Ptab(1))xlp=Ptab(1)
      if (xlp.gt.Ptab(itable))xlp=Ptab(itable)
      if (xlt.lt.Ttab(1)) xlt=Ttab(1)
      
      it1=int((xlt-Ttab(1))/(Ttab(itable)-Ttab(1))*
     &     dble(itable-1)+1.)
      it2=it1+1
      
      if (it2.gt.itable) it2=it1
      
      ip1=int((xlp-Ptab(1))/(Ptab(itable)-Ptab(1))*
     &     dble(itable-1)+1.)
      ip2=ip1+1

      if (ip2.gt.itable) ip2=ip1
      

      l=1
      do it=it1,it2
         j=1
         do ip=ip1,ip2
            o1(j)=ROStable(ip,it)
            j=j+1
c	    P
         enddo 
         
         if (ip2.ne.ip1) then
            o2(l)=o1(1)+(o1(2)-o1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
         else
            o2(l)=o1(1)
         endif
         
         l=l+1
c	 T
      enddo 
      
      if (it2.ne.it1) then
         Ors=o2(1)+(o2(2)-o2(1))*(xlt-Ttab(it1))/
     &        (Ttab(it2)-Ttab(it1))
      else
         Ors=o2(1)
      endif
      
      Ors=10.**Ors

      if (xlt.ge.6.) then
         Ors=0.3429
      endif

      RETURN
      END
 
      
      SUBROUTINE DALESSIO(TK,Pcgs,Ors,Opl,Oab,Oas)

      implicit real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"

C     By courtesy of Paola D'Alessio (2002)  
c     limits in T&P: 10**0.5 to 1.e7 K, 1.e-12 to 1.e9 din/cm^2
c     tables in log T, log P, 
      real*8 o1(2),o2(2),p1(2),p2(2),g1(2),g2(2),s1(2),s2(2)
      real*8 TK,Pcgs,Ors,Opl,Oab,Oas

      Ors=0.0
      Opl=0.0
      Oab=0.0
      Oas=0.0
            
      xlt=dlog10(TK)
      xlp=dlog10(Pcgs)

c     Out of bounds
      
      if (xlt.gt.Ttab(itable)) xlt=Ttab(itable)
      if (xlp.lt.Ptab(1))xlp=Ptab(1)
      if (xlp.gt.Ptab(itable))xlp=Ptab(itable)
      if (xlt.lt.Ttab(1)) xlt=Ttab(1)
      
      it1=dint((xlt-Ttab(1))/(Ttab(itable)-Ttab(1))*
     &     dfloat(itable-1)+1.)
      it2=it1+1
      
      if (it2.gt.itable) it2=it1
      
      ip1=dint((xlp-Ptab(1))/(Ptab(itable)-Ptab(1))*
     &     dfloat(itable-1)+1.)
      ip2=ip1+1

      if (ip2.gt.itable) ip2=ip1
      

      l=1
      do it=it1,it2
         j=1
         do ip=ip1,ip2
            o1(j)=ROStable(ip,it)
            p1(j)=PLKtable(ip,it)
            g1(j)=IRRtable(ip,it)
            s1(j)=SCAtable(ip,it)
            j=j+1
c	    P
         enddo 
         
         if (ip2.ne.ip1) then
            o2(l)=o1(1)+(o1(2)-o1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
            p2(l)=p1(1)+(p1(2)-p1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
            g2(l)=g1(1)+(g1(2)-g1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
            s2(l)=s1(1)+(s1(2)-s1(1))*(xlp-Ptab(ip1))/(Ptab(ip2)
     &           -Ptab(ip1))
         else
            o2(l)=o1(1)
            p2(l)=p1(1)
            g2(l)=g1(1)
            s2(l)=s1(1)
         endif
         
         l=l+1
c	 T
      enddo 
      
      if (it2.ne.it1) then
         Ors=o2(1)+(o2(2)-o2(1))*(xlt-Ttab(it1))/
     &        (Ttab(it2)-Ttab(it1))
         Opl=p2(1)+(p2(2)-p2(1))*(xlt-Ttab(it1))/
     &        (Ttab(it2)-Ttab(it1))
         Oab=g2(1)+(g2(2)-g2(1))*(xlt-Ttab(it1))/
     &        (Ttab(it2)-Ttab(it1))
         Oas=s2(1)+(s2(2)-s2(1))*(xlt-Ttab(it1))/
     &        (Ttab(it2)-Ttab(it1))
      else
         Ors=o2(1)
         Opl=p2(1)
         Oab=g2(1)
         Oas=s2(1)
      endif
      
      Ors=10.**Ors
      Opl=10.**Opl
      Oab=10.**Oab
      Oas=10.**Oas

      if (xlt.ge.6.) then
         Ors=0.3429
      endif

      RETURN
      END
      
      FUNCTION FLUXLMDF(T0,T1,Ors,Ors1,Rho0,Rho1,jpos,dim)
      
      IMPLICIT real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"
#include "units.h"

      real*8 T0,T1,Ors,Ors1,Rho0,Rho1,y,beta,ddim,dtau
      integer jpos,dim

      if (dim.eq.1) then
         ddim=rof3n
      else if (dim.eq.2) then
         ddim=zof3n
      else if (dim.eq.3) then
         ddim=rhf(jpos)*dtheta
      endif

      dtau = (ors+ors1)*sconv*(rho0+rho1)*ddim*quarter
! if optical depth limiting is required, uncomment
!      if (dtau > 10.) dtau = 10.
      dtau = dtau*four

C     this y-beta thing is the flux limiter
      
      y=32d0*abs(T0-T1)/
     &     ( dtau*(T0+T1) )

      beta=(two + y) / (six + three*y + y**2)
c      beta= 1.0/3.0  
c ... Turn the flux limiter off

      Fluxlmdf=-( eight*sigma*beta*(T0+T1)**3*(T0-T1) ) /
     &     ( dtau*(Tconv**4) )      
      
      RETURN
      END
c
      SUBROUTINE HEIGHTS ! THIS IS KEPT HERE FOR BONUS POINTS.
          ! IF YOU FIND IT, YOU GET THE ARCHAEOLOGY AWARD.
      
      IMPLICIT real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"
#include "units.h"
      real*8 limiter
      limiter = den*phylim 
      do j=2,jmax
         do k=2,kmax
            if ((rho(j,k,1).ge.limiter).and.
     &           (rho(j,k+1,1).lt.limiter))print*,j,'klim ',k
            if ((tau(j,k,1,1).ge.(twothree)).and.
     &           (tau(j,k+1,1,1).lt.(twothree)))print*,j,' kross ',k
            if ((tau(j,k,1,2).ge.(twothree)).and.
     &           (tau(j,k+1,1,2).lt.(twothree)))print*,j,' kplck ',k
            if ((tau(j,k,1,3).ge.(twothree)).and.
     &           (tau(j,k+1,1,3).lt.(twothree)))print*,j,'kirr ',k
         enddo
      enddo
C      do k=2,kmax
C         do j=2,jmax
C            if ((rho(j,k,1).ge.limiter).and.
C     &           (rho(j-1,k,1).lt.limiter))print*,j,'jlim ',k 
C         enddo
C      enddo

      RETURN
      END      

      subroutine TempFind() ! ACB find the temperature based on the energy/particle 
                            ! temperature relation. The temperature is interpolated
                            ! from the temperature table temptable, which is 
                            ! calculated in setup based on X, Y, Z composition
                            ! and other assumptions.
         implicit none
#include "hydroparam.h"
#include "globals.h"
#include "units.h"
        
         integer J,K,L,I 
         real*8 limiter
         real*8 eng,dummy

         limiter = den*phylim 

!$OMP DO SCHEDULE(STATIC) REDUCTION(+:eflufftot)
         do L = 1, LMAX
           do K = 1, KMAX2
             do J = 1, JMAX2

               eng = zero
            
               if (rho(J,K,L) >= limiter) then

                 eng = eps(J,K,L)/rho(J,K,L)*engconv
                 
                 if (eng < engtable(1)) then

                   tempk(J,K,L) = temptable(1)/engtable(1)*eng
                   eps(J,K,L) = engtable(1)*rho(J,K,L)/engconv

                   I = 2

                   gamma1(J,K,L) = gammatable(I-1) +
     & (gammatable(I)-gammatable(I-1))/(temptable(I)-temptable(I-1))
     & * (tempk(J,K,L)-temptable(I-1))

                 else

                   I = 2
                   do while ( engtable(I) < eng .and. I < TTABLE )
                    I = I+1
                   enddo

                   tempk(J,K,L) = temptable(I-1) +
     & (temptable(I)-temptable(I-1))/(engtable(I)-engtable(I-1))*
     & (eng - engtable(I-1))

                   gamma1(J,K,L) = gammatable(I-1) +
     & (gammatable(I)-gammatable(I-1))/(temptable(I)-temptable(I-1))
     & * (tempk(J,K,L)-temptable(I-1))


                 endif

               else

                 I = 2
                 do while ( temptable(I) < tbgrnd .and. I < TTABLE )
                  I = I+1
                 enddo

                 tempk(J,K,L) = tbgrnd

                 dummy   = engtable(I-1) + (engtable(I)-engtable(I-1))
     &                   / (temptable(I)-temptable(I-1))
     &                   * (tbgrnd-temptable(I-1))

                 dummy   = dummy*rho(J,K,L)*rhoconv/pconv

                 eflufftot = eflufftot + (eps(J,K,L)-dummy)*rhf(J)
     &                     * rof3n*dtheta*zof3n
                 eps(J,K,L) = dummy

                 gamma1(J,K,L) = gammatable(I-1) +
     & (gammatable(I)-gammatable(I-1))/(temptable(I)-temptable(I-1))
     & * (tempk(J,K,L)-temptable(I-1))


               endif

             enddo
           enddo
         enddo
!$OMP END DO

      return

      end subroutine TempFind

      subroutine TempFindSpec(eps,rho,temp)! like TempFind, but for a single cell. 

         implicit none

#include "hydroparam.h"

         integer I
         real*8 eng,eps,rho,temp,kk
         real*8 Msyscgs,PKcgs,Tconv,Sconv,Dconv,Pconv,sigma,rhoconv,
     &       engconv,bkmpcode
          COMMON /CONVERT/
     &     Msyscgs,PKcgs,Tconv,Sconv,
     &     Dconv,Pconv,sigma,rhoconv,engconv,bkmpcode
         real*8 :: temptable,engtable,gammatable,muc
         common /engtables/temptable(TTABLE),engtable(TTABLE),
     &           gammatable(TTABLE),muc

         eng = (eps)/rho*engconv

         if (eng < engtable(1)) then

            temp = temptable(1)/engtable(1)*eng
            eps = engtable(1)*rho/engconv+kk*rho*rho

         else

            I = 2
            do while ( engtable(I) < eng .and. I < TTABLE )
              I = I+1
            enddo

            temp = temptable(I-1) +
     & (temptable(I)-temptable(I-1))/(engtable(I)-engtable(I-1))*
     & (eng - engtable(I-1))

         endif

      return

      end subroutine TempFindSpec
