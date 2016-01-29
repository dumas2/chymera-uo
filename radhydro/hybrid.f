      subroutine Hybrid()
      use eos, only: get_gamma,get_gamma_from_tk
      implicit none
#include "hydroparam.h"
#include "globals.h"
#include "units.h"

      real(KIND=8) :: dumb ,oldcrate,eps1,tirr,dumbexp
      real(KIND=8) :: limiter,aconst,bconst,cconst, mu,mmw,gam
      real(KIND=8) :: FluxLmDf,sum,crate,ctime,maxt,taufld
      
      integer :: J,K,L,I,LM,LP,ITER,JL,KL,LL
      integer, parameter::maxiter=8

      REAL*8 Oross(jmax2,kmax2,lmax),Oplck(jmax2,kmax2,lmax)
     &     ,Oabs(jmax2,kmax2,lmax),Otot(jmax2,kmax2,lmax),absfrac,
     &      oplck_env(JMAX2,KMAX2,LMAX),tacc,dlimit
     

      COMMON /COOLINGSHARED/Oross,Oplck,Oabs,Otot,sum,oplck_env

      logical :: use_acc_luminosity=.false.

      mu = one/sqrt(three)
      limiter = den*phylim
      LO=0
      CCOUNTER = 0
      HCOUNTER = 0
      ctime=zero
      iter=0
      taufld=one

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(tacc,dlimit)

      if (mdot < zero .and.use_acc_luminosity) then
              
      tacc=SQRT(SQRT(((Mstar+tmassacc)*(-mdot))   !referenced at 1AU
     &   /((ten+six)*pi*sigma)))
!$OMP MASTER
      print*, 'Lacc (Lsun)', (mstar+tmassacc)*(-mdot)/(rstar*two)
     &*engconv*msuncgs/5d6/3.8d33,
     & tmassacc, -mdot, tacc
!$OMP END MASTER
      else
        tacc=zero
      endif 


!$OMP DO SCHEDULE(STATIC)
      do L =  1, LMAX 
       do K = 1, KMAX2
         do J = 1, JMAX2
           tau(J,K,L,4) = zero
           l_tau_z(J,K,L)=zero
         enddo
       enddo
       do J = 1, JMAX2
          ! THIS IS USERDEFINED USER SPECIFIED
         init_int_in(J,L) = (((tenvk/tconv/(rhf(J)**(.75))+tbgrnd)**4
     & + ((tacc/tconv/(rhf(J)**(.75)))**4)))/ pi*sigma

!         init_int_in(J,L) = (tenvk/tconv)**4
                 
         KFITA(J,L) = 2
       enddo
      enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC) ! set optical depths of positive plane in l_tau
      do L = 1, LMAX 
        do K = KMAX1, 2, -1
          do J = JCOOL, JMAX

            ! interpolate between rosseland and planck means.

            if(rho(J,K,L)>limiter)then
             dtau_z(J,K,L) = rho(J,K,L)*sconv*zof3n *oross(J,K,L)/mu
            else
             dtau_z(J,K,L)=zero
            endif

            l_tau_z(J,K,L) = l_tau_z(J,K+1,L) + dtau_z(J,K,L)
          enddo
        enddo
      enddo
!$OMP ENDDO
!$OMP END PARALLEL

      oldcrate=zero
      ctime=zero
      do while (ctime<delt)
        crate=delt
        maxt=zero
        JL=JMIN;KL=2;ll=2
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(dumb,cconst,bconst,aconst,I,     &
!$OMP& LM,LP,J,K,L,tirr,dlimit) REDUCTION(+:CCOUNTER,HCOUNTER,LO)       &
!$OMP& REDUCTION(min:crate) REDUCTION(max:maxt)                          


!$OMP DO SCHEDULE(STATIC)
      do L =  1, LMAX 
        do K = 1, KMAX2
         do J = 1, JMAX2
           tempk(J,K,L)=tempk(J,K,L)/tconv
           tau(J,K,L,1)=tau(J,K,L,1)/mu
           divflux(J,K,L) = zero
           temporary(J,K,L) = zero
           radflux(J,K,L,:) = zero
           intensity_z(J,K,L)=zero
           intensity_in_z(J,K,L)=zero
           int_temp(J,K,L) = zero
           sfunc(J,K,L) = sigma*tempk(J,K,L)**4/pi
         enddo
       enddo
       do J = JCOOL, JMAX
         sfunc(J,1,L) = sfunc(J,2,L)
       enddo
      enddo
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC)! set optical depths of positive plane in l_tau
      do L = 1, LMAX
        do K = KMAX, 3, -1
          do J = JCOOL, JMAX

            if ( rho(J,K,L) >= limiter 
     &          .and. rho(J,K+1,L) >= limiter 
     &          .and. rho(J,K-1,L) >= limiter ) then
       
              dsdt(J,K,L) =  two*sigma*tempk(J,K,L)**3 * 
     &            (tempk(J,K-1,L)-tempk(J,K+1,L))/(pi*dtau_z(J,K,L) )

            else
              dsdt(J,K,L) = zero
            endif

          enddo
        enddo

        do J = JCOOL, JMAX
          if ( rho(J,2,L) >= limiter 
     &          .and. rho(J,3,L) >= limiter 
     &          .and. rho(J,4,L) >= limiter ) then

            dsdt(J,2,L) = sigma*tempk(J,2,L)**3 * 
     &        ( eight*tempk(J,2,L)-seven*tempk(J,3,L)-tempk(J,4,L)) 
     &        / (three*pi*dtau_z(J,2,L))


          else

            dsdt(J,2,L) = zero

          endif

          dsdt(J,1,L)=-dsdt(J,2,L)

        enddo
      enddo
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX 
        do K = KMAX, 2, -1
          do J = JCOOL, JMAX 

            if ( rho(J,K,L) >= limiter 
     &         .and. rho(J,K+1,L) >= limiter 
     &         .and. rho(J,K-1,L) >= limiter 
     &         .and. dsdt(J,K,L) /= zero) then

              ddsdtt(J,K,L) = ( dsdt(J,K-1,L)-dsdt(J,K+1,L) ) 
     &                      / (two*dtau_z(J,K,L) )


! if the second derivative becomes important, kill it. This is a stability 
! consideration
 
            if(abs(ddsdtt(J,K,L)) >= two*abs(dsdt(J,K,L)))then
               ddsdtt(J,K,L) = zero
            endif

                     
            else
              ddsdtt(J,K,L) = zero
            endif

          enddo
        enddo
  
        do J = JCOOL, JMAX

          if ( rho(J,2,L) >= limiter 
     &          .and. rho(J,3,L) >= limiter 
     &          .and. rho(J,4,L) >= limiter 
     &         .and. dsdt(J,2,L) /= zero) then

            ddsdtt(J,2,L) = half*ddsdtt(J,2,L) 
     &        + (dsdt(J,2,L)-dsdt(J,3,L))/(two*dtau_z(J,2,L))


! if the second derivative becomes important, kill it. This is a stability
! consideration

        if(abs(ddsdtt(J,2,L)) >= two*abs(dsdt(J,2,L)))then
               ddsdtt(J,2,L) = zero
            endif


          else
   
            ddsdtt(J,2,L) = zero

          endif

        enddo
      enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX
        LM = modulo(L-2,LMAX)+1
        LP = modulo(L,LMAX)+1
        do K = 2, KMAX
          do J = JCOOL, JMAX                                           
            if ( tau(J,K,L,1) >= taufld .and. rho(J,K,L)>=limiter ) then
            KFITA(J,L) = K

                radflux(J,K,L,1) = FluxLmDf(tempk(J,K,L)*tconv,
     &                  tempk(J-1,K,L)*tconv,
     &                  dtau_z(J,K,L)*mu/(zof3n*sconv*rho(J,K,L)), 
     &                  dtau_z(J-1,K,L)*mu/(zof3n*sconv*rho(J-1,K,L)),
     &                  rho(J,K,L),rho(J-1,K,L),J,1)

                radflux(J,K,L,3) = FluxLmDf(tempk(J,K,LP)*tconv,
     &                  tempk(J,K,L)*tconv,
     &                  dtau_z(J,K,LP)*mu/(zof3n*sconv*rho(J,K,LP)), 
     &                  dtau_z(J,K,L)*mu/(zof3n*sconv*rho(J,K,L)),
     &                  rho(J,K,LP),rho(J,K,L),J,3)

            endif
          enddo
        enddo
      enddo
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX 
        do K = KMAX1, 2, -1
          do J = JCOOL, JMAX  
   
              I = K+1

              aconst  = l_tau_z(J,I-1,L)**2
     &                - (two*l_tau_z(J,I-1,L)-two)
              aconst  = aconst - (l_tau_z(J,I,L)**2
     &                - (two*l_tau_z(J,I,L)-two)) 
     &                * exp (-dtau_z(J,I-1,L))

              aconst  = aconst * half*ddsdtt(J,I-1,L)

              bconst  =  l_tau_z(J,I-1,L) - one 
     &                - ( l_tau_z(J,I,L) - one ) 
     &                * exp (-dtau_z(J,I-1,L))

              bconst  = bconst * (dsdt(J,I-1,L) 
     &                - ddsdtt(J,I-1,L)*tau(J,I-1,L,1) )

              cconst  = ( sfunc(J,I-1,L) 
     &                - tau(J,I-1,L,1)*dsdt(J,I-1,L) 
     &                + half*tau(J,I-1,L,1)**2*ddsdtt(J,I-1,L) )
     &                * ( one - exp (-dtau_z(J,I-1,L))) 

              dumb    = cconst + bconst + aconst

              if ( dumb < zero ) dumb = sfunc(J,I-1,L) 
     &                  * ( one - exp ( -dtau_z(J,I-1,L) ) )

            intensity_in_z(J,K,L) = dumb + intensity_in_z(J,I,L)
     &                            * exp(-dtau_z(J,K,L))

               
          enddo
        enddo
        do K = KMAX1,2,-1
          do J = JCOOL, JMAX
            intensity_in_z(J,K,L) = intensity_in_z(J,K,L) 
     &              + init_int_in(J,L)*exp(-l_tau_z(J,K,L))
!     &              + init_int_in(J,L)*exp(-tau(J,K,L,4))

          enddo
        enddo
      enddo
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC) ! calculate intensity for each cell
      do L = 1, LMAX
        do J = JCOOL, JMAX 
          intensity_z(J,2,L) = intensity_in_z(J,2,L)
        enddo
        do K = 3 , KMAX1
          do J = JCOOL, JMAX

            I = K-1

              aconst  = l_tau_z(J,I+1,L)**2
     &                + (two*l_tau_z(J,I+1,L)+two)
              aconst  = aconst - (l_tau_z(J,I,L)**2
     &                + (two*l_tau_z(J,I,L)+two)) 
     &                * exp (-dtau_z(J,I,L))
              aconst  = aconst * half*ddsdtt(J,I,L)

              bconst  =  l_tau_z(J,I+1,L) +  one
     &                - ( l_tau_z(J,I,L) + one ) 
     &                * exp (-dtau_z(J,I,L))
              bconst  = bconst * (dsdt(J,I,L) 
     &                - ddsdtt(J,I,L)*tau(J,I,L,1) )

              cconst  = ( sfunc(J,I,L) - tau(J,I,L,1)*dsdt(J,I,L) 
     &                + half*tau(J,I,L,1)**2*ddsdtt(J,I,L) )
     &                * ( one - exp ( -dtau_z(J,I,L)))

              dumb    = aconst + bconst + cconst

              if ( dumb < zero ) dumb = sfunc(J,I,L) 
     &                   * ( one - exp ( -dtau_z(J,I,L) ) )

            intensity_z(J,K,L) = intensity_z(J,I,L)
     &                         * exp(-dtau_z(J,I,L)) + dumb

           enddo
         enddo
         do K = 2, KMAX1
           do J = JCOOL, JMAX
             int_temp(J,K,L) = two*pi*mu*(intensity_z(J,K,L) 
     &                       - intensity_in_z(J,K,L))
           enddo

        enddo
      enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX
        do J = JCOOL, JMAX
          teffk(J,L)= tconv * ( SQRT(SQRT( pi
     & * intensity_z(J,KMAX1,L) / sigma )) )
        enddo
      enddo
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX 
        LM = modulo(L-2,LMAX)+1
        LP = modulo(L,LMAX)+1
        do K = 2, KMAX 
          do J = JCOOL, JMAX
            if ( tau(J,K,L,1) >= taufld  ) then
            if (rho(J,K,L) .ge. limiter) then ! must be optically thick regime. 
                                               ! heat using divF
                       ! set BC for r and theta

                if ( tau(J+1,K,L,1) < taufld ) then

                     radflux(J+1,K,L,1)=abs(int_temp(J,KFITA(J,L)+1,L))
         
                endif

                if ( tau(J-1,K,L,1) < taufld ) then

                     radflux(J,K,L,1)=-abs(int_temp(J,KFITA(J,L)+1,L))

                endif

                if ( tau(J,K,LM,1) < taufld ) then

                     radflux(J,K,LM,3)=-abs(int_temp(J,KFITA(J,L)+1,L))
        
                endif

                if ( tau(J,K,LP,1) < taufld ) then

                     radflux(J,K,L,3)=abs(int_temp(J,KFITA(J,L)+1,L))

                endif

            endif
            endif ! end cooling if structure
            enddo
          enddo
        enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      do L = 1, LMAX 
        LM = modulo(L-2,LMAX)+1
        LP = modulo(L,LMAX)+1
        do K = 2, KMAX 
          do J = JCOOL, JMAX
            if ( tau(J,K,L,1) >= taufld  ) then
            if (rho(J,K,L) .ge. limiter) then ! must be optically thick regime. 
                                               ! heat using divF
 
                temporary(J,K,L) = ((r(J+1)*radflux(J+1,K,L,1))
     &                        - (r(J)*radflux(J,K,L,1)))
     &                        / (rhf(J)*rof3n)        
     &                        + (radflux(J,K,L,3)-radflux(J,K,LM,3))
     &                        / (rhf(J)*dtheta)

            endif
            endif
          enddo
        enddo
      enddo
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC) private(dlimit)  ! calculate divflux
      do L = 1, LMAX 
        LM = modulo(L-2,LMAX)+1
        LP = modulo(L,LMAX)+1
        do K = 2, KMAX
          do J = JCOOL, JMAX
               
           if (rho(J,K,L) > limiter )then
	
             divflux(J,K,L) =                                          
     &         (int_temp(J,K+1,L)-int_temp(J,K,L))/zof3n      
     &       + temporary(J,K,L)
           else 
            divflux(J,K,L)=zero
          endif

          enddo
        enddo
      enddo
!$OMP END DO


C$OMP DO SCHEDULE(STATIC) 
        do L = 1, LMAX
          do K = 2, KMAX
            do J = JCOOL, JMAX

           if (abs(divflux(J,K,L)) > zero )then
	       if (rho(J,K,L) > limiter)then
                dumb=abs(eps(J,K,L)/(divflux(J,K,L)*ten))
                if (dumb<crate)then 
                  JL=J;KL=K;LL=L
                  crate=dumb
                endif
               endif
	      endif    
              maxt=max(tempk(J,K,L),maxt)

           enddo
         enddo
       enddo
!$OMP END DO
!$OMP END PARALLEL 
      if(iter<maxiter-1) then 
       crate=min(crate,(delt-ctime))
      else
       crate=delt-ctime
      endif 
       oldcrate=crate
      
!$OMP PARALLEL REDUCTION(+:totdflux) private(tacc)

      if (mdot < zero .and.use_acc_luminosity) then
              
       tacc=SQRT(SQRT(((Mstar+tmassacc)*(-mdot))   !referenced at 1AU
     &   /((ten+six)*pi*sigma)))
      else
        tacc=zero
      endif 

!$OMP DO SCHEDULE(STATIC) PRIVATE(tirr,mmw,gam)
      do L = 1,LMAX 
       do K = 2,KMAX
         do J = JCOOL,JMAX
            !USER SPECIFIED
	        tirr=((tenvk/rhf(J)**(.75)+tbgrnd)**4
     & +(tacc/rhf(J)**(.75))**4)**.25d0
              if (rho(J,K,L) > limiter)then
	        if (iter==maxiter-1) then 
		  if (divflux(J,K,L)<-0.1*eps(J,K,L)/crate) then
		    divflux(J,K,L)=-0.1*eps(J,K,L)/crate
		  endif
	        endif
                
 	           eps(J,K,L)=eps(J,K,L)-divflux(J,K,L)*crate  
	           !call tempfindspec(eps(J,K,L),rho(J,K,L),tempk(J,K,L))
                   call get_gamma(eps(J,K,L),rho(J,K,L),tempk(J,K,L)
     & ,mmw,gam)
                   if(tempk(J,K,L)<tirr)then
                      tempk(J,K,L)=tirr
                      !call EngFind(eps(J,K,L),rho(J,K,L),tirr)
                      call get_gamma_from_tk(eps(J,K,L),rho(J,K,L),
     &  tirr,mmw,gam)
                   endif
	      else
                tempk(J,K,L)=tbgrnd
                !!call EngFind(eps(J,K,L),rho(J,K,L),tbgrnd)
                call get_gamma_from_tk(eps(J,K,L),rho(J,K,L),
     &  tbgrnd,mmw,gam)
 
	
      	      endif
	      
	      totdflux=totdflux+divflux(J,K,L)*crate
     & *rhf(J)*rof3n*zof3n*dtheta
            !endif
         enddo
       enddo
       do K=1,KMAX2
         do J=1,JMAX2
           tau(J,K,L,1) = tau(J,K,L,1)*mu
         enddo
       enddo 
      enddo
!$OMP END DO 
!$OMP END PARALLEL

       ctime=ctime+crate
       iter=iter+1
!       print*, ctime,crate,delt
!       print *, "TEMPK ",tempk(JL,KL,LL),maxt,JL,KL,LL
!       print *, " "
!       print *, int_temp(JL,KL,LL),     
!     &         int_temp(JL,KL+1,LL),radflux(JL,KL,LL,1),
!     &         radflux(JL,KL,LL,3),tau(JL,KL,LL,1)
!       print *, " "
!       print *, " "

      enddo ! end while loop

      print *, "Maximum iteration count in hybrid, ", iter
      return


!!!$OMP DO SCHEDULE(STATIC)
!!      do j = jcool,jmax
!!        write(65,*) j,divflux(j,2,1)
!!!      enddo
!!!!!$OMP END DO


      end subroutine Hybrid
  
