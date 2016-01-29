      subroutine Hybrid()
      implicit none
#include "hydroparam.h"
#include "globals.h"
#include "units.h"

      real(KIND=8) :: dumb ,oldcrate
      real(KIND=8) :: limiter,aconst,bconst,cconst, mu
      real(KIND=8) :: FluxLmDf,sum,crate,ctime,maxt,taufld

      integer :: J,K,L,I,LM,LP,ITER,JL,KL,LL

      REAL*8 Oross(jmax2,kmax2,lmax),Oplck(jmax2,kmax2,lmax)
     &     ,Oabs(jmax2,kmax2,lmax),Otot(jmax2,kmax2,lmax),absfrac,
     &      oplck_env(JMAX2,KMAX2,LMAX)

      COMMON /COOLINGSHARED/Oross,Oplck,Oabs,Otot,sum,oplck_env

      mu = one/sqrt(three)
      efl = zero
      limiter = den*phylim 
      LO=0
      CCOUNTER = 0
      HCOUNTER = 0
      ctime=zero
      iter=0
      taufld=one

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(LM,dumb,cconst,bconst,aconst,I,  &
!$OMP& LP,J,K,L)                                                        &
!$OMP&  SHARED(ccounter,hcounter,LO,delt,dtheta,rof3n,limiter,zof3n,    &
!$OMP& sconv,sigma,pi,mu,tconv) REDUCTION(+:efl)
!$OMP DO SCHEDULE(STATIC)
      do L =  1, LMAX 
       do K = 1, KMAX2
         do J = 1, JMAX2
           tau(J,K,L,4) = zero
         enddo
       enddo
       do J = 1, JMAX2
         init_int_in(J,L) = (tenvk/tconv/rhf(J)+ten/tconv)**4
!        init_int_in(J,L) =
!     &      (6d2/tconv/sqrt(rhf(J)+(rhf(J)/six)**2))**4
     &                    / pi*sigma
         KFITA(J,L) = 2
       enddo
      enddo
!$OMP END DO

!$OMP DO SCHEDULE(STATIC) ! set optical depths of positive plane in l_tau
      do L = 1, LMAX 
        do K = KMAX1, 2, -1
          do J = JCOOL, JMAX

            ! interpolate between rosseland and planck means.

            dtau_z(J,K,L) = rho(J,K,L)*sconv*zof3n *oross(J,K,L)/mu

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
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(LM,dumb,cconst,bconst,aconst,I,  &
!$OMP& LP,J,K,L)                                                        &
!$OMP&  SHARED(ccounter,hcounter,LO,delt,dtheta,rof3n,limiter,zof3n,    &
!$OMP& sconv,sigma,pi,mu,tconv) REDUCTION(+:efl)
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
!$OMP END DO NOWAIT

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
!$OMP END DO NOWAIT

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
!$OMP END DO NOWAIT

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
     &              + init_int_in(J,L)*exp(-tau(J,K,L,1))
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
!$OMP END DO NOWAIT

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
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)  ! calculate divflux
      do L = 1, LMAX 
        LM = modulo(L-2,LMAX)+1
        LP = modulo(L,LMAX)+1
        do K = 2, KMAX
          do J = JCOOL, JMAX
               
          if (rho(J,K,L) >= limiter ) then
             divflux(J,K,L) =                                          
     &         (int_temp(J,K+1,L)-int_temp(J,K,L))/zof3n      
     &       + temporary(J,K,L)
          endif

          enddo
        enddo
      enddo
!$OMP END DO NOWAIT

      JL=JMIN;KL=2;ll=2
C$OMP DO SCHEDULE(STATIC) REDUCTION(+:CCOUNTER,HCOUNTER,LO)
!$OMP& REDUCTION(min:crate) REDUCTION(max:maxt) ! calculate divflux
        do L = 1, LMAX
          do K = 2, KMAX1
            do J = JCOOL, JMAX

              if(abs(divflux(J,K,L))>zero.and.eps(J,K,L)
     &            >engtable(1)/engconv*rho(J,K,L)*1.1)then
               dumb=abs(eps(J,K,L)/(divflux(J,K,L)*four))
               !crate=min(crate,abs(eps(J,K,L)/(divflux(J,K,L)*four)))
               if (dumb<crate)then 
                 JL=J;KL=K;LL=L
                 crate=dumb
               endif
              endif    
              maxt=max(tempk(J,K,L),maxt)

           enddo
         enddo
       enddo
!$OMP END DO
!$OMP MASTER
      crate=min(crate,(delt-ctime))
      oldcrate=crate
!$OMP END MASTER
!$OMP DO SCHEDULE(STATIC)
      do L = 1,LMAX 
       do K = 1,KMAX2
         do J = 1,JMAX2
           eps(J,K,L)=max(eps(J,K,L)-divflux(J,K,L)*crate,
     &                   rho(J,K,L)*engtable(1)/engconv)
           call TempFindSpec(eps(J,K,L),rho(J,K,L),tempk(J,K,L))
           tau(J,K,L,1) = tau(J,K,L,1)*mu
         enddo
       enddo
      enddo
!$OMP END DO 
!$OMP END PARALLEL

      ctime=ctime+crate
      iter=iter+1
!      print*, ctime,crate,delt,maxt,tempk(JL,KL,LL),iter,JL,KL,LL
!      print *, " "
!      print *, tempk(JL,KL,LL),     
!     &         int_temp(JL,KL+1,LL),radflux(JL,KL,LL,1),
!     &         radflux(JL,KL,LL,3),tau(JL,KL,LL,1)
!      print *, " "
!      print *, " "
      enddo

      print *, "Maximum iteration count in hybrid, ", iter
      return

      end subroutine Hybrid
  
