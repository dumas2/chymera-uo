      SUBROUTINE SOURCE
      use eos
      IMPLICIT real*8 (a-h,o-z)      

#include "hydroparam.h"
#include "globals.h"
#include "units.h"

      integer jstart
      REAL*8 VT(JMAX2,KMAX2,LMAX),dume,limiter,surflux,gam,mu
      COMMON /SOURCESHARED/VT,surflux

      real*8 epred,tpred,ppred,ppred2
#if BINARY>0
      real*8 thisTime,rBin,thetaBin,omegaBin,massDisk,muMass
#endif

      DR=ROF3N
      DZ=ZOF3N

      HCOUNT = 0
      HHIT   = 0 !these are used to tally heating limiter
      
c  Should be ok to run the following loops over jmax1,kmax1 instead of
c  pot3jmax1,pot3kmax1 since we are only interested in sourcing the
c  regions where there is mass; the potential on the outer part of
c  the big grid does not really matter.

      limiter = den*phylim 

      if(jmin.gt.2) then
         jstart=jmin
      else   
         jstart=2
      endif 

#if BINARY>0
      call getBinaryPosition(thisTime,rBin,thetaBin,omegaBin,
     &   massDisk,muMass)
#endif

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k,LP,LM,l)                  &
!$OMP&  SHARED(dtheta,dz,delt,DR,jstart)
      DO L=1,LMAX
        LM=L-1
        LP=L+1
        IF(L.EQ.1) LM=LMAX
        IF(L.EQ.LMAX) LP=1
        DO K=2,KMAX1
 
C......SOURCE S
          DO J=jstart+1,JMAX1
            S(J,K,L)=S(J,K,L)+DELT*
     &          (R(J)*(RHO(J,K,L)+RHO(J-1,K,L))*
     &          ((OMEGA(J,K,L)+OMEGA(J-1,K,L))**2)*quarter*half
     &          -(RHO(J,K,L)+RHO(J-1,K,L))*(PHI(J,K,L)-PHI(J-1,K,L))
     &          /(two*DR)-(P(J,K,L)-P(J-1,K,L))/DR)
#if BINARY>0
     &    +omegaBin*omega(J,K,L)*rhf(J)*two*rho(J,K,L)*delt
#endif

          END DO

C......SOURCE T.
          DO J=jstart,JMAX1
            T(J,K,L)=T(J,K,L)-DELT*
     &        ( half*(RHO(J,K,L)+RHO(J,K-1,L))*(PHI(J,K,L)-PHI(J,K-1,L))
     &        /DZ+(P(J,K,L)-P(J,K-1,L))/DZ)
          END DO

C......SOURCE A.
          DO J=jstart,JMAX1
            A(J,K,L)=A(J,K,L)-DELT*(RHO(J,K,L)*(PHI(J,K,LP)-
     &           PHI(J,K,LM))+P(J,K,LP)-P(J,K,LM))/(two*DTHETA)
#if BINARY>0
     &           -rho(J,K,L)*omegaBin*(u(J+1,K,L)+u(J,K,L))*delt*rhf(J)
#endif
          END DO

        END DO
      END DO
!$OMP END PARALLEL DO

      CALL RADTRANS

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(dume,ppred,tpred,epred,j,k,LP,LM,&
!$OMP& l,mu,gam)                                                        &
!$OMP&  SHARED(etotfl,efl,rholmt,epslmt,pconv,muc,rhoconv,HHIT,         &
!$OMP& HCOUNT,dtheta,dz,delt,DR,jstart,zof3n,rof3n) REDUCTION(+:surflux,&
!$OMP& totirr,totdflux,totheat,totcool)
      CALL VELOCITY
      CALL VLIMIT
      CALL AVISC
  
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            LM=L-1
            LP=L+1
            IF(L.EQ.1) LM=LMAX
            IF(L.EQ.LMAX) LP=1
            DO  K=2,KMAX1

               DO  J=jstart+1,JMAX1
                  S(J,K,L)=S(J,K,L)-DELT*(RHO(J,K,L)*QRR(J,K,L)
     &                 -RHO(J-1,K,L)*QRR(J-1,K,L))/(DR)
               END DO

               DO  J=jstart,JMAX1
                  T(J,K,L)=T(J,K,L)-DELT*(RHO(J,K,L)*QZZ(J,K,L)
     &                 -RHO(J,K-1,L)*QZZ(J,K-1,L))/(DZ)
               END DO
         
               DO  J=jstart,JMAX1
                  A(J,K,L)=A(J,K,L)-DELT*(RHO(J,K,LP)*QTT(J,K,LP)
     &                 -RHO(J,K,LM)*QTT(J,K,LM))/(two*DTHETA) 
               END DO

               DO J=jstart,JMAX1
                  VT(J,K,L)=half*RHF(J)*(OMEGA(J,K,L)+OMEGA(J,K,LM))
               END DO

            END DO
         END DO
!$OMP END DO

C.......SOURCE EPS. (Includes artificial viscosity, L cooling, divflux,
c...... and irradiation.)

!$OMP DO SCHEDULE(STATIC) REDUCTION(+:HHIT,HCOUNT)
         DO L=1,LMAX
            LM=L-1
            LP=L+1
            IF(L.EQ.1) LM=LMAX
            IF(L.EQ.LMAX) LP=1
            DO K=2,KMAX
               DO J=jstart,JMAX

C    Artificial viscosity

                  Hgamma(J,K,L)=-RHO(j,k,l)*
     &                 (QRR(J,K,L)*(U(J+1,K,L)-U(J,K,L))/DR
     &                 +QZZ(J,K,L)*(W(J,K+1,L)-W(J,K,L))/DZ
     &                 +QTT(J,K,L)*(VT(J,K,LP)-VT(J,K,L))
     &                 /(DTHETA*RHF(J)))

C    Introduce a limiter just like in Igamma

                  if ((Hgamma(j,k,l).gt.zero)) then

                        HCOUNT = HCOUNT + 1

                     if((eps(j,k,l)/Hgamma(j,k,l)).lt.(theatlmt
     &                 *delt)) then

                          Hgamma(j,k,l)=eps(j,k,l)/(theatlmt*delt)
                  
                          HHIT = HHIT + 1

                     endif

                  endif

C    Depending on the type of cooling...


                  if ((ictype.eq.1).or.(ictype.eq.2)) then
                     Divflux(j,k,l)=zero
                  else if (ictype.eq.3) then
                     Lambda(j,k,l)=zero
                  else if (ictype.eq.4) then
                     if (tau(j,k,l,1).lt.(twothree)) then 
                        Divflux(j,k,l)=zero
                     else
                        Lambda(j,k,l)=zero
                    endif
                  elseif(ictype.eq.0)then
                     lambda(J,K,L)=zero
                     divflux(J,K,L)=zero
                  endif

                  epred = eps(J,K,L) - p(J,K,L) * (
     &                 (r(J+1)*u(J+1,K,L)-r(J)*u(J,K,L))/(rhf(J)*rof3n)
     &               + (vt(J,K,LP)-vt(J,K,L))/(dtheta*rhf(J))
     &               + (w(J,K+1,L)-w(J,K,L))/zof3n)*delt
     &               +
     &                 ( Hgamma(j,k,l)  - Lambda(j,k,l)
     &               -   Divflux(j,k,l) ) * delt

                  !print *, eps(J,K,L),p(J,K,L)
                  !print *, r(J+1),rhf(J),rof3n,zof3n,delt,dtheta
                  !print *, vt(J,K,LP),w(J,K+1,L),u(J+1,K,L)


                  !call TempFindSpec(epred,rho(J,K,L),tpred) ! replaced by full table
c                  call get_gamma(epred,rho(J,K,L),tpred,mu,gam)

                  !print *, epred,rho(J,K,L),tpred,mu,gam
                  !stop

                  ppred2= rho(J,K,L)*rhoconv
     &                  * bkmpcgs/(mu*pconv)*tpred

                  ppred = epred*(gamma-one)
                  eps(J,K,L) = eps(J,K,L) - half *(p(J,K,L)+ppred) * (
     &                 (r(J+1)*u(J+1,K,L)-r(J)*u(J,K,L))/(rhf(J)*rof3n)
     &               + (vt(J,K,LP)-vt(J,K,L))/(dtheta*rhf(J))
     &               + (w(J,K+1,L)-w(J,K,L))/zof3n)*delt
     &               +  
     &                 ( Hgamma(j,k,l) + Igamma(j,k,l) - Lambda(j,k,l) 
     &               -   Divflux(j,k,l) ) * delt               

                  eps(J,K,L)=max(eps(J,K,L),epslmt)

                  if(rho(j,k,l).le.rholmt) eps(j,k,l)=epslmt

               enddo
            enddo
         enddo
!$OMP END DO         

! master switch added by aaron to tally limited heating cells.

!                  print *, 'ppred =',ppred
!                  print *, 'ppred2 =',ppred2
!                  print *, 'rhoconv =',rhoconv
!                  print *, 'bkmpcgs =',bkmpcgs
!                  print *, 'mu =',mu
!                  print *, 'pconv =',pconv
!                  print *, 'tpred =',tpred
!$OMP MASTER
       if(HCOUNT==0)then
         dume=zero
       else
         dume = dble(HHIT)/dble(HCOUNT)*ten*ten
       endif

       if (dume > one)print *, " PERCENT TOTAL LIMITED HEATING CELLS ", 
     & dume

      etotfl = etotfl + efl 
      efl = zero
!$OMP END MASTER


C....Calculate total heating and total cooling 

!$OMP DO SCHEDULE(STATIC)
      do l=1,lmax
         do k=2,kmax
           do j=jcool,jmax
               totcool=totcool+Lambda(j,k,l)*0.5*dtheta*dz
     &              *(r(j+1)**2-r(j)**2)*delt
               totheat=totheat+Hgamma(j,k,l)*0.5*dtheta*dz
     &              *(r(j+1)**2-r(j)**2)*delt
               totdflux=totdflux+Divflux(j,k,l)*0.5*dtheta*dz
     &              *(r(j+1)**2-r(j)**2)*delt
               totirr=totirr+Igamma(j,k,l)*0.5*dtheta*dz
     &              *(r(j+1)**2-r(j)**2)*delt
               if((tau(j,k,l,1).ge.(twothree)).and.(tau(j,k+1,l,1).lt.(
     &              twothree))) surflux=surflux+Radflux(j,k,l,2)*half
     &              *dtheta*(r(j+1)**2-r(j)**2)*delt
               
            enddo
         enddo
      enddo
!$OMP END DO NOWAIT


C....Set quantities around the z-axis.
      
      if (jmin.gt.2) then

!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            DO K=2,KMAX2
               do j=1,jmin1
                  S(j,K,L)  = zero
                  T(j,K,L)  = zero
                  A(j,K,L)  = zero
                  U(j,k,l)  = zero
                  W(j,k,l)  = zero
                  EPS(j,K,L)= epslmt
               enddo
            END DO
         END DO
!$OMP END DO NOWAIT

      else   
         
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            LP=1+MOD(L-1+LMAX/2,LMAX)
            DO K=2,KMAX2
               S(1,K,L)  = -S(3,K,LP)
               T(1,K,L)  = T(2,K,LP)
               A(1,K,L)  = A(2,K,LP)
               EPS(1,K,L)= EPS(2,K,LP)
            END DO
         END DO
!$OMP END DO NOWAIT

      endif


!$OMP DO SCHEDULE(STATIC) 
      DO L=1,LMAX

C.....Set S and U on the z-axis.
        DO K=2,KMAX
          U(2,K,L)=zero
          S(2,K,L)=zero
        END DO

C....Set quantities on the bottom of the grid.
        DO J=1,JMAX2
          S(J,1,L)  = S(J,2,L)
          T(J,1,L)  = -T(J,3,L)
          A(J,1,L)  = A(J,2,L)
          EPS(J,1,L)= EPS(J,2,L)
          W(J,2,L)  = zero
          T(J,2,L)  = zero
        END DO

      END DO
!$OMP END DO 

      CALL CLEANUP
!$OMP END PARALLEL

      RETURN
      END
