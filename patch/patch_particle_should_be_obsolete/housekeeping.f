C***********************************************************************
C   LIMIT THE VELOCITIES
      SUBROUTINE VLIMIT
      IMPLICIT real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"
      integer jstart
      VLIM=two*SOUND
      if (jmin.gt.2) then
        jstart=jmin
      else
        jstart=2
      endif
!$OMP DO SCHEDULE(STATIC)
        DO L=1,lmax,1
          DO K=2,kmax1,1
            DO J=jstart,jmax1,1
            if (ABS(u(J,K,L)).gt.VLIM) then
              u(J,K,L)=SIGN(VLIM,u(J,K,L))
              s(J,K,L)=u(J,K,L)*half*(rho(J,K,L)+rho(J-1,K,L))
            endif
            if (ABS(w(J,K,L)).gt.VLIM) then
              w(J,K,L)=SIGN(VLIM,w(J,K,L))
              t(J,K,L)=w(J,K,L)*half*(rho(J,K,L)+rho(J,K-1,L))
            endif
            if (ABS(omega(J,K,L)*rhf(J)).gt.VLIM) then
              omega(J,K,L)=VLIM*omega(J,K,L)/(ABS(omega(J,K,L)*rhf(J)))
              JN(J,K,L)=omega(J,K,L)*rhf(J)**2
              a(J,K,L)=rho(J,K,L)*JN(J,K,L)
            endif
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
      RETURN
      END


C***********************************************************************
      SUBROUTINE DELTA(TINFO)
      IMPLICIT real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"
      logical TINFO
      common /TIMEST/INDX,ISOADI,ALLOW,dmax,CHGMAX
      real*8 factor,amin
      integer J1,K1,L1
      integer jstart
      real*8 SP(7),RD(jmax1),ZD(kmax1)
      real*8 amin_cap1
      integer j1_cap1
      integer k1_cap1
      integer l1_cap1
      real*8 SP_CAP7
      real*8 SP_CAP6
      real*8 SP_CAP5
      real*8 SP_CAP4
      real*8 SP_CAP3
      real*8 SP_CAP2
      real*8 SP_CAP1
C  FACTOR is fraction of courant time being used.
C  FACTORQ is the fraction of courant time for the av terms.
C  ALLOW is maximum allowable fractional change in maximum density,
C  during one time step.
      ALLOW=0.15d0
      DMAXO=dmax
      factor = six/ten
      amin=.1d+14
      if (jmin.gt.2) then
        jstart=jmin
      else
        jstart=2
      endif
      J1=0;K1=0;L1=0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(SP_CAP7,l1_cap1,k1_cap1,j1_cap1, &
!$OMP& amin_cap1,DELT3,DELT2,DELT1,speed7,speed6,speed5,SPEED4,SPEED3,  &
!$OMP& SPEED2,SPEED1,l,k,j,SP_CAP1,SP_CAP2,SP_CAP3,SP_CAP4,SP_CAP5,     &
!$OMP& SP_CAP6)                                                         &
!$OMP&  SHARED(dmax,DTHETA,muc,jstart,l1,k1,j1,amin)
!$OMP DO SCHEDULE(STATIC)
        DO j=jstart,jmax1,1
        RD(j)=r(j+1)-r(j)
        ENDDO
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
        DO k=2,kmax1,1
        ZD(k)=z(k+1)-z(k)
        ENDDO
!$OMP END DO
      amin_cap1=amin
      j1_cap1=j1
      k1_cap1=k1
      l1_cap1=l1
      SP_CAP7=SP(7)
      SP_CAP6=SP(6)
      SP_CAP5=SP(5)
      SP_CAP4=SP(4)
      SP_CAP3=SP(3)
      SP_CAP2=SP(2)
      SP_CAP1=SP(1)
!$OMP DO SCHEDULE(STATIC)
        DO l=1,lmax,1
          DO k=2,kmax1,1
            DO j=jstart,jmax1,1
            SPEED1=sqrt(gamma1(j,k,l)*p(j,k,l)/(rho(j,k,l)))
            SPEED2=half*ABS(u(j+1,k,l)+u(j,k,l))
            SPEED3=half*ABS(w(j,k+1,l)+w(j,k,l))
            SPEED4=ABS(JN(j,k,l)/rhf(j))
            speed5=four*sqrt(qrr(j,k,l))
            speed6=four*sqrt(qzz(j,k,l))
            speed7=four*sqrt(qtt(j,k,l))
            DELT1=RD(j)/(SPEED1+SPEED2+speed5)
            DELT2=ZD(k)/(SPEED1+SPEED3+speed6)
            DELT3=rhf(j)*DTHETA/(SPEED1+SPEED4+speed7)
            DELT1=MIN(DELT1,DELT2)
            DELT1=MIN(DELT1,DELT3)
            if (DELT1.LT.amin_cap1) then
              amin_cap1=DELT1
              j1_cap1=j
              k1_cap1=k
              l1_cap1=l
              SP_CAP1=SPEED1
              SP_CAP2=SPEED2
              SP_CAP3=SPEED3
              SP_CAP4=SPEED4
              SP_CAP5=speed5
              SP_CAP6=speed6
              SP_CAP7=speed7
            endif
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO 
!$OMP CRITICAL
      if (amin_cap1.LT.amin) then
        amin=amin_cap1
        j1=j1_cap1
        k1=k1_cap1
        l1=l1_cap1
        SP(7)=SP_CAP7
        SP(6)=SP_CAP6
        SP(5)=SP_CAP5
        SP(4)=SP_CAP4
        SP(3)=SP_CAP3
        SP(2)=SP_CAP2
        SP(1)=SP_CAP1
      endif
!$OMP END CRITICAL
! Now chop this Courant-determined time by density change criterion. We
! scrupulously avoid possible division by 0 when calculating F1 and F2,
! and possibility of rounding errors yielding  DELT(new) > DELT(old).
      if (jmin.gt.2) then
!$OMP MASTER
        dmax=rho(jmin,2,1)
!$OMP END MASTER
      endif
!$OMP DO SCHEDULE(STATIC) REDUCTION(MAX:dmax)
        DO l=1,lmax,1
          DO k=2,kmax,1
            DO j=jstart,jmax,1
            dmax=MAX(dmax,rho(j,k,l))
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
      delt = amin*factor
      CHGMAXN=ABS(1.0-DMAXO/dmax)
      if ((INDX.LE.1).AND.(DMAXO.LT.den)) then
        F1=one
        F2=one
      else
        F1=ALLOW/MAX(ALLOW,CHGMAXN)
        F2=ALLOW/MAX(ALLOW,CHGMAX)
        CHGMAX=CHGMAXN
      endif
      delt=MIN(delt,delt*F1*F2)
      if (TINFO) then
        WRITE(3,102)delt,j1,k1,l1,SP,ALLOW,F1,F2,CHGMAX
      endif
c     IF(AMINQ.LT.AMIN) WRITE(3,100) DELT,J1,K1,L1,SP,ALLOW,F1,F2,CHGMAX
  100 FORMAT('TINFQ:',2X,1PE10.2,3I5,3X,1P7E10.2,3X,1P4E10.2)
! 101 FORMAT('TINFO:',2X,1PE10.2,3I5,3X,1P7E10.2,3X,1P4E10.2)
  102 FORMAT('TINFO:',/,'  DELT      ',1PE12.3,/,'  J1,K1,L1  ',3I5,/,  &
     &  '  SP(1:7)   ',1P7E12.3,/,'  ALLOW     ',1PE12.3,/,             &
     &  '  F1,F2     ',1P2E12.3,/,'  CHGMAX    ',1PE12.3)
      RETURN
      END


C***********************************************************************
      SUBROUTINE CLEANUP
      IMPLICIT real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"
      integer jstart
      dr=rof3n
      if (jmin.gt.2) then
        jstart=1
      else
        jstart=2
      endif
C.... SET QUANTITIES AT THE GRID EDGE.
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:tmassadd)
        do l=1,lmax,1
          do k=2,kmax2,1
            do j=jmax1,jmax2,1
            s(j,k,l)=abs(s(j-1,k,l))
            t(j,k,l)=t(j-1,k,l)
            a(j,k,l)=a(j-1,k,l)
            rho(j,k,l)=rho(j-1,k,l)
#if PASSIVE>0
            do I=1,PAS
              PASSFLUX(J,K,L,I)=PASSFLUX(J-1,K,L,I)
            enddo
#endif
            eps(j,k,l)=eps(j-1,k,l)
            enddo
          enddo
          do k=kmax1,kmax2,1
            do j=jstart,jmax2,1
            s(j,k,l)=s(j,k-1,l)
            t(j,k,l)=abs(t(j,k-1,l))
            a(j,k,l)=a(j,k-1,l)
            rho(j,k,l)=rho(j,k-1,l)
#if PASSIVE>0
            do I=1,PAS
              PASSFLUX(J,K,L,I)=PASSFLUX(J,K-1,L,I)
            enddo
#endif
            eps(j,k,l)=eps(j,k-1,l)
            enddo
          enddo
C....Limit the lowest density.
          do k=2,kmax2,1
            do j=2,jmax2,1
            if (((((rho(j,k,l).lt.rholmt).and.(j.le.jmax)).and.(k.le.   &
     &      kmax)).and.(j.ge.jmin)).and.(k.ge.2)) then
              tmassadd=tmassadd+(rholmt-rho(j,k,l))*dtheta*zof3n*(r(j+1)&
     &        **2-r(j)**2)
            endif
            rho(j,k,l)=MAX(rho(j,k,l),rholmt)
            eps(j,k,l)=MAX(eps(j,k,l),epslmt)
#if PASSIVE>0
            do I=1,PAS
              PASSFLUX(J,K,L,I)=MAX(PASSFLUX(J,K,L,I),rholmt)
            enddo
#endif
            enddo
          enddo
        enddo
!$OMP END DO
C.....Set quantities around the z- axis.
      if (jmin.gt.2) then
C.....Use linear extrapolation for the ghost zones in the inner bdry (jm
C.....S(jmin)=Sjmin+1). Restrict calculated values if the inner value te
C.....acm2001
!$OMP DO SCHEDULE(STATIC)
          do l=1,lmax,1
            do k=2,kmax2,1
              do j=1,jmin2-1,1
              s(1,k,l)=0.0
              t(1,k,l)=0.0
              a(1,k,l)=0.0
              rho(1,k,l)=rholmt
              eps(1,k,l)=epslmt
              enddo
            s(jmin,k,l)=2*s(jmin+1,k,l)-s(jmin+2,k,l)
            if (abs(s(jmin,k,l)).gt.abs(s(jmin+1,k,l))) then
              s(jmin,k,l)=s(jmin+1,k,l)
            endif
            s(jmin,k,l)=-abs(s(jmin,k,l))
              do j=jmin1,jmin2,-1
              rho(j,k,l)=-((rho(jmin+1,k,l)-rho(jmin,k,l))*(rhf(jmin)-  &
     &        rhf(j))/dr)+rho(jmin,k,l)
              if (rho(j,k,l).gt.rho(j+1,k,l)) then
                rho(j,k,l)=rho(j+1,k,l)
              endif
              if (rho(j,k,l).lt.rholmt) then
                rho(j,k,l)=rholmt
              endif
              eps(j,k,l)=-((eps(jmin+1,k,l)-eps(jmin,k,l))*(rhf(jmin)-  &
     &        rhf(j))/dr)+eps(jmin,k,l)
              if (eps(j,k,l).gt.eps(j+1,k,l)) then
                eps(j,k,l)=eps(j+1,k,l)
              endif
              if ((eps(j,k,l).lt.epslmt).or.(rho(j,k,l).le.rholmt)) then
                eps(j,k,l)=epslmt
              endif
              a(j,k,l)=-((a(jmin+1,k,l)-a(jmin,k,l))*(rhf(jmin)-rhf(j))/&
     &        dr)+a(jmin,k,l)
              if (a(j,k,l).gt.a(j+1,k,l)) then
                a(j,k,l)=a(j+1,k,l)
              endif
              if ((a(j,k,l).lt.0.0).or.(rho(j,k,l).le.rholmt)) then
                a(j,k,l)=0.0
              endif
              t(j,k,l)=-((t(jmin+1,k,l)-t(jmin,k,l))*(rhf(jmin)-rhf(j))/&
     &        dr)+t(jmin,k,l)
              if (abs(t(j,k,l)).gt.abs(t(j+1,k,l))) then
                t(j,k,l)=t(j+1,k,l)
              endif
              if (rho(j,k,l).le.rholmt) then
                t(j,k,l)=0.0
              endif
              s(j,k,l)=-((s(jmin+2,k,l)-s(jmin+1,k,l))*(r(jmin+1)-r(j))/&
     &        dr)+s(jmin+1,k,l)
              if (abs(s(j,k,l)).gt.abs(s(j+1,k,l))) then
                s(j,k,l)=s(j+1,k,l)
              endif
              s(j,k,l)=-abs(s(j,k,l))
              enddo
            enddo
          enddo
!$OMP END DO nowait
      else
!$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            LP=1+MOD(L-1+LMAX/2,LMAX)
            DO K=2,KMAX2
               S(1,K,L)   = -S(3,K,LP)
               T(1,K,L)   = T(2,K,LP)
               A(1,K,L)   = A(2,K,LP)
               RHO(1,K,L) = RHO(2,K,LP)
               EPS(1,K,L) = EPS(2,K,LP)
#if PASSIVE>0
            do I=1,PAS
              PASSFLUX(1,K,L,I)=PASSFLUX(2,K,LP,I)
            enddo
#endif
            ENDDO
         ENDDO
!$OMP END DO

      endif
C.....Set quantities below the equatorial plane.
!$OMP DO SCHEDULE(STATIC)
        do l=1,lmax,1
          do j=1,jmax2,1
          s(j,1,l)=s(j,2,l)
          t(j,1,l)=-t(j,3,l)
          a(j,1,l)=a(j,2,l)
          rho(j,1,l)=rho(j,2,l)
          eps(j,1,l)=eps(j,2,l)
#if PASSIVE>0
            do I=1,PAS
              PASSFLUX(J,1,L,I)=PASSFLUX(J,2,L,I)
            enddo
#endif
          enddo
        enddo
!$OMP END DO
      RETURN
      END


C***********************************************************************
      SUBROUTINE CENTMASS(DISP)
      IMPLICIT real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"
      real*8 THETA(lmax+1)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(RM,K,J,L)                        &
!$OMP&  SHARED(CM,Y,X,DTHETA)
!$OMP DO SCHEDULE(STATIC)
        DO L=1,lmax+1,1
        THETA(L)=DTHETA*(L-1)
        ENDDO
!$OMP END DO
!$OMP MASTER
      X=0.0
      Y=0.0
      CM=0.0
!$OMP END MASTER
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:CM,Y,X)
        DO J=2,jmax,1
          DO K=2,kmax,1
            DO L=1,lmax,1
            RM=rho(J,K,L)*(r(J+1)**3-r(J)**3)
            X=X+(SIN(THETA(L+1))-SIN(THETA(L)))*RM
            Y=Y+(COS(THETA(L))-COS(THETA(L+1)))*RM
            CM=CM+rho(J,K,L)*(r(J+1)**2-r(J)**2)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
      X=2.*X/(DTHETA*CM*3.)
      Y=2.*Y/(DTHETA*CM*3.)
      DISP=SQRT(X**2+Y**2)
      IF ((X.EQ.0.0).AND.(Y.GE.0.0)) THEN
        ALFA=PI/2.0
      ELSE
        IF ((X.EQ.0.0).AND.(Y.LT.0.0)) THEN
          ALFA=3*PI/2.0
        ELSE
          ALFA=ATAN(Y/X)
          IF (X.LT.0.0) THEN
            ALFA=PI+ALFA
          ENDIF
          IF ((X.GT.0.0).AND.(Y.LE.0.0)) THEN
            ALFA=2*PI+ALFA
          ENDIF
        ENDIF
      ENDIF
      ALF=ALFA*180/PI
      WRITE(3,40)TIME,X,Y,DISP,ALF
   40 FORMAT(' CENTER OF MASS:   TIME='1PE12.4,'  X=',1PE12.4,'     Y=',&
     &1PE12.4,'   DISPLACEMENT=',1PE12.4,'   ANGLE=',0PF7.2)
      RETURN
      END


C***********************************************************************
      SUBROUTINE CLEARUP
      IMPLICIT real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"
      integer jstart
      if (jmin.gt.2) then
        jstart=jmin2
      else
        jstart=2
      endif
!$OMP DO SCHEDULE(STATIC)
        DO L=1,lmax,1
          DO K=1,kmax2,1
            DO J=jstart,jmax2,1
            rho(J,K,L)=ABS(rho(J,K,L))
            eps(J,K,L)=ABS(eps(J,K,L))
#if PASSIVE>0
            do I=1,PAS
              PASSFLUX(J,K,L,I)=ABS(PASSFLUX(J,K,L,I))
            enddo
#endif
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO nowait
      RETURN
      END


