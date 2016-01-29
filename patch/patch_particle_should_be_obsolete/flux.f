
c*******************************************************************************
c.....Translate cvmgp and cvmgz into equivalent if statements
      FUNCTION SLOPE(f,b,c)
      implicit real*8 (a-h,o-z)
      real*8 f,b,c,tmp,tmm
      tmp=0.0
      if ((b-c)*(f-b).ge.0.0) tmp=1.0
      tmm=f-c
      if (f-c.eq.0.0) tmm=1.0
      slope=2.*(b-c)*(f-b)*tmp/tmm
      return
      end


c*******************************************************************************
c.........and van leer interpolation
      FUNCTION VLI(qf,qb,slopf,slopb,vel,dt,dx)
      implicit real*8 (a-h,o-z)
c     real*8 qf,qb,slopf,slopb,vel,dt,dx,tmn
      if(vel.lt.0.0) then
        vli=qf-(dx+vel*dt)*slopf/(2.*dx)
      else
        vli=qb+(dx-vel*dt)*slopb/(2.*dx)
      endif
      return
      end


c*******************************************************************************

#if PASSIVE>0
      SUBROUTINE FLUX(SS,TT,AA,RRHO,EEPS,PPASSFLUX)
#else
      SUBROUTINE FLUX(SS,TT,AA,RRHO,EEPS)
#endif
      IMPLICIT real*8 (a-h,o-z)      
#include "hydroparam.h"
#include "globals.h"
      
#if PASSIVE>0
      REAL*8 SS(JMAX2,KMAX2,LMAX),
     &       TT(JMAX2,KMAX2,LMAX),
     &       AA(JMAX2,KMAX2,LMAX),
     &       RRHO(JMAX2,KMAX2,LMAX),
     &       EEPS(JMAX2,KMAX2,LMAX),
     &       PPASSFLUX(JMAX2,KMAX2,LMAX,PAS)
#else
      REAL*8 SS(JMAX2,KMAX2,LMAX),
     &       TT(JMAX2,KMAX2,LMAX),
     &       AA(JMAX2,KMAX2,LMAX),
     &       RRHO(JMAX2,KMAX2,LMAX),
     &       EEPS(JMAX2,KMAX2,LMAX)
#endif


C The following arrays are local to subroutine flux, and should not need 
C to be placed in common. But the OpenMP version of the code may fail unless 
C they are.
      REAL*8 
     &       VINV(JMAX2),
     &       VINVH(JMAX2),
     &       RD(JMAX2),
     &       RHD(JMAX2),
     &       VR(JMAX2,KMAX2,LMAX),
     &       VZ(JMAX2,KMAX2,LMAX),
     &       VT(JMAX2,KMAX2,LMAX),
     &       VLOP(JMAX2,KMAX2,LMAX),
     &       FR(JMAX2,KMAX2,LMAX),
     &       FZ(JMAX2,KMAX2,LMAX),
     &       FT(JMAX2,KMAX2,LMAX)

      COMMON  /FLUXPRIVATE/VINV,VINVH,RD,RHD
      COMMON  /FLUXSHARED/VR,VZ,VT,VLOP,FR,FZ,FT
      
      integer jstart,LP,LM,I,J,K,L
      
C The quantities are centered as follows:
C    S(J,K,L) are centered at r-directioned cell surface centers, i.e.
C        at [R(J), ZHF(K), (L-1/2)*DTHETA].
C    T(J,K,L) are centered at z-directioned cell surface centers, i.e.
C        at [RHF(J), Z(K), (L-1/2)*DTHETA].
C    A(J,K,L) AND RHO(J,K,L) are centered at cell centers, i.e. at
C        [RHF(J), ZHF(K), (L-1/2)*DTHETA].

      DR=ROF3N
      DZ=ZOF3N

C Verify boundary conditions are met at axis and equatorial plane

      if (jmin.gt.2) then        

!$OMP DO SCHEDULE(STATIC)
         do l=1,lmax
            do k=1,kmax1
               do j=1,jmin2-1
                  fz(j,k,l)=0.d0
                  fr(j,k,l)=0.d0
               enddo
            enddo
         enddo
!$OMP END DO NOWAIT

      endif

!$OMP DO SCHEDULE(STATIC)
      do  l=1,lmax
         do j=2,jmax1
            fz(j,2,l)=0.d0
         end do
         do k=2,kmax1
            fr(2,k,l)=0.d0
         end do
      end do
!$OMP END DO NOWAIT
      
    
C Set up a few arrays to be used in this routine.
C rd(j) is located at rhf(j), rhd(j) is located at r(j).
C Vinv and vinvh are inverses of control volumes at rhf(j) and r(j).

!$OMP DO SCHEDULE(STATIC)
      DO J=3,JMAX2
        RD(J)  = R(J)**2-R(J-1)**2
        RHD(J) = RHF(J)**2-RHF(J-1)**2
        VINV(J)  = two/(DTHETA*DZ*RD(J))
        VINVH(J) = two/(DTHETA*DZ*RHD(J))
      END DO
!$OMP END DO     

C######FLUX S
C...........VELOCITIES TO FLUX S(J,K,L) IN THE R-DIRECTION AT RHF(J-1)
C...........VELOCITIES TO FLUX S(J,K,L) IN THE Z-DIRECTION AT Z(K)
C...........VELOCITIES TO FLUX S(J,K,L) IN THE THETA-DIRECTION AT L.
C...........   (AREA WEIGHTED AVERAGE)



!$OMP DO SCHEDULE(STATIC)
      DO 100 L=1,LMAX
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=3
        endif  
        DO K=2,KMAX2
          DO J=jstart,JMAX2
            VR(J,K,L)=half*(U(J-1,K,L)+U(J,K,L))
            VZ(J,K,L)=half*(W(J-1,K,L)+W(J,K,L))
            RDM=R(J)**2-RHF(J-1)**2
            RDP=RHF(J)**2-R(J)**2
            VT(J,K,L)=((OMEGA(J,K,L)+OMEGA(J,K,LM))*RHF(J)*RDP+
     &           (OMEGA(J-1,K,L)+OMEGA(J-1,K,LM))*RHF(J-1)*RDM)/
     $           (two*RHD(J)) 
          END DO
        END DO

         
C............VAN LEER SLOPES FOR S(J,K,L) IN THE R-DIRECTION

        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif  
        
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(S(J+1,K,L),S(J,K,L),S(J-1,K,L))
          END DO
        END DO
         
C............FLUXES FOR S(J,K,L) IN THE R-DIRECTION AT RHF(J-1)

        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=3
        endif  
        
        DO K=2,KMAX1
c          FR(3,K,L)=0.d0
          DO J=jstart,JMAX1
            FR(J,K,L)=VLI(S(J,K,L),S(J-1,K,L),VLOP(J,K,L),VLOP(J-1,K,L),
     $           VR(J,K,L),DELT,DR)*VR(J,K,L)*(RHF(J-1)-VR(J,K,L)
     $           *0.5*DELT)*DZ*DTHETA
          END DO
        END DO     

C.............VAN LEER SLOPES FOR S(J,K,L) IN THE Z-DIRECTION

        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=3
        endif  

        DO J=jstart,JMAX1
          VLOP(J,1,L)=SLOPE(S(J,2,L),S(J,1,L),S(J,3,L))
        END DO
        DO K=3,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(S(J,K+1,L),S(J,K,L),S(J,K-1,L))
          END DO
        END DO

C............FLUXES FOR S(J,K,L) IN THE Z-DIRECTION AT Z(K)

        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=3
        endif 

        DO K=3,KMAX1
          DO J=jstart,JMAX1
            FZ(J,K,L)=VLI(S(J,K,L),S(J,K-1,L),VLOP(J,K,L),VLOP(J,K-1,L),
     $           VZ(J,K,L),DELT,DZ)*VZ(J,K,L)*RHD(J)*0.5*DTHETA
          END DO
        END DO
  100 CONTINUE
!$OMP END DO NOWAIT

C............VAN LEER SLOPES FOR S(J,K,L) IN THE THETA-DIRECTION
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1 
        IF(L.EQ.LMAX) LP=1 
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=3
        endif         
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(S(J,K,LP),S(J,K,L),S(J,K,LM))
          END DO
        END DO
      END DO
C$OMP END DO 
      
C............FLUXES FOR S(J,K,L) IN THE L DIRECTION AT L.
C            parallelizable once all vlop(j,k,l) calculated.

C$OMP DO SCHEDULE(STATIC)
      do l=1,lmax
         if(l.eq.1) then
            lm = lmax
         else
            lm = l-1
         end if
         if (jmin.gt.2) then
            jstart=jmin
         else
            jstart=3
         endif 
         do k=2,kmax1
            do j=jstart,jmax1
               RDTHETA=R(J)*DTHETA
               FT(J,K,L)=VLI(S(J,K,L),S(J,K,LM),VLOP(J,K,L),VLOP(J,K,LM)
     $              ,VT(J,K,L),DELT,RDTHETA)*VT(J,K,L)*DZ*DR
               
            end do
         end do
      end do
C$OMP END DO

      
C.............COMPUTE NEW S(J,K,L) BY FLUXING
c     This loop parallelizable only after all FT(j,k,l) have been
c     calculated.
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
         LP=L+1
         IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=3
        endif 
         DO K=2,KMAX
            DO J=jstart,JMAX
               SS(J,K,L)=SS(J,K,L)
     &              -(FR(J+1,K,L)-FR(J,K,L)+FZ(J,K+1,L)
     &              -FZ(J,K,L)+FT(J,K,LP)-FT(J,K,L))*VINVH(J)*DELT
            END DO
         END DO
      END DO
C$OMP END DO NOWAIT
      
      
C#######FLUX T
C...........VELOCITIES TO FLUX T(J,K,L) IN THE R-DIRECTION AT R(J)
C...........VELOCITIES TO FLUX T(J,K,L) IN THE Z-DIRECTION AT ZHF(K-1)
C...........VELOCITIES TO FLUX T(J,K,L) IN THE THETA-DIRECTION AT L.
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=2,KMAX2
          DO J=jstart,JMAX2
            VR(J,K,L)=0.5*(U(J,K-1,L)+U(J,K,L))
            VZ(J,K,L)=0.5*(W(J,K,L)+W(J,K-1,L))
            VT(J,K,L)=0.25*(OMEGA(J,K,L)+OMEGA(J,K,LM)+
     &           OMEGA(J,K-1,L)+OMEGA(J,K-1,LM))*RHF(J)
          END DO
        END DO
      END DO
C$OMP END DO NOWAIT
      
C...........VAN LEER SLOPE FOR T IN THE R-DIRECTION
C$OMP DO SCHEDULE(STATIC)
      DO 210 L=1,LMAX
         if (jmin.gt.2) then
            jstart=jmin1
         else
            jstart=2
         endif 
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(T(J+1,K,L),T(J,K,L),T(J-1,K,L))
          END DO
        END DO 

C...........FLUXES FOR T(J,K,L) IN THE R-DIRECTION AT R(J)

        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=3
        endif 
        DO K=3,KMAX1
          DO J=jstart,JMAX1
            FR(J,K,L)=VLI(T(J,K,L),T(J-1,K,L),VLOP(J,K,L),VLOP(J-1,K,L),
     &           VR(J,K,L),DELT,DR)*VR(J,K,L)*(R(J)-VR(J,K,L)*0.5
     &           *DELT)*DZ*DTHETA
          END DO
        END DO
         
C...........VAN LEER SLOPE FOR T IN THE Z-DIRECTION

        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO J=jstart,JMAX1
          XT= -T(J,4,L)
          VLOP(J,1,L)=SLOPE(T(J,2,L),T(J,1,L),XT)
        END DO
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(T(J,K+1,L),T(J,K,L),T(J,K-1,L))
          END DO
        END DO
         
C...........FLUXES FOR T(J,K,L) IN THE Z-DIRECTION AT ZHF(K-1)

        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=3,KMAX1
          DO J=jstart,JMAX1
            FZ(J,K,L)=VLI(T(J,K,L),T(J,K-1,L),VLOP(J,K,L),VLOP(J,K-1,L),
     $           VZ(J,K,L),DELT,DZ)*VZ(J,K,L)*RD(J+1)*0.5*DTHETA
          END DO
        END DO
  210 CONTINUE
C$OMP END DO NOWAIT
      
C...........VAN LEER SLOPE FOR T IN THE THETA-DIRECTION
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(T(J,K,LP),T(J,K,L),T(J,K,LM))
          END DO
        END DO
      END DO
C$OMP END DO
      
C...........FLUXES FOR T(J,K,L) IN THE THETA-DIRECTION
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=3,KMAX1
          DO J=jstart,JMAX1
            RDTHETA=RHF(J)*DTHETA
            FT(J,K,L)=VLI(T(J,K,L),T(J,K,LM),VLOP(J,K,L),VLOP(J,K,LM),
     $           VT(J,K,L),DELT,RDTHETA)*VT(J,K,L)*DZ*DR
          END DO
        END DO
      END DO
C$OMP END DO
      
C...........COMPUTE NEW T(J,K,L) BY FLUXING
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1
        IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=3,KMAX
          DO J=jstart,JMAX
            TT(J,K,L)=TT(J,K,L)-(FR(J+1,K,L)-FR(J,K,L)+FZ(J,K+1,L)
     &           -FZ(J,K,L)+FT(J,K,LP)-FT(J,K,L))*VINV(J+1)*DELT
          END DO
        END DO
      END DO
C$OMP END DO nowait
      
      
C#######FLUX A
C...........VELOCITIES TO FLUX A(J,K,L) IN THE R-DIRECTION AT R(J)
C...........VELOCITIES TO FLUX A(J,K,L) IN THE Z-DIRECTION AT Z(K)
C...........VELOCITIES TO FLUX A(J,K,L) IN THE THETA-DIRECTION AT L.
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=2,KMAX2
          DO J=jstart,JMAX2
            VR(J,K,L)=U(J,K,L)
            VZ(J,K,L)=W(J,K,L)
            VT(J,K,L)=0.5*(OMEGA(J,K,L)+OMEGA(J,K,LM))*RHF(J)
          END DO
        END DO
      END DO 
C$OMP END DO NOWAIT

      
C...........VAN LEER SLOPE FOR A IN THE R-DIRECTION
C$OMP DO SCHEDULE(STATIC)

      DO 310 L=1,LMAX
        if (jmin.gt.2) then
           jstart=jmin1
        else
            jstart=2
        endif 
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(A(J+1,K,L),A(J,K,L),A(J-1,K,L))
          END DO

C...........FLUXES FOR A(J,K,L) IN THE R-DIRECTION AT R(J)

          if (jmin.gt.2) then
             jstart=jmin
          else
             jstart=3
          endif 
          DO J=jstart,JMAX1
            FR(J,K,L)=VLI(A(J,K,L),A(J-1,K,L),VLOP(J,K,L),VLOP(J-1,K,L),
     $           VR(J,K,L),DELT,DR)*VR(J,K,L)*(R(J)-VR(J,K,L)*0.5
     $           *DELT)*DZ*DTHETA
          END DO
        END DO
         

C...........VAN LEER SLOPE FOR A IN THE Z-DIRECTION

        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO J=jstart,JMAX1
          VLOP(J,1,L)=SLOPE(A(J,2,L),A(J,1,L),A(J,3,L))
        end do
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(A(J,K+1,L),A(J,K,L),A(J,K-1,L))
          END DO
        END DO

C...........FLUXES FOR A(J,K,L) IN THE Z-DIRECTION AT Z(K)

        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif         
        DO K=3,KMAX1
          DO J=jstart,JMAX1
            FZ(J,K,L)=VLI(A(J,K,L),A(J,K-1,L),VLOP(J,K,L),VLOP(J,K-1,L),
     $           VZ(J,K,L),DELT,DZ)*VZ(J,K,L)*RD(J+1)*0.5*DTHETA
          END DO
        END DO
  310 END DO
C$OMP END DO
      
C...........VAN LEER SLOPE FOR A IN THE THETA-DIRECTION
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(A(J,K,LP),A(J,K,L),A(J,K,LM))
          END DO
        END DO
      END DO
C$OMP END DO

C...........FLUXES FOR A(J,K,L) IN THE THETA-DIRECTION AT L.
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            RDTHETA=RHF(J)*DTHETA
            FT(J,K,L)=VLI(A(J,K,L),A(J,K,LM),VLOP(J,K,L),VLOP(J,K,LM),
     $           VT(J,K,L),DELT,RDTHETA)*VT(J,K,L)*DZ*DR
          END DO
        END DO
      END DO
C$OMP END DO
      
C...........COMPUTE NEW A BY FLUXING
c     parallelizable once all FT(j,k,l) have been calc'ed
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1
        IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=2,KMAX
          DO J=jstart,JMAX
            AA(J,K,L)=AA(J,K,L)-(FR(J+1,K,L)-FR(J,K,L)+FZ(J,K+1,L)
     &           -FZ(J,K,L)+FT(J,K,LP)-FT(J,K,L))*VINV(J+1)*DELT
          END DO
        END DO
      END DO
C$OMP END DO NOWAIT
      
#if PASSIVE>0
! BEGIN PASSIVE FLUX ARRAYS IF DEFINED       
C######FLUX PASSIVE SCALAR
      do I=1,PAS
C$OMP DO SCHEDULE(STATIC)
      DO 410 L=1,LMAX

C............VAN LEER SLOPE FOR PASSFLUX IN THE R-DIRECTION
         
        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(PASSFLUX(J+1,K,L,I),PASSFLUX(J,K,L,I),    &
     &                  PASSFLUX(J-1,K,L,I))
          END DO

C............FLUXES FOR PASSFLUX(J,K,L) IN THE R-DIRECTION AT R(J)
          FR(2,K,L)=0.D0
          if (jmin.gt.2) then
             jstart=jmin
          else
             jstart=3
          endif 
          DO J=jstart,JMAX1 
            FR(J,K,L)=VLI(PASSFLUX(J,K,L,I),PASSFLUX(J-1,K,L,I),        &
     &           VLOP(J,K,L),
     &           VLOP(J-1,K,L),VR(J,K,L),DELT,DR)*VR(J,K,L)*
     &           (R(J)-VR(J,K,L)*0.5*DELT)*DZ*DTHETA
          END DO
        END DO
         
C............VAN LEER SLOPE FOR RHO IN THE Z-DIRECTION

        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO J=jstart,JMAX1
          VLOP(J,1,L)=SLOPE(PASSFLUX(J,2,L,I),PASSFLUX(J,1,L,I),        &
     &                PASSFLUX(J,3,L,I))
        END DO
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(PASSFLUX(J,K+1,L,I),PASSFLUX(J,K,L,I),    &
     &                 PASSFLUX(J,K-1,L,I))
          END DO
        END DO

C............FLUXES FOR PASSFLUX(J,K,L) IN THE Z-DIRECTION AT Z(K)

        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO J=jstart,JMAX1
          FZ(J,2,L)=0.D0
        END DO
        DO K=3,KMAX1
          DO J=jstart,JMAX1
            FZ(J,K,L)=VLI(PASSFLUX(J,K,L,I),PASSFLUX(J,K-1,L,I),        &
     &           VLOP(J,K,L),
     &           VLOP(J,K-1,L),VZ(J,K,L),DELT,DZ)*VZ(J,K,L)
     &           *RD(J+1)*0.5*DTHETA
          END DO
        END DO
  410 CONTINUE
C$OMP END DO NOWAIT

C............VAN LEER SLOPE passflux RHO IN THE THETA-DIRECTION
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif         
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(PASSFLUX(J,K,LP,I),PASSFLUX(J,K,L,I),     &
     &                  passflux(J,K,LM,I))
          END DO
        END DO
      END DO
C$OMP END DO
      
C............FLUXES FOR passflux(J,K,L) IN THE THETA-DIRECTION AT L
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif         
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            RDTHETA=RHF(J)*DTHETA
            FT(J,K,L)=VLI(passflux(J,K,L,I),passflux(J,K,LM,I),         &
     &            VLOP(J,K,L),VLOP(J,K
     $           ,LM),VT(J,K,L),DELT,RDTHETA)*VT(J,K,L)*DZ*DR
          END DO
        END DO
      END DO
C$OMP END DO
      
C............COMPUTE NEW passflux BY FLUXING
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1
        IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=2,KMAX
          DO J=jstart,JMAX
            ppassflux(J,K,L,I)=ppassflux(J,K,L,I)-(FR(J+1,K,L)-FR(J,K,L)&
     &           +FZ(J,K+1,L)
     $           -FZ(J,K,L)+FT(J,K,LP)-FT(J,K,L))*VINV(J+1)*DELT
          END DO
        END DO
      END DO
C$OMP END DO NOWAIT
      enddo
#endif
      
C######FLUX RHO
C$OMP DO SCHEDULE(STATIC)
      DO 460 L=1,LMAX

C............VAN LEER SLOPE FOR RHO IN THE R-DIRECTION
         
        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(RHO(J+1,K,L),RHO(J,K,L),RHO(J-1,K,L))
          END DO

C............FLUXES FOR RHO(J,K,L) IN THE R-DIRECTION AT R(J)
          FR(2,K,L)=0.D0
          if (jmin.gt.2) then
             jstart=jmin
          else
             jstart=3
          endif 
          DO J=jstart,JMAX1 
            FR(J,K,L)=VLI(RHO(J,K,L),RHO(J-1,K,L),VLOP(J,K,L),
     &           VLOP(J-1,K,L),VR(J,K,L),DELT,DR)*VR(J,K,L)*
     &           (R(J)-VR(J,K,L)*0.5*DELT)*DZ*DTHETA
          END DO
        END DO
         
C............VAN LEER SLOPE FOR RHO IN THE Z-DIRECTION

        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO J=jstart,JMAX1
          VLOP(J,1,L)=SLOPE(RHO(J,2,L),RHO(J,1,L),RHO(J,3,L))
        END DO
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(RHO(J,K+1,L),RHO(J,K,L),RHO(J,K-1,L))
          END DO
        END DO

C............FLUXES FOR RHO(J,K,L) IN THE Z-DIRECTION AT Z(K)

        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO J=jstart,JMAX1
          FZ(J,2,L)=0.D0
        END DO
        DO K=3,KMAX1
          DO J=jstart,JMAX1
            FZ(J,K,L)=VLI(RHO(J,K,L),RHO(J,K-1,L),VLOP(J,K,L),
     &           VLOP(J,K-1,L),VZ(J,K,L),DELT,DZ)*VZ(J,K,L)
     &           *RD(J+1)*0.5*DTHETA
          END DO
        END DO
  460 CONTINUE
C$OMP END DO NOWAIT

C............VAN LEER SLOPE FOR RHO IN THE THETA-DIRECTION
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif         
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(RHO(J,K,LP),RHO(J,K,L),RHO(J,K,LM))
          END DO
        END DO
      END DO
C$OMP END DO
      
C............FLUXES FOR RHO(J,K,L) IN THE THETA-DIRECTION AT L
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif         
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            RDTHETA=RHF(J)*DTHETA
            FT(J,K,L)=VLI(RHO(J,K,L),RHO(J,K,LM),VLOP(J,K,L),VLOP(J,K
     $           ,LM),VT(J,K,L),DELT,RDTHETA)*VT(J,K,L)*DZ*DR
          END DO
        END DO
      END DO
C$OMP END DO
      
C............COMPUTE NEW RHO BY FLUXING
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1
        IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=2,KMAX
          DO J=jstart,JMAX
            RRHO(J,K,L)=RRHO(J,K,L)-(FR(J+1,K,L)-FR(J,K,L)+FZ(J,K+1,L)
     $           -FZ(J,K,L)+FT(J,K,LP)-FT(J,K,L))*VINV(J+1)*DELT
          END DO
        END DO
      END DO
C$OMP END DO NOWAIT
      
C######FLUX EPS
C$OMP DO SCHEDULE(STATIC)
      DO 510 L=1,LMAX
        if (jmin.gt.2) then
           jstart=jmin2
        else
           jstart=1
        endif  
         
C............VAN LEER SLOPE FOR EPS IN THE R-DIRECTION

        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(eps(J+1,K,L),eps(J,K,L),eps(J-1,K,L))
          END DO
            
C............FLUXES FOR EPS(J,K,L) IN THE R-DIRECTION AT R(J)

          if (jmin.gt.2) then
             jstart=jmin
          else
             jstart=3
          endif 
          DO J=jstart,JMAX1
            FR(J,K,L)=VLI(eps(J,K,L),eps(J-1,K,L),VLOP(J,K,L),
     &            VLOP(J-1,K,L),VR(J,K,L),DELT,DR)*VR(J,K,L)*     
     &            (R(J)-VR(J,K,L)*0.5*DELT)*DZ*DTHETA
          END DO
        END DO
         
C............VAN LEER SLOPE FOR EPS IN THE Z-DIRECTION

        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO J=jstart,JMAX1
          VLOP(J,1,L)=SLOPE(eps(J,2,L),eps(J,1,L),eps(J,3,L))
        END DO
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(eps(J,K+1,L),eps(J,K,L),eps(J,K-1,L))
          END DO
        END DO
         
C............FLUXES FOR EPS(J,K,L) IN THE Z-DIRECTION AT Z(K)

        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=3,KMAX1
          DO J=jstart,JMAX1
            FZ(J,K,L)=VLI(eps(J,K,L),eps(J,K-1,L),VLOP(J,K,L),
     &            VLOP(J,K-1,L),VZ(J,K,L),DELT,DZ)*VZ(J,K,L)*
     &            RD(J+1)*0.5*DTHETA
          END DO
        END DO
  510 CONTINUE
C$OMP END DO NOWAIT

C............VAN LEER SLOPE FOR EPS IN THE THETA-DIRECTION
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin1
        else
           jstart=2
        endif 
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            VLOP(J,K,L)=SLOPE(eps(J,K,LP),eps(J,K,L),eps(J,K,LM))
          END DO
        END DO
      END DO
C$OMP END DO

C............FLUXES FOR EPS(J,K,L) IN THE THETA-DIRECTION AT L
C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LM=L-1
        IF(L.EQ.1) LM=LMAX
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=2,KMAX1
          DO J=jstart,JMAX1
            RDTHETA=RHF(J)*DTHETA
            FT(J,K,L)=VLI(eps(J,K,L),eps(J,K,LM),VLOP(J,K,L),VLOP(J
     $           ,K,LM),VT(J,K,L),DELT,RDTHETA)*VT(J,K,L)*DZ*DR
          END DO
        END DO
      END DO
C$OMP END DO


C............COMPUTE NEW EPS BY FLUXING

C$OMP DO SCHEDULE(STATIC)
      DO L=1,LMAX
        LP=L+1
        IF(L.EQ.LMAX) LP=1
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=2
        endif 
        DO K=2,KMAX
          DO J=jstart,JMAX
            IF(RRHO(J,K,L).LE.RHOLMT) THEN 
              eeps(J,K,L)=epslmt
            ELSE 
              eeps(J,K,L)=eeps(J,K,L)-(FR(J+1,K,L)-FR(J,K,L)+FZ(J
     $             ,K+1,L)-FZ(J,K,L)+FT(J,K,LP)-FT(J,K,L))*VINV(J+1)
     $             *DELT
            END IF
          END DO
        END DO
      END DO
C$OMP END DO NOWAIT

C-----------------------------------------------------------------------
C######FINISH FLUXING ALL THE QUANTITIES, NOW CLEAN UP.

      
C.... SET QUANTITIES AT THE GRID EDGE.

C$OMP DO SCHEDULE(STATIC) REDUCTION(+:tmassadd)
      DO L=1,LMAX
        DO K=2,KMAX2
          DO J=JMAX1,JMAX2
            SS(J,K,L)   = ABS(SS(J-1,K,L))
            TT(J,K,L)   = TT(J-1,K,L)
            AA(J,K,L)   = AA(J-1,K,L)
            RRHO(J,K,L) = RRHO(J-1,K,L)
            EEPS(J,K,L) = EEPS(J-1,K,L)
#if PASSIVE>0
            do I=1,PAS
              PPASSFLUX(J,K,L,I)=PPASSFLUX(J-1,K,L,I)
            enddo
#endif
          END DO
        END DO
        if (jmin.gt.2) then
           jstart=jmin
        else
           jstart=1
        endif 
        DO K=KMAX1,KMAX2
          DO J=jstart,JMAX2
            SS(J,K,L)   = SS(J,K-1,L)
            TT(J,K,L)   = ABS(TT(J,K-1,L))
            AA(J,K,L)   = AA(J,K-1,L)
            RRHO(J,K,L) = RRHO(J,K-1,L)
            EEPS(J,K,L) = EEPS(J,K-1,L)
#if PASSIVE>0
            do I=1,PAS
              PPASSFLUX(J,K,L,I)=PPASSFLUX(J,K-1,L,I)
            enddo
#endif
          END DO
        END DO

C.....THE LOWEST DENSITY
        DO K=2,KMAX2
          DO J=2,JMAX2
             if ((rrho(j,k,l).lt.rholmt).and.(j.le.jmax).and.(k.le
     &            .kmax).and.(j.ge.jmin).and.(k.ge.2)) then 
                tmassadd=tmassadd+(rholmt-rrho(j,k,l))
     &               *dtheta*zof3n*(r(j+1)**2-r(j)**2)
             endif   
            RRHO(J,K,L) = MAX(RRHO(J,K,L),RHOLMT)
            EEPS(J,K,L) = MAX(EEPS(J,K,L),EPSLMT)
#if PASSIVE>0
            do I=1,PAS
              PPASSFLUX(J,K,L,I)=MAX(PPASSFLUX(J,K,L,I),RHOLMT)
            enddo
#endif
          END DO
        END DO
      END DO
C$OMP END DO
      

C.....Set quantities around the z- axis.
C.....For inner outflow boundary, extrapolate values of ghost zones. 
C.....acm2000

      if (jmin.gt.2) then

C$OMP DO SCHEDULE(STATIC)
         do l=1,lmax
            do  K=2,KMAX2
               do j=1,jmin2-1
                  TT(j,K,l)=0.0
                  SS(j,K,l)=0.0
                  RRHO(j,K,l)=rholmt
                  EEPS(j,K,l)=epslmt
                  AA(j,K,l)=0.0
                  U(j,K,l)=0.d0
                  W(j,k,l)=0.0
                  OMEGA(j,k,l)=0.0
               enddo
               
               SS(jmin,k,l)=2*ss(jmin+1,k,l)-ss(jmin+2,k,l)
               if (abs(ss(jmin,k,l)).gt.abs(ss(jmin+1,k,l))) then
                  ss(jmin,k,l)=ss(jmin+1,k,l)
               endif
               ss(jmin,k,l)=-abs(ss(jmin,k,l))
               
               do j=jmin1,jmin2,-1
                  
                  RRHO(j,k,l)=-((rrho(jmin+1,k,l)-rrho(jmin,k,l))
     &                 *(rhf(jmin)-rhf(j))/dr) + rrho(jmin,k,l)
                  if (RRHO(j,k,l).gt.RRHO(j+1,k,l)) then 
                     RRHO(j,k,l)=RRHO(j+1,k,l)
                  endif
                  if (RRHO(j,k,l).lt.rholmt) then 
                     RRHO(j,k,l)=rholmt
                  endif
                  
                  EEPS(j,k,l)=-((eeps(jmin+1,k,l)-eeps(jmin,k,l))
     &                 *(rhf(jmin)-rhf(j))/dr) + eeps(jmin,k,l)
                  if (EEPS(j,k,l).gt.EEPS(j+1,k,l)) then 
                     EEPS(j,k,l)=EEPS(j+1,k,l)
                  endif
                  if ((EEPS(j,k,l).lt.epslmt).or.(RRHO(j,k,l).le.rholmt)
     &                 )then 
                     EEPS(j,k,l)=epslmt
                  endif
                  
                  AA(j,k,l)=-((aa(jmin+1,k,l)-aa(jmin,k,l))*(rhf(jmin)
     &                 -rhf(j))/dr) + aa(jmin,k,l)
                  if (AA(j,k,l).gt.AA(j+1,k,l)) then 
                     AA(j,k,l)=AA(j+1,k,l)
                  endif
                  if ((AA(j,k,l).lt.0.0).or.(RRHO(j,k,l).le.rholmt))
     &                 then 
                     AA(j,k,l)=0.0
                  endif
                  
                  TT(j,k,l)=-((tt(jmin+1,k,l)-tt(jmin,k,l))*(rhf(jmin)
     &                 -rhf(j))/dr) + tt(jmin,k,l)
                  if (abs(TT(j,k,l)).gt.abs(TT(j+1,k,l))) then 
                     TT(j,k,l)=TT(j+1,k,l)
                  endif
                  if (RRHO(j,k,l).le.rholmt) then 
                     TT(j,k,l)=0.0
                  endif
                  
                  SS(j,k,l)=-((ss(jmin+2,k,l)-ss(jmin+1,k,l))*(r(jmin+1)
     &                 -r(j))/dr)+ ss(jmin+1,k,l)
                  if (abs(ss(j,k,l)).gt.abs(ss(j+1,k,l))) then 
                     ss(j,k,l)=ss(j+1,k,l) 
                  endif
                  SS(j,k,l)=-abs(ss(j,k,l))
                  
               enddo   
            enddo
         enddo
C$OMP END DO NOWAIT

      else   

C$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            LP=1+MOD(L-1+LMAX/2,LMAX)
            DO K=2,KMAX2
               SS(1,K,L)    = -SS(3,K,LP)
               TT(1,K,L)    = TT(2,K,LP)
               AA(1,K,L)    = AA(2,K,LP)
               RRHO(1,K,L)  = RRHO(2,K,LP)
               EEPS(1,K,L)  = EEPS(2,K,LP)
#if PASSIVE>0
            do I=1,PAS
              PPASSFLUX(1,K,L,I)=PPASSFLUX(2,K,LP,I)
            enddo
#endif
            END DO
         END DO
C$OMP END  DO
      
      
C.....Set S on the z-axis.
C.....Set velocities of neighbouring cells.
C$OMP DO SCHEDULE(STATIC)
         DO L=1,LMAX
            DO K=2,KMAX1
               VR(3,K,L)=2*SS(3,K,L)/(RRHO(3,K,L)+RRHO(2,K,L))
               VT(3,K,L)=AA(2,K,L)/(RRHO(2,K,L)*RHF(2))
            END DO
            DO K=2,KMAX2
               U(2,K,L)=0.d0
               SS(2,K,L)=0.d0
            END DO
         ENDDO   
C$OMP END DO NOWAIT
         
      endif
        
C.....Set quantities below equatorial plane.

C$OMP DO SCHEDULE(STATIC)
        DO L=1,LMAX
           if (jmin.gt.2) then
              jstart=jmin2
           else
              jstart=1
           endif 
           DO J=jstart,JMAX2
              SS(J,1,L)=SS(J,2,L)
              W(J,2,L)=0.D0
              TT(J,2,l)=0.D0
              TT(J,1,L)=-TT(J,3,L)
              AA(J,1,L)=AA(J,2,L)
              RRHO(J,1,L)=RRHO(J,2,L)
              EEPS(J,1,L)=EEPS(J,2,L)
#if PASSIVE>0
            do I=1,PAS
              PPASSFLUX(J,1,L,I)=PPASSFLUX(J,2,L,I)
            enddo
#endif
           END DO
        END DO
C$OMP END DO NOWAIT

      RETURN
      END
      
