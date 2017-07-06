!=============================================================================
! Name        : pot3MPI.F90
! Version     : 1.0
! Description : Solves the Poisson equation in two dimensions.
!
! Method      :
!
!   Driver for the poisson solver, performs fourier transform, then
!   solves via multigrid method. 
!==============================================================================


      SUBROUTINE POT3MPI(NPOINT,IPRINT,jmax,kmax,lmax,rho,dx,phi)
      USE defines_mod
      use io, only : writeData
      use MPI_F08  , only : MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD

      implicit none

      real, intent(INOUT):: rho(-1:jmax+1,-1:kmax+1,1:lmax)
      real, intent(inout):: phi(-1:jmax+1,-1:kmax+1,1:lmax)
      integer, intent(in):: jmax,kmax,lmax,NPOINT,IPRINT
      real, allocatable  :: A1(:),B1(:),coef(:,:,:,:)
      real, allocatable  :: Y(:,:)
      integer            :: N,j,k,l, rank, numRanks
      real, intent(in)   :: dx
      real               :: sf1,sf2,x1,philag,pi
!      integer jmax,kmax,lmax
!     COMMON /COEFS/COEF(pot3JMAX2,pot3KMAX2,LMAX2,2)

call MPI_Comm_size(MPI_COMM_WORLD, numRanks)
call MPI_Comm_rank(MPI_COMM_WORLD, rank)


! A1 and B1  should be dimensioned lmax.  They're used in fft.
allocate(A1(1:lmax),B1(1:lmax))
allocate(coef(jmax,kmax,lmax,2),Y(-1:jmax+1,-1:kmax+1))
! These next arrays are used in blktri. The size of wfw(max) is given by
! max. ge. (2*(kmax3+2)*(log2(kmax3+1)-1) + 6*jmax3+2).
print *,"rank",rank,"has kmax=",kmax
print *, """I'm inside pot3"", said rank", rank
      N=LMAX/2
      IF(LMAX.EQ.1) N=1
      SF1=0.5/FLOAT(N)
      SF2=0.5

      DO L=1,LMAX                                                     
        A1(L)=0.0
        B1(L)=0.0
      ENDDO
pi = acos(-1.d0)
! For convenience, put all rho's into phi array -- leave boundary phi's
!- alone.
      DO L=1,LMAX                                                     

        DO K=-1,KMAX+1
          DO J=-1,JMAX+1
        if (rank==numRanks-1) then
           if((k.ne.kmax).and.(j.ne.jmax)) then
            PHI(J,K,L)=4.0d0*pi*RHO(J,K,L) ! not boundaries
           else
            PHI(J,K,L)=RHO(J,K,L) ! boundaries
           end if
        else
           if(j.ne.jmax) then
            PHI(J,K,L)=4.0d0*pi*RHO(J,K,L) ! not boundaries
           else
            PHI(J,K,L)=RHO(J,K,L) ! boundaries
           end if
        end if 
          ENDDO
        ENDDO

      ENDDO

!call writeData(-1,jmax,kmax,rho(:,:,1),'rho')
!call writeData(-1,jmax,kmax,phi(:,:,1),'phi0')
!-----------------------------------------------------------------------
!
!  If lmax>1 then for all j,k, calculate fourier transform of rho's
!  and boundary phi's.

      IF(LMAX.gt.1) THEN
         
        DO K=-1,KMAX+1
          DO J=-1,JMAX
            DO L=1,N
              A1(L) = PHI(J,K,2*L-1)
              B1(L) = PHI(J,K,2*L)
            ENDDO 
            CALL FFT(A1,B1,N,N,N,1)
            CALL REALTR(A1,B1,N,1)
!     Cosine coefficients are in a1? Sine coef's are in b1.
!     Put cosine coef's in phi(j,k,1:n+1)?  Put sine coef's in      
!     phi(j,k,n+2:lmax).  Normalize computed values with 0.5/n.     
!     Compute amplitudes of modes.
!            X1=ABS(A1(1))
!            IF(X1.NE.0.0) THEN
!              DO L=2,N
!                COEF(J,K,L-1,1) = 2.0*SQRT( A1(L)**2 + B1(L)**2 )/X1
!              ENDDO
!              COEF(J,K,N,1) = ABS(A1(N+1))/X1
!            ENDIF

!     Compute phase angles of modes. The following code simply calculates
!     philag=arctan(b/a) in the range -pi/2 <= philag < 3*pi/2. It sets
!     philag=4*pi if a and b are both 0.
!            DO L=2,N
!              IF(A1(L).GT.0.0) THEN
!                PHILAG = ATAN(B1(L)/A1(L))
!              ELSE IF(A1(L).LT.0.0) THEN
!                PHILAG = ATAN(B1(L)/A1(L))+PI
!              ELSE
!                IF(B1(L).GT.0.0) THEN
!                  PHILAG = 0.5*PI
!                ELSE IF(B1(L).LT.0.0) THEN
!                  PHILAG = -0.5*PI
!                ELSE
!                  PHILAG = 4.0*PI
!                ENDIF
!              ENDIF
!              COEF(J,K,L-1,2)=PHILAG*180./PI
           ! ENDDO
!            COEF(J,K,N,2)=A1(1)

!     Copy fourier transform back into phi array.
            DO L=1,N+1
              PHI(J,K,L)=A1(L)*SF1
              IF(L.LT.N) THEN
                PHI(J,K,L+N+1)=B1(L+1)*SF1
              ENDIF
!             write(100+L,*) j,k,L, A1(L)*SF1, B1(L)*SF1
            ENDDO
            DO L=N+1,LMAX
              A1(L) = 0.0
              B1(L) = 0.0
            ENDDO
          ENDDO
        ENDDO                                                           
         
      ENDIF

!    call writeData(-1,jmax,kmax,phi(:,:,1),'fftrho1')
!    call writeData(-1,jmax,kmax,phi(:,:,2),'fftrho2')
!-----------------------------------------------------------------------
!  Now solve the transformed problem. The fourier transform in the L
!  direction has turned the 3-D Poisson problem into a set of LMAX
!  uncoupled 2-D problems which may be solved concurrently.
!     
! Call the mutigrid potential solver
do l = 1,lmax
 Call poissMultigrid(phi(:,:,l),Y,jmax,kmax,lmax,l,dx,dx)
!     Solution on 2-D grid is complete.
!     Put transformed phi's from y into phi array.
 PHI(:,:,l)=Y(:,:)
end do

!-----------------------------------------------------------------------
!
!     All transformed phi's have been calculated.  Now obtain real
!     phi's by inverse fourier transform.

      IF(LMAX.NE.1) THEN
        DO K=1,KMAX                                                  
          DO J=1,JMAX
            DO L=1,LMAX
              A1(L)=0.0
              B1(L)=0.0
            ENDDO                                                           
            DO L=1,N+1
              A1(L)=PHI(J,K,L)
              IF(L.LT.N) THEN
                B1(L+1)=PHI(J,K,L+N+1)
              ENDIF
            ENDDO                                                 
            CALL REALTR(A1,B1,N,-1)
            CALL FFT(A1,B1,N,N,N,-1)
!     Now real phi's are in first n places of both a1 and b1.  With
!     increasing l, phi's alternate between a1 and b1.  Put phi's into
!     array phi, normalizing them by 0.5.
!     
            DO L=1,N
              PHI(J,K,2*L-1) = SF2*A1(L)
              PHI(J,K,2*L)   = SF2*B1(L)
            ENDDO
          ENDDO
        ENDDO
      ENDIF                                                        

!     Calculation of phi's is now finished. 

!    call writeData(0,jmax,kmax,phi(:,:,2),'phir2')
!    call writeData(0,jmax,kmax,phi(:,:,5),'phir5')
!    call writeData(0,jmax,kmax,phi(:,:,16),'phir16')
      DO L=1,N                                                       
        DO K=1,KMAX 
          PHI(0,K,L)   = PHI(1,K,N+L)
          PHI(0,K,N+L) = PHI(1,K,L)
        ENDDO                                                           
       if(rank==0) then
        DO J=1,JMAX+1                                                  
          PHI(J,0,L)   = PHI(J,1,L)
          PHI(J,0,N+L) = PHI(J,1,N+L)
        ENDDO
       end if
      ENDDO                                                           

      print *, "calling zaxphi"
      CALL ZAXPHI(NPOINT,IPRINT)

    call writeData(1,jmax,kmax,phi(:,:,1),'phir1')

      RETURN
      END                                                                


!*******************************************************************************

      SUBROUTINE ZAXPHI(NPOINT,IPRINT)
      USE defines_mod
      IMPLICIT real*8 (a-h,o-z)      

      INCLUDE "hydroparam.h"
      INCLUDE "globals.h"

      COMMON /INSIDE/TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT
      COMMON /COEFS/COEF(POT3JMAX2,POT3KMAX2,LMAX2,2)

      DIMENSION X(10),PHI2(10),PH(LMAX),PHAC(KMAX2)
!     PHAC IS DIMENSIONED KMAX2? PH IS DIMENSIONED LMAX/2.
! But it would not hurt to dimension PH to LMAX.

 

! With X and PHI2 dimensioned 10, 10-point interpolation is the
! maximum you can do.
      IF(NPOINT.GT.10) THEN
        WRITE(3,10000) NPOINT
10000 FORMAT(///,' YOU ARE LIMITED TO A 10-POINT INTERPOLATION HERE IN', &
     &       ' SUBROUTINE ZAXPHI.',/,' BUT YOU PUT NPOINT =',I4,///)
        NPOINT=10
      ENDIF
      XPOINT=0.0
      LM=LMAX/2
      N=NPOINT
      ICHK=MOD(N,2)
      N=N-ICHK
      NUM=N/2
      ISPECL=NUM+NUM+1
      JSP=NUM+1
      DO J=2,JSP
        I1=NUM+(J-1)
        I2=NUM-(J-2)
        X(I1)=RHF(J)
        X(I2)=-RHF(J)
      ENDDO
      IF(ICHK.EQ.1) X(ISPECL)=RHF(NUM+2)

      DO 50 K=2,KMAX1
        PHMAX=0.0
!       For given k, calculate phi on z-axis at all angles.  store
!       results in PH(L) array and in PHI(JMAX2,K,L).
        DO 20 L=1,LM
          LP=L+LM
          DO J=2,JSP
            I1=NUM+(J-1)
            I2=NUM-(J-2)
            PHI2(I1)=PHI(J,K,L)
            PHI2(I2)=PHI(J,K,LP)
          ENDDO
          IF(ICHK.EQ.1) PHI2(ISPECL)=PHI(NUM+2,K,L)
          IM=NPOINT-1
          DO I=1,IM
            P1=PHI2(I)
            XINV=1.0/X(I)
            XR=XPOINT*XINV
            IST=I+1
            DO J=IST,NPOINT
              XRATIO=X(J)*XINV
              P2=PHI2(J)
              PHI2(J)=(P1*(XRATIO-XR)+P2*(XR-1.0))/(XRATIO-1.0)
            ENDDO
          ENDDO
          PH(L)=PHI2(NPOINT)
          IF(ABS(PH(L)).GT.ABS(PHMAX))PHMAX=PH(L)
          PHI(JMAX2,K,L)=PH(L)
          PHI(JMAX2,K,LP)=PH(L)
   20   CONTINUE
!       At this k, find maximum deviation in ph(l)'s? Put result in phac(k).
        ERR=0.0
        DO L=1,LM
          ER=1.0-PH(L)/PHMAX
          IF(ABS(ER).GT.ABS(ERR))ERR=ER
        ENDDO
        PHAC(K)=ERR
   50 CONTINUE

! Find largest of all deviations in PHI's on z-axis and put value in PHICHK
! and k-value in KLOCAT.
      ERR=0.0
      DO K=2,KMAX1
        XX=ABS(PHAC(K))
        IF(XX.GT.ABS(ERR)) THEN
          ERR=PHAC(K)
          KLOCAT=K
        ENDIF
      ENDDO
      PHICHK=ERR
 
! Write the amplitude and phase to the results file.
      IF(IPRINT.EQ.1)THEN
        WRITE(3,10100) 'PHAC(1:KMAX1):'
        WRITE(3,10040)(PHAC(K),K=2,KMAX1)
        WRITE(3,10100) 'Fourier amplitudes and phases:'
10040   FORMAT(1P11E11.3)
10100   FORMAT(a)
        IF(LMAX.GT.8) THEN
          MR=1
          MP=8
!         LQ=LM/8
!         DO K=1,LQ
            WRITE(3,10110)(J,((COEF(J,2,M,I),M=MR,MP),I=1,2),J=2,JMAX)
            MR=MR+8
            MP=MP+8
!         ENDDO
        END IF
      END IF
10110 FORMAT(20X,I5,' COEFS=',1P8E12.3,/,25X,' PHASE=',0P8F12.2)


      RETURN
      END
