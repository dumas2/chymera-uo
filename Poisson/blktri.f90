SUBROUTINE BLKTRI (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,IERROR,W)

!*******************************************************************************
!***BEGIN PROLOGUE  BLKTRI
!***PURPOSE  Solve a block tridiagonal system of linear equations
!            (usually resulting from the discretization of separable
!            two-dimensional elliptic equations).
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B4B
!***TYPE      SINGLE PRECISION (BLKTRI-S, CBLKTR-C)
!***KEYWORDS  ELLIPTIC PDE, FISHPACK, TRIDIAGONAL LINEAR SYSTEM
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine BLKTRI Solves a System of Linear Equations of the Form
!
!          AN(J)*X(I,J-1) + AM(I)*X(I-1,J) + (BN(J)+BM(I))*X(I,J)
!
!          + CN(J)*X(I,J+1) + CM(I)*X(I+1,J) = Y(I,J)
!
!               for I = 1,2,...,M  and  J = 1,2,...,N.
!
!     I+1 and I-1 are evaluated modulo M and J+1 and J-1 modulo N, i.e.,
!
!          X(I,0) = X(I,N),  X(I,N+1) = X(I,1),
!          X(0,J) = X(M,J),  X(M+1,J) = X(1,J).
!
!     These equations usually result from the discretization of
!     separable elliptic equations.  Boundary conditions may be
!     Dirichlet, Neumann, or Periodic.
!
!
!     * * * * * * * * * *     ON INPUT     * * * * * * * * * *
!
!     IFLG
!       = 0  Initialization only.  Certain quantities that depend on NP,
!            N, AN, BN, and CN are computed and stored in the work
!            array  W.
!       = 1  The quantities that were computed in the initialization are
!            used to obtain the solution X(I,J).
!
!       NOTE   A call with IFLG=0 takes approximately one half the time
!              as a call with IFLG = 1  .  However, the
!              initialization does not have to be repeated unless NP, N,
!              AN, BN, or CN change.
!
!     NP
!       = 0  If AN(1) and CN(N) are not zero, which corresponds to
!            periodic boundary conditions.
!       = 1  If AN(1) and CN(N) are zero.
!
!     N
!       The number of unknowns in the J-direction. N must be greater
!       than 4. The operation count is proportional to MNlog2(N), hence
!       N should be selected less than or equal to M.
!
!     AN,BN,CN
!       One-dimensional arrays of length N that specify the coefficients
!       in the linear equations given above.
!
!     MP
!       = 0  If AM(1) and CM(M) are not zero, which corresponds to
!            periodic boundary conditions.
!       = 1  If AM(1) = CM(M) = 0  .
!
!     M
!       The number of unknowns in the I-direction. M must be greater
!       than 4.
!
!     AM,BM,CM
!       One-dimensional arrays of length M that specify the coefficients
!       in the linear equations given above.
!
!     IDIMY
!       The row (or first) dimension of the two-dimensional array Y as
!       it appears in the program calling BLKTRI.  This parameter is
!       used to specify the variable dimension of Y.  IDIMY must be at
!       least M.
!
!     Y
!       A two-dimensional array that specifies the values of the right
!       side of the linear system of equations given above.  Y must be
!       dimensioned at least M*N.
!
!     W
!       A one-dimensional array that must be provided by the user for
!       work space.
!             If NP=1 define K=INT(log2(N))+1 and set L=2**(K+1) then
!                     W must have dimension (K-2)*L+K+5+MAX(2N,6M)
!
!             If NP=0 define K=INT(log2(N-1))+1 and set L=2**(K+1) then
!                     W must have dimension (K-2)*L+K+5+2N+MAX(2N,6M)
!
!       **IMPORTANT** For purposes of checking, the required dimension
!                     of W is computed by BLKTRI and stored in W(1)
!                     in floating point format.
!
!     * * * * * * * * * *     On Output     * * * * * * * * * *
!
!     Y
!       Contains the solution X.
!
!     IERROR
!       An error flag that indicates invalid input parameters.  Except
!       for number zero, a solution is not attempted.
!
!       = 0  No error.
!       = 1  M is less than 5.
!       = 2  N is less than 5.
!       = 3  IDIMY is less than M.
!       = 4  BLKTRI failed while computing results that depend on the
!            coefficient arrays AN, BN, CN.  Check these arrays.
!       = 5  AN(J)*CN(J-1) is less than 0 for some J. Possible reasons
!            for this condition are
!            1. The arrays AN and CN are not correct.
!            2. Too large a grid spacing was used in the discretization
!               of the elliptic equation.
!            3. The linear equations resulted from a partial
!               differential equation which was not elliptic.
!
!     W
!       Contains intermediate values that must not be destroyed if
!       BLKTRI will be called again with IFLG=1.  W(1) contains the
!       number of locations required by W in floating point format.
!
! *Long Description:
!
!     * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!     Dimension of   AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N)
!     Arguments      W(See argument list)
!
!     Latest         June 1979
!     Revision
!
!     Required       BLKTRI,BLKTRI,PROD,PRODP,CPROD,CPRODP,COMPB,INDXA,
!     Subprograms    INDXB,INDXC,PPADD,PSGF,PPSGF,PPSPF,BSRH,TEVLS,
!                    R1MACH
!
!     Special        The Algorithm may fail if ABS(BM(I)+BN(J)) is less
!     Conditions     than ABS(AM(I))+ABS(AN(J))+ABS(CM(I))+ABS(CN(J))
!                    for some I and J. The Algorithm will also fail if
!                    AN(J)*CN(J-1) is less than zero for some J.
!                    See the description of the output parameter IERROR.
!
!     Common         CBLKT
!     Blocks
!
!     I/O            None
!
!     Precision      Single
!
!     Specialist     Paul Swarztrauber
!
!     Language       FORTRAN
!
!     History        Version 1 September 1973
!                    Version 2 April     1976
!                    Version 3 June      1979
!
!     Algorithm      Generalized Cyclic Reduction (See Reference below)
!
!     Space
!     Required       Control Data 7600
!
!     Portability    American National Standards Institute Fortran.
!                    The machine accuracy is set using function R1MACH.
!
!     Required       None
!     Resident
!     Routines
!
!     References     Swarztrauber,P. and R. Sweet, 'Efficient FORTRAN
!                    Subprograms For The Solution Of Elliptic Equations'
!                    NCAR TN/IA-109, July, 1975, 138 PP.
!
!                    Swarztrauber P. ,'A Direct Method For The Discrete
!                    Solution Of Separable Elliptic Equations', S.I.A.M.
!                    J. Numer. Anal.,11(1974) PP. 1136-1150.
!
!***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
!                 subprograms for the solution of elliptic equations,
!                 NCAR TN/IA-109, July 1975, 138 pp.
!               P. N. Swarztrauber, A direct method for the discrete
!                 solution of separable elliptic equations, SIAM Journal
!                 on Numerical Analysis 11, (1974), pp. 1136-1150.
!***ROUTINES CALLED  BLKTR1, COMPB, CPROD, CPRODP, PROD, PRODP
!***COMMON BLOCKS    CBLKT
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BLKTRI
!*******************************************************************************

      IMPLICIT real*8 (a-h,o-z)                                         

      DIMENSION AN(*),BN(*),CN(*),AM(*),BM(*),CM(*),Y(IDIMY,*),W(*)     
      COMMON /CBLKT/ NPP,Kx,EPSLON,CNV,Lx,NCMPLX,IK,IZ,DUM(1)
      COMMON /BLOKJ1/IWCN,IW1,IW2,IW3,IWD,IWW,IWU
!$OMP THREADPRIVATE(/CBLKT/,/BLOKJ1/)

! Initialization.
      IF(IFLG.EQ.0) THEN
        IERROR = 1                                                        
        IF(M.LT.2) RETURN
        NH = N                                                            
        NPP = NP                                                          
        IF(NPP.NE.0) NH = NH+1
        IK = 2                                                            
        Kx = 0                                                             
   25   IK = IK+IK                                                        
        Kx = Kx+1                                                           
        IF(NH/IK) 35,35,25                                                 
   35   IERROR = 2                                                        
        IWCN = 2*(N+1)*(Kx-1)-N+3                                          
        IW1 = IWCN+N                                                      
        IWAH = IW1                                                        
        IWBH = IWAH+N                                                     
        IF (Kx.LT.2) RETURN
        NCK = 2**Kx-1                                                      
        IF (N.NE.NCK) RETURN
        IERROR = 5                                                        
        IF (IDIMY.LT.M) RETURN
        IERROR = 0                                                        
        IW2 = IW1+M                                                       
        IW3 = IW2+M                                                       
        IWD = IW3+M                                                       
        IWW = IWD+M                                                       
        IWU = IWW+M                                                       
        W(IWCN) = CN(N)                                                   
        I = IWCN                                                          
        DO  65 J=2,N                                                      
           I = I+1                                                        
           W(I) = CN(J-1)                                                 
   65   CONTINUE                                                          
        CALL COMPB(N,IERROR,AN,BN,W(IWCN),W,W(IWAH),W(IWBH))             

! Call the cyclic reduction solver.
      ELSE
        CALL BLKTR1(N,AN,BN,W(IWCN),M,AM,BM,CM,IDIMY,Y,W,W(IW1),W(IW2),  &
     &              W(IW3),W(IWD),W(IWW),W(IWU))                         

      ENDIF                                                          

      RETURN                                                            
END SUBROUTINE BLKTRI

!*******************************************************************************

SUBROUTINE BLKTR1 (N,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,B,W1,W2,W3,WD,WW,WU)
      IMPLICIT REAL*8 (A-H,O-Z)                                         

      DIMENSION AN(*),BN(*),CN(*),AM(*),BM(*),CM(*),B(*),W1(*),W2(*),   &
     &          W3(*),WD(*),WW(*),WU(*),Y(IDIMY,*)

      COMMON /CBLKT/ NPP,Kx,EPSLON,CNV,Lx,NCMPLX,IK,IZ,DUM(1)
!$OMP THREADPRIVATE(/CBLKT/)

      NH = 2**Kx                                                         
      IH1 = NH+NH                                                       
      IF (NPP.EQ.0) THEN
        KDO = Kx                                                           
      ELSE
        KDO = Kx-1                                                         
      ENDIF
      DO  40 Lx=1,KDO                                                    
        IR = Lx-1                                                       
        IZ = 2**IR                                                     
        ISGN = (-1)**IR                                                
        MSGN = -ISGN                                                   
        IH2 = 2**(Kx-IR+1)                                              
        LM = (IR-2)*IH1+IH2+IH2+1                                      
        LZ = (IR-1)*IH1+IH2+1                                          
        IIM = IZ-1                                                     
        IIZ = IIM+IIM+1                                                
        JM1 = IIM+IIM+LM                                               
        I1 = IZ+IZ                                                     
        CALL PRDCT (IIZ,B(LZ),IIM,B(LM),IIM,B(JM1),0,DUM,Y(1,IZ),W3,M, &
     &              AM,BM,CM,WD,WW,WU,ISGN)                            
        IF = 2**(Kx-IR-1)                                               
        DO  35 JJ=1,IF                                                 
          I = JJ*I1                                                   
          I6 = I                                                      
          I7 = I-IZ                                                   
          I9 = I+IZ                                                   
          J2 = JJ+JJ                                                  
          J4 = J2+J2                                                  
          JM1 = (J4-2)*IIM+LM                                         
          JP1 = J4*IIM+LM                                             
          JP2 = J2*IIZ+LZ                                             
          JP3 = (J4+2)*IIM+LM                                         
          IF (JJ.LT.IF)  GOTO 25
          IF (NPP.NE.0) GOTO 35
   20     JP1 = LM                                                    
          JP2 = LZ                                                    
          JP3 = IIM+IIM+LM                                            
          I6 = 0                                                      
          I9 = IZ                                                     
   25     CALL PRDCT(IIM,B(JM1), 0,  DUM,   0,  DUM, IZ,AN(I7+1),W3,W1,&
     &                M,AM,BM,CM,WD,WW,WU,MSGN)                            
          CALL PRDCT(IIZ,B(JP2),IIM,B(JP1),IIM,B(JP3),0,DUM,Y(1,I9),W3,&
     &                M,AM,BM,CM,WD,WW,WU,ISGN)                    
          CALL PRDCT(IIM,B(JP1), 0,  DUM,   0,  DUM, IZ,CN(I6+1),W3,W2,&
     &                M,AM,BM,CM,WD,WW,WU,MSGN)                            
          DO  J=1,M                                                
             Y(J,I) = W1(J)+W2(J)-Y(J,I)                              
          ENDDO                                                    
   35   CONTINUE
   40 CONTINUE

   60 DO 130 LL=1,Kx                                                     
        Lx = Kx-LL+1                                                     
        IR = Lx-1                                                       
        ISGN = (-1)**IR                                                
        MSGN = -ISGN                                                   
        IZ = 2**IR                                                     
        IIM = IZ-1                                                     
        IIZ = IIM+IIM+1                                                
        IH2 = 2**(Kx-IR+1)                                              
        LM = (IR-2)*IH1+IH2+IH2+1                                      
        LZ = (IR-1)*IH1+IH2+1                                          
        IF = 2**(Kx-IR)-1                                               
        DO 125 JJ=1,IF,2                                               
          I = JJ*IZ                                                   
          I5 = I-IZ                                                   
          I6 = I+IZ                                                   
          I7 = I5                                                     
          J2 = JJ+JJ                                                  
          JM1 = (J2-2)*IIM+LM                                         
          JZ = (JJ-1)*IIZ+LZ                                          
          JP1 = J2*IIM+LM                                             
!         write(6,9000) ll,jj,if,IIZ,JZ,IIM,JM1,JP1,I                  !DKB debug
!9000     format('LL=',i3,'  JJ=',i3,'  IF=',i3,'  IIZ=',i3,           !DKB debug
!    &         '  JZ=',i3,'  IIM=',i3,'  JM1=',i3,'  JP1=',i3,         !DKB debug
!    &         '  I=',i3)                                              !DKB debug
          IF (JJ.GT.1)  GOTO 85                                       
          IF (NPP.NE.0) GOTO 75                                        
          I7 = N                                                      
          GO TO  85                                                   
   75     DO J=1,M                                                
             W1(J) = 0.                                               
          ENDDO
          GO TO  90                                                   
   85     CALL PRDCT(IIM,B(JM1),0,DUM,0,DUM,IZ,AN(I5+1),Y(1,I7),W1,  &
     &                M,AM,BM,CM,WD,WW,WU,MSGN)                       
   90     IF (JJ.LT.IF .OR. NPP.EQ.0) GOTO 110
            DO J=1,M                                                
               W2(J) = 0.                                               
            ENDDO
          GO TO 115                                                   
  110       CALL PRDCT(IIM,B(JP1),0,DUM,0,DUM,IZ,CN(I+1),Y(1,I6),W2,&
     &                  M,AM,BM,CM,WD,WW,WU,MSGN)                         
  115     DO J=1,M                                                
             W1(J) = Y(J,I)-W1(J)-W2(J)                               
          ENDDO
          CALL PRDCT(IIZ,B(JZ),IIM,B(JM1),IIM,B(JP1),0,DUM,W1,Y(1,I),&
     &                M,AM,BM,CM,WD,WW,WU,ISGN)                
  125   CONTINUE                                                       
  130 CONTINUE

      RETURN                                                            
END SUBROUTINE BLKTR1

!*******************************************************************************
      
SUBROUTINE PRDCT(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,U,IS)
      IMPLICIT real*8 (a-h,o-z)                                         

      DIMENSION A(*),B(*),C(*),X(*),Y(*),D(*),W(*),BD(*),BM1(*),BM2(*), &
     &          AA(*),U(*)

      IF (ND.LE.0)  GOTO 10
      IF (IS.LE.0)  GOTO 20
   10 DO J=1,M                                                      
        W(J) = X(J)                                                    
        Y(J) = X(J)                                                         
      ENDDO
      GO TO  30                                                         
   20 DO J=1,M                                                      
        W(J) = -X(J)                                                   
        Y(J)=W(J)                                                         
      ENDDO
   30 MM = M-1                                                          
      ID = ND                                                           
      IBR = 0                                                           
      M1 = NM1                                                          
      M2 = NM2                                                          
      IA = NA                                                           

   35 CONTINUE
      IF(IA.GT.0) THEN
        RT = AA(IA)                                                       
        IA = IA-1                                                         
        DO J=1,M                                                      
           Y(J) = RT*W(J)                                                 
        ENDDO
      ENDIF

   50 IF (ID.LE.0) GOTO 150
      RT = BD(ID)                                                       
      ID = ID-1                                                         
      IF (ID .EQ. 0) IBR = 1                                            
      D(M) = A(M)/(B(M)-RT)                                             
      W(M) = Y(M)/(B(M)-RT)                                             
      DO J=2,MM                                                     
        K = M-J                                                        
        DEN = B(K+1)-RT-C(K+1)*D(K+2)                                  
        D(K+1) = A(K+1)/DEN                                            
        W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN                            
      ENDDO
      DEN = B(1)-RT-C(1)*D(2)                                           

      IF (DEN.EQ.0) THEN
        W(1) = 1.                                                         
      ELSE
        W(1) = (Y(1)-C(1)*W(2))/DEN                                       
      ENDIF
      DO J=2,M                                                      
        W(J) = W(J)-D(J)*W(J-1)                                        
      ENDDO
      IF (NA)  90, 90, 35                                               

   80 DO J=1,M                                                      
        Y(J) = W(J)                                                    
      ENDDO
      IBR = 1                                                           
      GO TO  35                                                         

   90 IF (M1.GT.0)  GOTO 100                                               
      IF (M2)  80, 80,125                                               

  100 IF (M2.LE.0) GOTO 110
      IF (ABS(BM1(M1)).LE.ABS(BM2(M2))) GOTO 125
  110 IF (IBR.GT.0) GOTO 120                                              
      IF (ABS(BM1(M1)-BD(ID)).LT.ABS(BM1(M1)-RT))  GOTO 80
  120 RT = RT-BM1(M1)                                                   
      M1 = M1-1                                                         
      GO TO 140                                                         
  125 IF (IBR.GT.0) GOTO 135                                              
      IF (ABS(BM2(M2)-BD(ID)).LT.ABS(BM2(M2)-RT))  GOTO 80
  135 RT = RT-BM2(M2)                                                   
      M2 = M2-1                                                         
  140 DO J=1,M                                                      
         Y(J) = Y(J)+RT*W(J)                                            
      ENDDO
      GO TO  35                                                         

  150 RETURN                                                            
END SUBROUTINE PRDCT

!*******************************************************************************
      
SUBROUTINE COMPB (N,IERROR,AN,BN,CN,B,AH,BH)                      
      IMPLICIT real*8 (a-h,o-z)                                         

      DIMENSION AN(*),BN(*),CN(*),B(*),AH(*),BH(*)                      
      COMMON /CBLKT/ NPP,Kx,EPSLON,CNV,Lx,NCMPLX,IK,IZ,DUM(1)
!$OMP THREADPRIVATE(/CBLKT/)

      EPSLON = 1.                                                          
    5 EPSLON = EPSLON/10.                                                     
      DIF = 1.+EPSLON                                                      
      DIFH = DIF                                                        
      IF (DIFH.GT.1.)  GOTO 5                                          

   10 EPSLON = 100.*EPSLON                                                    
      BNORM = ABS(BN(1))                                                
      DO J=2,N                                                      
         BNORM = MAX(BNORM,ABS(BN(J)))                                
      ENDDO
      CNV = EPSLON*BNORM                                                   
      DO  45 J=1,N                                                      
        IF (J.NE.1)  GOTO 25
        IF (NPP.NE.0)  GOTO 45
   25   ARG = AN(J)*CN(J)                                              
        IF (ARG.LT.0) GOTO 140
        CHLD = SQRT(ARG)                                               
        IF (CN(J).LT.0)  GOTO 40
        B(J) = CHLD                                                    
        GO TO  45                                                      
   40   B(J) = -CHLD                                                   
   45 CONTINUE                                                          
      IF = 2**Kx                                                         
      LH = IF                                                           
      LH1 = IF+1                                                        
      I1 = 1                                                            

      DO 105 Lx=1,Kx                                                      
        I1 = I1+I1                                                     
        NN = I1+I1-1                                                   
        IFL = IF-I1+1                                                  
        DO 100 I=1,IFL,I1                                              
          IF (I.LT.IFL)  GOTO 75
          IF (NPP.EQ.0)  GOTO 60                                        
          LH = LH+NN                                                  
          GO TO 100                                                   
   60     LS = 0                                                      
          DO J=I,IF                                               
            LS = LS+1                                                
            BH(LS) = BN(J)                                           
            AH(LS) = B(J)                                            
          ENDDO
          I1M = I1-1                                                  
          DO J=1,I1M                                              
            LS = LS+1                                                
            BH(LS) = BN(J)                                           
            AH(LS) = B(J)                                            
          ENDDO
          GO TO  85                                                   
   75     LS = 0                                                      
          JFL = I+NN-1                                                
          DO J=I,JFL                                              
            LS = LS+1                                                
            BH(LS) = BN(J)                                           
            AH(LS) = B(J)                                            
          ENDDO
   85     CALL TQLRAT (NN,BH,AH,IERROR)                               
          IF (IERROR.NE.0) GOTO 140
          DO J=1,NN                                               
            LH = LH+1                                                
            B(LH) = -BH(J)                                           
          ENDDO
  100   CONTINUE                                                       
        LH1 = LH+1                                                     
  105 CONTINUE                                                          

      DO J=1,N                                                      
        B(J) = -BN(J)                                                  
      ENDDO
      RETURN                                                            
  140 IERROR = 4                                                        
      RETURN                                                            
END SUBROUTINE COMPB

!*******************************************************************************
      
SUBROUTINE TQLRAT (N,D,E2,IERR)                                   
      IMPLICIT real*8 (a-h,o-z)                                         

      REAL*8 D(N),E2(N),B,C,F,G,H,P,R,S,MACHEP                          
!     COMMON /CBLKT/ NPP,Kx,MACHEP,CNV,LDZ,NCMPLX,IK,IZ,DUM(1)
      COMMON /CBLKT/ NPP,Kx,MACHEP,CNV,Lx,NCMPLX,IK,IZ,DUM(1)
!$OMP THREADPRIVATE(/CBLKT/)

      IERR = 0                                                          
      IF (N .EQ. 1) GO TO  70                                           
      DO I=2,N                                                      
        E2(I-1) = E2(I)*E2(I)                                          
      ENDDO
      F = 0.0                                                           
      B = 0.0                                                           
      E2(N) = 0.0                                                       
      DO  60 L=1,N                                                      
        J = 0                                                          
        H = MACHEP*(ABS(D(L))+SQRT(E2(L)))                             
        IF (B .GT. H) GO TO  10                                        
        B = H                                                          
        C = B*B                                                        
   10   DO M=L,N                                                   
          IF (E2(M) .LE. C) GO TO  20                                 
        ENDDO
   20   IF (M .EQ. L) GO TO  40                                        
   25   IF (J .EQ. 30) GO TO  65                                       
        J = J+1                                                        
        L1 = L+1                                                       
        S = SQRT(E2(L))                                                
        G = D(L)                                                       
        P = (D(L1)-G)/(2.0*S)                                          
        R = SQRT(P*P+1.0)                                              
        D(L) = S/(P+SIGN(R,P))                                         
        H = G-D(L)                                                     
        DO  I=L1,N                                                  
          D(I) = D(I)-H                                               
        ENDDO
        F = F+H                                                        
        G = D(M)                                                       
        IF (G .EQ. 0.0) G = B                                          
        H = G                                                          
        S = 0.0                                                        
        MML = M-L                                                      
        DO II=1,MML                                                
          I = M-II                                                    
          P = G*H                                                     
          R = P+E2(I)                                                 
          E2(I+1) = S*R                                               
          S = E2(I)/R                                                 
          D(I+1) = H+S*(H+D(I))                                       
          G = D(I)-E2(I)/G                                            
          IF (G .EQ. 0.0) G = B                                       
          H = G*P/R                                                   
        ENDDO
        E2(L) = S*G                                                    
        D(L) = H                                                       
        IF (H .EQ. 0.0) GO TO  40                                      
        IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO  40                        
        E2(L) = H*E2(L)                                                
        IF (E2(L) .NE. 0.0) GO TO  25                                  
   40   P = D(L)+F                                                     
        IF (L .EQ. 1) GO TO  50                                        
        ABSP = ABS(P)                                                  
        DO  II=2,L                                                  
          I = L+2-II                                                  
          IF (ABSP .GE. ABS(D(I-1))) GO TO  55                        
          D(I) = D(I-1)                                               
        ENDDO
   50   I = 1                                                          
   55   D(I) = P                                                       
   60 CONTINUE                                                          
      GO TO  70                                                         
   65 IERR = L                                                          
   70 RETURN                                                            

END SUBROUTINE TQLRAT
