C***********************************************************************
      SUBROUTINE SETBDY(NCALL,ISYM)                                     
C                                                                       
C  This routine simply initializes various arrays before a call to   
C  BDYGEN is made.
C  If NCALL = 0, then innitialize all arrays in common block /BDY/.  
C           .GT.0, grid has moved, so re-evaluate RBDY, JPOS and KPOS
C                                                                       
C                                                                       
C  Dimension COSM and SINM (LMAX,10).                                
C  The dimension of RBDY,JPOS and KPOS depends on symmetries  
C  used (see also last dimension of BDYTRM in routine BDYGEN):       
C  If ABS(ISYM) = 1 or 8, dimension them (2*POT3JMAX + POT3KMAX -1).
C  = 2,3, or 9, dimension them (POT3JMAX + POT3KMAX -1).        
C                                                                       
      IMPLICIT real*8 (a-h,o-z)                                         

#include "hydroparam.h"
#include "globals.h"

C     
      PARAMETER(JKM1=2*POT3JMAX+POT3KMAX-1)
      COMMON /BDY/ COSM(LMAX,10),SINM(LMAX,10),RBDY(JKM1),BDYCHK
      COMMON /BDY1/ JPOS(JKM1),KPOS(JKM1),JKMAX 
C$OMP THREADPRIVATE(/BDY/,/BDY1/)
C                                                                       
C                                                                       
      ISYMA = IABS(ISYM)
      IF(NCALL.EQ.0) THEN
        BDYCHK = 1.0
        DO L=1,LMAX
          PSI = (L-0.5)*DTHETA
          DO M=1,10
            COSM(L,M)=COS(M*PSI)
            SINM(L,M)=SIN(M*PSI)
          ENDDO
        ENDDO
      ENDIF

      DO JK=1,JKM1                                           
        RBDY(JK)=0.0                                                 
        JPOS(JK)=0                                                   
        KPOS(JK)=0                                                   
      ENDDO

C  Fill RBDY with the distances from the origin to the boundary
C  points around a vertical slice. (Imagine slicing a cake.) J runs
C  along the cylindrical radius, k runs vertically. First we do
C  the points along the top edge of the slice:
      I=0
      ZK2=ZHF(POT3KMAX1)**2
      DO J=2,POT3JMAX1
        I=I+1
        RBDY(I)=SQRT(ZK2+RHF(J)**2)
        JPOS(I)=J
        KPOS(I)=POT3KMAX1
      ENDDO

C  Then along the bottom edge, unless equatorial plane symmytry
C  is assumed:
      IF(ISYMA.EQ.1.OR.ISYMA.EQ.8) THEN
        ZK2=ZHF(1)**2
        DO J=2,POT3JMAX1
          I=I+1
          RBDY(I)=SQRT(ZK2+RHF(J)**2)
          JPOS(I)=J
          KPOS(I)=1
        ENDDO
      ENDIF

C  Finally, along the outside vertical edge:
      RJ2=RHF(POT3JMAX1)**2
      DO K=2,POT3KMAX
        I=I+1
        RBDY(I)=SQRT(RJ2+ZHF(K)**2)
        JPOS(I)=POT3JMAX1
        KPOS(I)=K
      ENDDO
      JKMAX=I

      RETURN
      END



C**********************************************************************

      SUBROUTINE BDYGEN(MAXTRM,ISYM,REDGE)                              
Cxxx NB:  as presently written the code must be linked with:
C      f77 -mp boundary.o -Wl,Xlocal,terms_
      
C                                                                       
C...   CALL PARAMETERS  ...                                             
C  MAXTRM = maximum l to be used in spherical harmonic expansion.       
C           for greatest efficiency, use 4, 8, or 10.
C           10 is maximum allowable value.
C  ISYM = negative, all mass totally within grid boundary r's.
C       = positive, use general expansion since some mass outside.
C  ABS(ISYM) = 1,  no symmetries.
C            = 2,  sym. thru equatorial plane, full 2-pi.
C            = 3,  pi-symmetry and sym. thru equatorial plane.
C            = 8,  2-D with no symmetries.
C            = 9,  2-D with sym. thru equatorial plane.
C  REDGE = 0.0,  mass can be anywhere in grid.
C        .GT.0.0, mass entirely within sphere of radius = redge.
C                                                                       
C                                                                       
      IMPLICIT real*8 (a-h,o-z)
#include "hydroparam.h"
#include "globals.h"

      COMMON /INSIDE/ TMASS,ENEW,ELOST,EDIF,PHICHK,KLOCAT              

      parameter(igrid=0)
C igrid never gets changed from 0 anyhow.
C                                                                       
      PARAMETER(JKM1=2*POT3JMAX+POT3KMAX-1)
      COMMON /BDY/ COSM(LMAX,10),SINM(LMAX,10),RBDY(JKM1),BDYCHK
      COMMON /BDY1/ JPOS(JKM1),KPOS(JKM1),JKMAX
C$OMP THREADPRIVATE(/BDY/,/BDY1/)

C     
C                                                                       
C                                                                       
C The dimensions of TERM,Q,CM,SM,TRMTOT,MASSIN, MASSOUT and MRING should
C never be altered. They allow for expansions through L=10, M=10.
C                                                                       
C The last dimension of array BDYTRM, however, should be set equal to
C the size of arrays RBDY, JPOS, and KPOS.  (Since BDYTRM is so large
C it can, if need be, be reduced further if MAXTRM.LT.10.  The third
C dimension of BDYTRM must be .GE.NTRM, where NTRM is determined in
C loop 5 of this subroutine.)
      DOUBLE PRECISION   TERM(66),Q(10),CM(10),SM(10)                            
      DOUBLE PRECISION   TRMTOT(2,66),BDYTRM(2,2,66,JKM1)
      DOUBLE PRECISION   MASSIN(2,66,JKM1),MASSOUT(2,66,JKM1)
      DOUBLE PRECISION   MRING(2,66)
C These arrays are used only in this subroutine, and thus should not need
C to be put in a common block. But the OpenMP version of the code crashes
C with SIGSEGV if they are not. So I have put the private arrays into 
C threadprivate common BDYGENPRIVATE, and the shared arrays in common block
C BDYGENSHARED.
      COMMON /BDYGENPRIVATE/TERM,Q,CM,SM,TRMTOT,BDYTRM,MRING
      COMMON /BDYGENSHARED/MASSIN,MASSOUT
C$OMP THREADPRIVATE(/BDYGENPRIVATE/)

C
C
C
C
C
C    In a spherical harmonic expansion, traditionally,                
C
C (1)  YLM(THETA,PSI) = SQRT((2L+1)/4PI*(L-M)!/(L+M)!)*PLM(X)*EXP(I*M*PSI),
C
C      YL,-M(THETA,PSI) = (-1)**M*COMPLEX CONJUGATE(YLM),            
C
C           where:  X=COS(THETA)   and   I=SQRT(-1).                 
C
C    The expansion of boundary potentials in terms of spherical
C    harmonics leads to products of YLM and complex conjugate(YLM), so
C    a more useful definition of YLM is
C
C (2)  YLM(THETA,PSI) = SQRT((2L+1)/4PI)*SQRT(DLM)*BLM(THETA,PSI)    
C                       *EXP(I*M*PSI).
C
C    In this expression for YLM, PLM(X) has been factored into a leading
C    numerical coefficient 'COEF' and a theta-dependent expression BLM
C    then,
C
C      DLM = ((L-M)!/(L+M)!)*COEF**2.                
C
C    In practice, DLM always consists of an odd integer divided by 2**N
C    The power N varies with L and M, but for terms through L=10, N varies
C    from 0 to 18.  The following parameter statement contains exact
C    values of 2**(-N) for N = 5 through 18; e.g., TW8 = 2**(-8).     
C    These terms will be used to calculate appropriate DLM'S.         
C
C
      real*8 TW5,TW6,TW7,TW8,TW9,TW10,TW11,TW12,TW13,TW14,TW15,TW16,TW17
     $     ,TW18 
      PARAMETER (TW5=0.03125, TW6 =0.015625,       TW7 =7.8125E-3,
     1   TW8 =3.90625E-3,     TW9 =1.953125E-3,    TW10=9.765625E-04,
     2   TW11=4.8828125E-4,   TW12=2.44140625E-4,  TW13=1.220703125E-4,
     3   TW14=6.103515625E-5, TW15=3.051757813E-5, TW16=1.525878907E-5,
     4   TW17=7.629394535E-6, TW18=3.814697266E-6 )
C -------  For test & debug  ------------------------------------------
C     DATA TW5,TW6,TW7,TW8,TW9,TW10,TW11,TW12,TW13,TW14,TW15,TW16,TW17, 
C    1 TW18/0.03125,0.015625,7.8125E-3,3.90625E-3,1.953125E-3,          
C    2 9.765625E-04,4.8828125E-4,2.44140625E-4,1.220703125E-4,          
C    3 6.103515625E-5,3.051757813E-5,1.525878907E-5,7.629394535E-6,     
C    4 3.814697266E-6/                                                  
C ---------------------------------------------------------------------
C
C
C The following statement functions are the BLM's used in the
C definition of YLM in equation (2), above.  Variables that are
C multiplied by numerical coefficients are always even powers of
C COS(THETA); a variable preceding a parenthetical expression is always
C COS(THETA) to the first power; a variable trailing a parenthetical
C expression is always some power (either even or odd) of SIN(THETA)
C
C
      
      B20(AA) = 3.0*AA - 1.0                                              
      B30(AA,D) = D*(5.0*AA - 3.0)                                        
      B31(AA,D) = (5.0*AA - 1.0)*D                                        
      B40(AA,D) = 35.0*AA - 30.0*D + 3.0                                  
      B41(AA,D,E) = D*(7.0*AA - 3.0)*E                                    
      B42(AA,D) = (7.0*AA - 1.0)*D                                        
      B50(AA,D,E) = E*(63.0*AA - 70.0*D + 15.0)                           
      B51(AA,D,E) = (21.0*AA - 14.0*D + 1.0)*E                            
      B52(AA,D,E) = D*(3.0*AA - 1.0)*E                                    
      B53(AA,D) = (9.0*AA - 1.0)*D                                        
      B60(AA,D,E) = 231.0*AA - 315.0*D + 105.0*E - 5.0                    
      B61(AA,D,E,F) = E*(33.0*AA - 30.0*D + 5.0)*F                        
      B62(AA,D,E) = (33.0*AA - 18.0*D + 1.0)*E                            
      B63(AA,D,E) = D*(11.0*AA - 3.0)*E                                   
      B64(AA,D) = (11.0*AA - 1.0)*D                                       
      B70(AA,D,E,F) = F*(429.0*AA - 693.0*D + 315.0*E - 35.0)             
      B71(AA,D,E,F) = (429.0*AA - 495.0*D + 135.0*E -5.0)*F               
      B72(AA,D,E,F) = E*(143.0*AA - 110.0*D + 15.0)*F                     
      B73(AA,D,E) = (143.0*AA - 66.0*D + 3.0)*E                           
      B74(AA,D,E) = D*(13.0*AA - 3.0)*E                                   
      B75(AA,D) = (13.0*AA - 1.0)*D                                       
      B80(AA,D,E,F) = 6435.0*AA - 12012.0*D + 6930.0*E - 1260.0*F + 35.0  
      B81(AA,D,E,F,GG) = F*(715.0*AA - 1001.0*D +385.0*E - 35.0)*GG       
      B82(AA,D,E,F) = (143.0*AA - 143.0*D + 33.0*E -1.0)*F                
      B83(AA,D,E,F) = E*(39.0*AA - 26.0*D + 3.0)*F                        
      B84(AA,D,E) = (65.0*AA - 26.0*D + 1.0)*E                            
      B85(AA,D,E) = D*(5.0*AA - 1.0)*E                                    
      B86(AA,D) = (15.0*AA - 1.0)*D                                       
      B90(AA,D,E,F,GG) = GG*(12155.0*AA - 25740.0*D + 18018.0*E           
     &                 -4620.0*F + 315.0)                               
      B91(AA,D,E,F,GG) = (2431.0*AA - 4004.0*D + 2002.0*E - 308.0*F
     &                    + 7.0)*GG
      B92(AA,D,E,F,GG) = F*(221.0*AA - 273.0*D + 91.0*E - 7.0)*GG         
      B93(AA,D,E,F) = (221.0*AA - 195.0*D + 39.0*E - 1.0)*F               
      B94(AA,D,E,F) = E*(17.0*AA - 10.0*D + 1.0)*F                        
      B95(AA,D,E) = (85.0*AA - 30.0*D + 1.0)*E                            
      B96(AA,D,E) = D*(17.0*AA - 3.0)*E                                   
      B97(AA,D) = (17.0*AA - 1.0)*D                                       
      B100(AA,D,E,F,GG) = 46189.0*AA - 109395.0*D + 90090.0*E
     &                  - 30030.0*F + 3465.0*GG - 63.0                              
      B101(AA,D,E,F,GG,X) = GG*(4199.0*AA - 7956.0*D + 4914.0*E
     &                    - 1092.0*F + 63.0)*X                                     
      B102(AA,D,E,F,GG) = (4199.0*AA - 6188.0*D + 2730.0*E - 364.0*F      
     &                  + 7.0)*GG                                       
      B103(AA,D,E,F,GG) = F*(323.0*AA - 357.0*D + 105.0*E - 7.0)*GG       
      B104(AA,D,E,F) = (323.0*AA - 255.0*D + 45.0*E - 1.0)*F              
      B105(AA,D,E,F) = E*(323.0*AA - 170.0*D + 15.0)*F                    
      B106(AA,D,E) = (323.0*AA - 102.0*D + 3.0)*E                         
      B107(AA,D,E) = D*(19.0*AA - 3.0)*E                                  
      B108(AA,D) = (19.0*AA - 1.0)*D                                      
C                                                                       
C                                                                       
C                                                                       
C                                                                       
C      The potential at any point (rb,thetab,psib) is given by:         
C                                                                       
C         -GRAV*SUM OVER ALL L,M (M.NE.0) OF                            
C                                                                       
C         2.0*COS(M*PSIB)*DLM*BLM(THETAB)                               
C            *(RB**-(L+1)*MASSIN1(L,M) + RB**L*MASSOUT1(L,M))           
C        +2.0*SIN(M*PSIB)*DLM*BLM(THETAB)                               
C            *(RB**-(L+1)*MASSIN2(L,M) + RB**L*MASSOUT2(L,M))           
C                                                                       
C         PLUS -GRAV*SUM OVER ALL L,M=0 OF                              
C                                                                       
C         DL0*BL0(THETAB)*(RB**-(L+1)*MASSIN1(L,0) + RB**L*MASSOUT1(L,0)
C                                                                       
C      Given that the point (rb,thetab,psib) is at a spherical radius   
C      rspher, the terms massin and massout (each a function of l and m)
C      are integrals over the mass distribution inside and outside,     
C      respectively, of rspher.   Letting term(r,theta,psi) =           
C      blm(theta)*rho(r,theta,psi)*dvolume(r,theta),                    
C                                                                       
C         MASSIN1(L,M) = sum over all theta,psi,r inside rspher of      
C                        R**L*COS(M*PSI)*TERM(R,THETA,PSI),             
C                                                                       
C         MASSIN2(L,M) = sum over all theta,psi,r inside rspher of      
C                        R**L*SIN(M*PSI)*TERM(R,THETA,PSI),             
C                                                                       
C         MASSOUT1(L,M) = sum over all theta,psi,r outside rspher of    
C                         R**-(L+1)*COS(M*PSI)*TERM(R,THETA,PSI),       
C                                                                       
C         MASSOUT2(L,M) = sum over all theta,psi,r outside rspher of    
C                         R**-(L+1)*SIN(M*PSI)*TERM(R,THETA,PSI).       
C                                                                       
C                                                                       
C  A Note About Indexing
C  ---------------------
C  There can be some confusion about indexing since (J,K,L) is used
C  elsewhere in the program to denote a point in the cylindrical grid,
C  while (L,M) is typically used to index spherical harmonics. To
C  clear things up, we will use double letters (JJ,KK,LL) in this
C  subroutine to denote cylindrical grid points, and retain (L,M)
C  exclusively for the spherical harmonics. JK is used to index points
C  on the boundary of a jk-plane.
C                                                                       
C                                                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
C                                                                   C   
C                                                                   C   
C                        NOW BEGIN PROGRAM.                         C   
C                                                                   C   
C                                                                   C   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC   
C                                                                       
C                                                                       

      IF(BDYCHK.NE.1.0)GO TO 990                                        
      IF(IGRID.GT.0) CALL SETBDY(IGRID,ISYM)                             

C$OMP   PARALLEL DEFAULT(PRIVATE)
C$OMP&  COPYIN(/BDY/,/BDY1/)
C$OMP&  SHARED(MAXTRM,ISYM,REDGE,
C$OMP&         TMASS,MASSIN,MASSOUT,
#if PARTICLE>0
C$OMP&         PHI,RHOTOT,r,rhf,z,zhf,rof3n,zof3n,dtheta,pi
#else
C$OMP&         PHI,RHO,r,rhf,z,zhf,rof3n,zof3n,dtheta,pi
#endif
C$OMP&         ,grav,bgden,JREQ,KZPOL)

      IF(ISYM.GT.0) THEN
        JMX = JKMAX                                                       
        ISMAX = 2                                                         
      ELSE
        JMX = 1                                                           
        ISMAX = 1                                                         
      ENDIF
      DO JK = 1,JMX                                                  
        DO IS = 1,ISMAX                                                
          DO I = 1,66                                                
            BDYTRM(IS,1,I,JK) = 0.0                                             
            BDYTRM(IS,2,I,JK) = 0.0                                             
          ENDDO
        ENDDO
      ENDDO

C$OMP SINGLE
      TMASS=0.0                                                         
C$OMP END SINGLE
      ISYMA = ABS(ISYM)                                                
      IF(ISYMA.EQ.2.OR.ISYMA.EQ.9) THEN
        F1 = 0.0                                
        F2 = 2.0
      ELSE IF(ISYMA.EQ.3) THEN
        F1 = 0.0                                              
        F2 = 4.0                                              
      ELSE
        F1 = 1.0                                                          
        F2 = 1.0                                                          
      ENDIF
      NTRM = 0
    5 DO L = 0,MAXTRM                                                
        NTRM = NTRM + (L+1)
      ENDDO
C                                                                       
C-----------------------------------------------------------------------
C
C...     BEGIN INTEGRALS OVER MASS DISTRIBUTION.     ...                
C
C$OMP BARRIER      
C$OMP DO  REDUCTION(+:TMASS)
      DO 300 KK = 2,KMAX
        ZZ = ZHF(KK)
        Z2 = ZZ*ZZ
        ZDEL = Z(KK+1) - Z(KK)
        DO 299 JJ = 2,JMAX

          RR = RHF(JJ)
          R2 = RR*RR
          RSPHER = SQRT(R2 + Z2)
      
C  If we are outside of redge, we skip directly to the end of the loop.
C  Thus, we only need to run the loop over the small grid, since we know
C  there will never be any mass outside that grid.
          IF(REDGE.GT.0.0.AND.RSPHER.GT.REDGE) GO TO 299
C                                                                       
C                                                                       
C  STEP ONE:  FOR THIS R,THETA, CALCULATE THE PRODUCTS BLM(THETA)*R**L
C             FOR ALL L,M.  STORE PRODUCTS IN ARRAY TERM.
C
           C1 = ZZ/RSPHER
           S1 = RR/RSPHER
           C2 = C1*C1
           C4 = C2*C2
           S2 = S1*S1
           S3 = S1*S2
           S4 = S2*S2
           IF(MAXTRM.LE.4)GO TO 602
           C6 = C2*C4
           C8 = C4*C4
           S5 = S1*S4
           S6 = S2*S4
           S7 = S1*S6
           S8 = S4*S4
           IF(MAXTRM.LE.8)GO TO 602
           C10 = C2*C8
           S9 = S1*S8
           S10 = S2*S8

  602      CONTINUE
           DO I = 1,NTRM                                              
             TERM(I) = 0.0                                                
           ENDDO
           DO L = 1,MAXTRM                                             
             Q(L) = RSPHER**L                                           
           ENDDO
C...   L.EQ.0 THRU 4                                               
           term(1) = 1.0                                                    
           term(2) = C1*Q(1)*F1                                              
           term(4) = B20(C2)*Q(2)                                           
           term(7) = B30(C2,C1)*Q(3)*F1                                      
           term(11) = B40(C4,C2) *Q(4)                                       
           IF(ISYMA.GE.8) GO TO 611                                     
           term(6) = S2*Q(2)                                                
           term(13) = B42(C2,S2) *Q(4)                                       
           term(15) = S4*Q(4)                                                
           IF(ISYMA.EQ.3) GO TO 611                                     
           term(3) = S1*Q(1)                                                 
           term(8) = B31(C2,S1) *Q(3)                                        
           term(10) = S3*Q(3)                                                
           IF(ISYMA.EQ.2) GO TO 611                                     
           term(5) = C1*S1*Q(2)                                               
           term(9) = C1*S2*Q(3)                                              
           term(12) = B41(C2,C1,S1) *Q(4)                                      
           term(14) = C1*S3*Q(4)                                              
  611      IF(MAXTRM.LE.4) GO TO 620                                    
C...   L.EQ.5 THRU 8                                               
           term(16) = B50(C4,C2,C1) *Q(5)*F1                                  
           term(22) = B60(C6,C4,C2) *Q(6)                                    
           term(29) = B70(C6,C4,C2,C1) *Q(7)*F1                               
           term(37) = B80(C8,C6,C4,C2) *Q(8)                                 
           IF(ISYMA.GE.8) GO TO 613                                     
           QH = Q(6)                                                    
           term(24) = B62(C4,C2,S2) *QH                                      
           term(26) = B64(C2,S4) *QH                                         
           term(28) = S6*QH                                                  
           QH = Q(8)                                                    
           term(39) = B82(C6,C4,C2,S2) *QH                                   
           term(41) = B84(C4,C2,S4) *QH                                      
           term(43) = B86(C2,S6) *QH                                         
           term(45) = S8*QH                                                  
           IF(ISYMA.EQ.3) GO TO 613                                     
           QH = Q(5)                                                    
           term(17) = B51(C4,C2,S1) *QH                                       
           term(19) = B53(C2,S3) *QH                                         
           term(21) = S5*QH                                                  
           QH = Q(7)                                                    
           term(30) = B71(C6,C4,C2,S1) *QH                                    
           term(32) = B73(C4,C2,S3) *QH                                      
           term(34) = B75(C2,S5) *QH                                         
           term(36) = S7*QH                                                  
           IF(ISYMA.EQ.2) GO TO 613                                     
           term(18) = B52(C2,C1,S2) *Q(5)                                     
           term(20) = C1*S4*Q(5)                                              
           QH = Q(6)                                                    
           term(23) = B61(C4,C2,C1,S1) *QH                                     
           term(25) = B63(C2,C1,S3) *QH                                       
           term(27) = C1*S5*QH                                                
           QH = Q(7)                                                    
           term(31) = B72(C4,C2,C1,S2) *QH                                    
           term(33) = B74(C2,C1,S4) *QH                                       
           term(35) = C1*S6*QH                                                
           QH = Q(8)                                                    
           term(38) = B81(C6,C4,C2,C1,S1) *QH                                  
           term(40) = B83(C4,C2,C1,S3) *QH                                    
           term(42) = B85(C2,C1,S5)*QH                                        
           term(44) = C1*S7*QH                                                
  613      IF(MAXTRM.LE.8)GO TO 620                                     
C...   L.EQ.9 THRU 10                                              
           term(46) = B90(C8,C6,C4,C2,C1) *Q(9)*F1                            
           term(56) = B100(C10,C8,C6,C4,C2) *Q(10)                          
           IF(ISYMA.GE.8)GO TO 617                                      
           QH = Q(10)                                                   
           term(58) = B102(C8,C6,C4,C2,S2) *QH                              
           term(60) = B104(C6,C4,C2,S4) *QH                                 
           term(62) = B106(C4,C2,S6) *QH                                    
           term(64) = B108(C2,S8) *QH                                       
           term(66) = S10*QH                                               
           IF(ISYMA.EQ.3) GO TO 617                                     
           QH = Q(9)                                                    
           term(47) = B91(C8,C6,C4,C2,S1)*QH                                  
           term(49) = B93(C6,C4,C2,S3)*QH                                    
           term(51) = B95(C4,C2,S5) *QH                                      
           term(53) = B97(C2,S7) *QH                                         
           term(55) = S9*QH                                                  
           IF(ISYMA.EQ.2) GO TO 617                                     
           term(48) = B92(C6,C4,C2,C1,S2) *QH                                 
           term(50) = B94(C4,C2,C1,S4) *QH                                    
           term(52) = B96(C2,C1,S6) *QH                                       
           term(54) = C1*S8*QH                                                
           QH = Q(10)                                                   
           term(57) = B101(C8,C6,C4,C2,C1,S1) *QH                             
           term(59) = B103(C6,C4,C2,C1,S3) *QH                               
           term(61) = B105(C4,C2,C1,S5) *QH                                  
           term(63) = B107(C2,C1,S7) *QH                                     
           term(65) = C1*S9*QH                                               
  617      CONTINUE                                                     
  620      CONTINUE                                                     

 
C
C
C  STEP TWO:
C    Calculate MASSIN1 and MASSIN2.  These MASSIN terms will depend on
C    L and M, and in general on the spherical radius of the boundary
C    zone being considered.  For each L,M term (I=1,NTRM), and for each
C    boundary cell (JK=1,JKMAX), MASSIN1 is stored in BDYTRM(1,1,I,JK)
C    and MASSIN2 is stored in BDYTRM(1,2,I,JK).            
C
 
C  Depending on chosen symmetry, F2 = 1.0, 2.0 or 4.0 to account for all mass.
           DV = 0.5*DTHETA*ZDEL*(R(JJ+1)**2 - R(JJ)**2)*F2                 
           DO M = 1,MAXTRM
             CM(M) = 0.0
             SM(M) = 0.0
           ENDDO

C  Calculate all the mass ring integrals for the ring that goes through
C  point (JJ,KK,0). First do the actual integrations over LL. H is the mass
C  in a single grid cell, with any symmetry accounted for in the elementary
C  volume DV.
           SUMASS = 0.0
  210      DO LL = 1,LMAX
#if PARTICLE>0
             H = RHOTOT(JJ,KK,LL)*DV
#else
             H = RHO(JJ,KK,LL)*DV
#endif
             SUMASS = SUMASS + H
             DO M = 1,MAXTRM
               CM(M) = CM(M) + COSM(LL,M)*H
               SM(M) = SM(M) + SINM(LL,M)*H
             ENDDO
           ENDDO

C  Then finish up by multiplying each integral by the corresponding TERM(I):
           I = 0                                                     
  225      DO L = 0,MAXTRM
             I = I+1                                              
             MRING(1,I) = TERM(I)*SUMASS                                
             MRING(2,I) = 0.0                                              
C            write(*,'(4i4,1p2e15.6)') jj,kk,l,0,mring(1,i),mring(2,i)  !DKB debug
             DO M = 1,L
               I = I+1                                              
               MRING(1,I) = TERM(I)*CM(M)                                 
               MRING(2,I) = TERM(I)*SM(M)                                
C              write(*,'(4i4,1p2e15.6)') jj,kk,l,m,mring(1,i),mring(2,i)  !DKB debug
             ENDDO
           ENDDO
           TMASS=TMASS+MRING(1,1)
           
C  Now accumulate the mass ring integrals in BDYTRM(1,1,I,JK) and
C  BDYTRM(1,2,I,JK) for all JK such that RSPHER<=RBDY(JK).
           IF(ISYM.LT.0) GO TO 1250                                     
           ICHK = 0                                                     
  235      DO JK = 1,JKMAX
             IF(RSPHER.LE.RBDY(JK)) THEN
               DO I = 1,NTRM                                          
                 BDYTRM(1,1,I,JK) = BDYTRM(1,1,I,JK) + MRING(1,I)                  
                 BDYTRM(1,2,I,JK) = BDYTRM(1,2,I,JK) + MRING(2,I)                  
               ENDDO
             ELSE
               ICHK = ICHK + 1                                              
             ENDIF
          ENDDO
          
C  If ICHK=0 this mass ring is inside all RBDY(JK), so  we don't need to
C  update any of the MASSOUT integrals.
           IF(ICHK.EQ.0)GO TO 299                                       

C                                                                       
C                                                                       
C  STEP THREE:
C    Calculate MASSOUT1 and MASSOUT2.  For each L,M term (I=1,NTRM) and
C    for each boundary cell (JK=1,JKMAX) MASSOUT1 is stored in
C    BDYTRM(1,2,I,JK) and MASSOUT2 is stored in BDYTRM(2,2,I,JK).
C  NOTE:  For any particular ring of mass at spherical radius r, its
C    contribution to MASSOUT is just R**-(2L+1) times its contribution
C    to MASSIN.
C                                                                       
           I = 0                                                     
           DO L = 0,MAXTRM                                         
             QH = RSPHER**(-(2*L+1))                                           
             DO M = 0,L                                           
               I = I + 1                                              
               MRING(1,I) = MRING(1,I)*QH                               
               MRING(2,I) = MRING(2,I)*QH                               
             ENDDO
           ENDDO
           
C  Accumulate the mass ring integrals in BDYTRM(2,1,I,JK) and
C  BDYTRM(2,2,I,JK) for all JK such that RSPHER>RBDY(JK).
  246      DO JK = 1,JKMAX                                          
             IF(RSPHER.GT.RBDY(JK)) THEN
               DO I = 1,NTRM
                 BDYTRM(2,1,I,JK) = BDYTRM(2,1,I,JK) + MRING(1,I)
                 BDYTRM(2,2,I,JK) = BDYTRM(2,2,I,JK) + MRING(2,I)
               ENDDO
             ENDIF
           ENDDO

           GO TO 299
C  If ISYM<0 then all mass is inside grid boundary r's. In this case
C  BDYTRM(1,*,*,JK) is the same for all JK and we can save work by just
C  accumulating BDYTRM(1,*,*,1). All BDYTRM(2,*,*,*) are 0.
 1250      CONTINUE
           DO I = 1,NTRM
             BDYTRM(1,1,I,1) = BDYTRM(1,1,I,1) + MRING(1,I)
             BDYTRM(1,2,I,1) = BDYTRM(1,2,I,1) + MRING(2,I)
           ENDDO

  299   CONTINUE
  300 CONTINUE
C$OMP END DO

C
C-----------------------------------------------------------------------
C  Combine the mass integrals from all processes.
C
C$OMP DO
      DO JK = 1,JKMAX
        DO I = 1,NTRM
          MASSIN(1,I,JK)  = 0.0
          MASSIN(2,I,JK)  = 0.0
          MASSOUT(1,I,JK) = 0.0
          MASSOUT(2,I,JK) = 0.0
        ENDDO
      ENDDO
C$OMP END DO

C ... perform parallel reduction
C$OMP CRITICAL (COMBINE_MASS_INTEGRALS)
      DO JK = 1,JKMAX
        DO I = 1,NTRM
          MASSIN(1,I,JK) = MASSIN(1,I,JK) + BDYTRM(1,1,I,JK)
          MASSIN(2,I,JK) = MASSIN(2,I,JK) + BDYTRM(1,2,I,JK)
          MASSOUT(1,I,JK) = MASSOUT(1,I,JK) + BDYTRM(2,1,I,JK)
          MASSOUT(2,I,JK) = MASSOUT(2,I,JK) + BDYTRM(2,2,I,JK)
        ENDDO
      ENDDO
C$OMP END CRITICAL (COMBINE_MASS_INTEGRALS)
      
C$OMP BARRIER

C
C      FINISHED CALCULATING MASSIN1, MASSIN2, MASSOUT1, AND MASSOUT2.   
C
C-----------------------------------------------------------------------
C
C
C  CALCULATE POTENTIAL FOR EACH BOUNDARY CELL.
C
C

C$OMP DO PRIVATE(KK)
      DO 400 JK = 1,JKMAX                                               
        JJ = JPOS(JK)                                                      
        KK = KPOS(JK)                                                      
        RSPHER = RBDY(JK)                                                 
C                                                                       
C                                                                       
C      STEP FOUR:  CALCULATE THE PRODUCT DLM*BLM AND STORE RESULTS IN   
C                  ARRAY TERM.                                          
C                                                                       
C                                                                       
        C1 = ZHF(KK)/RSPHER                                            
        S1 = RHF(JJ)/RSPHER                                            
        C2 = C1*C1                                                     
        C4 = C2*C2                                                   
        S2 = S1*S1                                                     
        S3 = S1*S2                                                    
        S4 = S2*S2                                                   
        IF(MAXTRM.LE.4)GO TO 304                                     
        C6 = C2*C4                                                   
        C8 = C4*C4                                                   
        S5 = S1*S4                                                    
        S6 = S2*S4                                                   
        S7 = S1*S6                                                    
        S8 = S4*S4                                                   
        IF(MAXTRM.LE.8)GO TO 304                                     
        C10 = C2*C8                                                  
        S9 = S1*S8                                                    
        S10 = S2*S8                                                  
  304   CONTINUE                                                     
        DO I = 1,66                                              
          TERM(I) = 0.0
        ENDDO
              
C...  L.EQ.0 THRU 4                                                     
        term(1) = 1.0                                                    
        term(2) = C1*F1                                                   
        term(4) = 0.25*B20(C2)                                           
        term(7) = 0.25*B30(C2,C1)*F1                                      
        term(11) = TW6 *B40(C4,C2)                                        
        IF(ISYMA.GE.8)GO TO 307                                      
        term(6) = 0.375*S2                                               
        term(13) = 5.0*TW5 *B42(C2,S2)                                    
        term(15) = 35.0*TW7 *S4                                           
        IF(ISYMA.EQ.3)GO TO 307                                      
        term(3) = 0.5*S1                                                  
        term(8) = 0.1875*B31(C2,S1)                                       
        term(10) = 0.3125*S3                                              
        IF(ISYMA.EQ.2) GO TO 307                                     
        term(5) = 1.5*C1*S1                                                
        term(9) = 1.875*C1*S2                                             
        term(12) = 0.3125 *B41(C2,C1,S1)                                    
        term(14) = 2.1875 *C1*S3                                           
  307   IF(MAXTRM.LE.4)GO TO 310                                     
C...  L.EQ.5 THRU 8                                               
        term(16) = TW6 *B50(C4,C2,C1) *F1                                  
        term(22) = TW8 *B60(C6,C4,C2)                                     
        term(29) = TW8 *B70(C6,C4,C2,C1) *F1                               
        term(37) = TW14 *B80(C8,C6,C4,C2)                                 
        IF(ISYMA.GE.8)GO TO 308                                      
        term(24) = 105.0*TW10 *B62(C4,C2,S2)                              
        term(26) = 63.0*TW9 *B64(C2,S4)                                   
        term(28) = 231.0*TW10 *S6                                         
        term(39) = 315.0*TW12 *B82(C6,C4,C2,S2)                           
        term(41) = 693.0*TW13 *B84(C4,C2,S4)                              
        term(43) = 429.0*TW12 *B86(C2,S6)                                 
        term(45) = 6435.0*TW15 *S8                                        
        IF(ISYMA.EQ.3)GO TO 308                                      
        term(17) = 15.0*TW7 *B51(C4,C2,S1)                                 
        term(19) = 35.0*TW8 *B53(C2,S3)                                   
        term(21) = 63.0*TW8 *S5                                           
        term(30) = 7.0*TW11 *B71(C6,C4,C2,S1)                              
        term(32) = 21.0*TW11 *B73(C4,C2,S3)                               
        term(34) = 231.0*TW11 *B75(C2,S5)                                 
        term(36) = 429.0*TW11 *S7                                         
        IF(ISYMA.EQ.2)GO TO 308                                      
        term(18) = 105.0*TW5*B52(C2,C1,S2)                                 
        term(20) = 315.0*TW7 *C1*S4                                        
        term(23) = 21.0*TW7 *B61(C4,C2,C1,S1)                               
        term(25) = 105.0*TW8 *B63(C2,C1,S3)                                
        term(27) = 693.0*TW8*C1*S5                                         
        term(31) = 21.0*TW10 *B72(C4,C2,C1,S2)                             
        term(33) = 231.0*TW9 *B74(C2,C1,S4)                                
        term(35) = 3003.0*TW10 *C1*S6                                      
        term(38) = 9.0*TW11 *B81(C6,C4,C2,C1,S1)                            
        term(40) = 1155.0*TW11 *B83(C4,C2,C1,S3)                           
        term(42) = 9009.0*TW11 *B85(C2,C1,S5)                              
        term(44) = 6435.0*TW11 *C1*S7                                      
  308   IF(MAXTRM.LE.8)GO TO 310                                     
C...  L.EQ.9 THRU 10                                              
        term(46) = TW14 *B90(C8,C6,C4,C2,C1)*F1                            
        term(56) = TW16 *B100(C10,C8,C6,C4,C2)                           
        IF(ISYMA.GE.8)GO TO 309                                      
        term(58) = 165.0*TW17 *B102(C8,C6,C4,C2,S2)                      
        term(60) = 2145.0*TW15 *B104(C6,C4,C2,S4)                        
        term(62) = 2145.0*TW18 *B106(C4,C2,S6)                           
        term(64) = 12155.0*TW17 *B108(C2,S8)                             
        term(66) = 46189.0*TW18 *S10                                    
        IF(ISYMA.EQ.3) GO TO 309                                     
        term(47) = 45.0*TW15 *B91(C8,C6,C4,C2,S1)                          
        term(49) = 1155.0*TW14 *B93(C6,C4,C2,S3)                          
        term(51) = 1287.0*TW14 *B95(C4,C2,S5)                             
        term(53) = 6435.0*TW16 *B97(C2,S7)                                
        term(55) = 12155.0*TW16*S9                                        
        IF(ISYMA.EQ.2)GO TO 309                                      
        term(48) = 495.0*TW12 *B92(C6,C4,C2,C1,S2)                         
        term(50) = 45045.0*TW13 *B94(C4,C2,C1,S4)                          
        term(52) = 2145.0*TW12 *B96(C2,C1,S6)                              
        term(54) = 109395.0*TW15 *C1*S8                                    
        term(57) = 55.0*TW15 *B101(C8,C6,C4,C2,C1,S1)                      
        term(59) = 2145.0*TW14 *B103(C6,C4,C2,C1,S3)                      
        term(61) = 429.0*TW14 *B105(C4,C2,C1,S5)                          
        term(63) = 36465.0*TW16 *B107(C2,C1,S7)                           
        term(65) = 230945.0*TW16 *C1*S9                                   
  309   CONTINUE                                                     
  310   CONTINUE                                                     
C
C
C  STEP FIVE:
C    Combine MASSIN and MASSOUT terms, keeping cosine dependent terms
C    TRMTOT(1,I) separate from sine dependent terms TRMTOT(2,I).
C
C
        IF(ISYM.GE.0) THEN
          I = 0
          DO L = 0,MAXTRM
            RIN = 1.0/RSPHER**(L+1)
            ROUT = RSPHER**L
            DO M = 0,L
              I = I + 1
              B = TERM(I)*RIN
              H = TERM(I)*ROUT
              TRMTOT(1,I) = B*MASSIN(1,I,JK) + H*MASSOUT(1,I,JK)
              TRMTOT(2,I) = B*MASSIN(2,I,JK) + H*MASSOUT(2,I,JK)
            ENDDO
          ENDDO
        ELSE
          I = 0
          DO L = 0,MAXTRM
            RIN = 1.0/RSPHER**(L+1)
            DO M = 0,L
              I = I+1
              B = TERM(I)*RIN
              TRMTOT(1,I) = B*MASSIN(1,I,1)
              TRMTOT(2,I) = B*MASSIN(2,I,1)
            ENDDO
          ENDDO
        END IF
                 
C
C
C  STEP SIX:
C    Sum over all L,M terms, taking into account the PSIB angle dependence.
C
C
        DO LL = 1,LMAX                                            
          DO M = 1,MAXTRM                                          
            CM(M) = COSM(LL,M)                                            
            SM(M) = SINM(LL,M)                                            
          ENDDO
          I = 0                                                     
          SUM = 0.0                                                    
          DO L = 0,MAXTRM
            I = I+1
            SUM = SUM + TRMTOT(1,I)
            DO M = 1,L                                           
              I = I+1
              SUM = SUM + 2.0*(CM(M)*TRMTOT(1,I) + SM(M)*TRMTOT(2,I))
            END DO
          END DO
          PHI(JJ,KK,LL) = -GRAV*SUM                                      
       ENDDO
       
       
  400 CONTINUE                                                          

C$OMP END DO

C$OMP END PARALLEL

c     write(22,*)                                       !dkb debug
c     do k=kmax1,1,-1                                   !dkb debug
c       write(22,9000) k,(phi(jmax1,k,l),l=1,lmax)      !dkb debug
c     enddo                                             !dkb debug
c9000 format(i4,1p12e14.5,/,(4x,1p12e14.5))             !dkb debug

c     do jk = 1,jkmax                                               !DKB debug
c       jj = jpos(jk)                                               !DKB debug
c       kk = kpos(jk)                                               !DKB debug
c       rspher = rbdy(jk)                                           !DKB debug
c       write(6,9020) jk,jj,kk,rspher,phi(jj,kk,2)                  !DKB debug
c     enddo                                                         !DKB debug
c9020 format('JK=',i3,'  JJ=',i3,'  KK=',i3,'  RSPHER=',1pe13.5,    !DKB debug
c    &       '  PHI=',1pe13.5)                                      !DKB debug

C                                                                       
C                                                                       
C-----------------------------------------------------------------------
C
C
C      FINISHED HARD WORK.  NOW TIDY UP BOUNDARY SYMMETRIES.            
C                                                                       
C                                                                       
      ISYMA = IABS(ISYM)                                                
      IF(ISYMA.EQ.1)GO TO 550                                           
C  For 2-D or pi-symmetry problems, set boundary conditions on z-axis.
        IF(ISYMA.NE.3.AND.ISYMA.NE.8.AND.ISYMA.NE.9)GO TO 525             
          DO LL=1,LMAX                                                   
            PHI(1,pot3KMAX1,LL) = PHI(2,pot3KMAX1,LL)
          ENDDO
          IF(ISYMA.NE.8)GO TO 525                                           
          DO LL=1,LMAX                                                   
            PHI(1,1,LL) = PHI(2,1,LL)                                           
          ENDDO
          GO TO 580                                                         
C                                                                       
C  For symmetry through equatorial plane, set boundry condition through plane. 
  525   IF(ISYMA.NE.2.AND.ISYMA.NE.3.AND.ISYMA.NE.9)GO TO 550             
        DO LL=1,LMAX                                                   
          PHI(pot3JMAX1,1,LL) = PHI(pot3JMAX1,2,LL)
        ENDDO
        IF(ISYMA.NE.2)GO TO 580 
C                                                                       
C  If no symmetry through z-axis, equate correct phi's at J=1.               
  550 CONTINUE
        LHAF = LMAX/2
        DO LL=1,LHAF
          LP=LL+LHAF
          PHI(1,pot3KMAX1,LL) = PHI(2,pot3KMAX1,LP)
          PHI(1,pot3KMAX1,LP) = PHI(2,pot3KMAX1,LL)
          IF(ISYMA.EQ.1) THEN
            PHI(1,1,LL) = PHI(2,1,LP)
            PHI(1,1,LP) = PHI(2,1,LL)
          ENDIF
        ENDDO
 580  CONTINUE


      RETURN                                                            
C                                                                       
C                                                                       
C                                                                       
  990 WRITE(3,101)BDYCHK                                                
  101 FORMAT(//,' STOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOP',
     1      5X,'BDYCHK =',1PE10.2,/,5X,'SETBDY HAS NOT BEEN CALLED.',/, 
     2      5X,'IT MUST BE CALLED AT LEAST ONCE IN ORDER TO INITIALIZE',
     3      5X,'THE ARRAYS IN COMMON BLOCK /BDY/',//,                   
     4      ' STOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOPSTOP',//) 
      STOP                                                              
      END                                                               
 
 
C***********************************************************************
       
      subroutine debug
      
      return
      
      end
 
C***********************************************************************
       
      SUBROUTINE MASS_REDUCE(JKMAX,NTRM,MASSIN,MASSOUT,BDYTRM)
      INTEGER :: JKMAX,NTRM
      DOUBLE PRECISION :: MASSIN(2,66,JKMAX),MASSOUT(2,66,JKMAX)
      DOUBLE PRECISION, TARGET :: BDYTRM(2,2,66,JKMAX)
c
c ... This code works only up to 512 threads.
c     It doesn't do error checking (maybe should).
      INTEGER, PARAMETER :: MAXTHREADS = 512
      TYPE VTYPE
         DOUBLE PRECISION, POINTER :: PBDY(:,:,:,:)
      END TYPE VTYPE
      TYPE (VTYPE), SAVE :: VA(MAXTHREADS)
C
!$      INTEGER, EXTERNAL  :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
C
      INTEGER :: MYID,NTHR,I,N,JK
C
!$      MYID = OMP_GET_THREAD_NUM() + 1
!$      NTHR = OMP_GET_NUM_THREADS()
      VA(MYID)%PBDY => BDYTRM
!$OMP BARRIER

! ... reduction is performed here
!$OMP DO
      DO JK = 1,JKMAX
        DO N = 1, NTHR
          DO I = 1,NTRM
            MASSIN(1,I,JK) = MASSIN(1,I,JK) + VA(N)%PBDY(1,1,I,JK)
            MASSIN(2,I,JK) = MASSIN(2,I,JK) + VA(N)%PBDY(1,2,I,JK)
            MASSOUT(1,I,JK) = MASSOUT(1,I,JK) + VA(N)%PBDY(2,1,I,JK)
            MASSOUT(2,I,JK) = MASSOUT(2,I,JK) + VA(N)%PBDY(2,2,I,JK)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
       print *, "Aaaah, this shouldn't be happening"
! ... deassociate the pointers
      NULLIFY(VA(MYID)%PBDY)

      RETURN
      END
