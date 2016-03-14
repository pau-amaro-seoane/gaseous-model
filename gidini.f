      SUBROUTINE GIDINI(III)
C***************************************************************
C     
C     Routine called by Henyey for each Grid Point
C     
C     Calls Equations, boundary conditions and
C     manages diagnostic printout (LS(25)=TRUE)
C--------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C     
      INCLUDE 'compar.f'
      COMMON/EBTEST/EBALA,EBEQ,EBEQ2,DEMASS,DEENER,DEMASS2
C***************************************************************
C     Initialize Auxiliary Control Variables
C---------------------------------------------------------------
      IV=ITER
      I=IABS(III)
      L1=I.EQ.1
      L2=I.EQ.2
      LNJ=I.EQ.NJ
      IF(L1)EBEQ=0.D0
      IF(L1)EBEQ2=0.D0
      IF(L1)DEMASS=0.D0
      IF(L1)DEENER=0.D0
      IF(L1)DEMASS2=0.D0
      F2=0.D0
      FNJ=0.D0
      IF(L2)F2=1.D0
      IF(LNJ)FNJ=1.D0
      IM=I-1
      IP=I+1
      IF (LNJ) IP=NJ
      IF (L1) IM=1
      LSTEVO = LS(21).AND.(ITER.EQ.1).AND.L2
C************************************************************
C     If LS(25)=TRUE print out is directed by values in IC-Vector
C---------------------------------------------------------------
      IF(LS(25))LDIAG=.TRUE.
      LS(25)=.FALSE.
C     
      IF(III.EQ.1.AND.IMOD.EQ.0)JZ=0
      IF(LDIAG.AND.(IMOD.EQ.IC(JZ+1)).AND.
     *     (ITER.EQ.IC(JZ+2)).AND.(IABS(III).EQ.IC(JZ+3)))THEN
         LS(25)=.TRUE.
      END IF
C*************************************************************
C     Set FX(J) auxiliary vectors for switching off equations
C     Set HX(J) for hybrid implicit-explicit scheme
C-------------------------------------------------------------
      DO 100 J=1,NG
         FX(J)=0.D0
         IF(LEQ(J))FX(J)=1.D0
         HX(J)=XIMPL*X(J,I)+(1.D0-XIMPL)*VX(J,I)
         HXM(J)=XIMPL*X(J,IM)+(1.D0-XIMPL)*VX(J,IM)
         HXP(J)=XIMPL*X(J,IP)+(1.D0-XIMPL)*VX(J,IP)
 100  CONTINUE
C     
C----------------------------------------------------------------------
C     Set non-used equations trivial (d/dt=0)
C-----------------------------------------------------------------------
      DO 7 KK=1,NG
         G(KK)=(X(KK,I)-VX(KK,I))*BREM
         DO 6 JJ=1,NG
            D(KK,JJ)=0.D0
            E(KK,JJ)=0.D0
            C(KK,JJ)=0.D0
 6       CONTINUE
C-------Divide by XIMPL due to modification of derivs. later
         D(KK,KK)=BREM/XIMPL
 7    CONTINUE
C***********************************************************************
C     Inner Boundary Condition
C-----------------------------------------------------------------------
      IF(L1)THEN
         CALL INNER
         GOTO 540
      END IF
C***********************************************************************
C     Initialize Auxiliary Quantities 
C-----------------------------------------------------------------------
      CALL GINIT
C**********************************************************************
C     Mass Equation (Poisson Equation)
C----------------------------------------------------------------------
      IF(LEQ(IMR))CALL IMREQ
 20   CONTINUE
C**********************************************************************
C     Equation of Continuity 
C----------------------------------------------------------------------
      CALL I00EQ
*     
C***********************************************************************
C     Equation of linear momentum density Cloud/Star/Gas Components ( /S/G)
C-----------------------------------------------------------------------
*     One of the options 30 or 31 means to call outer boundary 
*     routine explicitly instead of using i10eq/i30eq/i12eq at R(NJ)
*     
      LOUTER = LS(30).OR.LS(31)
      IF(.NOT.LNJ.OR..NOT.LOUTER)CALL I10EQ
C     
C***********************************************************************
C     Second Order Moment Equations
C     LS(13)=TRUE:    Isotropic Version
C-----------------------------------------------------------------------
      IF(LS(13))THEN
         CALL I2IEQ
      ELSE
         CALL I20EQ
         CALL I02EQ
      END IF
C***********************************************************************
C     Binary Energy generation 
C-----------------------------------------------------------------------
      IF(TIME.GT.TBINI)THEN
         IF(LS(1))CALL FORMB2
         IF(LS(2))CALL FORMB3
      END IF
C     
      IF(LS(3))THEN
         IF(I.EQ.2)DEBIN=0.D0
         CALL FORMBS
      END IF
C*************************************************************************
C     Equipartition and Dynamical Friction Terms
C-------------------------------------------------------------------------
      IF(LS(12))CALL RXTERM
C     
C*************************************************************************
C     Mass Loss due to stellar evolution
C-------------------------------------------------------------------------
      IF(LSTEVO)CALL STEVOL
C**************************************************************************
C     Star-accreting black hole terms
C--------------------------------------------------------------------------
      IF(.NOT.LNJ.AND.LS(19))THEN
*     Energy diffusion at innermost shell or grav. radiation
*     inside core radius or redistribution of energy in core if XBIN>0
         IF(L2.OR.LS(20).OR.XBIN.NE.0.D0)CALL EDIFF
         IF(LS(14))CALL LOSSC
      END IF
C**************************************************************************
C     Third Order Equations (Heat Transfer Equations)
C     (Only taken if LS(9)=.FALSE. and not at bouter boundary)
C--------------------------------------------------------------------------
      IF(LS(9))GOTO 540
C--------------------------------------------------------------------------
C     LS(22)=.TRUE.: Two-flux closure see Louis and Spurzem 1991
C     LS(22)=.FALSE.:One-flux    "          "         "      "
C--------------------------------------------------------------------------
      IF(.NOT.LNJ.OR..NOT.LOUTER)THEN
         IF(LS(22))THEN
            CALL I30EQ
            CALL I12EQ
         ELSE
            CALL I30EQ1
         END IF
      END IF
C     
 540  CONTINUE
C------------------------------------------------------------------------
C     Outer Boundary Conditions only for optiosn 30 or 31 taken
C     
      IF(LNJ.AND.LOUTER)CALL OUTER
C***********************************************************************
C     Turn Equation into trivial equation d/dt = 0 if switch LEQ(J) is true.
C     Multiply XIMPL for hybrid implicit explicit scheme at all derivatives
C     Note: All D/Dt-Derivatives are to be divided by XIMPL before!
C----------------------------------------------------------------------
      DO 2000 K=1,NG
         DO 2002 J=1,NG
C     
            IF(LEQ(J))THEN
C     
               IF(LEQ(K))THEN
                  D(K,J)=D(K,J)*XIMPL
                  E(K,J)=E(K,J)*XIMPL
                  C(K,J)=C(K,J)*XIMPL
               END IF
C     
            ELSE
C     
               D(K,J)=0.D0
               E(K,J)=0.D0
               C(K,J)=0.D0
C     
            END IF
C     
 2002    CONTINUE
C     
         IF(.NOT.LEQ(K))THEN
C     
            DO 2001 J=1,NG
               D(K,J)=0.D0
               E(K,J)=0.D0
               C(K,J)=0.D0
 2001       CONTINUE
C     
            G(K)=(X(K,I)-VX(K,I))*BREM
            D(K,K)=BREM
         END IF
C     
 2000 CONTINUE
C     
C*********************************************************************
C     Diagnostic Printout follows
C------------------------------------------------------------------------
      IF(LDIAG.AND.(IMOD.EQ.IC(JZ+1)).AND.
     *     (ITER.EQ.IC(JZ+2)).AND.(I.EQ.IC(JZ+3)))THEN
         PRINT*,' HERE IS PRINTING INFORMATION OF SUBROUTINE',
     *        ' GIDINI -- IMOD=',IMOD,' ITER=',ITER,' GRID=',I
C     
         WRITE(6,1007)(X(K,I),K=1,NG)
 1007    FORMAT(' VARIABLES X(1-NG)=',1P,10(1X,D10.2),/,
     *        10(1X,D10.2))
         WRITE(6,1006)(G(K),K=1,NG)
 1006    FORMAT(' J=1-NG, G(J)=',1P,10(1X,D10.2),/,10(1X,D10.2))
         DO 1108 K=1,NG
            WRITE(6,1008)K,(D(K,J),J=1,NG)
 1108    continue
         DO 1109 K=1,NG
            WRITE(6,1009)K,(C(K,J),J=1,NG)
 1109    CONTINUE
         DO 1110 K=1,NG
            WRITE(6,1010)K,(E(K,J),J=1,NG)
 1110    CONTINUE
 1008    FORMAT(' K=,D(',I2,',J)=',1P,10D9.2,(/,10X,10D9.2))
 1009    FORMAT(' K=,C(',I2,',J)=',1P,10D9.2,(/,10X,10D9.2))
 1010    FORMAT(' K=,E(',I2,',J)=',1P,10D9.2,(/,10X,10D9.2))
      END IF
      IF(LDIAG.AND.(IMOD.EQ.IC(JZ+1)).AND.
     *     (ITER.EQ.IC(JZ+2)).AND.(IABS(III).EQ.IC(JZ+3)))THEN
         IF(IC(JZ+4).LT.0.AND.IC(JZ+5).EQ.0)STOP
         JZ=JZ+3
      END IF
C     
      RETURN
      END
