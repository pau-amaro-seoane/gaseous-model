      SUBROUTINE HENYEY
********************************************************************************
*  Modernized version of 3-Pt-Henyey from H.W.Yorke                            *
*       and matrix inversion routine (Gauss elim. method) GIRL                 *
*  Rainer Spurzem, August 1992                                                 *
********************************************************************************
*  Logarithmic version of the equations, i.e. convergence for I=1,NPOS         *
*  means ABS(X-VX) < EPS     NOT(!!)    ABS(X-VX)/X < EPS                      *
*  for positive quantities; all vector quantities unchanged                    *
********************************************************************************
*  Note that substitution of GIRL-Routine by routines of special               *
*  library packages (ESSL, CRAY-SCILIB, e.g.) may cause erroneous              *
*  results due to limited accuracy of those programs for numerical             *
*  problems with bad condition numbers.                                        *
********************************************************************************
*                                                                              *
* Variables used:  (P) = Parameter contained in Include-File compar            *
*                  (C) = Common-Variable  "        "     "     "               *
*                  (L) = Local Variable                                        *
********************************************************************************
*    NG: Number of equations (P) 
*    NJ: Number of radial grid points (P)
*    NGO=2*(NG/2+1), NJO=2*(NJ/2+1) ensure efficient alignment of COMMONs on
*                                   double word computers (both (P))
*    (in the following comments assumed NGO=NG, NJO=NJ)
*
* NCOMP: Number of dynamical components in multi-comp. calc. (P) 
*  NPOS: Number of positive variables (P)
*
*  ITER: (Return value) Total number of iterations 
*        (in case of failure: return with ITER > 100) (C)
*   EPS: Convergence criterion for sum of squared corrections
*        of positive quantities  (C)
*        (Note: EPS is also used to define minimum absolute values 
*        W(Ixx) (C) of quantities which may become zero to avoid 
*        divergence of relative correction value)
*  FVEC: Convergence criterion for sum of squared corrections
*        of vector quantities (if=0 no check at all) (C)
* ITVEC: (Return value) Number of additional iterations needed
*        to reach convergence of vector quantities, if any. (C) 
* ITMIN: Minimum number of iterations in HENYEY (C)
* ITMAX: Maximum number of iterations in HENYEY (C)
*  FHEN: Factor applied at corrections for ITER<3 (C)
* FACDX: Modified FHEN-Factor applied at corrections (L)
*  FACC: Final correction value (L)
*    K2: Counter of mesh number in working loop (L)
*    K1: K2-1   (L)
*  IMOD: Model counter for whole run defined outside (C)
*
* Failure conditions:
* -------------------
* ITER > ITMAX: to be handled outside, expected action: reduce
*        timestep and retry.
* IFUT > 3 : if positive variable should be corrected to negative
*        value this is inhibited and IFUT incremented by one.
* For IFUT>3 an error message appears and return with ITER>100.
* (Note: For first model (IMOD=0) infinite number of retries are
*        allowed)
*
* Output options chosen from outside the program:
* -----------------------------------------------
* LWRIT1=IPR1.EQ.1: if true output of largest corrections per iteration
*                   and largest changes finally (L)
* LWRIT2=IPR2.EQ.1: if true output of LWRIT1 and square deviations 
*                   and max. gi per iteration (L)
*  (IPR1, IPR2: input coming from outside (C))
*
* Internally chosen Output options:
* ---------------------------------
* For first model of Run (IMOD=0, IMOD: (C))  LWRIT1=TRUE
* For ITER>ITMAX/2 additional output of square dev's (LSQUAR: (L))
*
* Global diagnostic output option:
* --------------------------------
*        (LS(25)=TRUE, defined outside (C))
* Option for checking matrix A or matrices B,C,D,E,F (select by ICTR) 
* by printout at certain triples of values (IMOD,ITER,K2), which are
* selected by the contents of the vector IC:
*
*    IC: vector for diag. print selection read outside (C)
*  ICTR: if= 1: print of matrices B,C,D,E,F
*        if=-1: print of matrix A       (C) 
*    JZ: variable controlling the current output (L)
* (for the output itself the subroutine CONTRO is called)
*        
*
* Henyey-specific quantities (dimension in brackets):
* ---------------------------------------------------
*
*     JG(NG): internal use, return with mesh number of max. changes  (C)
*     DG(NG): internal use, return value of max. gi (C)
*     GG(NG): internal use as sum of square corrections (L)
*     CS(17): returned as max. rel. change of positive quantities (C)
*     CS(18): returned as max. rel. change of vector quantities (C)
*   X(NG,NJ): dependent actual variables corrected in each iteration (C)
* XXY(NG*NJ): equivalenced to X to get one-index vector in HENYEY (L)
*  VX(NG,NJ): values of X at last model, not changed by HENYEY (C)
* VXY(NG*NJ): equivalenced to VX to get one-index vector in HENYEY (L)
* VVX(NG,NJ): values of X at model before last model, not used here
*
* A(NG*(2*NG+1)): matrix passed to subroutine GIRL; subroutine GIRL
*       inverts the first (NG*NG)-block of A and multiplies
*       the result with the remaining blocks, returning these
*       results in place of the remaining blocks. Returned
*       first block of A is not used. In any comments to this method
*       and in the inline comment statements the significant returned 
*       blocks of A are usually denoted as matrix Xjp of 
*       length NG*NG and as vector yjp, length NG, j=K2 actual mesh.
*  
* HY(NG*(NG+1)*NJ): Linearly arranged vector collecting the 
*       significant returned A-blocks (Xjp,yjp) from GIRL 
*       in ascending order from j=1 to j=NJ. For j=NJ
*       there are redundant non-used places in HY.
*  IH:  Integer flag used to memorize present position in vector HY
*
*      G(NG): value of equations (C)
*   C(NG,NG): being at mesh K2 C(J,K) is derivative of equation G(J) 
*             w.r.to variable X(K,K2-1) (C)
*  CC(NG*NG): equivalenced to C to get one-index vector in HENYEY (L)
*   D,DD: as B, but derivative w.r.to X(K,K2)
*   E,EE: as B, but derivative w.r.to X(K,K2+1)
*  (Note: CC,DD are also used internally for other quantities!)
*
*  Note: G,C,D,E are expected to be calculated in subr. GIDINI
*
* PROGRAMS USED: GIRL               H.W.YORKE 
* PROGRAMS USED: GIDINI ff          R.SPURZEM
*
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-K,M-N), LOGICAL(L)
*
         INCLUDE 'compar.f'
*
      DIMENSION CC(NGMAX*NGO),DD(NGMAX*NGO),EE(NGMAX*NGO)
      DIMENSION VVXY(NGMAX*NJO),VXY(NGMAX*NJO),XXY(NGMAX*NJO)
      DIMENSION GG(NGO)
      EQUIVALENCE (CC(1),C(1,1)),(DD(1),D(1,1)),(EE(1),E(1,1)),
     *    (VVXY(1),VVX(1,1)),(VXY(1),VX(1,1)),(XXY(1),X(1,1))
*
*          Defining convenient short notation
         NG1=NG+1
         NPOS1=NPOS+1
         NN=NG*NG
         NNG=NN+NG
         NNG1=NNG+1
*
*          Setting output options
         LWRIT1=IPR1.EQ.1
         LWRIT2=IPR2.EQ.1
         IF(IMOD.LE.1)LWRIT1=.TRUE.
         IF(LS(25))LWRIT2=.TRUE.
         IF(LWRIT2)LWRIT1=.TRUE.
         LSQUAR=LWRIT2
*
*          Initializing counters
         ITER=0
         ITVEC=0
         IFUT=0
         JZ=0
*
*          Iteration Loop
  500    ITER=ITER+1
*
*          Apply at first iterations smaller corrections if FHEN.NE.1.D0
      FACDX=1.D0
      IF(FHEN.LT.1.D0)THEN
      FACDX=FHEN*2.D0**ITER/2.D0
      IF(FACDX.GT.1.D0)FACDX=1.D0
      END IF
*
         IH=0
*
      IF(ITER.GT.ITMAX) RETURN
*         More Output for high iteration number
      IF(ITER.GT.ITMAX/2)THEN
      LSQUAR=.TRUE.
      LWRIT1=.TRUE.
      END IF
*
*         Start with inner boundary 
         K2=1
      CALL GIDINI(K2)
*
*         Optional Control Writing
      LCONTR=ICTR.EQ.1.AND.IMOD.EQ.IC(JZ+1).AND.
     *   ITER.EQ.IC(JZ+2).AND.K2.EQ.IC(JZ+3)
      IF(LCONTR)THEN
      CALL CONTRO(K2,2)
*
      IF(IC(JZ+4).LT.0.AND.IC(JZ+5).EQ.0)STOP
      JZ=JZ+3
      END IF
*         End Control Writing
*
*         Prepare Sections of matrix A = (Dj,-Ej,-gj)
      I2=NN
*
      DO 12 I1=1,NN
         I2=I2+1
         A(I1)= DD(I1)
         A(I2)=-EE(I1)
  12  CONTINUE
*
      DO 14 II=1,NG
         I2=I2+1
         A(I2)=-G(II)
         JG(II)=K2
         DG(II)=G(II)
  14  CONTINUE
*
      CALL GIRL(A,NG,NG1)
*
*         Optional Control Writing
      LCONTR=ICTR.EQ.-1.AND.IMOD.EQ.IC(JZ+1).AND.
     *   ITER.EQ.IC(JZ+2).AND.K2.EQ.IC(JZ+3)
*
      IF(NG1.LT.0)THEN
      PRINT*,' ERROR EXIT IN SUBROUTINE GIRL'
      CALL CONTRO(K2,2)
      CALL FLUSH(6)
      STOP
      END IF
*
      IF(LCONTR)THEN
      CALL CONTRO(K2,1)
*
      IF(IC(JZ+4).LT.0.AND.IC(JZ+5).EQ.0)STOP
      JZ=JZ+3
      END IF
*         End Control Writing
*
*         Store returned A-values in HY
      I2=NN
      DO 18 II=1,NNG
         IH=IH+1
         I2=I2+1
         HY(IH)=A(I2)
 18   CONTINUE
* 
*         Second Mesh Point
         K2=2
         K1=1
*
      CALL GIDINI(K2)
*
*         Optional Control Writing
      LCONTR=ICTR.EQ.1.AND.IMOD.EQ.IC(JZ+1).AND.
     *   ITER.EQ.IC(JZ+2).AND.K2.EQ.IC(JZ+3)
      IF(LCONTR)THEN
      CALL CONTRO(K2,2)
*
      IF(IC(JZ+4).LT.0.AND.IC(JZ+5).EQ.0)STOP
      JZ=JZ+3
      END IF
*         End Control Writing
*
*         Prepare Sections of matrix A
*         A=(  Dj+Cj*Xj,  -Ej,  -gj-Cj*yj  )
*
      DO 22 II=1,NG
         I1=II-NG
         IXJ=IH-NNG
*
        DO 23 J=1,NG
           I1=I1+NG
           I3=II-NG
           A(I1)=DD(I1)
           A(I1+NN)=-EE(I1)
*
          DO 231 K=1,NG
             I3=I3+NG
             IXJ=IXJ+1
             A(I1)=A(I1)+CC(I3)*HY(IXJ)
 231      CONTINUE
*
 23     CONTINUE
*
 22   CONTINUE
*
      I2=2*NN
      DO 24 II=1,NG
         I2=I2+1
         I3=II-NG
         IYJ=IH-NG
         A(I2)=-G(II)
*
        DO 241 K=1,NG
           I3=I3+NG
           IYJ=IYJ+1
           A(I2)=A(I2)-CC(I3)*HY(IYJ)
 241    CONTINUE
*
 24   CONTINUE
*         Matrix A ready
*
*        Check for Max G(I)
      DO 27 II=1,NG
      IF(DABS(G(II)).LE.DABS(DG(II))) GO TO 27
         JG(II)=K2
         DG(II)=G(II)
 27    CONTINUE
*
*         Loop for mesh points
  501    CONTINUE
*
      CALL GIRL(A,NG,NG1)
*
*         Optional Control Writing
      LCONTR=ICTR.EQ.-1.AND.IMOD.EQ.IC(JZ+1).AND.
     *   ITER.EQ.IC(JZ+2).AND.K2.EQ.IC(JZ+3)
*
      IF(NG1.LT.0)THEN
      PRINT*,' ERROR EXIT IN SUBROUTINE GIRL'
      CALL CONTRO(K2,2)
      CALL FLUSH(6)
      STOP
      END IF
*
      IF(LCONTR)THEN
      CALL CONTRO(K2,1)
*
      IF(IC(JZ+4).LT.0.AND.IC(JZ+5).EQ.0)STOP
      JZ=JZ+3
      END IF
*         End Control Writing
*
*         Store returned A-values in HY    
         I2=NN
      DO 28 II=1,NNG
         IH=IH+1
         I2=I2+1
         HY(IH)=A(I2)
  28  CONTINUE
*
      K1=K2
      K2=K2+1
*
      IF(K2.GT.NJ) GOTO 502
*
      CALL GIDINI(K2)
*
*         Prepare Matrix A
*    A=( Dj+Cj*Xj,  -Ej,  -gj-Cj*yj  )
*
      DO 42 II=1,NG
         I1=II-NG
         IXJ=IH-NNG
*
        DO 43 J=1,NG
           I1=I1+NG
           I3=II-NG
           A(I1)= DD(I1)
           A(I1+NN)=-EE(I1)
*
          DO 431 K=1,NG
             I3=I3+NG
             IXJ=IXJ+1
             A(I1)=A(I1)+CC(I3)*HY(IXJ)
 431      CONTINUE
*
 43     CONTINUE
*
 42   CONTINUE
*
         I2=2*NN
      DO 46 II=1,NG
         I2=I2+1
         I3=II-NG
         IYJ=IH-NG
         A(I2)=-G(II)
*
        DO 461 K=1,NG
           I3=I3+NG
           IYJ=IYJ+1
           A(I2)=A(I2)-CC(I3)*HY(IYJ)
 461    CONTINUE
*
 46   CONTINUE
*        Matrix A ready
*
*        Check for Max G(I)
      DO 47 II=1,NG
      IF(DABS(G(II)).LE.DABS(DG(II))) GO TO 47
         JG(II)=K2
         DG(II)=G(II)
 47    CONTINUE
*
*        Back to next mesh point
*
      GO TO 501
*
  502 CONTINUE
*
      IF(LWRIT2)WRITE(6,102)ITER,(JG(II),DG(II),II=1,NG)
*
      IH=IH+NNG-NG
      I1=NG*(K2-1)
*
 503  K2=K2-1
*
      I1=I1-NG
C**Min values for computation of relative corrections*****
         IGR=I1/NG+1
         IF(IGR.LT.2)IGR=2
         W(IMR)=0.D0
         RHOS=0.D0
         DO 239 J=1,NCOMP
 239     RHOS=RHOS+DEXP(VXY(I1+I00+J))
         TDY=1.D0/DSQRT(4.D0*PI*GRAV*RHOS)
         DO 240 J=1,NCOMP
         W(I10+J)=EPS*R(IGR)/TDY
         W(I30+J)=W(I10+J)
         W(I12+J)=W(I30+J)
         W(I00+J)=0.D0
         W(I20+J)=0.D0
         W(I02+J)=0.D0
 240     CONTINUE
C*******************************************************************
*         DXj = Xjp*DXjp + yjp
*         Temporarily: Dxjp  = DDj
*                      Dxj   = CCj
********************************************************************
*         Note: Take Care that by Gidini ff it is set
*               EE=0  (for I=NJ)    which is
*         implicitly assumed in Do-Loop 53 at I=NJ
*********************************************************************
         IH=IH-NNG
*
      DO 54 II=1,NG
         CC(II)=HY(IH+II)
         IYJ=IH+II-NG
         IXJ=IH+II-NNG
*
        DO 53 K=1,NG
           IXJ=IXJ+NG
         CC(II)=CC(II)+HY(IXJ)*DD(K)
  53    CONTINUE
*
  54  CONTINUE
*
        DO 56 II=1,NG
         FACC=CC(II)*FACDX
         XXY(II+I1)=XXY(II+I1)+FACC
         EE(II)=DD(II)
         DD(II)=CC(II)
C*****Normalization of Corrections: W contains lowest values
C*****for limiting the relative corrections; if W is zero the
C*****relative correction is restricted if the variable is
C*****exactly zero. For logarithmic quantities (I<NPOS) no
C*****restriction - absolute correction is taken
          IF(II.LE.NPOS)THEN
          DX=FACC
          ELSE
          XEPS=EPS*DABS(XXY(II+I1))+W(II)
          IF(XXY(II+I1).EQ.0.D0)XEPS=EPS+W(II)
          DX=FACC/(DABS(XXY(II+I1))+XEPS)
          END IF
*
      IF(K2.EQ.NJ)THEN
          GG(II)=DX*DX
          JG(II)=K2
          G(II)=DX
      ELSE
          GG(II)=GG(II)+ DX*DX
       IF(DABS(DX).GT.DABS(G(II)))THEN 
          JG(II)=K2
          G(II)=DX
       END IF
      END IF
  56  CONTINUE
*
      IF(K2.GT.1) GO TO 503
C+
      IF(LWRIT1) WRITE(6,103) (JG(II),G(II),II=1,NG)
      IF(ITER.LT.ITMIN) GO TO 500
          DO 59 II=1,NG
          GG(II)=DSQRT(GG(II))/NJ
 59       CONTINUE
          IF (LSQUAR) WRITE(6,105) (GG(II),II=1,NG)
         LITER=.FALSE.
         LVEC=.FALSE.
      DO 62 II=1,NPOS
 62     LITER=LITER .OR. GG(II).GT.EPS
      IF(LITER) GO TO 500
*
*        Check convergence of vector quantities
      LVEC2=FVEC.GT.0.D0
      IF(LVEC2)THEN
      LVEC=.FALSE.
      DO 63 II=NPOS+1,NG
 63   LVEC=LVEC.OR.(GG(II).GT.FVEC)
      IF(LVEC) ITVEC=ITVEC+1
      IF(LVEC) GOTO  500
      END IF
 66   CONTINUE
*
*        Calculate relative changes after final iteration
 64       CONTINUE
      DO 72 II=1,NG
      IF(II.GT.NPOS)XEPS=EPS*VXY(II)
      IF(VXY(II).EQ.0.D0)XEPS=EPS
         JG(II)=1
      IF(II.LE.NPOS)DD(II)=XXY(II)-VXY(II)
      IF(II.GT.NPOS)DD(II)=(XXY(II)-VXY(II))/(VXY(II)+XEPS)
  72  CONTINUE
*
         I1=0
      DO 76 K=2,NJ
         I1=I1+NG
      DO 74 II=1,NG
      IF(II.GT.NPOS.AND.K.EQ.NJ)GOTO 74
      IF(II.GT.NPOS)XEPS=EPS*VXY(II+I1)
      IF(VXY(II+I1).EQ.0.D0)XEPS=EPS
      XXX=XXY(II+I1)-VXY(II+I1)
      IF(II.GT.NPOS)XXX=XXX/(VXY(II+I1)+XEPS)
      IF(DABS(XXX).LE.DABS(DD(II))) GO TO 74
         JG(II)=K
         DD(II)=XXX
  74  CONTINUE
  76  CONTINUE
*
*         Compute max. relative changes of positive quantities
      DMAX=DABS(DD(1))
      DO 82 II=2,NPOS
       IF(DABS(DD(II)).LT.DMAX) GOTO 82
       DMAX=DABS(DD(II))
 82    CONTINUE
      CS(17)=DMAX
*
*         Compute max. relative changes of vector quantities
      DMAX=DABS(DD(NPOS1))
       DO 83 II=NPOS1,NG
       IF(DABS(DD(II)).LT.DMAX) GOTO 83
       DMAX=DABS(DD(II))
 83    CONTINUE
      CS(18)=DMAX
*
*         Final printout
      IF(LWRIT1) WRITE(6,104) ITER,(JG(II),DD(II),II=1,NG)
      RETURN
  102 FORMAT(' ITER=',I3,1X,'LARGEST GI*DT:',1X,1P,7(1X,I4,D10.2),
     1     4(/,8(1X,I4,D10.2)))
  103 FORMAT(1X,'(3L)LARGEST CORRECTIONS:',1X,1P,
     1  7(1X,I4,D10.2),4(/,8(1X,I4,D10.2)))
  105     FORMAT(3X,'SQUARE DEVIATION:',1P,
     1  10(1X,D9.2),3(/,12(1X,D9.2)))
  104 FORMAT(1X,'LARGEST CHANGES: ITER=',I3,1P,7(1X,I4,D10.2),
     1     4(/,8(1X,I4,D10.2)))
          END
*
**************************************************************************
*
*         Routine for diagnostic output used if LS(25)=.TRUE.
      SUBROUTINE CONTRO(K2,ICON)
*
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-K,M-N), LOGICAL(L)
         INCLUDE 'compar.f'
*
      DIMENSION CC(NGO*NGO),DD(NGO*NGO),EE(NGO*NGO)
      DIMENSION VVXY(NGO*NJO),VXY(NGO*NJO),XXY(NGO*NJO)
      DIMENSION GG(NGO)
      EQUIVALENCE (CC(1),C(1,1)),(DD(1),D(1,1)),(EE(1),E(1,1)),
     *    (VVXY(1),VVX(1,1)),(VXY(1),VX(1,1)),(XXY(1),X(1,1))
*
*        Control output for check of output in A-matrix (ICON=1)
*                             and of B,C,D,E,F matrices (ICON=2)
*
      IF(ICON.EQ.1)THEN
*
      NG1=NG+1
      IX=-NG
       WRITE(6,112) K2
 112  FORMAT(' STUETZSTELLE NR.',I3)
*
      DO 143  IY=1,NG1+NG
      IX=IX+NG
      IF(IY.EQ.1) WRITE(6,113)
      IF(IY.EQ.NG1) WRITE(6,115)
      IF(IY.EQ.NG1+NG)WRITE(6,116)
*
      WRITE(6,111)(A(IX+K),K=1,NG)
*
 143      CONTINUE
*
       END IF
*
 111  FORMAT(1X,1P,10D10.2,/,6D10.2)
 113  FORMAT(' MATRIXELEMENTE JUNK')
 115  FORMAT(' MATRIXELEMENTE X(I+1)')
 116  FORMAT('         VEKTOR Y(I+1)')
*
      IF(ICON.EQ.2)THEN
*
        WRITE(6,119) K2,ITER,IMOD,(G(II),II=1,NG)
 119  FORMAT(' STUETZST.NR.',I3,' ITER=',I3,' IMOD=',I8,
     *     ' G(II)=',1P,10D10.2,/,6D10.2)
*
      DO 142 IY=1,3
      IX=-NG
      IF(IY.EQ.1)WRITE(6,137)
      IF(IY.EQ.2)WRITE(6,138)
      IF(IY.EQ.3)WRITE(6,139)
*
      DO 141 IZ=1,NG
      IX=IX+NG
      IF(IY.EQ.1)WRITE(6,135) (CC(IX+II),II=1,NG)
      IF(IY.EQ.2)WRITE(6,135) (DD(IX+II),II=1,NG)
      IF(IY.EQ.3)WRITE(6,135) (EE(IX+II),II=1,NG)
*
 141  CONTINUE
*
 142  CONTINUE
*    set printout parameters
       LS(25) = .TRUE.
       IC(1) = IMOD
       IC(2) = ITER
       IC(3) = K2
       CALL GIDINI(K2)
      END IF
*
 137  FORMAT(' MATRIXELEMENTE C')
 138  FORMAT(' MATRIXELEMENTE D')
 139  FORMAT(' MATRIXELEMENTE E')
 135  FORMAT(1X,1P,10D10.2,/,6D10.2)
*
      RETURN
      END
