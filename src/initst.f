      SUBROUTINE INITST
C***********************************************
C     INITIALIZE STARTING MODEL FOR STELLAR COMPONENT
C     CALCULATING VELOCITY DISPERSION VARIABLES
C****************************************************************
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
         INCLUDE 'compar.f'
      DIMENSION XTMP(NCOMPO)
C       Allow large number of pressure integrations
      JENMAX = 1000
      FGRA=0.D0
      IF(LS(19).AND.LS(16))FGRA=1.D0
C
      RHOCT = 0.D0
      RHTOUT = 0.D0
      RHXOUT = 0.D0
      DO 5555 K=1,NCOMP
      RHXOUT = RHXOUT + DEXP(XOUT(K))
      RHTOUT = RHTOUT + DEXP(X(I00+K,NJ))
 5555 RHOCT = RHOCT + DEXP(X(I00+K,1))
*        Note for King model this loop is only passed once for total pressure
      DO 1000 IST=1,NCOMP
C--------Exact Calculation of pressure for non-loaded Plummer's model
      LEXACT=.NOT.LS(19).AND.NOL(IST).LT.0
      IATTE = 0
*
 2727 CONTINUE
*
      JEN = 0
      NJ1 = NJ - 1
*
      IF(.NOT.LEXACT)THEN
*
*       Prepare some data for the case of King model start
      IF(NOL(IST).GT.100)THEN
      W0 = DBLE(NOL(IST)-100)
      SIGC2=DEXP(X(I20+IST,1)-X(I00+IST,1))
      RKING=DSQRT(9.D0*SIGC2/4.D0/PI/GRAV/RHOCT)
*       Only determine PROUT at first attempt
      IF(IATTE.EQ.0)THEN
      IF(LS(27))THEN
      PROUT(IST)=DLOG(RHTOUT)+DLOG(DABS(PHTID))
      ELSE
      PROUT(IST)=-3.D1
      END IF
      END IF
*
      PRINT*,' Start with King model W0=',W0,' PHTID=',PHTID
      PRINT*,' SIGC should be=',DSQRT(SIGC2),' RKING=',RKING
      PRINT*,' Start with PROUT=',PROUT(IST),' SIGOUT=',
     *    DEXP(C12*(PROUT(IST)-DLOG(RHXOUT)))
*
      ELSE
*       Outer Pressures for non-King-models
      PROUT(IST)=XOUT(IST)-DLOG(6.D0)-
     *     C12*DLOG((R(NJ)*STEU1)**2+RSC(IST)**2)
      PTOUT(IST)=PROUT(IST)
*
      END IF
*
*       Main iteration loop to solve Jeans equation
      LSIGMA = .FALSE.
*
 27   CONTINUE
      AMAX=0.D0
C 
      DO 16 J=NJ,2,-1
*
      IF(JEN.EQ.0)THEN
      VX(I20+IST,J)=-3.D1
      ELSE
      VX(I20+IST,J)=X(I20+IST,J)
      END IF
*
      JP=J+1
      JM=J-1
      XLR=DLOG(R(J))
      IF(J.EQ.2)THEN
      RAV=C12*DLOG(R(J)*R(J)/STEU3)
      ELSE
      RAV=C12*DLOG(R(J)*R(JM))
      END IF
*
      IF(J.EQ.NJ)THEN
      RAVP=C12*DLOG(R(J)*R(J)*STEU1)
*
*       In case of King model use total density
      IF(NOL(IST).GT.100)THEN
      XAV = C12*(DLOG(RHTOUT)+DLOG(RHXOUT))
      ELSE
      XAV=C12*(X(I00+IST,J)+XOUT(IST))
      END IF
*
      IF(JEN.EQ.0)THEN
      PMEAN=PROUT(IST)
      ELSE
      PMEAN=C12*(X(I20+IST,J)+PROUT(IST))
      END IF
*
      ELSE
*
      RAVP=C12*DLOG(R(JP)*R(J))
*
*       In case of King model use total density
      IF(NOL(IST).GT.100)THEN
      RHOTI = 0.D0
      RHOTP = 0.D0
      DO 6666 K=1,NCOMP
      RHOTI = RHOTI + DEXP(X(I00+K,J))
 6666 RHOTP = RHOTP + DEXP(X(I00+K,JP))
      XAV = C12*(DLOG(RHOTI)+DLOG(RHOTP)) 
      ELSE
      XAV=C12*(X(I00+IST,J)+X(I00+IST,JP))
      END IF
*
      IF(JEN.EQ.0)THEN
      PMEAN=X(I20+IST,JP)
      ELSE
      PMEAN=C12*(X(I20+IST,J)+X(I20+IST,JP))
      END IF
*
      END IF
      DERAV=RAVP-RAV
      GRA=GRAV*DERAV*(DEXP(X(IMR,J)+XAV-XLR-PMEAN)+
     *     (XMRBIN(J)+FGRA*XMHOLE)*DEXP(XAV-XLR-PMEAN))
      IF(J.EQ.NJ)X(I20+IST,J)=PROUT(IST)+GRA
      IF(J.LT.NJ)X(I20+IST,J)=X(I20+IST,JP)+GRA
*
*       Take care of clean core vel. dispersion for King models
       IF(NOL(IST).GT.100)THEN
        SIGI2 = DEXP(X(I20+IST,J)-X(I00+IST,J))
        SIGP2 = DEXP(X(I20+IST,JP)-X(I00+IST,JP))
        IF(SIGI2.LT.SIGP2.AND.R(J).LT.RKING)THEN
*        PRINT*,' Sig inversion J=',J,
*    *           ' logrho=',X(I00+IST,J),X(I00+IST,JP)
         LSIGMA =.TRUE.
         X(I20+IST,J)=X(I20+IST,JP)-X(I00+IST,JP)+X(I00+IST,J)
        END IF
       END IF
*
      X(I02+IST,J)=X(I20+IST,J)
*
      XCORR=X(I20+IST,J)-VX(I20+IST,J)
      IF(DABS(XCORR).GT.AMAX)AMAX=XCORR
 16   CONTINUE
*      If central object present let velocity dispersion increase as 1/r
      IF(LS(16))THEN
      X(I20+IST,1)=X(I20+IST,2)+X(I00+IST,1)-X(I00+IST,2)+DLOG(STEU3)
      X(I02+IST,1)=X(I02+IST,2)+X(I00+IST,1)-X(I00+IST,2)+DLOG(STEU3)
      ELSE
      X(I20+IST,1)=X(I20+IST,2)
      X(I02+IST,1)=X(I02+IST,2)
      END IF
*
      IF(JEN.LT.JENMAX)THEN
      JEN=JEN+1
      IF(LS(25))THEN
      PRINT*,' Initial Pressure Iteration '
      PRINT*,' ITER=',JEN-1,' AMAX=',AMAX
      PRINT*,' LGPCEN=',X(I20+IST,1),' LGPOUT=',X(I20+IST,NJ)
      END IF
      IF(DABS(AMAX).GT.EPS) GOTO 27
      END IF
C       Retry with higher outer pressure in case of failure
      PRINT*,' End of initial Pressure Iteration with ITER=',JEN-1,
     *       ' AMAX=',AMAX,' Comp.=',IST
      IF(LSIGMA)PRINT*,' WARNING There was a sigma inversion!! '
      IF(LSIGMA.AND.JEN.EQ.JENMAX)THEN
      PROUT(IST) = PROUT(IST) + 1.D0
      PTOUT(IST) = PROUT(IST)
      IATTE = IATTE + 1
      PRINT*,' Next attempt with PROUT=',PROUT(IST)
      IF(IATTE.LT.30)GOTO 2727
      END IF
C
      ELSE
C       Exact calculation of Sigma for Plummers model
      PRINT*,'*INITIAL MODEL: SIGMA OF PLUMMERS MODEL EXACTLY TAKEN'
      DO 890 I=2,NJ
      X(I20+IST,I)=X(I00+IST,I)-DLOG(6.D0)-
     *                    C12*DLOG(R(I)**2+RSC(IST)**2)
 890  X(I02+IST,I)=X(I20+IST,I)
      X(I20+IST,1)=X(I20+IST,2)
      X(I02+IST,1)=X(I02+IST,2)
C
      END IF
*       In case of King model distribute total into partial pressures
      IF(NOL(IST).GT.100)THEN
*
      DO 7776 I=1,NJ
      RHOT = 0.D0
      DO 7777 K=1,NCOMP
 7777 RHOT = RHOT + DEXP(X(I00+K,I))
*
      DO 7778 K=NCOMP,1,-1
      X(I20+K,I)=X(I20+1,I) + X(I00+K,I) - DLOG(RHOT)
 7778 X(I02+K,I)=X(I20+K,I)
 7776 CONTINUE
*
      RHOT = 0.D0
      DO 7779 K=1,NCOMP
 7779 RHOT = RHOT + DEXP(XOUT(K))
*
      DO 7770 K=NCOMP,1,-1
      PROUT(K)=PROUT(1) + XOUT(K) - DLOG(RHOT)
 7770 PTOUT(K)=PROUT(K)
*
      GOTO 1234
      END IF
*
 1000 CONTINUE
 1234 CONTINUE
C
C  ENERGY FLUXES
C
      RFAC=R(3)/R(2)
      DO 600 I=2,NJ
*
      IM=I-1
      XMTNEW=0.D0
      DO 35 J=1,NCOMP
      IF(I.EQ.2)XLMIM=DLOG(PI43)+X(I00+J,1)+3.D0*DLOG(R(2)/RFAC)
*
      IF(I.EQ.2)THEN
         XLMIM=DLOG(PI43)+X(I00+J,1)+3.D0*DLOG(R(2)/RFAC)
         DELTR=DLOG(RFAC)
         RAVER=C12*DLOG(R(I)*R(I)/RFAC)
         RM=R(I)/RFAC
         ELSE
         XLMIM=DLOG(XTMP(J))
         DELTR=DLOG(R(I)/R(IM))
         RAVER=C12*DLOG(R(I)*R(IM))
         RM=R(IM)
         END IF
*
         IMASS=0
         AMAX=0.D0
 133     IMASS=IMASS+1
         IF(IMASS.EQ.1)THEN
         EFAC=DEXP(3.D0*DLOG(RM)-XLMIM)
         VMASS=0.D0
         ELSE
         XMAV=C12*(XLMI+XLMIM)
         EFAC=DEXP(3.D0*RAVER-XMAV)
         VMASS=XLMI
         END IF
*
         XLMI=XLMIM+PI4*DEXP(X(I00+J,I))*DELTR*EFAC
         XMCORR=XLMI-VMASS
         IF(IMASS.GT.100)THEN
         PRINT*,' Failure of initial mass iteration at:'
         PRINT*,' J=',J,' I=',I
         XMCORR=0.D0
         END IF
         IF(DABS(XMCORR).GT.EPS)GOTO 133
*
         XLMIM=XLMI
*
         XTMP(J)=DEXP(XLMI)
*
 35   CONTINUE
*
      DO 1001 J=1,NCOMP
      CLAM=PI4*GRAV/CS(21)
      CRX=-DLOG(CS(60))-DLOG(XMIND(J))-DLOG(XCOUL(J))
C
      IP=I+1
      IM=I-1
      XLR=DLOG(R(I))
      IF(I.EQ.2)THEN
      RAV=C12*DLOG(R(I)*R(I)/STEU3)
      ELSE
      RAV=C12*DLOG(R(I)*R(IM))
      END IF
      IF(I.EQ.NJ)THEN
      RAVP=C12*DLOG(R(I)*R(I)*STEU1)
      ELSE
      RAVP=C12*DLOG(R(IP)*R(I))
      END IF
C
      DERAV=RAVP-RAV
C
      IF(I.EQ.NJ)THEN
      XAVR(J)=C12*(XOUT(J)+X(I00+J,I))
      PAVR(J)=C12*(PROUT(J)+X(I20+J,I))
      PAVT(J)=PAVR(J)
      DSIGR=(PROUT(J)-X(I20+J,I)-XOUT(J)+X(I00+J,I))/DERAV
      DSIGT=(PTOUT(J)-X(I02+J,I)-XOUT(J)+X(I00+J,I))/DERAV
      ELSE
      XAVR(J)=C12*(X(I00+J,IP)+X(I00+J,I))
      PAVR(J)=C12*(X(I20+J,IP)+X(I20+J,I))
      PAVT(J)=C12*(X(I02+J,IP)+X(I02+J,I))
      DSIGR=(X(I20+J,IP)-X(I20+J,I)-X(I00+J,IP)+X(I00+J,I))/DERAV
      DSIGT=(X(I02+J,IP)-X(I02+J,I)-X(I00+J,IP)+X(I00+J,I))/DERAV
      END IF
C
      XLAIR=CLAM*DEXP(XLR+C12*PAVR(J)-C12*XAVR(J)+CRX)
      XLAIT=CLAM*DEXP(XLR+C12*PAVT(J)-C12*XAVR(J)+CRX)
C
      CORR1=1.D0+XMHOLE/DEXP(X(IMR,I))
*
      X(I30+J,I)=-DSIGR/XLAIR/CORR1
*      LS(22) Two Flux model
*             Different Lambda for energy fluxes if XLFAC.NE.0.D0
      IF(LS(22))THEN
      X(I30+J,I)=X(I30+J,I)*(1.D0+XLFAC)
      X(I12+J,I)=-DSIGT/XLAIT/CORR1/(1.D0+XLFAC)
      ELSE
      X(I12+J,I)=X(I30+J,I)
      END IF
 1001 CONTINUE
 600  CONTINUE
C
	  DO 1002 J=1,NCOMP
          DO 228 I=2,NJ
          IF(X(I12+J,I).EQ.0.D0)X(I12+J,I)=X(I12+J,I+1)
          IF(X(I30+J,I).EQ.0.D0)X(I30+J,I)=X(I30+J,I+1)
          IF(LS(9))THEN
          X(I12+J,I)=1.D-30
          X(I30+J,I)=1.D-30
          END IF
 228      CONTINUE
C
          X(I30+J,1)=X(I30+J,2)/STEU3
          X(I12+J,1)=X(I12+J,2)/STEU3
C
 1002     CONTINUE
C
  229     CONTINUE
C**************************************************************
C      INITIAL VELOCITY PROFILE
C**************************************************************
      CLFAC=CLICHT/299725.D0/CUNU
C
      DO 1003 IST=1,NCOMP
C
      DO 20 J=1,NJ
      X(I10+IST,J)=-2.D0*X(I30+IST,J)
 20   CONTINUE
C
 1003 CONTINUE
C
      RETURN
C
      END
