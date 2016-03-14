      SUBROUTINE FORMBS
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
      INCLUDE 'compar.f'
      COMMON/EBTEST/EBALA,EBEQ,EBEQ2,DEMASS,DEENER,DEMASS2
*
*       Energy generation by three body binary formation and hardening
*       Take Lee et al. formula
*       Possible radial/tangential difference in energy generation
*         by using factor AEXP
*
      C52=2.5D0
      DYB=0.0D0
      RHOT=0.D0
      CBIN=DLOG(18.D0*DSQRT(3.D0))+5.D0*DLOG(GRAV)
*
      DVOL=PI43*(R(I)**3-R(I-1)**3)
*     Mass Distribution Cutoff
      IF(VX(IMR,I)-VX(IMR,NJ).GT.-1.609D0)THEN
      CFACM=2.D0/(1.D0+(R(I)/RIR20)**4)
      CFACE=2.D0/(1.D0+(R(I)/RIR20)**4)
      CFACF=2.D0/(1.D0+(R(I)/RIR20)**4)
      ELSE
      RIR20=R(I)
      CFACM=1.D0
      CFACE=1.D0
      CFACF=1.D0
      END IF
*     Cutoff for Feedback Energy
*      IF(VX(IMR,I)-VX(IMR,NJ).GT.-0.693D0)THEN
*      CFACE=2.D0/(1.D0+(R(I)/RHALF)**2)
*      ELSE
*      RHALF=R(I)
*      CFACE=1.D0
*      END IF
*     Energy Distribution Cutoff
*     IF(R(I).GT.2.D0*RCORE)THEN
*     CFACE=2.D0/(1.D0+(R(I)/2.D0/RCORE)**2)
*     CFACF=CFACE
*     ELSE
*     CFACE=1.D0
*     CFACF=1.D0
*     END IF
*
      DO 901 J = 1,NCOMP
*
      IF(XFBINY(J)+XSBINY(J)+XMASSY(J).EQ.0.D0)GOTO 901
*
      RHOT = RHOT + DEXP(VX(I00+J,I))
* 
      AFCC = DEXP(VX(I02+J,I)-VX(I20+J,I))
      CBTERM = DLOG(XMIND(J))
      DBTERM = DEXP(CBTERM+C52*VX(I00+J,I)-C32*VX(I20+J,I))/
     *            (1.D0+2.D0*AFCC)**C32
      DYB = DYB + DBTERM
*
 901  CONTINUE
*
      DO 1010 J = 1,NCOMP
*
      IF(XFBINY(J)+XSBINY(J)+XMASSY(J).EQ.0.D0)GOTO 1010
*
*      Take corrected XBIN-value
      CBR=CBIN+AEXP*VX(I00+J,1)
      CBT=CBIN-AEXP*VX(I00+J,1)
*
* Radial energy generation (DYA determined in GINIT central sigma)
*
      C3 = DEXP(CBR+(1.D0-AEXP)*VX(I00+J,I)-HX(I20+J))*DYA*DYB**3/RHOT*
     *    XSBINY(J)*CFACE
      CRHO3=DEXP(CBR+(1.D0-AEXP)*VX(I00+J,I)-HX(I00+J))*DYA*DYB**3
     *    /RHOT**3*XMASSY(J)*C32*CFACM
      CFEED=DEXP(CBR+(1.D0-AEXP)*VX(I00+J,I)-HX(I20+J))*DYA*DYB**3/RHOT*
     *    XFBINY(J)*CFACF
*
*     DEBIN = DEBIN + C12*C3*DEXP(HX(I20+J))*DVOL/BREM
*        IF(I.EQ.20.AND.MOD(IMOD,10).EQ.0)THEN
*        PRINT*,' CRHO3=',CRHO3,' XMASSY=',XMASSY(J)
*        PRINT*,' CHCK I,DYA,DYB,RHOT=',I,DYA,DYB,RHOT
*        PRINT*,' DETERM=',C32*C3*DEXP(HX(I20+J))/XSBINY(J)*XBIN*
*    *     DVOL/BREM,
*    *     ' E**CBR=',C32*DEXP(CBR)
*        PRINT*,' I,RHO=',I,DEXP(X(I00+J,I)),' DRHO=',
*    *      CRHO3/BREM*DEXP(X(I00+J,I)),
*    *     ' CFACM=',CFACM,' CFACE=',CFACE
*        PRINT*,'----------------------------------'
*        END IF

*
      G(I20+J) = G(I20+J) - C3 - CFEED
*
      D(I20+J,I20+J) = D(I20+J,I20+J)+C3+CFEED
*
*       Mass loss due to binary events
*     CRHOB=0.D0
*     DO 1005 IB=IBIN,1,-1
*     IF(IEV(IB).EQ.2.AND.DMBIN(IB,I).GT.0)THEN
*     CRHOB = CRHOB - DMBIN(IB,I)/DVOL*BREM
*     PRINT*,' I,IB=',I,IB,' RHOB=',CRHOB/BREM,
*    *  ' RHO=',DEXP(HX(I00+J)),' DMBIN=',DMBIN(IB,I)
*     ELSE
*     GOTO 1005
*     END IF
*1005 CONTINUE
*
*     IF(CRHOB.EQ.0.D0)GOTO 1003
*
*     CRHO3 = CRHOB/DEXP(HX(I00+J))
*
      G(I00+J) = G(I00+J) - CRHO3
      G(I20+J) = G(I20+J) - CRHO3
      G(I02+J) = G(I02+J) - CRHO3
*
      D(I00+J,I00+J) = D(I00+J,I00+J)+CRHO3
      D(I20+J,I00+J) = D(I20+J,I00+J)+CRHO3
      D(I02+J,I00+J) = D(I02+J,I00+J)+CRHO3
*
      DEMASS = DEMASS + CRHO3*DEXP(VX(I00+J,I))*DVOL/BREM
      DEMASS2=DEMASS2 + CRHO3*DEXP(HX(I00+J))*DVOL/BREM
      DEENER=DEENER+CRHO3*C12*(DEXP(HX(I20+J))+2.D0*DEXP(HX(I02+J)))*
     *           DVOL/BREM
*        IF(IMOD.EQ.20.AND.CRHO3.NE.0.D0)THEN
*        DERHO=CRHO3*DEXP(X(I00+J,I))/BREM
*        DEM=PI43*DERHO*DVOL
*      PRINT*,' I,DERHO,DEM=',I,DERHO,DEM
*        PRINT*,' RHO,DEMS,DEMBIN=',DEXP(X(I00+J,I)),DEXP(X(IMR,I))-
*    *   DEXP(X(IMR,I-1)),XMRBIN(I)-XMRBIN(I-1)
*        PRINT*,' I,DEMASS,DEENER=',I,DEMASS,DEENER
*        END IF
*
*1003 CONTINUE
*        Tangential energy generation
*
      IF(.NOT.LS(13))THEN
*
      CT3 = DEXP(CBT+(1.D0+AEXP)*VX(I00+J,I)-HX(I02+J))*DYA*DYB**3
     *     /RHOT*XSBINY(J)*CFACE
      CTFEED=DEXP(CBT+(1.D0+AEXP)*VX(I00+J,I)-HX(I02+J))*DYA*DYB**3
     *     /RHOT*XFBINY(J)*CFACF
*     DEBIN = DEBIN + CT3*DEXP(HX(I02+J))*DVOL/BREM
*
      G(I02+J) = G(I02+J) - CT3 - CTFEED
*
      D(I02+J,I02+J) = D(I02+J,I02+J)+CT3+CTFEED
*
      END IF
*
*       Diagnostic output
*
      IF(LS(25))THEN
      PRINT*,' FORMBS: I,J=',I,J,' C3 (RAD) = ',C3,
     *' CT3 (TANG) = ',CT3,' C3*RHOT=',C3*RHOT
            PRINT*,' FORMBS: J=',J,' CFEED (RAD) = ',CFEED,
     *' CTFEED (TANG) = ',CTFEED
      PRINT*,' XFBINY,CFACF=',XFBINY(J),CFACF
      END IF
*
 1010 CONTINUE
*
      RETURN
*
      END
