          SUBROUTINE GINIT
C**************************************************************
C
C Routine called by GIDINI; computes useful  quantities
C         for use in the equation subroutines (..EQ)
C----------------------------------------------------------------
       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
         INCLUDE 'compar.f'
      DIMENSION XTMP(NCOMPO)
C
C**********************************************************
C
      IF(LNJ)IP=NJ
      IF(L2) IM=2
      RFAC=R(3)/R(2)
C
      XLR=DLOG(R(I))
      IF(L2)THEN
      RFAC=R(IP)/R(I)
      DER=DLOG(RFAC)
      RAV=C12*DLOG(R(I)*R(I)/RFAC)
      ELSE
      DER=DLOG(R(I)/R(IM))
      RAV=C12*DLOG(R(I)*R(IM))
      END IF
      RAVI=DEXP(RAV)
      IF(LNJ)THEN
      RFAC=R(I)/R(IM)
      DERP=DLOG(RFAC)
      RAVP=C12*DLOG(R(I)*R(I)*RFAC)
      ELSE
      DERP=DLOG(R(IP)/R(I))
      RAVP=C12*DLOG(R(IP)*R(I))
      END IF
      DERAV=RAVP-RAV
C
      IF(ITER.EQ.1)THEN
*
      XMTNEW=0.D0
      DO 35 J=1,NCOMP
      IF(I.EQ.2)XLMIM=DLOG(PI43)+VX(I00+J,1)+3.D0*DLOG(R(2)/RFAC)
*
      IF(I.EQ.2)THEN
         XLMIM=DLOG(PI43)+VX(I00+J,1)+3.D0*DLOG(R(2)/RFAC)
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
         XLMI=XLMIM+PI4*DEXP(VX(I00+J,I))*DELTR*EFAC
         XMCORR=XLMI-VMASS
         IF(IMASS.GT.100)THEN
         PRINT*,' Failure of initial mass iteration at:'
         PRINT*,' J=',J,' I=',I,' XLMI,XLMIM,DELTR,EFAC=',
     *   XLMI,XLMIM,DELTR,EFAC
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
      END IF
*
      CRX=-DLOG(CS(60))-C32*DLOG(3.D0)
C
      DO 40 J=1,NCOMP
C
      IF(LNJ)THEN
*       Iteration to get exact outer pressure value, fixed for ITER>1
      IF(ITER.EQ.1)THEN
*       Start values
        XOUT(J)=2.D0*VX(I00+J,I)-VX(I00+J,IM)
        PROUT(J)=2.D0*VX(I20+J,I)-VX(I20+J,IM)
        PTOUT(J)=PROUT(J)-VX(I20+J,I)+VX(I02+J,I)
*
        IPITER=0
*
 150    IPITER=IPITER+1
*
        XAVR(J)=C12*(VX(I00+J,I)+XOUT(J))
        PAVR(J)=C12*(VX(I20+J,I)+PROUT(J))
        PAVT(J)=C12*(VX(I02+J,I)+PTOUT(J))
        AFACAV=DEXP(PAVT(J)-PAVR(J))
        XMFAC=GRAV*DEXP(VX(IMR,I)+XAVR(J)-XLR-PAVR(J))
*
        VPR=PROUT(J)
        VPT=PTOUT(J)
        VRHO=XOUT(J)
        PROUT(J)=VX(I20+J,I)-DERAV*(XMFAC+2.D0*(1.D0-AFACAV))
        PTOUT(J)=PROUT(J)-VX(I20+J,I)+VX(I02+J,I)
        XOUT(J)=PROUT(J)-VX(I20+J,I)+VX(I00+J,I)+DLOG(RFAC)
*
          IF(IPITER.GT.100)THEN
          PRINT*,' Failure of outer P-iteration at IPITER=',IPITER
          PRINT*,' VPR,PR=',VPR,PROUT(J)
          STOP
          END IF
*
      DPR=DABS(VPR-PROUT(J))
      DPT=DABS(VPT-PTOUT(J))
      DRHO=DABS(VRHO-XOUT(J))
      IF(LS(25))PRINT*,' Outer P-Iter=',IPITER,' DPR=',DPR,
     *                 ' DPT=',DPT,' DRHO=',DRHO
      IF(DPR.GT.EPS.OR.DPT.GT.EPS.OR.DRHO.GT.EPS)GOTO 150
*
      IF(IPITER.GT.50.OR.LS(25))PRINT*,
     *   ' Pressure Iterations IPITER=',IPITER
*
      END IF
*
      HXP(I00+J)=XOUT(J)
      HXP(I20+J)=PROUT(J)
      HXP(I02+J)=PTOUT(J)
      HXP(I10+J)=0.D0
      HXP(I30+J)=0.D0
      HXP(I12+J)=0.D0
*
      END IF
C
      UAV(J)=C12*(HX(I10+J)+HXM(I10+J))
      UAVP(J)=C12*(HXP(I10+J)+HX(I10+J))
      VRAV(J)=C12*(HX(I30+J)+HXM(I30+J))
      VTAV(J)=C12*(HX(I12+J)+HXM(I12+J))
      PAVR(J)=C12*(HXP(I20+J)+HX(I20+J))
      PAVRM(J)=C12*(HX(I20+J)+HXM(I20+J))
      PAVT(J)=C12*(HXP(I02+J)+HX(I02+J))
      PAVTM(J)=C12*(HX(I02+J)+HXM(I02+J))
      XAVR(J)=C12*(HXP(I00+J)+HX(I00+J))
      XAVRM(J)=C12*(HX(I00+J)+HXM(I00+J))
C
      AFAC(J)=DEXP(HX(I02+J)-HX(I20+J))
      ATFAC(J)=1.D0/AFAC(J)
*
 40   CONTINUE
*
      IF(LS(2).OR.LS(3))THEN
      DYA = 0.D0
      RHOT = 0.D0
      END IF
*
      DO 50 J=1,NCOMP
*       Relaxation rate for multicomponent - if J=K normal relax. time
      DO 60 K=1,NCOMP
      TRX(J,K)=((1.D0+2.D0*AFAC(J)+(1.D0+2.D0*AFAC(K))*
     *       DEXP(HX(I20+K)-HX(I20+J)+HX(I00+J)-HX(I00+K)))/2.D0)**C32*
     * DEXP(C32*HX(I20+J)-25.D-1*HX(I00+J)+CRX
     *                   -DLOG(XMIND(J))-DLOG(XCOUL(J)))
 60   CONTINUE
*
*        Calculate mass averaged central vel. dispersion
*        for use in Lee et al. 3 body binary energy generation
*        (here 3-D - factor 3 is taken into acc. in FORMB3)
*
      IF(LS(2).OR.LS(3))THEN
*
      DYA = DYA + DEXP(VX(I20+J,1))+2.D0*DEXP(VX(I02+J,1))
      RHOT = RHOT + DEXP(VX(I00+J,1))
*
      END IF
*
*        Calculate escape masses and energies
      IF(LS(11))THEN
*      
      VSIGR = DEXP((VX(I20+J,I)-VX(I00+J,I))*C12)
      VSIGT = DEXP((VX(I02+J,I)-VX(I00+J,I))*C12)
      VESC = DSQRT(2.D0*(PHTID-PHI(I)-PHIBIN(I)))
      ROUT = RTIDAL
      VESCR = VESC
*       Select Apocentre Criterion
      VESCT = VESC
      IF(LS(29))VESCT = VESC*ROUT/DSQRT(ROUT**2-R(I)**2)
      XARG1 = C12*VESCR/VSIGR
      XARG2 = C12*VESCT/VSIGT
      XARG3 = VESCR/DSQRT(2.D0)/VSIGR
      XARG4 = VESCT/DSQRT(2.D0)/VSIGT
*       Keep old approximate values for comparison
      XF1 = DERF(XARG1)
      XF2 = 1.D0 - DEXP(-XARG2**2)
      XF3 = DERF(XARG3)
      XF4 = 1.D0 - DEXP(-XARG4**2)
      XFAC = 1.D0 - C12*(XF1*XF4 + XF3*XF2)
      XA1 = C12*XARG1/DSQRT(PI)*DEXP(-XARG1**2)
      XA2 = C12*XARG2**2*DEXP(-XARG2**2)
      XA3 = C12*XARG3/DSQRT(PI)*DEXP(-XARG3**2)
      XA4 = C12*XARG4**2*DEXP(-XARG4**2)
      YRFAC = XA1*XF4 + XA3*XF2
      YTFAC = XF1*XA4 + XF3*XA2
      XRFAC = XFAC + YRFAC
      XTFAC = XFAC + YTFAC
*       Get new expressions using Dawson's Integral Num.Rec. March 2002
      AA = XARG3
      BB = XARG4
      CALL DAFAC(AA,BB,XFAC,XRFAC,XTFAC)
*
      IF(XTFAC.LT.0.D0)THEN
      PRINT*,' I,A,B,Xe,Xr,Xt=',I,AA,BB,XFAC,XRFAC,XTFAC
      END IF
      TOUT = XALPHA*(ROUT-R(I))/VESC
      VSIG=DSQRT(C13*(VSIGR**2+2.D0*VSIGT**2))
      VTRX=VSIG**3/(VX(I00+J,I)*XMIND(J)*XCOUL(J)*CS(60))
*
      TIN = XBETA*TRX(J,J)
*     PF = 1.D0/(1.D0+TOUT/TIN)
      PF = 1.D0
      XCUT = 1.D0/(1.D0+(R(I)/R(NJ-2))**4)
*     XCUT = 1.D0
      DRHESC(J)=Y(J,I)*PF*XFAC/TOUT*XCUT
      DERESC(J)=Y(J,I)*PF*XRFAC/TOUT*XCUT
*
       IF(LS(13))THEN
       DETESC(J)=DERESC(J)
       ELSE
       DETESC(J)=Y(J,I)*PF*XTFAC/TOUT*XCUT
       END IF
*       New Fokker-Planck like treatment
      EE = C12*VSIGR**2 + VSIGT**2 - C12*VESC**2 + PHTID
      XXAV = DSQRT(GRAV*(XMTOT0+XMBTOT)/RTID0**3)/2.D0/PI
*
      FFAC = XXESC*DSQRT(1.D0-(EE/DABS(PHTID))**3)*XXAV
      IF(EE.GT.PHTID)THEN
*
      IF(FFAC/BREM.GT.2.D0*CORR)THEN
*     IF(FFAC/BREM.GT.CORR)THEN
*     PRINT*,' FF-CORR: IMOD,ITER,I=',IMOD,ITER,I
*     PRINT*,' FF-CORR: FFAC/BREM=',FFAC/BREM
      FFAC = 2.D0*CORR
*     FFAC = CORR
      END IF
*
      DRHESC(J) = FFAC
      DERESC(J) = FFAC
      DETESC(J) = FFAC
      END IF
*
*     LPR=(MOD(I,10).EQ.0.OR.I.GT.NJ-10).AND.MOD(IMOD,NCOL).EQ.0
*     IF(LPR)THEN
*     PRINT*,' I,EE,FFAC=',I,EE,FFAC
*     PRINT*,' I,DR,P,Pt=',I,DRHESC(J),DERESC(J),DETESC(J)
*     PRINT*,'            A,B,Xe,Xr,Xt=',AA,BB,XFAC,XRFAC,XTFAC
*     END IF
*
      END IF
*       End tidal mass and energy loss
*       Upstream/downstream flags
      IF(ITER.LE.5)THEN
      FP(J,I)=0.D0
      IF (HX(I10+J).GT.0.D0) FP(J,I)=1.D0
      IF (LNJ) FP(J,I)=1.D0
      FM(J,I)=1.D0-FP(J,I)
*
      FP(J,IM)=0.D0
      IF(HXM(I10+J).GT.0.D0) FP(J,IM)=1.D0
      FM(J,IM)=1.D0-FP(J,IM)
          IF(L2)FP(J,IM)=0.D0
      END IF
*
 50   CONTINUE
*
      IF(LS(2).OR.LS(3))DYA=DYA/RHOT
*
          RETURN
C
          END
