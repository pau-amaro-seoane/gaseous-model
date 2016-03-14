      SUBROUTINE INFORM
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
         INCLUDE 'compar.f'
          DIMENSION M(NJO),XLAGR(10),RLAGR(10)
          DIMENSION EPOT(NCOMPO),EKIN(NCOMPO),ETERM(NCOMPO)
          DIMENSION ETOT(NCOMPO),VIR(NCOMPO)
          DIMENSION E0TOT(NCOMPO),EBINR(NCOMPO),EBINT(NCOMPO)
          DIMENSION EEDIFF(NCOMPO),ELOSSC(NCOMPO),EGRAV(NCOMPO)
          DIMENSION EBRTOT(NCOMPO),EBTTOT(NCOMPO),EBCHCK(NCOMPO)
          DIMENSION EEDTOT(NCOMPO),ELCTOT(NCOMPO),EGRTOT(NCOMPO)
          DIMENSION EESTOT(NCOMPO),EEKICK(NCOMPO),XCNEW(NCOMPO)
*         DIMENSION EETEST(NCOMPO),EETTOT(NCOMPO)
      COMMON/EBTEST/EBALA,EBEQ,EBEQ2,DEMASS,DEENER,DEMASS2
        DATA XLAGR/0.01D0,0.02D0,0.05D0,0.1D0,0.2D0,0.5D0,0.75D0,
     *             0.0D0,0.0D0,0.0D0/
        DATA NLAGR/7/
C
C INFORMATION ABOUT GLOBAL QUANTITIES CALCULATED HERE
C E.G.    THE VARIOUS FORMS OF ENERGIES
*     Logarithmic Variables
*
         LWRITE=MOD(NMOD,NCOL).EQ.0
         TINY = 1.D-30
*
*       Save initial total mass
      IF(NMOD.EQ.0)THEN
      XMTOT0 = 0.D0
      DO 355 K = 1,NCOMP
      XMTOTI(K) = XMTOT(K)
 355  XMTOT0 = XMTOT0 + XMTOT(K)
      END IF
*
*       Change of Tidal Boundary
*
      IF(LS(27))THEN
*       Save initial total mass and tidal radius for later use
      IF(NMOD.EQ.0)THEN
      RTID0 = RTIDAL
      PRINT*,' Start with RTIDAL=',RTIDAL,' Mstar=',XMTOT0,
     *  ' Mbin=',XMBTOT
      END IF
*       Update tidal radius assuming constant mean density of cluster
      ITRT = 0
 391  CONTINUE
      ITRT = ITRT + 1
      RTOLD = RTIDAL
*       New Fokker-Planck like treatment
      J = 1
      SIGR = DEXP((X(I20+J,NJ)-X(I00+J,NJ))*C12)
      SIGT = DEXP((X(I02+J,NJ)-X(I00+J,NJ))*C12)
      VESC = DSQRT(2.D0*(PHTID-PHI(NJ)-PHIBIN(NJ)))
      EE = C12*SIGR**2 + SIGT**2 - C12*VESC**2 + PHTID
      XXAV = DSQRT(GRAV*(XMTOT0+XMBTOT)/RTID0**3)/2.D0/PI
      IF(EE.GT.PHTID)THEN
      FFAC = XXESC*DSQRT(1.D0-(EE/DABS(PHTID))**3)*XXAV
      ELSE
      FFAC=0.D0
      END IF
      DLOGRT = C13*(DLOG(DEXP(X(IMR,NJ)+XMRBIN(NJ)) -
     *              DLOG(XMTOT0+XMBTOT)))
      DLOGYT = C13*(DLOG(DEXP(X(IMR,NJ)+XMRBIN(NJ)) -
     *              DLOG(XMTOT0+XMBTOT)))-FFAC
      RTIDAL = RTID0*DEXP(DLOGRT)
      RTIDAY = RTID0*DEXP(DLOGYT)
         PHTID = -GRAV*(XMTOT0+XMBTOT)/RTID0*(RTIDAL/RTID0)**2
*
        DO 390 I = NJ,2,-1
*       IRR = 0
        IMEM = I
        IF(R(I).LE.RTIDAL)GOTO 401
 390    CONTINUE
*
 401    NJOLD = NJ
        NJ = IMEM
        PHTID = -GRAV*XMTOT0/RTID0*(RTIDAL/RTID0)**2
*       PHTID = -GRAV*DEXP(X(IMR,NJ))/RTIDAL
*
        LWRIT2 = .TRUE.
      IF(LS(23).AND.(ITER.LT.ITMAX/2).AND.(MOD(IMOD,NCOL).NE.0).
     * AND.(IMOD.GT.6))LWRIT2=.FALSE.
        IF(LWRIT2.OR.NJ.NE.NJOLD)THEN
          PRINT*,' Tidal EE,PHTID,XXX=',EE,PHTID,
     *     DSQRT(1.D0-(EE/DABS(PHTID))**3)
          PRINT*,' Tidal R(NJ),RTIDAY=',R(NJ),RTIDAY
          PRINT*,' Tidal R(NJ),RTIDAL,RTID0=',R(NJ),RTIDAL,RTID0
          PRINT*,' Tidal M(NJ),MTOT0=',DEXP(X(IMR,NJ)),XMTOT0
          PRINT*,' Tidal MB(NJ),MBTOT=',XMRBIN(NJ),XMBTOT
          PRINT*,' Tidal Cutoff IMOD,ITRT,IRR,NJ,NJOLD=',
     *      IMOD,ITRT,IRR,NJ,NJOLD
          PRINT*,' Tidal Cutoff PHINJ, PHTID=',PHI(NJ),PHTID
          PRINT*,' Tidal Cutoff M(NJ,NJ-1)=',DEXP(X(IMR,NJ)),
     *    DEXP(X(IMR,NJ-1))
        END IF
*     LCORR=DABS((RTOLD-RTIDAL)/RTIDAL).GT.EPS
*     IF(ITRT.LT.100.AND.LCORR)GOTO 391
      END IF
*          Calculation of gravitational potential
*          FP(), FM() used as temporary storage
*          and update of mass 
*
*          Start at outer boundary
         PHI(NJ)=-GRAV*DEXP(X(IMR,NJ))/R(NJ)
         XNTOT = 0.D0
*
         DO 340 I=1,NJ-1
 340     PHI(I)=0.D0
*
         DO 350 J=1,NCOMP
         FM(J,1)=0.D0
         XMTOT(J)=0.D0
         DO 360 I=2,NJ
         DER3=R(I)**3-R(I-1)**3
         FM(J,I)=FM(J,I-1)+DEXP(X(I00+J,I))*DER3
 360     CONTINUE
*
         DER3=R(NJ)**3-R(NJ-1)**3
         RAV=(R(NJ)+R(NJ-1))/2.D0
         FP(J,NJ)=DEXP(X(I00+J,NJ))*DER3/RAV
         DO 370 I=NJ-1,2,-1
         DER3=R(I)**3-R(I-1)**3
         RAV=(R(I)+R(I-1))/2.D0
         FP(J,I)=FP(J,I+1)+DEXP(X(I00+J,I))*DER3/RAV
 370     CONTINUE
*
         PHI(1)=PHI(1)-PI43*GRAV*FP(J,2)
         RFAC=R(3)/R(2)
         IF(LS(19))PHI(1)=PHI(1)-GRAV*XMHOLE/R(2)*RFAC
         DO 380 I=2,NJ-1
         RAV=(R(I)+R(I-1))/2.D0
         PHI(I)=PHI(I)-PI43*GRAV*(FM(J,I)/R(I)+FP(J,I+1))
         IF(LS(19))PHI(I)=PHI(I)-GRAV*XMHOLE/R(I)
 380     CONTINUE
*          and update of mass
         XMTOT(J) = FM(J,NJ)*PI43
         XNTOT = XNTOT + XMTOT(J)/XMIND(J)
*
 350     CONTINUE
*          update Coulomb Logarithm
         DO 385 J=1,NCOMP
         XCNEW(J) = DLOG(GAMMA*XNTOT)
 385     CONTINUE
*
*       Remove shell if E > 0 (unbound even in isolated system)
      IF(LS(27))THEN
        DO 395 I = 1,NJ
        LREMOVE = .TRUE.
        DO 1000 J=1,NCOMP
      SIGR = DEXP((X(I20+J,I)-X(I00+J,I))*C12)
      SIGT = DEXP((X(I02+J,I)-X(I00+J,I))*C12)
      EFAC = C12*(SIGR**2+2.D0*SIGT**2)+PHI(I)+PHIBIN(I)
*
      IF(NMOD.GT.0)LREMOVE = LREMOVE.AND.(EFAC.GT.0.D0)
*
 1000   CONTINUE
*
      IF(LREMOVE) THEN
      PRINT*,' Tidal: Shell ',I,' with E>0 removed '
      IRR = 1
      IMEM = I - 1
      GOTO 401
      END IF
*
  395   CONTINUE
      END IF
*
*        Calculate cumulative escape masses and energies
*
      IF(LS(3))THEN
*
      DO 404 J=1,NCOMP
      DMESC(J)=XMASS(J)
      TMESC(J)=XMASS(J)*BREM
 404  CONTINUE
*
      END IF
*
*       Mass loss by escapers
      IF(LS(11))THEN
*
*
      DO 505 J=1,NCOMP
      DEM=0.D0
      DEE=0.D0
      DO 506 I=2,NJ
      HX(I00+J) = DEXP(X(I00+J,I))
      SIGR = DEXP((X(I20+J,I)-X(I00+J,I))*C12)
      SIGT = DEXP((X(I02+J,I)-X(I00+J,I))*C12)
      ROUT = RTIDAL
      IF(PHTID.GT.PHI(I)+PHIBIN(I))THEN
      VESC = DSQRT(2.D0*(PHTID-PHI(I)-PHIBIN(I)))
      ELSE
      VESC = DSQRT(2.D0*(PHTID-PHI(I-1)-PHIBIN(I-1)))
      END IF
      VESCR = VESC
*       Select Apocentre Criterion
      VESCT = VESC
      IF(LS(29))VESCT = VESC*ROUT/DSQRT(ROUT**2-R(I)**2)
*       Old Expressions
      XARG1 = C12*VESCR/SIGR
      XARG2 = C12*VESCT/SIGT
      XARG3 = VESCR/DSQRT(2.D0)/SIGR
      XARG4 = VESCT/DSQRT(2.D0)/SIGT
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
      TOUT = XALPHA*(ROUT-R(I))/VESC
      SIG=DSQRT(C13*(SIGR**2+2.D0*SIGT**2))
      TRX(J,J)=SIG**3/(HX(I00+J)*XMIND(J)*XCOUL(J)*CS(60))
*     EFAC = C12*(XBRFAC*SIGR**2+2.D0*XBTFAC*SIGT**2)+PHI(I)+PHIBIN(I)
*     TIN = (EFAC - PHTID)/EFAC*TRX(J,J)
      TIN = XBETA*TRX(J,J)
*     PF = 1.D0/(1.D0+TOUT/TIN)
      PF = 1.D0
*
      XK=1.D0+TOUT/PF/TIN
*      At initialization take full loss cones
*     IF(NMOD.EQ.0)Y(J,I)=1.D0
*      At initialization take zero filling factor
*     IF(NMOD.EQ.0)Y(J,I)=0.D0
*      At initialization take steady state filling factor
      IF(NMOD.EQ.0)Y(J,I)=1.D0/(1.D0+PF*TIN/TOUT)
*
*     DRHOR = Y(J,I)*PF*XFAC/TOUT/BREM
*     IF(DRHOR.GT.CORR)THEN
*     PRINT*,' Warning J,I=',J,I,' DRHOR=',DRHOR
*     END IF
*
*     DERR = Y(J,I)*PF*(XFAC+YRFAC)/TOUT/BREM
*     IF(DERR.GT.CORR)THEN
*     PRINT*,' Warning J,I=',J,I,' DERR=',DERR
*     END IF
*
*     DERT = Y(J,I)*PF*(XFAC+YTFAC)/TOUT/BREM
*     IF(DERT.GT.CORR)THEN
*     PRINT*,' Warning J,I=',J,I,' DERT=',DERT
*     END IF
*       Determine Energy balance
      XCUT = 1.D0/(1.D0+(R(I)/R(NJ-2))**4)
*     XCUT = 1.D0
      DVOL=PI43*(R(I)**3-R(I-1)**3)
      DMLOC = Y(J,I)*PF*XFAC*DEXP(X(I00+J,I))/TOUT*DVOL/BREM*XCUT
      DELOC = Y(J,I)*PF*C12*(XRFAC*DEXP(X(I20+J,I))+
     *       2.D0*XTFAC*DEXP(X(I02+J,I)))/TOUT*DVOL/BREM*XCUT
      DEM=DEM + DMLOC
      DEE=DEE + DELOC
*       Energy change of self-gravitating system is 0.5*phi*dm !
      XEERR = XEERR - DELOC - C12*DMLOC*(PHI(I)+2.D0*PHIBIN(I))
      XMERR = XMERR - DMLOC
*       Determine new filling factor
          IF(NMOD.GT.0)THEN
          XEXP=DEXP(-PF*XK/TOUT/BREM)
          Y(J,I)=Y(J,I)*XEXP+TOUT/PF/TIN/XK*(1.D0-XEXP)
          IF(Y(J,I).GT.1.D0)Y(J,I)=1.D0
          IF(Y(J,I).LE.0.D0)Y(J,I)=1.D-30
          END IF
*     Y(J,I) = 6.D-3
*
 506  CONTINUE
      Y(J,1) = Y(J,2)
*
      TMESC(J)=DEM
      TEESC(J)=DEE
      DMESC(J)=DMESC(J)+TMESC(J)
      DEESC(J)=DEESC(J)+TEESC(J)
 505  CONTINUE
      END IF
*
      IF(LS(11).OR.LS(3))THEN
      IF(LWRITE)THEN
      IF(TIME.GT.0.D0)THEN
      TTT=TIME
      ELSE
      TTT=1.D30
      END IF
      IF(LS(11))THEN
      WRITE(6,281)(TMESC(J)/XMIND(J),J=1,NCOMP),(TMESC(J),J=1,NCOMP),
     *            (TEESC(J),J=1,NCOMP),(TMESC(J)*BREM,J=1,NCOMP),
     *            (TEESC(J)*BREM,J=1,NCOMP)
      WRITE(6,284)(DMESC(J)/XMIND(J),J=1,NCOMP),(DMESC(J),J=1,NCOMP),
     *            (DEESC(J),J=1,NCOMP),(DMESC(J)/TTT,J=1,NCOMP),
     *            (DEESC(J)/TTT,J=1,NCOMP)
      END IF
      END IF
 281  FORMAT(1X,' Esc-Info       DN,DM,DE,DMDT,DEDT=',1P,
     *      8(1X,D12.5),3(/,16X,8(1X,D12.5)))
 284  FORMAT(1X,' Esc-Info total DN,DM,DE,DMDT,DEDT=',1P,
     *      8(1X,D12.5),3(/,16X,8(1X,D12.5)))
*
      END IF
*
*      Start to calculate and eventually print out
*      global energy balance
*
*-------------------------------------------------------------------------
*    Notation:
*
*     Following variables contain space integrated values for
*     actual timestep and each component:
*
*     EPOT(J),EKIN(J),ETERM(J),ETOT(J)   actual energy values for comp. J
*     EEDIFF(J),ELOSSC(J) energies connection with accretion of black hole
*     EGRAV(J) energy lost by gravitational radiation near black hole
*     EBINR(J),EBINT(J) radial/tangential energies generated by binaries
*     DEESC(J): energy lost by escapers
*
*     Following variables contain sum over components of the above quantities:
*     TEPOT, TEKIN, TETERM, TETOT
*     TEDIFF, TLOSSC, TEGRAV, TEBINR, TEBINT, TEESCT
*
*     Following variables contain time integrated values for comp. J:
*     EEDTOT(J),ELCTOT(J),EGRTOT(J),EBRTOT(J),EBTTOT(J),EESTOT(J)
*
*     Following variables containg values summed over components and
*     integrated over time:
*     EEDT,ELCT,EEGT,EBTR,EBTT,EESC
*
*     T0POT, T0KIN, etc. same as TEPOT, TEKIN etc. but values
*      at the beginning of the run (stored for file in vector AEI)
*-------------------------------------------------------------------------
*
      IF(LWRITE)WRITE(6,29)
*
 29   FORMAT(1X,' GLOBAL INFO    NCOMP',
     *   '   TIME=    ETOT(J)= ETERM(J)=  EKIN(J)=  EPOT(J)=',
     *   '   VIR(J)=')
*
         TEKIN=0.D0
         TETERM=0.D0
         TEPOT=0.D0
         TETOT=0.D0
         TEBINR=0.D0
         TEBINT=0.D0
         TEESCT=0.D0
         TEDIFF=0.D0
         TLOSSC=0.D0
         TEGRAV=0.D0
         TETEST=0.D0
	 C43=4.D0/3.D0
	 C116=11.D0/6.D0
	 C76=7.D0/6.D0
         C52=2.5D0
         GF5=GRAV**5*27.D0*DSQRT(3.D0)
*
         DO 20 J=1,NCOMP
         EPOT(J)=0.D0
         EKIN(J)=0.D0
         ETERM(J)=0.D0
         EBINR(J)=0.D0
         EBINT(J)=0.D0
*        EETEST(J)=0.D0
         EBCHCK(J)=0.D0
         EEDIFF(J)=0.D0
         ELOSSC(J)=0.D0
         EGRAV(J)=0.D0
         XBINY(J)=1.D0
* 
        IF(LS(2))THEN
*
        DYA = 0.D0
        RHOT = 0.D0
*
        DO 70 K=1,NCOMP
        DYA = DYA + DEXP(X(I20+K,1))+2.D0*DEXP(X(I02+K,1))
  70    RHOT = RHOT + DEXP(X(I00+K,1))
*   Take into account binaries for core radius
        IF(LS(5))RHOT = RHOT + RHOBAV
        DYA = DYA/RHOT
      RCORE=DSQRT(DYA/(4.D0*PI*GRAV*RHOT))
*
        END IF
C
         DO 10 I=2,NJ
*
         IM=I-1
*
         HX(I00+J)=DEXP(X(I00+J,I))
         HX(I20+J)=DEXP(X(I20+J,I))
         HX(I02+J)=DEXP(X(I02+J,I))
	 RHOC=0.D0
	 DO 504 K=1,NCOMP
 504     RHOC=RHOC+DEXP(X(I00+K,1))
C
         DER3=PI43*(R(I)**3-R(IM)**3)
C
         PHIAV=C12*(PHI(I)+2.D0*PHIBIN(I)+PHI(IM)+2.D0*PHIBIN(IM))
         EPOT(J)=EPOT(J)+C12*HX(I00+J)*PHIAV*DER3
         EKIN(J)=EKIN(J)+(X(I10+J,I)**2/2.D0)*HX(I00+J)*DER3
         ETERM(J)=ETERM(J)+(2.D0*HX(I02+J)+HX(I20+J))/2.D0*DER3
*
         IF(TIME.GT.TBINI-1.D0/BREM)THEN
*
         IF(LS(1))THEN
*
         DYB=0.D0
         DO 507 K=1,NCOMP
         AFAC(K)=DEXP(X(I02+K,I)-X(I20+K,I))
         DYB=DYB+XMIND(K)**C43*DEXP(X(I00+K,I))**C116/
     *     DEXP(X(I20+K,I))**C76/(1.D0+2.D0*AFAC(K))**C76
 507     CONTINUE
*
         EBCHCK(J)=EBCHCK(J)+XBIN*HX(I00+J)*GF5/XMIND(J)*DYB**3*
     *      DER3/BREM
         EBINR(J)=EBINR(J)+C13*XBIN*HX(I00+J)**(1.D0-AEXP)*RHOC**AEXP*
     *      GF5/XMIND(J)*DYB**3*DER3/BREM
         EBINT(J)=EBINT(J)+C23*XBIN*HX(I00+J)**(1.D0+AEXP)*RHOC**AEXP*
     *      GF5/XMIND(J)*DYB**3*DER3/BREM
         END IF
*
         IF(LS(2))THEN
         DYB=0.D0
         RHOT=0.D0
*
         DO 517 K=1,NCOMP
         AFAC(K)=DEXP(X(I02+K,I)-X(I20+K,I))
         DYB=DYB+XMIND(K)*DEXP(X(I00+K,I))**C52/
     *   DEXP(X(I20+K,I))**C32/(1.D0+2.D0*AFAC(K))**C32
         RHOT=RHOT+DEXP(X(I00+K,I))
 517     CONTINUE
*
         EBCHCK(J)=EBCHCK(J)+XBIN*HX(I00+J)*GF5*
     *      DYA*DYB**3/RHOT*DER3/BREM
         EBINR(J)=EBINR(J)+C13*XBIN*HX(I00+J)**(1.D0-AEXP)*GF5*
     *      RHOC**AEXP*DYA*DYB**3/RHOT*DER3/BREM
         EBINT(J)=EBINT(J)+C23*XBIN*HX(I00+J)**(1.D0+AEXP)*GF5*
     *      RHOC**AEXP*DYA*DYB**3/RHOT*DER3/BREM
         END IF
*
         END IF
*
*        Energy Check for black hole related energy sources/sinks
*        First: energy diffusion (source)
         IF(LS(19))THEN
         SIGR2=HX(I20+J)/HX(I00+J)
         SIGT2=HX(I02+J)/HX(I00+J)
         SIG=DSQRT(C13*(SIGR2+2.D0*SIGT2))
         TRX(J,J)=SIG**3/(HX(I00+J)*XMIND(J)*XCOUL(J)*CS(60))
         XLPT=C23*GRAV*XMHOLE/RTIDE(J)
         IF(I.EQ.2.AND.AEXP.EQ.0.D0)EEDIFF(J)=EEDIFF(J)+
     *    C32*HX(I00+J)*(XLPT-C13*(SIGR2+SIGT2))/
     *     TRX(J,J)*DER3/BREM
*
*        Distribute energy of size of energy diffusion in core
         IF(AEXP.GT.0.D0)THEN
         IF(I.EQ.2)EBCHCK(J)=C32*HX(I00+J)*
     *    (XLPT-C13*(SIGR2+SIGT2))/TRX(J,J)*DER3/BREM
*
         IF(R(I).GT.RGRAV)EEDIFF(J)=EEDIFF(J)+HX(I00+J)**AEXP*DER3/BREM
         END IF
*
*        Second: loss cone accretion (sink)
         IF(LS(14))THEN
*     do not use XALPHA here because it is used for tidal field.
         YALPHA = 1.D0
         RLC=RTIDE(J)
         IF(R(I).LT.1.D1*R(2))RLC=RLC*R(I)/1.D1/R(2)
C
*        Prepare HX for use in routine LOSSCO
          HX(I00+J)=X(I00+J,I)
          HX(I20+J)=X(I20+J,I)
          HX(I02+J)=X(I02+J,I)
          HX(I10+J)=X(I10+J,I)
*
        CALL LOSSCO(XMHOLE,RLC,J)
C
          TCR=R(I)/DEXP((HX(I20+J)-HX(I00+J))/2.D0)
*
          ERARG=DSQRT(YALPHA*4.D0*DOMEGA*TRX(J,J)/TCR)
          IF(ERARG.LT.1.D1)THEN
          EARG=ERARG
          PF=DERF(EARG)
          ELSE
          PF=1.D0
          END IF
      GTERM=DEXP(HX(I00+J))*YALPHA*Y(J,I)*PF/TCR
      DGEPSR=-GTERM*EPSR*SIGR2
      DGEPST=-GTERM*EPST*SIGT2
      ELOSSC(J)=ELOSSC(J)+C12*(DGEPSR+2.D0*DGEPST)*DER3/BREM
          END IF
*
*     Third: gravitational radiation energy sink
          IF(LS(20))THEN
          HX(I20+J)=DEXP(X(I20+J,I))
          HX(I02+J)=DEXP(X(I02+J,I))
      TGRAV=CS(25)*(R(I)/XMIND(J))**2*R(I)**2/
     *           XMHOLE/(1.D0+XMHOLE/XMIND(J))
      EGRAV(J)=EGRAV(J)+C12*(HX(I20+J)+2.D0*HX(I02+J))/TGRAV*DER3/BREM
          END IF
*
          END IF
C
 10     CONTINUE
*
*        Correction of total core distributed energy diffusion energy
         IF(LS(19).AND.AEXP.GT.0.D0)THEN
*
         XBINY(J)=XBIN*EBCHCK(J)/EEDIFF(J)
*
*          Energy distributed in core must be removed near tidal radius
         EEDIFF(J)=(1.D0-XBIN)*EBCHCK(J)
*
         DO 112 I=2,NJ
         IM=I-1
         DER3=PI43*(R(I)**3-R(IM)**3)
         HX(I00+J)=DEXP(X(I00+J,I))
*
         IF(R(I).GT.RGRAV)EEDIFF(J)=EEDIFF(J)+
     *          XBINY(J)*HX(I00+J)**AEXP*DER3/BREM
 112     CONTINUE
*
         END IF
*
*        Correction of total binary energy generation
         IF((LS(1).OR.LS(2)).AND.TIME.GT.TBINI-1.D0/BREM)THEN
*
         XBINY(J)=XBIN*EBCHCK(J)/(EBINR(J)+EBINT(J))
*
         EBINR(J)=0.D0
         EBINT(J)=0.D0
         DO 111 I=2,NJ
*
         IM=I-1
*
         HX(I00+J)=DEXP(X(I00+J,I))
         HX(I20+J)=DEXP(X(I20+J,I))
         HX(I02+J)=DEXP(X(I02+J,I))
*
         DER3=PI43*(R(I)**3-R(IM)**3)
*
	 IF(LS(1))THEN
*
         DYB=0.D0
         DO 508 K=1,NCOMP
         AFAC(K)=DEXP(X(I02+K,I)-X(I20+K,I))
         DYB=DYB+XMIND(K)**C43*DEXP(X(I00+K,I))**C116/
     *   DEXP(X(I20+K,I))**C76/(1.D0+2.D0*AFAC(K))**C76
 508     CONTINUE
C
         EBINR(J)=EBINR(J)+C13*XBINY(J)*HX(I00+J)**(1.D0-AEXP)*GF5*
     *      RHOC**AEXP/XMIND(J)*DYB**3*DER3/BREM
         EBINT(J)=EBINT(J)+C23*XBINY(J)*HX(I00+J)**(1.D0+AEXP)*GF5*
     *      RHOC**AEXP/XMIND(J)*DYB**3*DER3/BREM
         END IF
*
	 IF(LS(2))THEN
*
         DYB=0.D0
         RHOT=0.D0
*
         DO 518 K=1,NCOMP
         AFAC(K)=DEXP(X(I02+K,I)-X(I20+K,I))
         DYB=DYB+XMIND(K)*DEXP(X(I00+K,I))**C52/
     *   DEXP(X(I20+K,I))**C32/(1.D0+2.D0*AFAC(K))**C32
         RHOT=RHOT+DEXP(X(I00+K,I))
 518     CONTINUE
*
         EBINR(J)=EBINR(J)+C13*XBINY(J)*HX(I00+J)**(1.D0-AEXP)*GF5*
     *      RHOC**AEXP*DYA*DYB**3/RHOT*DER3/BREM
         EBINT(J)=EBINT(J)+C23*XBINY(J)*HX(I00+J)**(1.D0+AEXP)*GF5*
     *      RHOC**AEXP*DYA*DYB**3/RHOT*DER3/BREM
         END IF
*
 111   CONTINUE
         END IF
*
*         Sum energy generation by stochastic binaries
*
         IF(LS(3))THEN
*
         EBINR(J)=C13*(XHEAT(J)+XFEED(J))
         EBINT(J)=C23*(XHEAT(J)+XFEED(J))
*
         END IF
*
         ETOT(J)=ETERM(J)+EKIN(J)+EPOT(J)
         VIR(J)=EPOT(J)/ETERM(J) + 2.D0
C---------Sum components-----------------------
       TEKIN=EKIN(J)+TEKIN
       TETERM=ETERM(J)+TETERM
       TEPOT=EPOT(J)+TEPOT
       TETOT=ETOT(J)+TETOT
       TEBINR=EBINR(J)+TEBINR
       TEBINT=EBINT(J)+TEBINT
*      TETEST=EETEST(J)+TETEST
       TEESCT=DEESC(J)+TEESCT
       TEDIFF=EEDIFF(J)+TEDIFF
       TLOSSC=ELOSSC(J)+TLOSSC
       TEGRAV=EGRAV(J)+TEGRAV
C-------------------------------------
      IF(LWRITE)
     * WRITE(6,30)J,TIME,ETOT(J),ETERM(J),EKIN(J),EPOT(J),VIR(J)
*
 30   FORMAT(1X,' For comp. No. ',4X,I2,1P,6D10.2)
C
 20      CONTINUE
*
*     IF(IMOD.GT.0)THEN
*     EBALA = EBALA + EBINR(1) + EBINT(1)
*     EBALA = EBALA + EETEST(1)
*     EBSUM = EBSUM + EBEQ
*     XMBALA = XMBALA + DEMASS
*
*     IF(MOD(NMOD,NCOL).EQ.0)THEN
*     WRITE(6,*)' This timestep bal/prog=',EBINR(1)+EBINT(1),EBEQ,EBEQ2
*     WRITE(6,*)' This timestep bal/prog=',EETEST(1),EBEQ
*     WRITE(6,*)' Total         bal/prog=',EBALA,EBSUM
*     WRITE(6,*)' Total             mass=',DEMASS,XMBALA,IBIN*2.D-4
*     CALL FLUSH(6)
*     END IF
*     END IF
C---------Sum total energy lost/gained by binaries 
C         and black hole (time integration)-------
      DO 701 J=1,NCOMP
      EBRTOT(J)=EBRTOT(J)+EBINR(J)
      EBTTOT(J)=EBTTOT(J)+EBINT(J)
      EEKICK(J)=EEKICK(J)+XKICK(J)
*     EETTOT(J)=EETTOT(J)+EETEST(J)
      EEDTOT(J)=EEDTOT(J)+EEDIFF(J)
      ELCTOT(J)=ELCTOT(J)+ELOSSC(J)
      EESTOT(J)=EESTOT(J)+DEESC(J)
      EGRTOT(J)=EGRTOT(J)+EGRAV(J)
 701  CONTINUE
C
C   Initial Values of Ekin, Etherm, Epot, Etot stored and saved;
C   If first run calls Inform (IMOD=0 and NMOD=0) the values are
C   saved in AEI-Vector, which is itself saved on file (PRINT)
C
C   If a restart calls Inform (IMOD=0 but NMOD.NE.0) the values
C   are restored from the AEI-vector
C---------------------------------------------------------------------
       LZERO = LS(5).AND.NMOD.EQ.1
       IF(IMOD.EQ.0.OR.LZERO)THEN
*       For primordial binaries wait one model to set energy t=0 values
          IF(NMOD.EQ.0.OR.LZERO)THEN
          T0KIN=TEKIN
          T0TERM=TETERM
          T0POT=TEPOT
          T0TOT=TETOT
          XMTOT0=DEXP(X(IMR,NJ))
          AEI(5)=TEKIN
          AEI(6)=TETERM
          AEI(7)=TEPOT
          AEI(8)=TETOT
          
*
          ISTA = 100
          DO 501 J=1,NCOMP
          EBINR(J)=0.D0
          EBINT(J)=0.D0
          EEDIFF(J)=0.D0
          ELOSSC(J)=0.D0
          E0TOT(J)=ETOT(J)
          EBRTOT(J)=0.D0
          EBTTOT(J)=0.D0
          EEKICK(J)=0.D0
          EEDTOT(J)=0.D0
          ELCTOT(J)=0.D0
          EESTOT(J)=0.D0
          EGRTOT(J)=0.D0
          AEI(ISTA+J)=ETOT(J)
 501      CONTINUE
C
          ELSE
C
          T0KIN=AEI(5)
          T0TERM=AEI(6)
          T0POT=AEI(7)
          T0TOT=AEI(8)
          ISTA=100
         DO 502 J=1,NCOMP
          E0TOT(J)=AEI(ISTA+J)
          EEKICK(J)=AEI(ISTA+NCOMP+J)
          EBRTOT(J)=AEI(ISTA+2*NCOMP+J)
          EBTTOT(J)=AEI(ISTA+3*NCOMP+J)
          EEDTOT(J)=AEI(ISTA+4*NCOMP+J)
          ELCTOT(J)=AEI(ISTA+5*NCOMP+J)
          EESTOT(J)=AEI(ISTA+6*NCOMP+J)
          EGRTOT(J)=AEI(ISTA+7*NCOMP+J)
          XINDHS(J)=AEI(ISTA+8*NCOMP+J)
          XINDHB(J)=AEI(ISTA+9*NCOMP+J)
 502     CONTINUE          
*
         END IF
*
        END IF
*
*         Store binary/hole total energy generation for save on file
      ISTA = 100 + NCOMP
      DO 503 J=1,NCOMP
        AEI(ISTA+J)=EEKICK(J)
        AEI(ISTA+NCOMP+J)=EBRTOT(J)
        AEI(ISTA+2*NCOMP+J)=EBTTOT(J)
        AEI(ISTA+3*NCOMP+J)=EEDTOT(J)
        AEI(ISTA+4*NCOMP+J)=ELCTOT(J)
        AEI(ISTA+5*NCOMP+J)=EESTOT(J)
        AEI(ISTA+6*NCOMP+J)=EGRTOT(J)
        AEI(ISTA+7*NCOMP+J)=XINDHS(J)
        AEI(ISTA+8*NCOMP+J)=XINDHB(J)
 503  CONTINUE
*
C---------Sum all components-----------------------------------
      EBTR=0.D0
      EBTT=0.D0
      EEKK=0.D0
*     EETT=0.D0
      EEDT=0.D0
      ELCT=0.D0
      EESC=0.D0
      EEGT=0.D0
      EERS=0.D0
      EESS=0.D0
      EESB=0.D0
      DO 602 J=1,NCOMP
      EBTR=EBTR+EBRTOT(J)
      EBTT=EBTT+EBTTOT(J)
      EEKK=EEKK+EEKICK(J)
      EESS=EESS+XINDHS(J)
      EESB=EESB+XINDHB(J)
*     EETT=EETT+EETTOT(J)
      EEDT=EEDT+EEDTOT(J)
      ELCT=ELCT+ELCTOT(J)
      EEGT=EEGT+EGRTOT(J)
      EERS=EERS+EREM(J)
*       Escaper energies calculated in subroutine PRINTS
      EESC=EESC+EESTOT(J)
 602  CONTINUE
C------------------------------------------------------------------
C       Output of distribution of energy forms and balance
      XXMASS = XXMASS + DEMASS
      XXENER = XXENER + DEENER
*
*     CALL CHECK
*
      IF(LWRITE)THEN
*     DETOT=TETOT-ETOLD
*     DEKIN=TEKIN-EKOLD
*     DETH=TETERM-ETHOLD
*     DEPOT=TEPOT-EPOLD
      WRITE(6,31)TIME,TETOT,TETERM,TEKIN,TEPOT
*     ETOLD=TETOT
*     EKOLD=TEKIN
*     ETHOLD=TETERM
*     EPOLD=TEPOT
*
 31   FORMAT(1X,' Sum components',6X,1P,5D10.2)
      WRITE(6,32)TIME,TETOT-T0TOT-EBTR-EBTT-EEDT-ELCT-EESC-EEGT-XEERR,
     *    TETERM-T0TERM,TEKIN-T0KIN,TEPOT-T0POT,
     *    DEXP(X(IMR,NJ))-XMTOT0-XMERR
      WRITE(6,33)TIME,EBTR+EBTT,EEKK,DPTOT,XEERR,EESS,EESB
      WRITE(6,331)TIME,EEDT+ELCT,EESC,EEGT
*     WRITE(6,*)' Direct computation XXMASS,XXENER=',XXMASS,XXENER
*     WRITE(6,*)' Indirect comp binsto XMERR,XEERR=',XMERR,XEERR
*
*     PRINT*,' ------Total energy changes------'
*     PRINT*,' Total: new,old,diff=',TETOT,ETOLD,DETOT
*     PRINT*,' Kin  : new,old,diff=',TEKIN,EKOLD,DEKIN
*     PRINT*,' Therm: new,old,diff=',TETERM,ETHOLD,DETH
*     PRINT*,' Pot  : new,old,diff=',TEPOT,EPOLD,DEPOT
*     PRINT*,'-----------------------------------'
      END IF
 32   FORMAT(1X,' Balance to T=0:',5X,1P,5D10.2,' energy   ',
     *    D10.2,' mass ')
 33   FORMAT(1X,' t,Ebin/kick/dptot/xeerr/iheats/iheatb',2X,1P,7D12.4)
 331  FORMAT(1X,' t,Ehole/esc/grav=',2X,1P,4D12.4)
C---------------------------------------------------------------------
*
      IF(LWRITE)THEN
*
      IF(LS(19))THEN
      WRITE(6,291)
*
 291  FORMAT(/,1X,' GL. INFO HOLE  NCOMP',
     *   ' EDIFF(J)= ELOSS(J)= EGRAV(J)=')
*
      DO 1001 J=1,NCOMP
      WRITE(6,301)J,EEDIFF(J),ELOSSC(J),EGRAV(J)
 1001 CONTINUE
      WRITE(6,321)EEDT,ELCT,EEGT
*
      END IF
*
      IF(LS(1).OR.LS(2).OR.LS(3).OR.LS(11))THEN
      WRITE(6,292)
*
 292  FORMAT(/,1X,' GL. INFO BINES NCOMP',
     *   ' EBINR(J)= EBINT(J)= EBTOT(J)=  EREM(J)=  EESC(J)=')
*
      DO 1002 J=1,NCOMP
      WRITE(6,301)J,EBINR(J),EBINT(J),EBINR(J)+EBINT(J),
*     WRITE(6,301)J,EETEST(J),0.D0,0.D0,
     *    EREM(J),DEESC(J)
 1002 CONTINUE
      WRITE(6,321)EBTR,EBTT,EBTR+EBTT,EERS,EESC
*     WRITE(6,321)EETT,0.D0,0.D0,EERS,EESC
*
      END IF
*
      END IF
*
 301  FORMAT(1X,' For comp. No.',5X,I2,1P,5D10.2)
 321  FORMAT(1X,' Time integral:',6X,1P,5D10.2)
*---------------------------------------------------------------------
*         CS(200) used to store relative energy error
      IF(NMOD.GT.0)THEN
      CS(200)=(TETOT-T0TOT-EBTR-EBTT-EEDT-ELCT-EESC-EEGT)/DABS(TETOT)
      IF(DABS(CS(200)).LT.EPS)CS(200)=0.D0
      ELSE
      CS(200)=0.D0
      END IF
*
*        Binary Energy in case of primordial binaries
*
      IF(LS(5))THEN
      EKINB = 0.D0
      EPOTB = 0.D0
*
      DO 1111 IB=1,IBIN
      XMBIN = BODY1(IB) + BODY2(IB)
      I = ISH(IB)
      EKINB = EKINB + XMBIN*(VR(IB)**2+VT(IB)**2)
      EPOTB = EPOTB + XMBIN*(PHIBIN(I) + PHI(I))
 1111 CONTINUE
      EKINB = C12*EKINB
      EPOTB = C12*EPOTB
*       To check total balance remember to take out factor 2 from PHIAV
*     WRITE(6,293)
*293  FORMAT(/,1X,' GLOBAL BINARY ENERGY',
*    *   '  EKINBIN= EPOTBIN= ETOTBIN= ')
*     WRITE(6,402)EKINB,EPOTB,EKINB+EPOTB
*     WRITE(6,403)EKINB+TETERM,EPOTB+TEPOT,EKINB+EPOTB+TETOT
*402  FORMAT(1X,' For binaries ',7X,1P,3D10.2)
*403  FORMAT(1X,' Summed up tot',7X,1P,3D10.2)
*
      END IF
*        Further lines only if output is due
      IF(.NOT.LWRITE)RETURN
*
*        Lagrangian radii
*
      IF(NMOD.EQ.0)THEN
      SIGC2 = DEXP(X(I20+1,1)-X(I00+1,1))
*   Take into account binaries for core radius
        IF(LS(5))RHOT = DEXP(X(I00+1,1)) + RHOBAV
      RCORE=DSQRT(3.D0*SIGC2/(4.D0*PI*GRAV*RHOT))
      END IF
*
      DO 700 J = 1,NLAGR
      IP = 1
 705  IP = IP+1
*
      IF(X(IMR,IP).LE.DLOG(XLAGR(J))+X(IMR,NJ)) GOTO 705
      I = IP-1
*        Linear interpolation in r**3
      RIP3 = R(IP)**3
      RI3 = R(I)**3
      ALPH = (RIP3 - RI3)/(DEXP(X(IMR,IP))-DEXP(X(IMR,I)))
      RLAGR(J) = RI3 + ALPH*(XLAGR(J)*DEXP(X(IMR,NJ)) - DEXP(X(IMR,I)))
      RLAGR(J) = RLAGR(J)**C13
*        Half mass average initial relaxation time
      IF(J.EQ.6.AND.IZ00.EQ.0)THEN
      IZ00=1
      SIG=0.D0
      RHO=0.D0
      XMC=0.D0
      XNTOT=0.D0
      DO 800 K=1,NCOMP
      RHO=RHO+DEXP(X(I00+K,I))
      SIG=SIG+DSQRT(C13*DEXP(X(I20+K,I)-X(I00+K,I))*
     *   (1.D0+2.D0*DEXP(X(I02+K,I)-X(I20+K,I))))
      XNTOT=XNTOT+XMTOT(K)/XMIND(K)
 800  XMC=XMC+XMIND(K)
      RHO=RHO/DBLE(NCOMP)
      SIG=SIG/DBLE(NCOMP)
      XMC=XMC/DBLE(NCOMP)
      TRH=SIG**3/(RHO*XMC*CS(60)*DLOG(GAMMA*XNTOT))
      TSPITZ=0.138D0/DLOG(GAMMA*XNTOT)*
     *        DSQRT(XNTOT*RLAGR(J)**3/XMC/GRAV)
      PRINT*,' TRH=',TRH,' TSPITZ=',TSPITZ
      END IF
*
 700  CONTINUE     
*
      RM2 = RLAGR(6)
*        End Lagrangian radii
C
      NINT=2*NWRITE
      IRUN=0
      DO 8880 I=1,NWRITE
      IRUN=IRUN+1
      M(IRUN)=I
 8880 CONTINUE
*
      DO 8890 I=NWRITE+1,NJ-NINT,NINT
      IRUN=IRUN+1
      M(IRUN)=I
 8890 CONTINUE
*
      DO 8990 I=NJ-NINT+1,NJ
      IRUN=IRUN+1
      M(IRUN)=I
 8990 CONTINUE
*
      MMAX=IRUN
C**********************************************************************
         DO 105 K=1,NCOMP
C
      SUM1=DEXP(X(I20+K,1)-X(I00+K,1))*
     *   (1.D0+2.D0*DEXP(X(I02+K,1)-X(I20+K,1)))
*   Take into account binaries for core radius
        IF(LS(5))RHOT = DEXP(X(I00+K,1)) + RHOBAV
      RCORE=DSQRT(SUM1/(4.D0*PI*GRAV*RHOT))
C
         WRITE(6,106)K
         DO 107 JJ=1,MMAX
         J=M(JJ)
         IMM=J-NINT
         IM=J-1
         IF(IM.EQ.0)IM=1
         I=J
         IP=J+1
         HX(I00+K)=DEXP(X(I00+K,I))
         RHOS=HX(I00+K)/UNRH
         XFN=HX(I00+K)/XMIND(K)*PI43*RIND(K)**3
         SIGR=DEXP((X(I20+K,I)-X(I00+K,I))/2.D0)
         SIGT=DEXP((X(I02+K,I)-X(I00+K,I))/2.D0)
         ANI=2.D0-2.D0*DEXP(X(I02+K,I)-X(I20+K,I))
         SIG=DSQRT(C13*(SIGR**2+2.D0*SIGT**2))
         TRX(K,K)=SIG**3/(HX(I00+K)*XMIND(K)*XCOUL(K)*CS(60))
         TDY=DSQRT(3.D0*PI/32.D0/GRAV/HX(I00+K)) 
      VSIG2=SIGR**2+2.D0*SIGT**2
      PMAX=DSQRT(PI*VSIG2/3.D0/GRAV/HX(I00+K))
      PMIN=2.D0*GRAV*XMIND(K)/VSIG2
      CLF=DLOG(PMAX/PMIN)
      CLOLD=DLOG(0.11D0*XMTOT(K)/XMIND(K))
      DLOGP=X(I20+K,IP)-X(I20+K,I)
      DLOGR=X(I00+K,IP)-X(I00+K,I)
      IF(DLOGP.EQ.0.D0)THEN
      QUANT=1.D60
      ELSE
      QUANT=DLOGR/DLOGP
      END IF
C
      DVOL=PI43*(R(I)**3-R(IM)**3)
C      Divide by sqrt(6)**3 to scale from rcore to r0
      XF=DSQRT(216.D0)
      PHII=PHI(I)+PHIBIN(I)
      PHI1=PHI(1)+PHIBIN(1)
      DEM=CS(22)*XMIND(K)*XMTOT(K)**2/RCORE**3/XF/BREM*
     *      DVOL*1.D1**(8.5D0*DLOG10(DABS(PHII/PHI1)))
      DEE=CS(24)*XMIND(K)*XMTOT(K)**2/RCORE**3/XF/BREM*
     *      DVOL*1.D1**(9.5D0*DLOG10(DABS(PHII/PHI1)))
C
         WRITE(6,108)J,RHOS,SIGR,ANI,CLF/CLOLD,XCNEW(K)/XCOUL(K),
     *         XMIND(K),TDY,TRX(K,K),DEM,DEE,PHI(I),PHIBIN(I)
  107    CONTINUE
C
 106  FORMAT(/,' ** Information: Component Number ',I3,' (stars)**',
     */,1X,' **I= ',' RHO=     ',' SIGMA=   ',' ANISO=   ',
     * ' CLF/OLD= ',' CLNEW/OLD',' XMIND=   ',' TDYN=    ',
     * ' TRX=     ',' DELTM=   ',' DELTE=   ',' PHI=     ',
     * ' PHIBIN=')
 108     FORMAT(1X,I5,1P,13(1X,D9.2))
C
 105     CONTINUE
C*********************************************************************
C        Escaping Stars information
         IF(LS(11))THEN
C
         DO 205 J=1,NCOMP
C
         WRITE(6,206)J
         DO 207 JJ=1,MMAX
         K=M(JJ)
         IMM=K-NINT
         IM=K-1
         IF(IM.EQ.0)IM=1
         I=K
         IP=K+1
         HX(I00+J)=DEXP(X(I00+J,I))
         HX(I20+J)=DEXP(X(I20+J,I))
         HX(I02+J)=DEXP(X(I02+J,I))
      VESC = DSQRT(2.D0*(PHTID-PHI(I)-PHIBIN(I)))
      ROUT = RTIDAL
      VESCR = VESC
*       Select Apocentre Criterion
      VESCT = VESC
      IF(LS(29))VESCT = VESC*ROUT/DSQRT(ROUT**2-R(I)**2)
      SIGR = DEXP((X(I20+J,I)-X(I00+J,I))*C12)
      SIGT = DEXP((X(I02+J,I)-X(I00+J,I))*C12)
      ANI=2.D0-2.D0*DEXP(X(I02+J,I)-X(I20+J,I))
      TOUT = XALPHA*(ROUT-R(I))/VESC
*
      XARG1 = C12*VESCR/SIGR
      XARG2 = C12*VESCT/SIGT
      XARG3 = VESCR/DSQRT(2.D0)/SIGR
      XARG4 = VESCT/DSQRT(2.D0)/SIGT
      XF1 = DERF(XARG1)
      XF2 = 1.D0 - DEXP(-XARG2**2)
      XF3 = DERF(XARG3)
      XF4 = 1.D0 - DEXP(-XARG4**2)
      XOLD = 1.D0 - C12*(XF1*XF4 + XF3*XF2)
      XA1 = C12*XARG1/DSQRT(PI)*DEXP(-XARG1**2)
      XA2 = C12*XARG2**2*DEXP(-XARG2**2)
      XA3 = C12*XARG3/DSQRT(PI)*DEXP(-XARG3**2)
      XA4 = C12*XARG4**2*DEXP(-XARG4**2)
      YRFAC = XA1*XF4 + XA3*XF2
      YTFAC = XF1*XA4 + XF3*XA2
      XROLD = XOLD + YRFAC
      XTOLD = XOLD + YTFAC
*       Get new expressions using Dawson's Integral Num.Rec. March 2002
      AA = XARG3
      BB = XARG4
      CALL DAFAC(AA,BB,XFAC,XRFAC,XTFAC)
*
      SIG=DSQRT(C13*(SIGR**2+2.D0*SIGT**2))
      TRX(J,J)=SIG**3/(HX(I00+J)*XMIND(J)*XCOUL(J)*CS(60))
*     EFAC = C12*(XBRFAC*SIGR**2+2.D0*XBTFAC*SIGT**2)+PHI(I)+PHIBIN(I)
*     TIN = (EFAC - PHTID)/EFAC*TRX(J,J)
      TIN = XBETA*TRX(J,J)
*     PF = 1.D0/(1.D0+TOUT/TIN)
      PF = 1.D0
      DVOL = PI43*(R(I)**3-R(I-1)**3)
      XCUT = 1.D0/(1.D0+(R(I)/R(NJ-2))**4)
*     XCUT = 1.D0
      ESCR = Y(J,I)*PF*XFAC/TOUT*XCUT
      DRH = ESCR/BREM
      TESC = 1.D0/(ESCR+TINY)
*
         WRITE(6,208)I,RHOS,SIGR,ANI,VESCR,VESCT,TOUT,TIN,TESC,
     *    TRX(J,J),Y(J,I),PF,XFAC,XOLD,XRFAC,XROLD,XTFAC,XTOLD,XCUT
  207    CONTINUE
C
 206  FORMAT(/,' ** Escaper Information: Component Number ',I3,
     */,1X,' **I= ',' RHO=     ',' SIGMA=   ',' ANISO=   ',
     * ' VESCR=   ',' VESCT=   ',' TOUT=    ',' TIN=     ',
     * ' TESC=    ',' TRX=     ',' Y=       ',' PF=      ',
     * ' XFAC=    ',' XOLD=    ',' XRFAC=   ',' XROLD=   ',
     * ' XTFAC=   ',' XTOLD=   ',' XCUT=    ')
 208     FORMAT(1X,I5,1P,18(1X,D9.2))
C
 205     CONTINUE
*
         END IF
C
                  RETURN
                  END
*
      SUBROUTINE CHECK
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
         INCLUDE 'compar.f'
          DIMENSION EPOT(NCOMPO),EKIN(NCOMPO),ETERM(NCOMPO)
          DIMENSION ETOT(NCOMPO),VETOT(NCOMPO)
          DIMENSION VEPOT(NCOMPO),VEKIN(NCOMPO),VETERM(NCOMPO)
          DIMENSION VPHI(NJO),VPHIB(NJO)
          DIMENSION XM(NJO),VXM(NJO),XE(NJO),VXE(NJO)
*
      COMMON/EBTEST/EBALA,EBEQ,EBEQ2,DEMASS,DEENER,DEMASS2
*          Calculation of gravitational potential
*          FP(), FM() used as temporary storage
*
*          Start at outer boundary
         PHI(NJ)=-GRAV*DEXP(X(IMR,NJ))/R(NJ)
*
         DO 340 I=1,NJ-1
 340     PHI(I)=0.D0
*
         DO 350 J=1,NCOMP
         FM(J,1)=0.D0
         DO 360 I=2,NJ
         DER3=R(I)**3-R(I-1)**3
         FM(J,I)=FM(J,I-1)+DEXP(X(I00+J,I))*DER3
 360     CONTINUE
*
         DER3=R(NJ)**3-R(NJ-1)**3
         RAV=(R(NJ)+R(NJ-1))/2.D0
         FP(J,NJ)=DEXP(X(I00+J,NJ))*DER3/RAV
         DO 370 I=NJ-1,2,-1
         DER3=R(I)**3-R(I-1)**3
         RAV=(R(I)+R(I-1))/2.D0
         FP(J,I)=FP(J,I+1)+DEXP(X(I00+J,I))*DER3/RAV
 370     CONTINUE
*
         PHI(1)=PHI(1)-PI43*GRAV*FP(J,2)
         RFAC=R(3)/R(2)
         IF(LS(19))PHI(1)=PHI(1)-GRAV*XMHOLE/R(2)*RFAC
         DO 380 I=2,NJ-1
         RAV=(R(I)+R(I-1))/2.D0
         PHI(I)=PHI(I)-PI43*GRAV*(FM(J,I)/R(I)+FP(J,I+1))
         IF(LS(19))PHI(I)=PHI(I)-GRAV*XMHOLE/R(I)
 380     CONTINUE
*
         EPOT(J)=0.D0
         ETERM(J)=0.D0
         EKIN(J)=0.D0
*
         DO 10 I=2,NJ
*
         IM=I-1
*
         HX(I00+J)=DEXP(X(I00+J,I))
         HX(I20+J)=DEXP(X(I20+J,I))
         HX(I02+J)=DEXP(X(I02+J,I))
         DER3=PI43*(R(I)**3-R(IM)**3)
C
         PHIAV=C12*(PHI(I)+2.D0*PHIBIN(I)+PHI(IM)+2.D0*PHIBIN(IM))
         EPOT(J)=EPOT(J)+C12*HX(I00+J)*PHIAV*DER3
         EKIN(J)=EKIN(J)+(X(I10+J,I)**2/2.D0)*HX(I00+J)*DER3
         ETERM(J)=ETERM(J)+(2.D0*HX(I02+J)+HX(I20+J))/2.D0*DER3
*
 10      CONTINUE
*
         ETOT(J)=EPOT(J)+EKIN(J)+ETERM(J)
*
 350     CONTINUE
*
*   Checks
*
         DO 650 J=1,NCOMP
*
         DT1 = 0.D0
         DT2 = 0.D0
         DT3 = 0.D0
         DT4 = 0.D0
         DT5 = 0.D0
         DT6 = 0.D0
         DTX1 = 0.D0
         DTX2 = 0.D0
         DTX3 = 0.D0
         DMM=0.D0
*
         DO 750 I=1,NJ
*
         IF(I.GT.1)THEN
         DVOL=PI43*(R(I)**3-R(I-1)**3)
         IM=I-1
         DER=R(I)-R(IM)
         ELSE
         DVOL=PI43*R(I)**3
         IM=I
         DER=1.D10
         END IF
*
         VRHO=DEXP(VX(I00+J,I))
         RHO=DEXP(X(I00+J,I))
         VSIG2=(DEXP(VX(I20+J,I))+2.D0*DEXP(VX(I02+J,I)))/VRHO
         SIG2=(DEXP(X(I20+J,I))+2.D0*DEXP(X(I02+J,I)))/RHO
         XM(I)=DVOL*DEXP(X(I00+J,I))
         XE(I)=SIG2
*
         DM=XM(I)-VXM(I)
         DE=XE(I)-VXE(I)
         DMM=DMM+DM
*
         DT1 = DT1 + C12*DM*VSIG2
         DT2 = DT2 + C12*DM*VPHI(I)
         DT3 = DT3 + C12*DM*VPHIB(I)
         DT4 = DT4 + C12*VXM(I)*DE
         DT5 = DT5 + C12*VXM(I)*(PHI(I)-VPHI(I))
         DT6 = DT6 + VXM(I)*(PHIBIN(I)-VPHIB(I))
*
         DPHIM=PHIBIN(IM)-VPHIB(IM)
         DPHI = PHIBIN(I)-VPHIB(I)
         IF(I.GT.1)THEN
         DMX = R(I)**2/GRAV*(DPHI-DPHIM)/DER
         ELSE
         DMX=0.D0
         END IF
*
         DTX1 = DTX1 + C12*DMX*C12*(VSIG2+SIG2)
         DTX2 = DTX2 + C12*DMX*C12*(PHI(I)+VPHI(I))
         DTX3 = DTX3 + DMX*(PHIBIN(I)+VPHIB(I))
         DTX4 = DTX4 + DMX*(PHIBIN(I)-VPHIB(I))
*
  750    CONTINUE
*
         DETOT=ETOT(J)-VETOT(J)
         DEPOT=EPOT(J)-VEPOT(J)
         DEKIN=EKIN(J)-VEKIN(J)
         DETERM=ETERM(J)-VETERM(J)
         DXERR=XEERR-VXERR
*
         IF(DABS(DT6).GT.0.D0)XXX = (DETOT-DT6)/DT6
         IF(IMOD.GT.0)THEN
         ACCDT = ACCDT + DT1+DT2+DT3+DT4+DT5+DT6
         ACCDE = ACCDE + DETOT
         ACCDT1 = ACCDT1 + DT1
         ACCDT2 = ACCDT2 + DT2
         ACCDT3 = ACCDT3 + DT3
         ACCDT4 = ACCDT4 + DT4
         ACCDT5 = ACCDT5 + DT5
         ACCDT6 = ACCDT6 + DT6
         ELSE
         ACCDT=0.D0
         ACCDE=0.D0
         ACCDT1=0.D0
         ACCDT2=0.D0
         ACCDT3=0.D0
         ACCDT4=0.D0
         ACCDT5=0.D0
         ACCDT6=0.D0
         END IF
*
         PRINT*,' ********************************************'
         PRINT*,' Routine Check IMOD=',IMOD,' TIME=',TIME
         PRINT*,' Kinetic DT1=',DT1,' DT4=',DT4,' Sum=',DT1+DT4
         PRINT*,' Potential DT2=',DT2,' DT3=',DT3,' DT5=',DT5,
     *      ' DT6=',DT6,' Sum=',DT2+DT3+DT5+DT6
         PRINT*,' Comb Single+Bin Pot (2+3),(5+6)=',DT2+DT3,DT5+DT6
         PRINT*,' Sum of all DT=',DT1+DT2+DT3+DT4+DT5+DT6
         PRINT*,' XEERR-Change=',DXERR,' DEMASS,DEMASS2,DMM=',
     *            DEMASS,DEMASS2,DMM
         PRINT*,' DETOT,POT,TERM,KIN=',DETOT,DEPOT,DETERM,DEKIN
         PRINT*,' ACCDT1 ',ACCDT1,' ACCDT2 ',ACCDT2,' ACCDT3 ',ACCDT3
         PRINT*,' ACCDT4 ',ACCDT4,' ACCDT5 ',ACCDT5,' ACCDT6 ',ACCDT6
         PRINT*,' ACCDT ',ACCDT
         PRINT*,' ********************************************'
*
*   Save for next time
*
         VETOT(J)=ETOT(J)
         VEPOT(J)=EPOT(J)
         VEKIN(J)=EKIN(J)
         VETERM(J)=ETERM(J)
         VXERR=XEERR
*
         DO 550 I=1,NJ
*
         VPHI(I)=PHI(I)
         VPHIB(I)=PHIBIN(I)
         VXM(I)=XM(I)
         VXE(I)=XE(I)
*
 550     CONTINUE
*
 650     CONTINUE
*
         RETURN
         END
