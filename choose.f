      SUBROUTINE CHOOSE
C******************Choose Timestep******************************
C
C     Routine called after each HENYEY-Iteration
C
C     Maintains:   i) timestep adjustment
C                 ii) time interpolation if wished by LS(6)
C                iii) movement of radial grids if wished by LS(26)
C                (note: iii) not tested in this program version)
C---------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
         INCLUDE 'compar.f'
C
          DTIME=1.D0/BREM
*
          IF(NMOD.EQ.0) GOTO 90
          IF(ITER.EQ.0)GOTO 4
*
          IM2=ITMIN+2
          IM3=ITMIN+3
C***************************************************************
C     TIMESTEP ADJUSTMENT governed by DMAX = max. rel. corr.
C                         of even quantities
C     If central black hole present: check also change of its mass
C
C     The timestep is adjusted such that DMAX comes close to
C                        CORR (Input parameter)
C                        Recommended: CORR=0.05
C***************************************************************
          DMAX=CS(17)
C     Check black hole correction
          IF(LS(19))THEN
              DMH=DABS(XMHOLE-VMHOLE)/XMHOLE
              IF(DMH.GT.DMAX)DMAX=DMH
          END IF
C
          LPOS=(DMAX.LT.CORR .AND. ITER.LT.IM2 )
C
          IF(DMAX.LT.CORR/1.D3)DTIME=DTIME*4.D0
          IF(DMAX.LT.CORR/5.D1)DTIME=DTIME*2.D0
          IF(DMAX.LT.CORR/1.D1)DTIME=DTIME*5.D0
c
          IF(LPOS) GOTO 2
 8        CONTINUE
C
          LNEG=((DMAX.GT.CORR).OR.(ITER.GE.IM2))
          IF(DMAX.GE.1.D1*CORR)DTIME=DTIME/1.D1
          IF((DMAX.LT.1.D1*CORR).AND.(DMAX.GE.5.D0*CORR))
     *   DTIME=DTIME/5.D0
          IF(LNEG) GOTO 6
C
          IF(ITER.LE.ITMIN+3) GO TO 4
C
  2   DTIME=DTIME*FACT
          GOTO 4
  6   DTIME=DTIME/FACT
C
  4   CONTINUE
C
C***********************************************************************
C          Time interpolation in x; take care of explicit/implicit part
C***********************************************************************
      IF(IMOD*ITER.EQ.0)THEN
      WRITE(6,991)(VX(I,1),I=1,NPOS)
      WRITE(6,992)(X(I,1),I=1,NPOS)
      END IF
C
      IF (LS(6))THEN
C
      XTFAC=(1.D0-XIMPL)/XIMPL
*
      DO 70 J=1,NJ
*      Do not extrapolate mass
      VX(1,J)=X(1,J)
      DO 71 K=2,NG
      XK=X(K,J)
      VXK=VX(K,J)
      VVXK=VVX(K,J)
      XXX=(XK-VXK)*BREM*DTIME
*     XXX=(XK-VXK)*(BREM*DTIME-XTFAC)+BREM*DTIME*XTFAC*(VXK-VVXK)
      VVX(K,J)=VXK
      VX(K,J)=XK
      X(K,J)=XK+XXX
 71   CONTINUE
 70       CONTINUE
C***********************************************************************
C          WITHOUT ANY EXTRAPOLATION IN TIME THE GLOBAL COSERVATION laws
C                         OF ENERGY AND MASS seem to BE BETTER FULFILLED
C  HOWEVER THE ITERATIONS NEEDED TO CONVERGE ARE LARGER (iter is larger)
C***********************************************************************
C
      IF(IMOD.EQ.0)THEN
      PRINT*,' AFTER EXTRAPOLATION:'
      WRITE(6,991)(VX(I,1),I=1,NPOS)
      WRITE(6,992)(X(I,1),I=1,NPOS)
      END IF
C*********************************************************************
C No extrapolation: only vx=x vvx=vx
C********************************************************************
      ELSE
          DO 101 J=1,NJ
          DO 101 KK=1,NG
          VVX(KK,J)=VX(KK,J)
 101      VX(KK,J)=X(KK,J)
      END IF
*
 90     CONTINUE
*         Compute minimum central relaxation time 
         I=3
         J=1
         SIG2=DEXP(X(I20+J,I)-X(I00+J,I))*
     *             (1.D0+2.D0*DEXP(X(I02+J,I)-X(I20+J,I)))
         TRXCEN=(SIG2/3.D0)**C32/
     *    CS(60)/XMIND(J)/XCOUL(J)/DEXP(X(I00+J,I))
*
         DO 567 J=2,NCOMP
         SIG2=DEXP(X(I20+J,I)-X(I00+J,I))*
     *             (1.D0+2.D0*DEXP(X(I02+J,I)-X(I20+J,I)))
         TRXC=(SIG2/3.D0)**C32/
     *    CS(60)/XMIND(J)/XCOUL(J)/DEXP(X(I00+J,I))
         IF(TRXC.LT.TRXCEN)TRXCEN=TRXC
 567     CONTINUE
*
*       Limit DT by central relaxation time and advance times
*
         IF(NMOD.GT.0)THEN
         IF(DTIME.GT.FTRX*TRXCEN)DTIME=FTRX*TRXCEN
*
*        limit DTIME to values < 1  *****  TEST *****
*         IF(DTIME.GT.5.D0)DTIME=5.D0 
*          IF(DTIME.GT.1.D0)DTIME=1.D0
*        IF(DTIME.GT.0.5D0)DTIME=0.5D0
*
         VBREM=BREM
         BREM=1.D0/DTIME
         TIME=TIME+DTIME         
         DTTRX=DTIME/TRXCEN
         TTRX=TTRX+DTTRX
         END IF
C
        VMHOLE=XMHOLE
C
        IF(IMOD.LT.0)GOTO 95
C*************Here radial grid movement if implemented************
        IF(LS(26))THEN
        DO 96 I=2,IS(12)
 96     R(I)=R(I)+CS(19)*R(I)/R(IS(12))/BREM
        DO 94 I=IS(12)+1,NJ
        DR=(DLOG(R(NJ))-DLOG(R(IS(12))))/DBLE(NJ-IS(12))
 94     R(I)=DEXP(DBLE(I-IS(12))*DR)*R(IS(12))
        END IF
C*****************************************************************
 95     CONTINUE
        IF(NMOD.EQ.0)GOTO 91
C
      IF((IMOD.EQ.0).OR.(IMOD.EQ.MODF))THEN
      WRITE(6,991)(VX(I,1),I=1,NPOS)
      WRITE(6,992)(X(I,1),I=1,NPOS)
      END IF
C       Limit Timestep for Escapers
      IF(LS(11))THEN
*
      DTC = 0.D0
*
      DO 2000 J=1,NCOMP
*
      DO 500 I=1,NJ
*
      VESC = DSQRT(2.D0*(PHTID-PHI(I)-PHIBIN(I)))
      ROUT = RTIDAL
      VESCR = VESC
*       Select Apocentre Criterion
      VESCT = VESC
      IF(LS(29))VESCT = VESC*ROUT/DSQRT(ROUT**2-R(I)**2)
      SIGR = DEXP((X(I20+J,I)-X(I00+J,I))*C12)
      SIGT = DEXP((X(I02+J,I)-X(I00+J,I))*C12)
      TOUT = XALPHA*(ROUT-R(I))/VESC
*       Old expressions
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
*       rho*XFAC is density of escaping stars
         SIG=DSQRT(C13*(SIGR**2+2.D0*SIGT**2))
         TRX(J,J)=SIG**3/(DEXP(X(I00+J,I))*XMIND(J)*XCOUL(J)*CS(60))
*        EFAC = C32*SIG**2 + PHI(I) + PHIBIN(I)
*        TIN = (EFAC - PHTID)/EFAC*TRX(J,J)
         TIN = XBETA*TRX(J,J)
*
*        PF = 1.D0/(1.D0+TOUT/TIN)
         PF = 1.D0
      XCUT = 1.D0/(1.D0+(R(I)/R(NJ-2))**4)
*     XCUT = 1.D0
*
      DRCORR = Y(J,I)*PF*XFAC/TOUT/BREM/CORR*XCUT
      DERCRR = Y(J,I)*PF*XRFAC/TOUT/BREM/CORR*XCUT
      DETCRR = Y(J,I)*PF*XTFAC/TOUT/BREM/CORR*XCUT
*
      IF(DRCORR.GT.DTC)DTC = DRCORR
      IF(DERCRR.GT.DTC)DTC = DERCRR
      IF(DETCRR.GT.DTC)DTC = DETCRR
*
 500  CONTINUE
 2000 CONTINUE
*
      IF(DTC.GT.BREM)THEN
      PRINT*,' Choose DTtid,DT=',1.D0/DTC,1.D0/BREM
      BREM = DTC
      END IF
*
      END IF
*
          RETURN
C***********************************************************************
 81   FORMAT(1X,' mistake IN .CHOOSE : variable ',I3,
     *   ' less than zero at grid no. ',I3)
C**********************************************************************
C     The following is executed only for a beginning run (NMOD=0)
C**********************************************************************
  91  CONTINUE
      DO 93 J=1,NJ
      DO 92 K=1,NG
         VX(K,J)=X(K,J)
  92     VVX(K,J)=VX(K,J)
  93  CONTINUE
C
      TIME=TIME+1.D0/BREM
C
      RETURN
C
 991  FORMAT(' CONTROL CHOOSE VX(I,1)=',1P,12E9.2,/,
     *      ('                     ...',1P,12E9.2,/))
 992  FORMAT(' CONTROL CHOOSE  X(I,1)=',1P,12E9.2,/,
     *      ('                     ...',1P,12E9.2,/))
C
          END
