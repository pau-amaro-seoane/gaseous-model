      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
      CHARACTER*30 TSTR
      DIMENSION RMIN1(100),RMIN2(100),RMIN3(100),RMIN4(100)
      DIMENSION RMAX1(100),RMAX2(100),RMAX3(100),RMAX4(100)
      DIMENSION SMIN1(100),SMIN2(100),SMIN3(100),SMIN4(100)
      DIMENSION SMAX1(100),SMAX2(100),SMAX3(100),SMAX4(100)
C
      DIMENSION TRMIN1(100),TRMIN2(100),TRMIN3(100),TRMIN4(100)
      DIMENSION TRMAX1(100),TRMAX2(100),TRMAX3(100),TRMAX4(100)
      DIMENSION TSMIN1(100),TSMIN2(100),TSMIN3(100),TSMIN4(100)
      DIMENSION TSMAX1(100),TSMAX2(100),TSMAX3(100),TSMAX4(100)
C
      DIMENSION XY(30,2),STAT(3)
      DIMENSION SLRMAX(4,2),SLRMIN(4,2),SLSMAX(4,2),SLSMIN(4,2)
C
      IREC=0
      IRMAX=0
      IRMIN=0
      ISMAX=0
      ISMIN=0
C
      PRINT*,' Attractor Search in Rho, Sigma, Xi Embedding Space ...'
C
      PRINT*,' Enter Plot Titles String:'
      READ(5,*,END=500,ERR=500)TSTR
 500  CONTINUE
      PRINT*,' Plot Titles: ',TSTR
      PRINT*,' Using fort.50 as input data...'
      PRINT*,' Step 1: Detrending Rho and Sigma '
      PRINT*,' Enter initial offset for linear regression:'
      READ*,IOFF
C
      IF(IOFF.GT.1000)THEN
      PRINT*,' Direct Input of Detrending Exponents:'
      PRINT*,' Enter Rho-Value, Sigma-Value:'
      READ*,ARHO,ASIG
      GOTO 5555
      END IF
C
 900  IREC=IREC+1
C
      VTIME=TIME
      VTTRX=TTRX
      IF(RHO.NE.VRHO)THEN
      VRHO=RHO
      V2DRHO=VDRHO
      VDRHO=DRHO
      END IF
      IF(SIG.NE.VSIG)THEN
      VSIG=SIG
      VDSIG=DSIG
      END IF
      VXI=XI
C
      READ(50,*,ERR=800,END=800)TIME,TTRX,RHO,SIG,XI
C
      DRHO=RHO-VRHO
      IF(DRHO.EQ.0.D0)THEN
      PRINT*,' WARNING -- WARNING -- WARNING '
      PRINT*,' At TIME=',TIME,' DRHO=',DRHO,' Extrema check ambiguous!'
      END IF
C
      IF(DRHO.EQ.0.D0)THEN
      DRHO=VDRHO
      VDRHO=V2DRHO
      END IF
      IF(VDRHO.EQ.0.D0)THEN
      VDRHO=V2DRHO
      END IF
C
      DSIG=SIG-VSIG
C        Density Extrema
      IF(DRHO*VDRHO.LT.0.D0)THEN
C        Density Maxima
      IF(DRHO.GT.0.D0)THEN
C      PRINT*,' Rho Maximum found at T=',1.D1**TIME
      IRMAX=IRMAX+1
C
      IF(MOD(IRMAX,4).EQ.1)THEN
      RMAX1(IRMAX/4)=RHO
      TRMAX1(IRMAX/4)=TIME
      END IF
      IF(MOD(IRMAX,4).EQ.2)THEN
      RMAX2(IRMAX/4)=RHO
      TRMAX2(IRMAX/4)=TIME
      END IF
      IF(MOD(IRMAX,4).EQ.3)THEN
      RMAX3(IRMAX/4)=RHO
      TRMAX3(IRMAX/4)=TIME
      END IF
      IF(MOD(IRMAX,4).EQ.0)THEN
      RMAX4(IRMAX/4)=RHO
      TRMAX4(IRMAX/4)=TIME
      END IF
C
      WRITE(11,*)TIME,RHO
      END IF
      IF(DRHO.LE.0.D0)THEN
C      PRINT*,' Rho Minimum found at T=',1.D1**TIME
      IRMIN=IRMIN+1
C
      IF(MOD(IRMIN,4).EQ.1)THEN
      RMIN1(IRMIN/4)=RHO
      TRMIN1(IRMIN/4)=TIME
      END IF
      IF(MOD(IRMIN,4).EQ.2)THEN
      RMIN2(IRMIN/4)=RHO
      TRMIN2(IRMIN/4)=TIME
      END IF
      IF(MOD(IRMIN,4).EQ.3)THEN
      RMIN3(IRMIN/4)=RHO
      TRMIN3(IRMIN/4)=TIME
      END IF
      IF(MOD(IRMIN,4).EQ.0)THEN
      RMIN4(IRMIN/4)=RHO
      TRMIN4(IRMIN/4)=TIME
      END IF
C
      WRITE(12,*)TIME,RHO
      END IF 
      END IF
C
      IF(DSIG*VDSIG.LT.0.D0)THEN
      IF(DSIG.GT.0.D0)THEN
C      PRINT*,' Sigma Maximum found at T=',1.D1**TIME
      ISMAX=ISMAX+1
C
      IF(MOD(ISMAX,4).EQ.1)THEN
      SMAX1(ISMAX/4)=SIG
      TSMAX1(ISMAX/4)=TIME
      END IF
      IF(MOD(ISMAX,4).EQ.2)THEN
      SMAX2(ISMAX/4)=SIG
      TSMAX2(ISMAX/4)=TIME
      END IF
      IF(MOD(ISMAX,4).EQ.3)THEN
      SMAX3(ISMAX/4)=SIG
      TSMAX3(ISMAX/4)=TIME
      END IF
      IF(MOD(ISMAX,4).EQ.0)THEN
      SMAX4(ISMAX/4)=SIG
      TSMAX4(ISMAX/4)=TIME
      END IF
C
      WRITE(13,*)TIME,SIG
      END IF
      IF(DSIG.LT.0.D0)THEN
C      PRINT*,' Sigma Minimum found at T=',1.D1**TIME
      ISMIN=ISMIN+1
C
      IF(MOD(ISMIN,4).EQ.1)THEN
      SMIN1(ISMIN/4)=SIG
      TSMIN1(ISMIN/4)=TIME
      END IF
      IF(MOD(ISMIN,4).EQ.2)THEN
      SMIN2(ISMIN/4)=SIG
      TSMIN2(ISMIN/4)=TIME
      END IF
      IF(MOD(ISMIN,4).EQ.3)THEN
      SMIN3(ISMIN/4)=SIG
      TSMIN3(ISMIN/4)=TIME
      END IF
      IF(MOD(ISMIN,4).EQ.0)THEN
      SMIN4(ISMIN/4)=SIG
      TSMIN4(ISMIN/4)=TIME
      END IF
C
      WRITE(14,*)TIME,SIG
      END IF
      END IF

C
      IF(MOD(IREC,500).EQ.0)PRINT*,IREC,' Records read...'
      GOTO 900
 800  CONTINUE
      PRINT*,IREC-1,' Records read - STOP '
      PRINT*,' Found were: '
      PRINT*,' =========== '
      PRINT*,IRMAX,' Maxima and',IRMIN,' Minima of rho '
      PRINT*,ISMAX,' Maxima and',ISMIN,' Minima of sigma '
C
      PRINT*,' Start Calculation of Slopes '
C
      NFIT=IRMAX/4-IOFF
      II=IRMAX/4
      PROD=TRMAX1(II)*TRMAX2(II)*TRMAX3(II)*TRMAX4(II)
      IF(PROD.EQ.0.D0)NFIT=NFIT-1
      PRINT*,' Rho max NFIT=',NFIT
      IF(NFIT.GT.30)PRINT*,' Warning DIM(XY) too small!!'
      DO 501 I=1,NFIT
      XY(I,1)=TRMAX1(I+IOFF)
 501  XY(I,2)=RMAX1(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLRMAX(1,1)=STAT(1)
      SLRMAX(1,2)=STAT(2)
      DO 502 I=1,NFIT
      XY(I,1)=TRMAX2(I+IOFF)
 502  XY(I,2)=RMAX2(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLRMAX(2,1)=STAT(1)
      SLRMAX(2,2)=STAT(2)
      DO 503 I=1,NFIT
      XY(I,1)=TRMAX3(I+IOFF)
 503  XY(I,2)=RMAX3(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLRMAX(3,1)=STAT(1)
      SLRMAX(3,2)=STAT(2)
      DO 504 I=1,NFIT
      XY(I,1)=TRMAX4(I+IOFF)
 504  XY(I,2)=RMAX4(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLRMAX(4,1)=STAT(1)
      SLRMAX(4,2)=STAT(2)
      PRINT*,' Density Maxima Slopes/S.Devs:'
      PRINT*,(SLRMAX(I,1),I=1,4)
      PRINT*,(SLRMAX(I,2),I=1,4)
      ARMAX=0.D0
      DO 505 I=1,4
 505  ARMAX=ARMAX+SLRMAX(I,1)
      ARMAX=ARMAX/4.D0
      PRINT*,' Average Density Maxima Slope=',ARMAX
C
      NFIT=IRMIN/4-IOFF
      II=IRMIN/4
      PROD=TRMIN1(II)*TRMIN2(II)*TRMIN3(II)*TRMIN4(II)
      IF(PROD.EQ.0.D0)NFIT=NFIT-1
      PRINT*,' Rhomin NFIT=',NFIT
      IF(NFIT.GT.30)PRINT*,' Warning DIM(XY) too small!!'
      DO 551 I=1,NFIT
      XY(I,1)=TRMIN1(I+IOFF)
 551  XY(I,2)=RMIN1(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLRMIN(1,1)=STAT(1)
      SLRMIN(1,2)=STAT(2)
      DO 552 I=1,NFIT
      XY(I,1)=TRMIN2(I+IOFF)
 552  XY(I,2)=RMIN2(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLRMIN(2,1)=STAT(1)
      SLRMIN(2,2)=STAT(2)
      DO 553 I=1,NFIT
      XY(I,1)=TRMIN3(I+IOFF)
 553  XY(I,2)=RMIN3(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLRMIN(3,1)=STAT(1)
      SLRMIN(3,2)=STAT(2)
      DO 554 I=1,NFIT
      XY(I,1)=TRMIN4(I+IOFF)
 554  XY(I,2)=RMIN4(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLRMIN(4,1)=STAT(1)
      SLRMIN(4,2)=STAT(2)
      PRINT*,' Density Minima Slopes/S.Devs:'
      PRINT*,(SLRMIN(I,1),I=1,4)
      PRINT*,(SLRMIN(I,2),I=1,4)
      ARMIN=0.D0
      DO 555 I=1,4
 555  ARMIN=ARMIN+SLRMIN(I,1)
      ARMIN=ARMIN/4.D0
      PRINT*,' Average Density Minima Slope=',ARMIN
C
      NFIT=ISMAX/4-IOFF
      II=ISMAX/4
      PROD=TSMAX1(II)*TSMAX2(II)*TSMAX3(II)*TSMAX4(II)
      IF(PROD.EQ.0.D0)NFIT=NFIT-1
      PRINT*,' Sigmax NFIT=',NFIT
      DO 601 I=1,NFIT
      XY(I,1)=TSMAX1(I+IOFF)
 601  XY(I,2)=SMAX1(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLSMAX(1,1)=STAT(1)
      SLSMAX(1,2)=STAT(2)
      DO 602 I=1,NFIT
      XY(I,1)=TSMAX2(I+IOFF)
 602  XY(I,2)=SMAX2(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLSMAX(2,1)=STAT(1)
      SLSMAX(2,2)=STAT(2)
      DO 603 I=1,NFIT
      XY(I,1)=TSMAX3(I+IOFF)
 603  XY(I,2)=SMAX3(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLSMAX(3,1)=STAT(1)
      SLSMAX(3,2)=STAT(2)
      DO 604 I=1,NFIT
      XY(I,1)=TSMAX4(I+IOFF)
 604  XY(I,2)=SMAX4(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLSMAX(4,1)=STAT(1)
      SLSMAX(4,2)=STAT(2)
      PRINT*,' Sigma Maxima Slopes/S.Devs:'
      PRINT*,(SLSMAX(I,1),I=1,4)
      PRINT*,(SLSMAX(I,2),I=1,4)
      ASMAX=0.D0
      DO 605 I=1,4
 605  ASMAX=ASMAX+SLSMAX(I,1)
      ASMAX=ASMAX/4.D0
      PRINT*,' Average Sigma Maxima Slope=',ASMAX
C
      NFIT=ISMIN/4-IOFF
      II=ISMIN/4
      PROD=TSMIN1(II)*TSMIN2(II)*TSMIN3(II)*TSMIN4(II)
      IF(PROD.EQ.0.D0)NFIT=NFIT-1
      PRINT*,' Sigmin NFIT=',NFIT
      DO 651 I=1,NFIT
      XY(I,1)=TSMIN1(I+IOFF)
 651  XY(I,2)=SMIN1(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLSMIN(1,1)=STAT(1)
      SLSMIN(1,2)=STAT(2)
      DO 652 I=1,NFIT
      XY(I,1)=TSMIN2(I+IOFF)
 652  XY(I,2)=SMIN2(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLSMIN(2,1)=STAT(1)
      SLSMIN(2,2)=STAT(2)
      DO 653 I=1,NFIT
      XY(I,1)=TSMIN3(I+IOFF)
 653  XY(I,2)=SMIN3(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLSMIN(3,1)=STAT(1)
      SLSMIN(3,2)=STAT(2)
      DO 654 I=1,NFIT
      XY(I,1)=TSMIN4(I+IOFF)
 654  XY(I,2)=SMIN4(I+IOFF)
      CALL SLOPE(XY,NFIT,STAT)
      SLSMIN(4,1)=STAT(1)
      SLSMIN(4,2)=STAT(2)
      PRINT*,' Sigma Minima Slopes/S.Devs:'
      PRINT*,(SLSMIN(I,1),I=1,4)
      PRINT*,(SLSMIN(I,2),I=1,4)
      ASMIN=0.D0
      DO 655 I=1,4
 655  ASMIN=ASMIN+SLSMIN(I,1)
      ASMIN=ASMIN/4.D0
      PRINT*,' Average Sigma Minima Slope=',ASMIN
C
      ARHO=(ARMAX+ARMIN)/2.D0
      ASIG=(ASMAX+ASMIN)/2.D0
C
      WRITE(21,*)'ERASE'
      WRITE(21,*)'TITLE Density Maxima ',TSTR
      WRITE(21,*)'DATA fort.11'
      WRITE(21,*)'XCOLUMN 1'
      WRITE(21,*)'YCOLUMN 2'
      WRITE(21,*)'LIMITS'
      WRITE(21,*)'BOX'
      WRITE(21,*)'XLABEL Log Nbody Time'
      WRITE(21,*)'YLABEL Max. \\gr\\dc '
      WRITE(21,*)'PTYPE 5 3'
      WRITE(21,*)'POINTS'
C
      WRITE(23,*)'ERASE'
      WRITE(23,*)'TITLE Sigma Maxima ',TSTR
      WRITE(23,*)'DATA fort.13'
      WRITE(23,*)'XCOLUMN 1'
      WRITE(23,*)'YCOLUMN 2'
      WRITE(23,*)'LIMITS'
      WRITE(23,*)'BOX'
      WRITE(23,*)'XLABEL Log Nbody Time'
      WRITE(23,*)'YLABEL Max. \\gs\\dc '
      WRITE(23,*)'PTYPE 5 3'
      WRITE(23,*)'POINTS'
C
      WRITE(22,*)'ERASE'
      WRITE(22,*)'TITLE Density Minima ',TSTR
      WRITE(22,*)'DATA fort.12'
      WRITE(22,*)'XCOLUMN 1'
      WRITE(22,*)'YCOLUMN 2'
      WRITE(22,*)'LIMITS'
      WRITE(22,*)'BOX'
      WRITE(22,*)'XLABEL Log Nbody Time'
      WRITE(22,*)'YLABEL Min. \\gr\\dc '
      WRITE(22,*)'PTYPE 5 3'
      WRITE(22,*)'POINTS'
C
      WRITE(24,*)'ERASE'
      WRITE(24,*)'TITLE Sigma Minima ',TSTR
      WRITE(24,*)'DATA fort.12'
      WRITE(24,*)'XCOLUMN 1'
      WRITE(24,*)'YCOLUMN 2'
      WRITE(24,*)'LIMITS'
      WRITE(24,*)'BOX'
      WRITE(24,*)'XLABEL Log Nbody Time'
      WRITE(24,*)'YLABEL Min. \\gr\\dc '
      WRITE(24,*)'PTYPE 5 3'
      WRITE(24,*)'POINTS'
C
      PRINT*,' Use input fort.21 in MONGO for plot of rho maxima '  
      PRINT*,' Use input fort.22 in MONGO for plot of rho minima '
      PRINT*,' Use input fort.23 in MONGO for plot of sigma maxima '
      PRINT*,' Use input fort.24 in MONGO for plot of sigma minima '
C
 5555 CONTINUE
C
      IREC=0
      REWIND 50
      XSIG=1.D0
 2222 IREC=IREC+1
C
      IF(RHO.NE.VRHO)THEN
      VRHO=RHO
      V2DRHO=VDRHO
      VDRHO=DRHO
      END IF
C
      READ(50,*,ERR=700,END=700)TIME,TTRX,RHO,SIG,XI
C
      DRHO=RHO-VRHO
C
      IF(IREC.GT.5.AND.DRHO*VDRHO.LT.0.D0)XSIG=-XSIG
C
      XI=XSIG*XI
C     
      WRITE(51,5005)TIME,TTRX,RHO-ARHO*TIME,SIG-ASIG*TIME,XI
 5005 FORMAT(1P,5(D14.7,1X))
      IF(MOD(IREC,500).EQ.0)PRINT*,IREC,' Records detrended written...'
      GOTO 2222
C
 700  PRINT*,IREC-1,' Records detrended written on fort.51 '
C
      WRITE(25,*)'ERASE'
      WRITE(25,*)'TITLE Detrended Density ',TSTR
      WRITE(25,*)'DATA fort.51'
      WRITE(25,*)'XCOLUMN 1'
      WRITE(25,*)'YCOLUMN 3'
      WRITE(25,*)'LIMITS'
      WRITE(25,*)'BOX'
      WRITE(25,*)'XLABEL Log Nbody Time'
      WRITE(25,*)'YLABEL \\gr\\dc '
      WRITE(25,*)'LTYPE 1'
      WRITE(25,*)'CONNECT'
C
      PRINT*,' Detrended Rho-Plot for MONGO: input fort.25 '
      STOP
      END
