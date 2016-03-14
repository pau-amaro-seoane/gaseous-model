          SUBROUTINE READPM
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
          CHARACTER*80 ZKOM
         INCLUDE 'compar.f'
         DIMENSION XX(NBINO),XY(NBINO),XZ(NBINO)
*
      READ(15,'(A80)',ERR=36,END=36) ZKOM
*
      WRITE(6,765)ZKOM
 765      FORMAT(1X,79('*'),/,A79,/,79('*'),/,1H1)
*
 36       CONTINUE
************************** OPTIONS READ ***********************
          DO 38 J=1,100
          LS(J)=.FALSE.
 38       CONTINUE
*
          DO 39 J=1,37,4
          READ(15,139)(LS(J+K),K=0,3)
 139      FORMAT(4(12X,L3))
 39       CONTINUE
************************** MODEL SPECIFICATIONS*************************
          READ(15,2200)MODA,MODM,MODF,MODFAC,MODFI,MAXMOD,NPR,NCOMP,NJ,
     *    NRIT,NCOL,NWRITE,ITIME,IUNIT,TSTEP,STEU1,STEU2,STEU3,NRAND
 2200     FORMAT(7X,I4,7X,I5,7X,I5,7X,I5,7X,I5,5X,I8,5X,I3,/,
     *           7X,I4,7X,I5,7X,I5,7X,I5,7X,I5,7X,I6,7X,I1,/,
     *           7X,1P,D9.2,7X,D9.2,7X,D9.2,7X,D9.2,7X,I8)
         WRITE(6,200)MODA,MODM,MODF,MODFAC,MODFI,MAXMOD,NPR,NCOMP,NJ,
     *    NRIT,NCOL,NWRITE,ITIME,IUNIT,TSTEP,STEU1,STEU2,STEU3,NRAND
 200      FORMAT(' MODA=',I4,' MODM=',I5,' MODF=',I5,' MODFAC=',
     *    I5,' MODFI=',I5,' MAX=',I8,' NPR=',I3,
     *    /,' NCOMP=',I4,' NJ=',I5,' NRIT=',I5,' NCOL=',I5,
     *    ' NWRITE=',I5,' ITIME=',I6,' IUNIT=',I1,/,
     *    ' TSTEP=',1P,D9.2,' STEU1=',D9.2,' STEU2=',D9.2,
     *    ' STEU3=',D9.2,' NRAND=',I8)
C
      READ(15,2201)ITMIN,ITMAX,EPS,FHEN,FACT,CORR,FAA,FVEC,
     *  IPR1,IPR2,ICTR
 2201  FORMAT(7X,I3,7X,I3,5X,1P,D9.2,6X,D9.2,6X,D9.2,/,
     *   6X,D9.2,5X,D9.2,7X,D9.2,6X,I3,6X,I3,6X,I3)
      WRITE(6,201)ITMIN,ITMAX,EPS,FHEN,FACT,CORR,FAA,FVEC,
     *  IPR1,IPR2,ICTR
 201  FORMAT(' ITMIN=',I3,' ITMAX=',I3,' EPS=',1P,D9.2,' FHEN=',
     *    D9.2,' FACT=',D9.2,/,' CORR=',D9.2,' FAA=',D9.2,
     *    '  FVEC=',D9.2,' IPR1=',I3,' IPR2=',I3,' ICTR=',I3)
*
*        Read in G, E5POLY, and MTOT in units of calculation
      READ(15,521)CG,CE5,CMT
 521  FORMAT(7X,1P,D12.5,7X,D12.5,7X,D12.5)
      WRITE(6,523)CG,CE5,CMT
 523  FORMAT(' Desired Units of calculation:',/,
     *       ' Grav. Constant =',1P,D12.5,' Energy of n=5 Polytr.=',
     *        D12.5,' Total Mass=',D12.5)
*
      READ(15,505)XLA,XTAN,XBIN,XIMPL
      WRITE(6,5050)XLA,XTAN,XBIN,XIMPL
 5050 FORMAT('  XLAM=',1P,D12.5,'  XTAN=',D12.5,'  XBIN=',D12.5,
     *       ' XIMPL=',D12.5)
 505  FORMAT(7X,1P,D12.5,7X,D12.5,7X,D12.5,7X,D12.5)
      READ(15,506)XMHOLE,ELL,THETA
      WRITE(6,5060)XMHOLE,ELL,THETA
 5060 FORMAT(' MHOLE=',1P,D12.5,' ELLIP=',D12.5,
     *       ' THETA=',D12.5)
 506  FORMAT(7X,1P,D12.5,7X,D12.5,7X,D12.5)
      READ(15,507)XESC,TBINI,AEXP,XLFAC
      WRITE(6,5070)XESC,TBINI,AEXP,XLFAC
 5070 FORMAT('  XESC=',1P,D12.5,' TBINI=',D12.5,'  AEXP=',D12.5,
     *       ' XLFAC=',D12.5)
 507  FORMAT(7X,1P,D12.5,7X,D12.5,7X,D12.5,7X,D12.5)
      READ(15,508)XTEQ,GAMMA,XBIN2B,FTRX
      WRITE(6,5080)XTEQ,GAMMA,XBIN2B,FTRX
 5080 FORMAT('  XTEQ=',1P,D12.5,' GAMMA=',D12.5,'XBIN2B=',D12.5,
     *       '  FTRX=',D12.5)
 508  FORMAT(7X,1P,D12.5,7X,D12.5,7X,D12.5,7X,D12.5)
*
      READ(15,509)ASPI,AHEG,BINFR,EBINI
      WRITE(6,5090)ASPI,AHEG,BINFR,EBINI
 5090 FORMAT('  ASPI=',1P,D12.5,'  AHEG=',D12.5,' BINFR=',D12.5,
     *       ' EBINI=',D12.5)
 509  FORMAT(7X,1P,D12.5,7X,D12.5,7X,D12.5,7X,D12.5)
      READ(15,510)EMIN,EMAX,XALPHA,XBETA
      WRITE(6,5100)EMIN,EMAX,XALPHA,XBETA
 5100 FORMAT('  EMIN=',1P,D12.5,'  EMAX=',D12.5,' ALPHA=',D12.5,
     *       '  BETA=',D12.5)
 510  FORMAT(7X,1P,D12.5,7X,D12.5,7X,D12.5,7X,D12.5)
      READ(15,511)TUPIN,SMBIN,XXESC,WASTE
      WRITE(6,5110)TUPIN,SMBIN,XXESC,WASTE
 5110 FORMAT(' TUPIN=',1P,D12.5,' SMBIN=',D12.5,' XXESC=',D12.5,
     *       '    CS=',D12.5)
 511  FORMAT(7X,1P,D12.5,7X,D12.5,7X,D12.5,7X,D12.5,10X)
      READ(15,512)WASTE,WASTE,WASTE,WASTE
      WRITE(6,5120)WASTE,WASTE,WASTE,WASTE
 5120 FORMAT('    CS=',1P,D12.5,'    CS=',D12.5,'    CS=',D12.5,
     *       '    CS=',D12.5)
 512  FORMAT(7X,1P,D12.5,7X,D12.5,7X,D12.5,7X,D12.5)
*
      DO 3310 I=1,NCOMP
      READ(15,3300)XMTOT(I),XMIND(I),XMREM(I),TAU(I),RTOT(I),
     *    RIND(I),TIND(I),RSC(I),EXP(I),UINIT(I),XBIND(I),NOL(I),
     *       (LEQ(K),K=NEQC*(I-1)+I00+1,NEQC*I+I00)
      WRITE(6,3311)I,XMTOT(I),XMIND(I),XMREM(I),TAU(I),RTOT(I),
     *    RIND(I),TIND(I),RSC(I),EXP(I),UINIT(I),
     *    XBIND(I),NOL(I),(LEQ(K),K=NEQC*(I-1)+I00+1,NEQC*I+I00)
 3310 CONTINUE
 3300 FORMAT(10X,1P,D12.5,6X,D12.5,7X,D12.5,5X,D12.5,/,
     *               10X,D12.5,6X,D12.5,7X,D12.5,5X,D12.5,/,
     *               10X,D12.5,6X,D12.5,7X,D12.5,5X,I4,/,
     *               10X,(6L3))
 3311 FORMAT(' COMPONENT DATA: ',/,I4,' MTOT=',1P,D12.5,' MIND=',D12.5, 
     *      '  MREM=',D12.5,' TAU=',D12.5,/,'     RTOT=',D12.5,
     *      ' RIND=',D12.5,'  TIND=',D12.5,' RSC=',D12.5,/,
     *      '      EXP=',D12.5,' UINI=',D12.5,' XBIND=',D12.5,' NOL=',I4,
     *    /,'      LEQ=',(6L3))
*
      READ(15,3312)XNTOTA,XIMF,LEQ(IMR)
 3312 FORMAT(7X,1P,D12.5,7X,D12.5,5X,L3)
      WRITE(6,3313)XNTOTA,XIMF,LEQ(IMR)
 3313 FORMAT(' XNTOTA=',1P,D12.5,' XIMF=',D12.5,' LEQ=',L3)
*
      I=1
      IF(NCOMP.GT.1)THEN
      I=2
      IF(XIMF.NE.0.D0)THEN
*
      IF(XIMF.LT.0.D0)THEN
*
      PRINT*,' This choice of XIMF=',XIMF,' enables the '
      PRINT*,' construction of a stellar mass spectrum with    '
      PRINT*,' MTOT=',XMTOT(1),' solar masses  and  '
      PRINT*,' stellar mass propto. MIND(I)**',XIMF
      PRINT*,' As mmin, mmax XMIND(1,2) are taken ',XMIND(1),XMIND(2)
      PRINT*,' So, all values of XMTOT(I) are not used - only XNTOTA!'
*
      XMSUM=0.D0
      XMEXP=-DABS(XIMF)
      XMTT=XMTOT(1)
      XMIN = XMIND(1)
      XMAX = XMIND(2)
*
      PRINT*,' mi determined in equal log steps '
      PRINT*,' mmin, mmax used are XMIND(1,2)=',XMIND(1),XMIND(2)
      PRINT*,' So, all other values of XMIND(I),I>2 are not used '
*
      DO 9 I=1,NCOMP
  9   XMIND(I)=XMIN*(XMAX/XMIN)**(DBLE(2*I-1)/DBLE(2*NCOMP))
      DO 11 I=1,NCOMP
      IF(I.EQ.1)THEN
      XLOW=XMIN
      ELSE
      XLOW=XMIN*(XMAX/XMIN)**(DBLE(2*I-2)/DBLE(2*NCOMP))
      END IF
      XUP=XMIN*(XMAX/XMIN)**(DBLE(2*I)/DBLE(2*NCOMP))
 11   XMSUM=XMSUM+(XUP-XLOW)*XMIND(I)**XMEXP*XMIND(I)
*
      DO 12 I=1,NCOMP
      IF(I.EQ.1)THEN
      XLOW=XMIN
      ELSE
      XLOW=XMIN*(XMAX/XMIN)**(DBLE(2*I-2)/DBLE(2*NCOMP))
      END IF
      XUP=XMIN*(XMAX/XMIN)**(DBLE(2*I)/DBLE(2*NCOMP))
 12   XMTOT(I)=(XUP-XLOW)*XMIND(I)**XMEXP*XMIND(I)/XMSUM*XMTT
*
      END IF
*
      IF(XIMF.GT.0.D0)THEN
*
      PRINT*,' mi taken from input data '
      PRINT*,' mmin, mmax used are XMIND(1,NCOMP)=',
     *   XMIND(1),XMIND(NCOMP)
*
      XMSUM=0.D0
      XMEXP=-DABS(XIMF)
*
      DO 13 I=1,NCOMP
      IF(I.EQ.1)THEN
      XLOW=XMIND(1)*DSQRT(XMIND(1)/XMIND(2))
      ELSE
      XLOW=DSQRT(XMIND(I)*XMIND(I-1))
      END IF
      IF(I.EQ.NCOMP)THEN
      XUP=XMIND(NCOMP)*DSQRT(XMIND(NCOMP)/XMIND(NCOMP-1))
      ELSE
      XUP=DSQRT(XMIND(I+1)*XMIND(I))
      END IF
 13   XMSUM=XMSUM+(XUP-XLOW)*XMIND(I)**XMEXP*XMIND(I)
*
      DO 14 I=1,NCOMP
      IF(I.EQ.1)THEN
      XLOW=XMIND(1)*DSQRT(XMIND(1)/XMIND(2))
      ELSE
      XLOW=DSQRT(XMIND(I)*XMIND(I-1))
      END IF
      IF(I.EQ.NCOMP)THEN
      XUP=XMIND(NCOMP)*DSQRT(XMIND(NCOMP)/XMIND(NCOMP-1))
      ELSE
      XUP=DSQRT(XMIND(I+1)*XMIND(I))
      END IF
 14   XMTOT(I)=(XUP-XLOW)*XMIND(I)**XMEXP*XMIND(I)/XMSUM
      END IF
      END IF
      END IF
*        Fix exact value of total particle number wanted 
*        (read into XNTOTA)
      XNTOT=0.D0
      DO 5235 J=1,NCOMP
 5235 XNTOT=XNTOT+XMTOT(J)/XMIND(J)
*
      IF(XNTOTA.GT.0.D0)THEN
      CFAC=XNTOT/XNTOTA
      ELSE
      CFAC=1.D0
      END IF
*
      DO 5236 J=1,NCOMP
 5236 XMTOT(J)=XMTOT(J)/CFAC
*
*        Derive Calcul. Units of Radius, Mass, and Time in CGS Units
*        Determine total mass in input units (solar masses)
*             and  average scale radius in input units (parsecs)
      CMTOT=0.D0
      CMRSC=0.D0
      DO 524 J=1,NCOMP
      CMRSC=CMRSC+RSC(J)
 524  CMTOT=CMTOT+XMTOT(J)
*        Transform to CGS units
      CMTOT=CMTOT*SUNM
      CMRSC=CMRSC/DBLE(NCOMP)*PC
*        Calculational Unit for mass
      CUNM=CMT*CMTOT
*        Calculational Unit for radius 
      CUNR=64.D0/3.D0/PI*CE5/CG/CMT**2*CMRSC
*        Calculational Unit for time
      CUNT=DSQRT(CUNR/GRAV/CUNM*CUNR**2)
*
      WRITE(6,522)CUNR,CUNM,CUNT
 522  FORMAT(' CUN: R=',1P,D12.5,' M=',D12.5,' T=',D12.5)
*
      READ(15,*)(IC(K),K=1,20)
      IF(LS(25))THEN
      WRITE(6,444)(IC(K),K=1,18)
 444  FORMAT(1X,79('*'),/,' THIS IS A DIAGNOSTIC RUN - PRINTING ',
     *   ' LARGE INFO AT:',/,' IMOD,ITER,GRID POINT AS FOLLOWS:',/,
     *   1X,6(3(1X,I3),'**'),/,1X,79('*'))
      IPR1=1
      END IF
*
*       Read tidal cross sections if 2B-binary heating is there
      IF(LS(1))THEN
      DO 600 J=1,NCOMP
      READ(9,*)(CR2B(J,K),ALP2B(J,K),K=1,NCOMP)
 600  CONTINUE
      END IF
*
C-------Construct other units
      CUNRH=CUNM/CUNR**3
      CUNU=CUNR/CUNT
      CUNL=CUNM/CUNR*CUNU**3
*        Transform gravitational constants
      GRAV=GRAV/(CUNR/CUNM/CUNT**2*CUNR**2)
*        Transform velocity of light
      CLICHT=CLICHT/CUNU
*        Factor used for relaxation time
      CS(60)=16.D0*GRAV*GRAV/9.D0*DSQRT(PI)
*        Lambda for heat flux equations
      CS(21)=XLA
*        Factors for escape rate
      CS(22)=0.1096D0/PI4*XESC
      CS(24)=0.01087D0/PI4*XESC
*
*         Factor for gravitational radiation time
      CS(25)=5.D0/256.D0*(CLICHT/GRAV)**3*CLICHT**2
C-------Transform Timestep to appropriate unit
      BREM=1.D0/(TSTEP*YEAR/CUNT)
C
C------Quantities inputted always in SUNM, etc.
C------Transformation to computational units by CUNR, etc.
      XMHOLE=XMHOLE*SUNM/CUNM
C
          IF(LS(19))WRITE(6,204)XMHOLE*CUNM/SUNM,DMTOT*CUNM/SUNM
 204      FORMAT(' START WITH BLACK HOLE MASS=',1P,E9.2,
     *  ' DM=',E9.2,' SOL.MASSES ')

      XNTOT=0.D0
      DO 102 I=1,NCOMP
      PRINT*,' I=',I,' NI=',XMTOT(I)/XMIND(I),' Mi=',XMTOT(I),' mi=',
     *     XMIND(I)
 102  XNTOT=XNTOT+XMTOT(I)/XMIND(I)
*
      PRINT*,' XNTOT=',XNTOT
*
      XCTOT=DLOG(GAMMA*XNTOT)
*
      DO 101 I=1,NCOMP
      TMESC(I)=0.D0
      TEESC(I)=0.D0
      XMTOT(I)=XMTOT(I)*SUNM/CUNM
      XMIND(I)=XMIND(I)*SUNM/CUNM
*
      IF(GAMMA.EQ.0.D0)THEN
      PRINT*,' Warning! Gamma not set in Input Data '
      PRINT*,' Standard Value of 0.11 for 1-component system taken ...'
      GAMMA=0.11D0
      END IF
*
      XCOUL(I)=DLOG(GAMMA*XNTOT)
      RTOT(I)=RTOT(I)*PC/CUNR
      RIND(I)=RIND(I)*SUNR/CUNR
      RSC(I)=RSC(I)*PC/CUNR
      UINIT(I)=UINIT(I)*CLFAC/CUNU
      RHOSTA(I)=XMIND(I)/PI43/RIND(I)**3
C
      IF(LS(5))XMTOT(I) = XMTOT(I)*(1.D0-BINFR)
*
 101  CONTINUE
*
*      Initialize stellar evolution
      IF(LS(21))CALL STEVOL
*
      IF(LS(5).AND..NOT.LS(10))THEN
*      Assume Mtot normalized to one
      XMBINT = BINFR
*      Assume binary mass as double single mass if not given
      IF(SMBIN.EQ.0.D0)SMBIN=2.D0
      IBIN = XNTOT*BINFR/SMBIN
*
      PRINT*,' Starting with ',IBIN,' Binaries and ',
     *    XMTOT(1)/XMIND(1),' Singles '
*
      READ(7,*)NBS
      READ(7,*)IDUM
      READ(7,*)DUMMY
*
      IF(IBIN.NE.NBS)THEN
      PRINT*,' Probably wrong binary input file NBS,IBIN=',NBS,IBIN
      STOP
      END IF
*
      XMBTOT=0.D0
      DO 103 IB=1,IBIN
      READ(7,*)BODY1(IB)
 103  XMBTOT=XMBTOT+SMBIN*BODY1(IB)
*
      DO 104 IB=1,IBIN
      READ(7,*)XX(IB),XY(IB),XZ(IB)
 104  RB(IB)=XX(IB)**2+XY(IB)**2+XZ(IB)**2
*
      DO 105 IB=1,IBIN
      READ(7,*)WX,WY,WZ
      XSIGN=1.D0
*
      VDOTR=WX*XX(IB)+WY*XY(IB)+WZ*XZ(IB)
      IF(VDOTR.LT.0.D0)XSIGN=-1.D0
      VR(IB)=(VDOTR)**2/RB(IB)
      VT(IB)=WX**2+WY**2+WZ**2-VR(IB)
      RB(IB)=DSQRT(RB(IB))
      VR(IB)=XSIGN*DSQRT(VR(IB))
      VT(IB)=DSQRT(VT(IB))
      BODY1(IB)=BODY1(IB)*XMBINT/XMBTOT
*       Planets if binary mass equals individual stellar mass
      IF(SMBIN.EQ.1.D0)THEN
      BODY2(IB)=1.D-10
      ELSE
      BODY2(IB)=SMBIN/2.D0*BODY1(IB)
      END IF
*      test scaling
*     VSCALE=DSQRT(1.D0/XNTOT/(BODY1(IB)+BODY2(IB)))
*     VR(IB)=VSCALE*VR(IB)
*     VT(IB)=VSCALE*VT(IB)
*     EB(IB)=EBINI
*     ECC(IB)=0.D0
*     SEMIA(IB) = C12*GRAV*BODY1(IB)*BODY2(IB)/EB(IB)
      SIZE1(IB)=RIND(1)
      SIZE2(IB)=RIND(1)
      NAMEB(IB)=IB
      ICO(IB)=1
 105  CONTINUE
*
      XMBTOT=0.D0
      DO 107 IB=1,IBIN
 107  XMBTOT=XMBTOT+(BODY1(IB)+BODY2(IB))
*
      END IF
          RETURN
          END
