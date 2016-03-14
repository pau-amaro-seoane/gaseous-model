      SUBROUTINE BINSTO
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
         INCLUDE 'compar.f'
      COMMON/EBTEST/EBALA,EBEQ,EBEQ2,DEMASS,DEENER
      DIMENSION EBCHCK(NCOMPO),XMCHCK(NCOMPO)
      DIMENSION EFCHCK(NCOMPO),EFREM(NCOMPO)
      DIMENSION EPOT(NCOMPO),ETERM(NCOMPO)
      DIMENSION IESC(NBINO),IESC2(NBINO)
      DIMENSION XKMASS(NCOMPO),XKMASSY(NCOMPO)
      DIMENSION XQMASS(NCOMPO),XQMASSY(NCOMPO)
      DIMENSION XFMASS(NCOMPO),XFMASSY(NCOMPO)
      DIMENSION TRXLOC(NJO),DTBACC(NBINO),TDYLOC(NJO)
      PARAMETER(NTAB=32)
      INTEGER IDUM2,K,IV,IY
      LOGICAL LPLAN
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
      COMMON/TIMER/TCOMP0
      REAL*4  TARRAY(2),ETIME
      DATA tttrip,ttquad /0.0d0,0.0d0/
*
      CHARACTER*9 ZS
      COMMON/INFOHELP/ZS(IVDIM)
*
      REAL RANF
      EXTERNAL RANF
*
*      Stochastic binary energy generation
*
      IEVENT = 0
*
      ISESC = 0
      IBESC = 0
      IQUAD = 0
      IBKICK = 0
      C34=0.75D0
      C43=4.D0/3.D0
      C83=8.D0/3.D0
      C103=10.D0/3.D0
      LPR=.FALSE.
      XMMIN = 1.D-13
*      Provide TOLD in case the time step has to be changed here.
*      Note that CHOOSE has been called before already.
      TOLD = TIME - 1.D0/BREM
*      Get all output here related to previous time (remember binsto
*      called after choose) - see update of TIME before RETURN
      TIME = TOLD
*
      ICORE = 0
      I2C = 0
      I3C = 0
      I4C = 0
      I5C = 0
      I10C = 0
      I20C = 0
      IHALF = 0
      IH100 = 0
      IH10 = 0
      ITRANS = 0
*
*       Initialize local relaxation time
      DO 10 I=2,NJ
      IF(X(IMR,I)-X(IMR,NJ).LT.-1.609D0)IR20=I
      IF(X(IMR,I)-X(IMR,NJ).LT.-0.693D0)IHALF=I
 10   CONTINUE
      IR20 = IR20 + 1
      IHALF = IHALF + 1
      RIR20 = R(IR20)
      RHALF = R(IHALF)
*       For many binaries take their density into account.
      RHOBAV = 0.D0
*
*       Initialize average stellar mass and central sigma for kT
      AVMASS = 0.D0
      SIGAV2 = 0.D0
      DO 20 K=1,NCOMP
      SIGAV2 = SIGAV2 + DEXP(X(I20+K,1)-X(I00+K,1))
 20   AVMASS = AVMASS + XMIND(K)
      SIGAV2 = SIGAV2/DBLE(NCOMP)
      AVMASS = AVMASS/DBLE(NCOMP)
*
      J = 1
      EKT=AVMASS*SIGAV2
      EKTPL = SIGAV2
*
      IF(LS(5))THEN
      ICNT = 0
*
      DO 13 I=2,NJ
*
      IF(RHOBIN(I).GT.0.D0)ICNT = ICNT + 1
*
      IF(ICNT.EQ.5)THEN
      VOL = PI43*R(I)**3
      RHOBAV = XMRBIN(I)/VOL
      GOTO 113
      END IF
*
 13   CONTINUE
 113  CONTINUE
      END IF
*       Use RHOBAV as approximate central binary density.
      DO 12 I=2,NJ
         IF(I.EQ.2)THEN
         RHOTOT = DEXP(X(I00+J,I)) + RHOBAV
         ELSE
         RHOTOT = DEXP(X(I00+J,I))
         END IF
*
         SIG2=C13*DEXP(X(I20+J,I)-X(I00+J,I))*
     *             (1.D0+2.D0*DEXP(X(I02+J,I)-X(I20+J,I)))
*                                                                             
      IF(I.EQ.2)RCORE=
     *   DSQRT(3.D0*SIG2/(4.D0*PI*GRAV*RHOTOT))
         TRXLOC(I)=SIG2**C32/
     *    CS(60)/XMIND(J)/XCOUL(J)/RHOTOT
       TDYLOC(I)=DSQRT(3.D0*PI/32.D0/GRAV/RHOTOT)
      IF(TRXLOC(I).LT.TDYLOC(I))ITRANS = I
      IF(R(I).LT.RCORE)ICORE = I
      IF(R(I).LT.2.D0*RCORE)I2C = I
      IF(R(I).LT.3.D0*RCORE)I3C = I
      IF(R(I).LT.4.D0*RCORE)I4C = I
      IF(R(I).LT.5.D0*RCORE)I5C = I
      IF(R(I).LT.10.D0*RCORE)I10C = I
      IF(R(I).LT.20.D0*RCORE)I20C = I
      IF(R(I).LT.RHALF/1.D2)IH100 = I
      IF(R(I).LT.RHALF/1.D1)IH10 = I
 12   CONTINUE
      TRXLOC(1) = TRXLOC(2)
      TDYLOC(1) = TDYLOC(2)
      IF(MOD(IMOD,NPR).EQ.0)
     * WRITE(6,666)TIME,DEXP(X(I00+J,1)),RHOBAV,
     * DEXP(X(I00+J,1)) + RHOBAV
 666  FORMAT(1X,' T=',1P,D13.5,' RHOSBT=',3D15.5)
*
      ITRANS = ITRANS + 1
      ICORE = ICORE + 1
      ITRANS = MAX(ITRANS,ICORE)
      I2C = I2C + 1
      I3C = I3C + 1
      I4C = I4C + 1
      I5C = I5C + 1
      I10C = I10C + 1
      I20C = I20C + 1
      IH10 = IH10 + 1
      IH100 = IH100 + 1
*
*        Start initializing or regaining the binary and events balance data
*
      IF(IMOD.EQ.0)THEN
*
      ITOT3 = 0
      ITOT4 = 0
      IDISS3 = 0
      IDISS4 = 0
*        Store binary balance data from binsto.f for file write
          IFR = IFRFIN
          ISTA = 100+10*NCOMP
          PRINT*,' Binsto: ISTA,IFRFIN=',ISTA,IFRFIN
*
*       keep track of numbers, masses, energies of
*       binary escapers due to relaxation
      XNRBES = VIT(IFR+1)
*       binary escapers due to 3b encounters
      XN3BES = VIT(IFR+2)
*       binary escapers due to 4b encounters
      XN4BES = VIT(IFR+3)
*       single escapers due to relaxation
      XNRSES = VIT(IFR+4)
*       single escapers due to 3b encounters
      XN3SES = VIT(IFR+5)
*       single escapers due to 4b encounters
      XN4SES = VIT(IFR+6)
*       binaries disrupted due to quad
      XN4BIS = VIT(IFR+7)
*       single stars created from 4b binaries
      XN4SBI = VIT(IFR+8)
*       single stars forming binaries
      XN3SBI = VIT(IFR+9)
*       binaries disrupted due to triple
      XN3BIS = VIT(IFR+10)
*
      IFR = IFR + 10
*       now masses ...
      XMRBES = VIT(IFR+1)
      XM3BES = VIT(IFR+2)
      XM4BES = VIT(IFR+3)
      XMRSES = VIT(IFR+4)
      XM3SES = VIT(IFR+5)
      XM4SES = VIT(IFR+6)
      XM4BIS = VIT(IFR+7)
      XM4SBI = VIT(IFR+8)
      XM3SBI = VIT(IFR+9)
      XM3BIS = VIT(IFR+10)
*
      IFR = IFR + 10
*       now (c.m., external) energies...
      XERBES = VIT(IFR+1)
      XE3BES = VIT(IFR+2)
      XE4BES = VIT(IFR+3)
      XERSES = VIT(IFR+4)
      XE3SES = VIT(IFR+5)
      XE4SES = VIT(IFR+6)
      XE4BIS = VIT(IFR+7)
      XE4SBI = VIT(IFR+8)
      XE3SBI = VIT(IFR+9)
      XE3BIS = VIT(IFR+10)
*
      IFR = IFR + 10
*       now binding (internal) energies
      XBRBES = VIT(IFR+1)
      XB3BES = VIT(IFR+2)
      XB4BES = VIT(IFR+3)
      XB3SBI = VIT(IFR+4)
      XB3BIS = VIT(IFR+5)
*
      END IF
*
*       End regaining the binary and events balance data
*
      IF(NMOD.EQ.0)THEN
*       Initialize the portable random number generator (range: 0 to 1).
      IDUM = -1
      IDUM2 = 123456789
      DO 11 KKK=1,NTAB
 11   IV(KKK) = 0
      IY = 0
*
      XRAN = RANF()
*       Skip the first random numbers (NRAND specified at input).
      DO 1 K = 1,NRAND
         XRAN = RANF()
    1 CONTINUE
*
      DO 850 I=1,NJ
 850  VMRBIN(I) = 0.D0
*
      IMCRIT = 0
      XMERR = 0.D0
      XEERR = 0.D0
      XMERR2 = 0.D0
      XMERR3 = 0.D0
      XMERR4 = 0.D0
      XEERR2 = 0.D0
      XEERR3 = 0.D0
      XEERR4 = 0.D0
*
      DO 2 K = 1,NCOMP
      XINDHS(K) = 0.D0
      XINDHB(K) = 0.D0
    2 EREM(K) = 0.D0
*       Reset binary number if there are no primordial binaries
      IF(.NOT.LS(5))THEN
      IBIN=0
      IBINT=0
      ELSE
      IBINT=IBIN
*       Save total initial binary mass
      XMBIN0=XMRBIN(NJ)
*
* Take care if planets should be initialized
      IF(DABS(XMBTOT-DBLE(IBIN)*AVMASS).LT.1.D-3)THEN
      PRINT*,' Binaries initialized as planetary systems!'
      LPLAN = .TRUE.
* Determine average planet mass
      AVPLMA = 0.D0
      DO 5 IB=1,IBIN
   5  AVPLMA = AVPLMA + MIN(BODY1(IB),BODY2(IB))
      AVPLMA = AVPLMA/DBLE(IBIN)
*
      PRINT*,' Binaries initialized as planetary systems!',
     & ' Average planet mass=',AVPLMA
*
      ELSE
      LPLAN = .FALSE.
      END IF
*
*       Initialize Binary Energy Distribution
      IF(LPLAN)THEN
      EBINI = EBINI*EKTPL
      EMIN = EMIN*EKTPL
      EMAX = EMAX*EKTPL
      ELSE
      EBINI = EBINI*EKT
      EMIN = EMIN*EKT
      EMAX = EMAX*EKT
      END IF
*
*       Input Parameters EMIN, EMAX
      IF(EBINI.GT.0.D0)THEN
*
      DO 112 IB=1,IBIN
      EB(IB) = EBINI
      ECC(IB)=0.D0
 112  CONTINUE
*
      PRINT*,' ************************************************ '
      IF(LPLAN)THEN
      PRINT*,' Initial Planets with constant v^2 of ',
     *  EBINI/EKTPL,' quasi-kt '
      ELSE
      PRINT*,' Initial Binaries with constant binding energy of ',
     *  EBINI/EKT,' kt '
      END IF
*
      ELSE
*
      PRINT*,' ************************************************ '
      IF(IBIN.GT.200)THEN
*
      PRINT*,' Random Binary Binding Energy Distribution'
      DLGE=DLOG(EMAX/EMIN)
      EEMIN = EMAX
      EEMAX = EMIN
*
      END IF
*
      IF(IBIN.LE.200)THEN
*
      PRINT*,' Equidistant Binary Binding Energy Distribution'
      DLGE=DLOG(EMAX/EMIN)/DBLE(IBIN-1)
      EEMIN = EMIN
      EEMAX = EMAX
*
      END IF
      PRINT*,' ************************************************ '
*
      DO 111 IB=1,IBIN
*       For prim. bin. number > 200 do stochastic distribution
*       Initialize relaxation intervals randomly
      XRAN = RANF()
      I = ISH(IB)
*
      DTBACC(IB) = XRAN*TRXLOC(I2C)/3.D0
*
      IF(IBIN.GT.200)THEN
      XRAN=RANF()
      EB(IB)=EMIN*DEXP(XRAN*DLGE)
      IF(EB(IB).GT.EEMAX)EEMAX=EB(IB)
      IF(EB(IB).LT.EEMIN)EEMIN=EB(IB)
      ELSE
*       For prim. bin. number <=200 use equidistant distribution
*       like in Heggie and Aarseth 92
      EB(IB) = EMIN*DEXP(DLGE*DBLE(IB-1))
      END IF
*
      ECC(IB)=0.D0
*
 111  CONTINUE
*
      IF(LS(3).AND.LS(5))THEN
*        For primordial binaries with 3-b integration thermal eccentricity distribution
*        (including as well 4-b integration)
*
      ECCAV = 0.D0
      ECCMIN = 1.D0
      ECCMAX = 0.D0
*
      DO 330 IB = 1,IBIN
*
      ICOUNT = 0
 333  CONTINUE
      ICOUNT = ICOUNT + 1
*
      QE = RANF()
*
      QV = 2.D0*RANF()
*
      IF(QV.GT.2.D0*QE) GOTO 333
*
      ECC(IB) = QE
      IF(QE.LT.ECCMIN)ECCMIN=QE
      IF(QE.GT.ECCMAX)ECCMAX=QE
      ECCAV = ECCAV + QE
*
 330  CONTINUE
*
      ECCAV = ECCAV/IBIN
*
      END IF
*        End thermal eccentricity distribution
*
      END IF
*
      IF(LPLAN)THEN
      DO 114 IB=1,IBIN
 114  EB(IB)=EB(IB)*AVPLMA
      EEMAX=EEMAX*AVPLMA
      EEMIN=EEMIN*AVPLMA
      END IF
*
      SMMAX = 0.D0
      SMMIN = 1.D30
      SMMAV = 0.D0
      DO 115 IB=1,IBIN
      SEMIA(IB) = C12*GRAV*BODY1(IB)*BODY2(IB)/EB(IB)
      XIN(IB) = 0.D0
      OMEG(IB) = 0.D0
      IF(SEMIA(IB).GT.SMMAX)SMMAX=SEMIA(IB)
      IF(SEMIA(IB).LT.SMMIN)SMMIN=SEMIA(IB)
      SMMAV = SMMAV + SEMIA(IB)
 115  CONTINUE
      SMMAV = SMMAV/DBLE(IBIN)
*
      IF(LPLAN)THEN
      PRINT*,' Initial Planets with distribution in log E'
      PRINT*,' Minimum, Maximum E=',EEMIN,EEMAX,' in quasi-kt=',
     *   EEMIN/EKTPL/AVPLMA,EEMAX/EKTPL/AVPLMA
      ELSE
      PRINT*,' Initial Binaries with distribution in log E'
      PRINT*,' Minimum, Maximum E=',EEMIN,EEMAX,' in kt=',
     *   EEMIN/EKT,EEMAX/EKT
      END IF
      RAU = CUNR/1.5D13
      PRINT*,' Minimum, Maximum, average ecc.=',ECCMIN,ECCMAX,ECCAV
      PRINT*,' Minimum, Maximum, average semi=',SMMIN,SMMAX,SMMAV
      PRINT*,' Minimum, Maximum, average semi/au=',SMMIN*RAU,
     *     SMMAX*RAU,SMMAV*RAU
*
      PRINT*,' kT=',EKT
*
      END IF
*        End binary initialization
*
      DO 2000 IB=1,IBIN
*      Write Binary data on separate file
      RMINB=R(IBSTA(IB)+1)
      RMAXB=R(IBSTA(IB)+ILEN(IB))
      K=ISH(IB)
      K1 = K - 1
      DPHDR=(PHI(K)+PHIBIN(K)-PHI(K1)-PHIBIN(K1))/(R(K)-R(K1))
      U1 = PHI(K1) + PHIBIN(K1) + DPHDR*(RB(IB)-R(K1))
      E1 = U1 + 0.5D0*(VR(IB)**2 + VT(IB)**2)
      A1 = RB(IB)*VT(IB)
      DTIME = 1.D0/BREM
*      Write binary data for initial time zero
      WRITE(61,2001)IB,NAMEB(IB),0.D0,RMINB,RMAXB,RB(IB),E1,A1,U1,RCORE,
     * EB(IB),EKT,ECC(IB),SEMIA(IB),BODY1(IB),BODY2(IB),SIZE1(IB),
     * SIZE2(IB),TORB(IB)*BREM
 2001 FORMAT(2I8,1P,18D12.4)
      CALL FLUSH(61)
      CALL FLUSH(6)
 2000 CONTINUE
*
      END IF
*
      DO 3 J=1,NCOMP
      XHEAT(J)=EREM(J)
      XFEED(J)=EFREM(J)
      XMASS(J)=0.D0
      XKMASS(J)=0.D0
      XKICK(J)=0.D0
      XMASSY(J)=0.D0
      XQMASS(J)=0.D0
      XQMASSY(J)=0.D0
      XFMASS(J)=0.D0
      XFMASSY(J)=0.D0
      XFBINY(J)=0.D0
    3 XSBINY(J)=0.D0
*
*       Check whether there are any binaries
       IF(IBIN.GT.0)THEN
*
*       Initialize feedback energy vector
*     DO 295 J=1,NCOMP
*295  XFEED(J) = 0.D0
*        First check hardening of binaries if present
*
       NESC = 0
       NESC2 = 0
*
       DO 300 IB=1,IBIN
*        If there was binary formation at last timestep reduce the flag
       IF(IEV(IB).EQ.2)IEV(IB)=1
*
       I=ISH(IB)
*        In case of multi-mass system pick component for relax
       J=ICO(IB)
*
*        Call Monte-Carlo Relaxation of Binary
*        Pick velocity of single star for encounter
      SIGR2 = DEXP(X(I20+J,I)-X(I00+J,I))
      SIGT2 = DEXP(X(I02+J,I)-X(I00+J,I))
      VESC = DSQRT(-2.D0*(PHI(I)-PHTID+PHIBIN(I)))
*
      ICOUNT = 0
*
 305  CONTINUE
*
      VRTEST = VESC*RANF()
      VTTEST = VESC*RANF()
      VPTEST = VESC*RANF()
*
      VR2 = VRTEST**2
      VT2 = VTTEST**2+VPTEST**2
      FV = DEXP(-C12*VR2/SIGR2 - C12*VT2/SIGT2)
*
      QV = RANF()
*
      IF(QV.GT.FV)GOTO 305
*
      VTOT=DSQRT(VR2+VT2)
      IF(VTOT.GT.VESC)GOTO 305
*
      VR2=DSQRT(VR2)
      VT2=DSQRT(VT2)
      IF(RANF().LT.0.5)VR2=-VR2
*
      XNPART=DEXP(X(I00+J,I))/XMIND(J)
      SM2=XMIND(J)
*
      VROLD=DABS(VR(IB))
      VTOLD=DABS(VT(IB))
      ROLD=RB(IB)
      I1=I-1
      U1=PHI(I1)-PHTID+PHIBIN(I1)+
     * (PHI(I)+PHIBIN(I)-PHI(I1)-PHIBIN(I1))/(R(I)-R(I1))*(RB(IB)-R(I1))
      
      ESOLD=XMIND(J)*(C12*(VR2*VR2 + VT2*VT2)+U1)
      RBOLD=RB(IB)
      IOLD = I
*
      IBD = 0
*        Update of TBIN until next relaxation
      IF(RB(IB).LT.R(ITRANS))THEN
      TBIN(IB) = TRXLOC(I2C)/1.D1
      ELSE IF(RB(IB).LT.R(IHALF))THEN
      TBIN(IB) = TRXLOC(I2C)/1.D1
      ELSE
      TBIN(IB) = TRXLOC(I2C)/1.D1
      END IF
*
      INAMEX(IB) = 0
      DTBACC(IB) = DTBACC(IB) + 1.D0/BREM
*        Apply relaxation in intervals of trx at rmin
      IF(DTBACC(IB).GT.TBIN(IB))THEN
         IORB = IORB + 1
         ISTA=IBSTA(IB)
         J=ICO(IB)
         IF(J.NE.1)THEN
         PRINT*,' J is in error IB=',IB,'J=',J
         STOP 'J is in error'
         END IF
         IF(ISTA.LT.1.OR.ISTA.GT.NJ)THEN
         PRINT*,' ISTA is in error IB=',IB,'ISTA=',ISTA
         STOP 'ISTA is in error'
         END IF
*       Force relaxation more often for eccentric binaries
      QQ = RBMAX(IB)/RBMIN(IB)
*
*     TBIN(IB) = TORB(IB)/QQ
*     TBIN(IB) = TORB(IB)/3.D0
*     IF(RB(IB).LT.15.D0*RCORE)THEN
*       TBIN(IB) = 10.D0*RB(IB)/SQRT(VR(IB)**2 + VT(IB)**2)
*     ELSE
*       TBIN(IB) = 20.D0*RCORE
*     END IF
*
*     if(rb(ib).lt.15.d0*rcore) then
*       tbin(ib) = trxloc(i)/1.d1
*     else
*       tbin(ib) = trxloc(i15c)/1.d1
*     endif
      DTBACC(IB) = 0.D0
*       Save times for printout
       TRXX = TRXLOC(I)
       TRXX20 = TRXLOC(I20C)
       TX2 = RB(IB)/RBMAX(IB)*TORB(IB)/3.D0
      IF(RB(IB).LT.15.D0*RCORE)THEN
       TX3 = 10.D0*RB(IB)/SQRT(VR(IB)**2 + VT(IB)**2)
      ELSE
       TX3 = 20.D0*RCORE
      END IF
      RBSMIN = RBMIN(IB)
      RBSMAX = RBMAX(IB)
      RBS = RB(IB)
*
*       Relaxation interval should not be smaller than dt
      IF(TBIN(IB).LT.1.D0/BREM)THEN
      TBIN(IB)=1.D0/BREM
      IBR = IBR + 1
      END IF
*
      LRELAX = .TRUE.
      INAMEX(IB) = 1
      IRX = IRX + 1
*
      ELSE
      LRELAX = .FALSE.
      END IF
*
      IB1 = IB
*
      IF(LRELAX)THEN
*
*       TEST PRINT
*
*      print*,'VRold,VTold,U1,ESold = ',Vr2,Vt2,u1,esold
*     vr2aa = vr2
*     vt2aa = vt2
*     u1aa = u1
*
      CALL RELAX(IB1,IBD,IPT1,IPT2,VR2,VT2,SM2,XNPART)
*       Reject if innermost rmin binary sinks deeper
       IF(IB1.LT.0.D0)THEN
       ICOUNT = ICOUNT + 1
*       Set DTBACC to large value to enforce relaxation next time
       DTBACC(IB) = 1.D30
       GOTO 305
       END IF
*
       IF(ICOUNT.GT.0)THEN
       PRINT*,' repeat rx TIME,IB=',TIME,IB,' icount ',ICOUNT
       END IF
*
*     WRITE(6,1111)NAMEB(IB),QQ,TBIN(IB),TRXX,TRXX20,TX2,TX3
*1111 FORMAT(' N(IB)=',I6,' QQ=',1PD9.2,
*    *  ' tbin,trxloc,trxloc20,t2-3=',5D9.2)
*     WRITE(6,1112)NAMEB(IB),QQ,TRXX/TBIN(IB),TRXX20/TBIN(IB),TX2/TBIN(IB),
*    *      TX3/TBIN(IB)
*1112 FORMAT(' N(IB)=',I6,' QQ=',1PD9.2,
*    *  ' (trxloc,trxloc20,t2-3)/tbin=',4D9.2)
*     WRITE(6,1113)NAMEB(IB),QQ,RBSMIN,RBS,RBSMAX,RCORE,20.*RCORE
*1113 FORMAT(' N(IB)=',I6,' QQ=',1PD9.2,' rmin,rb,rmax,rc,20*rc=',5D9.2)
*     WRITE(6,1114)NAMEB(IB),QQ,ITRANS,ICORE,I2C,I5C,I10C,I20C,
*    *   IH100,IH10,IHALF
*1114 FORMAT(' N(IB)=',I6,' RR=',1PD9.2,' ITRANS=',I4,' IC,2,5,10,20=',
*    *    5I4,' IH100,10,1=',3I4)
*     WRITE(6,1115)NAMEB(IB),QQ,TRXLOC(ICORE),TRXLOC(I2C),TRXLOC(I5C),
*    *   TRXLOC(I10C),TRXLOC(I20C)
*1115 FORMAT(' N(IB)=',I6,' RR=',1PD9.2,' TRXLOC(IC,2,5,10,20)=',
*    *  5D9.2)
*     WRITE(6,1116),NAMEB(IB),QQ,TDYLOC(ICORE),TDYLOC(I2C),TDYLOC(I5C),
*    *   TDYLOC(I10C),TDYLOC(I20C)
*1116 FORMAT(' N(IB)=',I6,' RR=',1PD9.2,' TDYLOC(IC,2,5,10,20)=',
*    *  5D9.2)
*
       END IF
*
*       Reset IEV after printout in relax
*
       IF(IEV(IB).NE.0)IEV(IB)=0
*
       I=ISH(IB)
       I1=I-1
*
       XMBIN=BODY1(IB)+BODY2(IB)
       U2=PHI(I1)-PHTID+PHIBIN(I1)+
     * (PHI(I)+PHIBIN(I)-PHI(I1)-PHIBIN(I1))/(R(I)-R(I1))*(RB(IB)-R(I1))
       ETBIN=XMBIN*(C12*(VR(IB)*VR(IB)+VT(IB)*VT(IB))+U2)
       ESNEW=XMIND(J)*(C12*(VR2*VR2 + VT2*VT2)+U1)
*
*       TEST PRINT
*     IF(LRELAX)THEN
*     nrelax = nrelax + 1
*     XX=(DBLE(NRELAX-1)*XX+(ESNEW-ESOLD)/DABS(ESNEW))/DBLE(NRELAX)
*     print*,' nrelax,xx,de/e=',nrelax,xx,(ESNEW-ESOLD)/DABS(ESNEW)
*     print*,' i,density = ',i,xnpart
*     print*,'VRold,VTold,U1,ESold = ',Vr2aa,Vt2aa,u1aa,esold
*     print*,'VRnew,VTnew,U2,ESnew= ',Vr2,Vt2,u2,esnew
*     END IF
*     
*       XHEAT(J) = XHEAT(J) + ESNEW - ESOLD
       XFEED(J) = XFEED(J) + ESNEW - ESOLD
*       Set flag to reduce timestep for large relax kick
*      IF(DABS((ESNEW-ESOLD)/ESOLD).GT.5.D-2)THEN
*      IRX2 = IRX2 + 1
*      PRINT*,' Strong Relax Event IRX=',IRX,IRX2,IBR,' DE=',
*    *  ESNEW - ESOLD,
*    *  ' Single E=',ESOLD,' DE/E=',DABS((ESNEW-ESOLD)/ESOLD)
*      IEVENT=2
*      ELSE
*      IF(LRELAX)PRINT*,' Weak Relax Event DE=',ESNEW - ESOLD,
*    *  ' Single E=',ESOLD,' DE/E=',DABS((ESNEW-ESOLD)/ESOLD)
*      END IF
*   test only
*      IF(IMOD.EQ.50.AND.IB.EQ.100)THEN
*      XMASS(J) = XMASS(J) - 1.D2*XMIND(J)
*      XKMASS(J) = XKMASS(J) - 1.D2*XMIND(J)
*      XMERR = XMERR - 1.D2*XMIND(J)
*      PRINT*,' J=',J,' XMASS=',XMASS(J),' Mass taken away'
*      END IF
*
       IF(ESNEW.GT.0.D0)THEN
*   Take out relaxation heating if single star escapes.
*       XHEAT(J) = XHEAT(J) - ESNEW + ESOLD
       XFEED(J) = XFEED(J) - ESNEW + ESOLD
*
       XMASS(J) = XMASS(J) - XMIND(J)
       XKMASS(J) = XKMASS(J) - XMIND(J)
       XNRSES = XNRSES + 1.D0
       XMRSES = XMRSES + XMIND(J)
       XERSES = XERSES + ESNEW
       ISESC = 1
*
       PRINT*,' Single Star Escaper due to relaxation IB,IBR=',IB,IBR,
     *   ' N=',NAMEB(IB),' XMASS, XKMASS=',XMASS(J),XKMASS(J)
       PRINT*,' Single Star Escaper due to relaxation IB=',IB,
     *   ' RBOLD=',RBOLD,' RB=',RB(IB),' ESOLD,N=',ESOLD,ESNEW
*
       END IF
*
*      IF(MOD(IMOD,NPR).EQ.0)THEN
*      IF(LRELAX)THEN
*      PRINT*,' Binary After Relax ',IB,' N=',NAMEB(IB),
*    *        ' After Relax VR,VT=',VR(IB),VT(IB)
*      PRINT*,' Binary After Relax ',IB,' N=',NAMEB(IB),
*    *        ' After Relax ROLD,RB=',ROLD,RB(IB)
*      PRINT*,' Binary After Relax ',IB,' N=',NAMEB(IB),
*    *        ' CM Energy=',ETBIN,' Phi=',U1
*      PRINT*,' Single After Relax Star Eold=',ESOLD,
*    *     ' ENEW=',ESNEW,' DE=',ESNEW-ESOLD
*      END IF
*
*      CALL FLUSH(6)
*      END IF
*
       SIG=DSQRT(C13*(DEXP(X(I20+J,I))+2.D0*DEXP(X(I02+J,I)))/
     *     DEXP(X(I00+J,I)))
*        Do not harden binary if it is escaping
       IF(IPT1.GT.0)THEN
       IF(IPT1.EQ.1)PRINT*,' Binary newpos1 relax escape ',IB,
     *    ' N=',NAMEB(IB),' EB=',EB(IB),' NESC=',NESC
       IF(IPT1.EQ.2)PRINT*,' Binary newpos2 relax escape ',IB,
     *    ' N=',NAMEB(IB),' EB=',EB(IB),' NESC=',NESC
       NESC = NESC + 1
       IESC(NESC) = IB
       NAMEB(IB) = - ABS(NAMEB(IB))
       XNRBES = XNRBES + 1.D0
       XMRBES = XMRBES + BODY1(IB) + BODY2(IB)
       XERBES = XERBES + ETBIN
       XBRBES = XBRBES + EB(IB)
       IPT1 = 0
       GOTO 300
       END IF
*
*       Use typical velocities of single stars for total probability
      VR1 = DSQRT(DEXP(X(I20+J,I) - X(I00+J,I)))
      VT1 = 2.D0*DSQRT(DEXP(X(I02+J,I) - X(I00+J,I)))
*
      TWOPI = 2.D0*PI
      FI = TWOPI*RANF()
      SINFI = SIN(FI)
      COSFI = COS(FI)
*
      VRELX = VR1-VR(IB)
      VRELY = VT1*SINFI
      VRELZ = VT1*COSFI-VT(IB)
      VRELT2 = VRELX**2+VRELY**2+VRELZ**2
      VRELT = DSQRT(VRELT2)
*
      BMMAX = MAX(BODY1(IB),BODY2(IB))
      XM = XMIND(J)
      XMTRIP = XM + BODY1(IB) + BODY2(IB)
       VC2 = BMMAX*XMTRIP/(XM*SEMIA(IB))
      VVC2 = DSQRT(VRELT2/VC2)
*
      IF(LPLAN) THEN
*    Use 100*a*(1+e) for planetary system run
         PMAX = 100.D0*SEMIA(IB)*(1.0D0 + ECC(IB))
      ELSE
         PMAX = SEMIA(IB)*(0.6D0*(1.0D0 + ECC(IB)) + 4.0D0/VVC2)
      END IF
*
      XNRHO = DEXP(X(I00+J,I))/XMIND(J)
*
      PBSNEW = XNRHO*PI*PMAX**2*VRELT/BREM
*
*        Probability of single star binary interaction
       PBSOLD=5.D0*ASPI*GRAV*GRAV*XMIND(J)**2*DEXP(X(I00+J,I))/BREM/
     *      4.D0/SIG/EB(IB)*PI*C12
       PBSTOT = PBSNEW
*       IF(EB(IB)/EKT.LT.3.D0)PBSTOT=0.D0
       XRAN = RANF()
*      test
*       PBSTOT = 0.D0
*   
*        If interaction takes place select energy change
*        by Monte-Carlo rejection technique (LS(4)=.FALSE.)
*        or calculate full three body interaction (LS(4)=.TRUE.)
*
*        Force interaction for test
*      LKICK=MOD(IMOD,5).EQ.0.AND.MOD(IB,100).EQ.0
*      IF(LKICK)PBSTOT=1.D0
*      IF(MOD(IMOD,5).EQ.0.AND.IMOD.GT.50)PBSTOT=1.D0
       IF(XRAN.LT.PBSTOT)THEN
*
       ITOT3 = ITOT3 + 1
       WRITE(6,556)TIME,ITOT3,IDISS3,PBSNEW,PBSOLD
 556   FORMAT(1X,10('-'),'START TRIPLE',5('-'),'T=',1P,D10.2,5('-'),
     *  ' ITOT3=',I8,' IDISS3=',I8,10('-'),/,' PBSNEW,OLD=',1P,2D12.4)
       VROLD=VR(IB)
       VTOLD=VT(IB)
*
      XMBIN = BODY1(IB) + BODY2(IB)
      XMTRIP = XMBIN + XMIND(J)
*
      ECMOLD = XMBIN*(C12*(VR(IB)**2+VT(IB)**2)+U2)
      EBOLD = EB(IB)
*
*       Pick star for encounter
       I=ISH(IB)
*        In case of multi-mass system pick component
       J=ICO(IB)
*        Marker for 3b interaction printout
       IEV(IB) = 1
*        Pick velocity of single star for encounter
*        needed for TRMAIN and/or to decide for single escaper due to kick
      SIGR2 = DEXP(X(I20+J,I)-X(I00+J,I))
      SIGT2 = DEXP(X(I02+J,I)-X(I00+J,I))
      VESC = DSQRT(-2.D0*(PHI(I)-PHTID+PHIBIN(I)))
*
 405  CONTINUE
*
      VRTEST = VESC*RANF()
      VTTEST = VESC*RANF()
      VPTEST = VESC*RANF()
      VR2 = VRTEST**2
      VT2 = VTTEST**2+VPTEST**2
      FV = DEXP(-C12*VR2/SIGR2 - C12*VT2/SIGT2)
*
      QV = RANF()
*
      IF(QV.GT.FV)GOTO 405
*
      VTOT=DSQRT(VR2+VT2)
      IF(VTOT.GT.VESC)GOTO 405
*
      ESOLD = C12*XMIND(J)*(VR2+VT2) + XMIND(J)*U2
*
      VR2=DSQRT(VR2)
      VT2=DSQRT(VT2)
      IF(RANF().LT.0.5)VR2=-VR2
*        Take cross sections or direct integration LS(4)
       IF(.NOT.LS(4))THEN
*
 1002  CONTINUE
       IREJ=0
       IACC=0
*        Calculate f(z)
       XRZ = PI/2.D0*RANF()
       FXR = DLOG10(DCOS(XRZ))
       FXR=1.D1**(6.D0*FXR)
*
       XRAN = RANF()
*        Accept interaction
       IF(XRAN.LT.FXR)THEN
*
          IACC=1
          DELTA = DTAN(XRZ)**2
*       Apply heating and reduce time step for large kick        
       IF(DELTA.GT.0.5D0)IEVENT = 2
*       Use old binary energy here, corrected Sep. 2000
       XKICK(J)=XKICK(J)+XMBIN/XMTRIP*DELTA*EBOLD
       XHEAT(J)=XHEAT(J)+XMBIN/XMTRIP*DELTA*EBOLD
*
          VB = DELTA*EB(IB)*XMIND(J)/(XMIND(J)+XMBIN)
          EB(IB)=EB(IB)*(1.D0+DELTA)
*       Calculate Kick on Binary
*
          CALL KICK(IB,VB)
*
       XM = XMIND(J)
       XMBIN = BODY1(IB) + BODY2(IB)
       XMTRIP = XM + XMBIN
       ETBIN=XMBIN*(C12*(VR(IB)*VR(IB)+VT(IB)*VT(IB))+U2)
       ESNEW = ESOLD + XMBIN/XMTRIP*DELTA*EBOLD
       ESNGL1 = 0.D0
       ESNGL2 = 0.D0
       ESNGL3 = 0.D0
*
       ELSE
          IREJ=IREJ+1
       END IF
*
       IF(IACC.EQ.0)GOTO 1002
*
       ELSE
*
*       Three body integration for hard binary encounters
*
      XR=RIND(J)
      XM=XMIND(J)
      EKT = C13*(SIGR2 + 2.D0*SIGT2)*XM
      XMOLD = XM
      VROLD = VR(IB)
      VTOLD = VT(IB)
      VR1 = VR2
      VT1 = VT2
*
       TCOMP1 = ETIME(TARRAY)
       CALL TRMAIN(IB,J,VR1,VT1,XM,XR,EKT,TIME)
       TCOMP2 = ETIME(TARRAY)
       tttrip = tttrip + TCOMP2-TCOMP1
       WRITE(6,1888)TCOMP2-TCOMP0,TCOMP2-TCOMP1,tttrip,ttquad
 1888  FORMAT(1X,' CPU tot 3bnow trip quad=',1P,4D13.5)
*
      XMBIN = BODY1(IB) + BODY2(IB)
      XMTRIP = XM + XMBIN
      ESB = C12*XMBIN*(VR(IB)**2+VT(IB)**2)
      ESS = C12*XM*(VR1**2+VT1**2)
      ESNEW = ESS + XM*U2
      ETBIN = ESB + XMBIN*U2
*
*       First treat case of binary dissolution
*       In case of planets (mass ratio < 1.e-3) use v^2 criterion
      BMRAT = MIN(BODY1(IB),BODY2(IB))/
     &              MAX(BODY1(IB),BODY2(IB))
      IF(BMRAT.LT.1.D-3)THEN
      BMKT = 0.1D0*BMRAT
      ELSE
      BMKT = 1.D0
      BMRAT = 1.D0
      END IF
*
      IF(EB(IB).LT.0.1D0*BMKT*EKT)THEN
*
        ESBOLD = C12*XMBIN*(VROLD**2 + VTOLD**2)
        DELTA = XMTRIP/XMBIN*(BMRAT*ESBOLD - EBOLD)/EBOLD
*
        IF(DELTA.GT.0.5D0)IEVENT = 2
*       Energy transfer to gas model is c.m. energy of dissolved bin. - EBOLD
        XKICK(J)=XKICK(J) + ESBOLD - EBOLD
        XHEAT(J)=XHEAT(J) + ESBOLD - EBOLD    
*
        XN3BIS = XN3BIS + 1.D0
        XM3BIS = XM3BIS + BODY1(IB) + BODY2(IB)
        XE3BIS = XE3BIS + ECMOLD
        XB3BIS = XB3BIS + EB(IB)
*
        IEV(IB) = -4
*
        IDISS3 = IDISS3 + 1
       WRITE(6,777)TIME,IB,NAMEB(IB),IDISS3,EB(IB),EB(IB)/BMRAT/EKT
 777   FORMAT(1X,' Dissolution by Triple T,IB2,N,IDISS3,EB,',
     *   'EB/BMRAT/kT=',1P,D13.5,3I6,2D13.5)
         XRAN = RANF()
*       If binary is dissolved distribute energy randomly between singles.
         ESB1 = XRAN*ESB
         ESB2 = (1.D0-XRAN)*ESB
         ETBIN = 0.D0
         ESNEW = 0.D0
         ESNGL1 = ESB1 + BODY1(IB)*U2
         ESNGL2 = ESB2 + BODY2(IB)*U2
         ESNGL3 = ESNEW
       ELSE
*
*        Now case of binary hardening
*
         DELTA = (EB(IB)-EBOLD)/EBOLD
*       Energy transfer to gas model is fraction of change of EB
         XKICK(J)=XKICK(J)+XMBIN/XMTRIP*(EB(IB)-EBOLD)
         XHEAT(J)=XHEAT(J)+XMBIN/XMTRIP*(EB(IB)-EBOLD)
         ESNGL1 = 0.D0
         ESNGL2 = 0.D0
         ESNGL3 = 0.D0
       END IF
*        End three-body integration
      END IF
*
*       Force Relaxation at next timestep
       DTBACC(IB) = 1.D30
*       Force mass update at next timestep
       IBKICK = 1
*
       WRITE(6,775)IB,NAMEB(IB),TIME,EB(IB)/EKT,DELTA
 775   FORMAT(1X,' Bin ',I6,' N=',I6,' T,EB/kT after kick=',1P,
     *  2D13.5,' DELTA=',D13.5)
       WRITE(6,776)TIME,IB,ETBIN,ESNEW,ESB,ESS
 776   FORMAT(' T,IB=',1P,D13.5,I6,' C.M.: Etot bin,s=',
     *    2D13.5,' Ekin bin,s=',2D13.5)
       CALL FLUSH(6)
*
*       First treat cases of dissolved binary and 1,2,3 stars may escape.
*
*       Add mass to single stars if they do not escape (triple dissolution)
       IF(ESNGL1.LT.0.D0)THEN
       XMASS(J) = XMASS(J) + BODY1(IB)
       XKMASS(J) = XKMASS(J) + BODY1(IB)
       ITRIP = 1
       WRITE(6,774)IB,TIME,ESNGL1,ESB1
 774   FORMAT(1X,' Single Star 1 from 3b dissolution of ',I6,
     *  ' T,Etot,Ekin=',1P,3D13.5)
      XN3SBI = XN3SBI + 1.D0
      XM3SBI = XM3SBI + BODY1(IB)
      XE3SBI = XE3SBI + ESB1
*       Delayed removal of binary only if both singles remain in system.
       END IF
*
       IF(ESNGL1.GT.0.D0)THEN
*       Take out kick heating in case of single star escaper
       XKICK(J)=XKICK(J) - ESB1
       XHEAT(J)=XHEAT(J) - ESB1
       WRITE(6,773)IB,TIME,ESNGL1,ESB1
 773   FORMAT(1X,' Single Star Escaper 1 from 3b dissolution of ',I6,
     *  ' T,Etot,Ekin=',1P,3D13.5)
*
       IEV(IB)=IEV(IB)+1
*
      XN4SES = XN4SES + 1.D0
      XM4SES = XM4SES + BODY1(IB)
      XE4SES = XE4SES + ESB1
       END IF
*
       IF(ESNGL2.LT.0.D0)THEN
       XMASS(J) = XMASS(J) + BODY2(IB)
       XKMASS(J) = XKMASS(J) + BODY2(IB)
       ITRIP = 1
       WRITE(6,772)IB,TIME,ESNGL2,ESB2
 772   FORMAT(1X,' Single Star 2 from 3b dissolution of ',I6,
     *  ' T,Etot,Ekin=',1P,3D13.5)
      XN3SBI = XN3SBI + 1.D0
      XM3SBI = XM3SBI + BODY2(IB)
      XE3SBI = XE3SBI + ESB2
*       Delayed removal of binary only if both singles remain in system.
       END IF
*
       IF(ESNGL2.GT.0.D0)THEN
*       Take out kick heating in case of single star escaper
       XKICK(J)=XKICK(J) - ESB2
       XHEAT(J)=XHEAT(J) - ESB2
       WRITE(6,771)IB,TIME,ESNGL2,ESB2
 771   FORMAT(1X,' Single Star Escaper 2 from 3b dissolution of ',I6,
     *  ' T,Etot,Ekin=',1P,3D13.5)
*
       IEV(IB)=IEV(IB)+1
*
      XN4SES = XN4SES + 1.D0
      XM4SES = XM4SES + BODY2(IB)
      XE4SES = XE4SES + ESB2
       END IF
*
       IF(ESNGL3.GT.0.D0)THEN
*
       XMASS(J) = XMASS(J) - XMOLD
       XKMASS(J) = XKMASS(J) - XMOLD
      XN4SES = XN4SES + 1.D0
      XM4SES = XM4SES + XM
      XE4SES = XE4SES + ESS
*       Take out kick heating in case of single star escaper
       XKICK(J)=XKICK(J) - ESS
       XHEAT(J)=XHEAT(J) - ESS
       WRITE(6,770)IB,TIME,ESNGL3,ESS
 770   FORMAT(1X,' Single Star Escaper 3 from 3b dissolution of ',I6,
     *  ' T,Etot,Ekin=',1P,3D13.5)
       END IF
*       Take care if an exchange has taken place
       IF(ESNGL3.LT.0.D0)THEN
       XMASS(J) = XMASS(J) + XM - XMOLD
       XKMASS(J) = XKMASS(J) + XM - XMOLD
       END IF
*
*       Delayed removal of binary only if both singles remain in system.
       IF(IEV(IB).EQ.-4)THEN
       NESC2 = NESC2 + 1
       IESC2(NESC2) = IB
*
       ELSE IF(IEV(IB).GT.-4.AND.IEV(IB).LT.-1) THEN
*       Immediate removal of binary if one single star escapes.
       NESC = NESC + 1
       IESC(NESC) = IB
       NAMEB(IB) = - ABS(NAMEB(IB))
*
       END IF
*
*       Now treat cases if the binary remains bound, but may escape (or single).
*
       IF(ESNEW.GT.0.D0)THEN
*
       XMASS(J) = XMASS(J) - XMOLD
       XKMASS(J) = XKMASS(J) - XMOLD
      XN3SES = XN3SES + 1.D0
      XM3SES = XM3SES + XM
      XE3SES = XE3SES + ESNEW
*       Take out kick heating in case of single star escaper
       XHEAT(J) = XHEAT(J) - XMBIN/XMTRIP*DELTA/(1.D0+DELTA)*EB(IB)
       XKICK(J) = XKICK(J) - XMBIN/XMTRIP*DELTA/(1.D0+DELTA)*EB(IB)
*
       ISESC = 1
*
       PRINT*,' Single Star Escaper related to Binary IB=',IB,
     *   ' N=',NAMEB(IB)
       END IF
*
       IF(ESNEW.LT.0.D0)THEN
*       Take into account exchange, but binary remains bound.    
       XMASS(J) = XMASS(J) + XM - XMOLD
       XKMASS(J) = XKMASS(J) + XM - XMOLD
*
       END IF
*
*      PRINT*,' Binary ',IB,' N=',NAMEB(IB),
*     *       ' After Kick VR,VT=',VR(IB),VT(IB)
       IF(ETBIN.GT.0.D0)THEN
       NESC = NESC + 1
       WRITE(6,769)IB,NAMEB(IB),TIME,EB(IB),NESC
 769   FORMAT(' Bin ',I6,' N=',I6,' triple kick escape T,EB,NESC=',
     *    1P,2D13.5,I6)
      XN3BES = XN3BES + 1.D0
      XM3BES = XM3BES + BODY1(IB) + BODY2(IB)
      XE3BES = XE3BES + ETBIN
      XB3BES = XB3BES + EB(IB)
       IESC(NESC) = IB
       NAMEB(IB) = - ABS(NAMEB(IB))
       GOTO 300
       END IF
*
       IF(IEVENT.EQ.0)IEVENT = 1 
*
*      ICASE = ICASE + 1
*      AVDE = (AVDE*DBLE(ICASE-1) + DELTA)/DBLE(ICASE)
*      PRINT*,' Binary # ',IB,' N=',NAMEB(IB),
*     *        ' hardened by fraction=',DELTA,' now EB=',
*     *   EB(IB),' Average now=',AVDE,' R=',RB(IB),' ICASE=',ICASE
*      PRINT*,' Binary Amount of Heating =',DELTA*EB(IB)/(1.D0+DELTA),
*     *   ' TIME=',TIME
*
*       Call newpos to determine new binary mass distribution via orbit
      DEFANG = 0.D0
      IPT1 = 0
      CALL NEWPOS(IB,IPT1,DEFANG)
*
*       Newpos could detect binary close to escape at r=rnj - should escape
       IF(IPT1.GT.0)THEN
*
       NESC = NESC + 1
       IF(IPT1.EQ.1)PRINT*,' Binary newpos1 kick escape ',IB,' N=',
     *        NAMEB(IB),' EB=',EB(IB),' NESC=',NESC
       IF(IPT1.EQ.2)PRINT*,' Binary newpos2 kick escape ',IB,' N=',
     *        NAMEB(IB),' EB=',EB(IB),' NESC=',NESC
      XN3BES = XN3BES + 1.D0
      XM3BES = XM3BES + BODY1(IB) + BODY2(IB)
      XE3BES = XE3BES + ETBIN
      XB3BES = XB3BES + EB(IB)
       IESC(NESC) = IB
       NAMEB(IB) = - ABS(NAMEB(IB))
       IPT1 = 0
       END IF
*
       PRINT*,'----------END TRIPLE-------------------------'
       END IF
*
  300  CONTINUE
       END IF
*
       IF(NESC.GT.0)THEN
       IBESC = 1
       PRINT*,' Removal of Binaries after TRIPLE ',
     *   ' NESC=',NESC,' list:',(IESC(III),III=1,NESC)
       CALL ESCAPE(NESC,IESC)
       END IF
*
*         Compute Binary Potential to save changes for energy balance
*         and for initialization in case of initial model
       IPOT = 0
       CALL PHIMAS(IPOT)
*
*        binary binary relaxation
*        only in single mass case
       XMBIN = 2.D0*BODY1(1) 
*        determine background binary density on radial mesh
       DO 340 I = 2,NJ
       I1 = I-1
       RI=R(I)
       RIM=R(I1)
       XM=XMRBIN(I)
       XMM=XMRBIN(I1)
       DM = XM - XMM
       DVOL=PI43*(RI**3-RIM**3)
*        set density to zero inside and outside binary mass distribution
       IF(DABS(XM).LT.XMMIN.OR.DABS(XM-IBIN*XMBIN).LT.XMMIN)THEN
           RHOBIN(I) = 0.D0
       ELSE
*        set density to zero if there is no binary in the middle
           IF(DABS(DM).LT.XMMIN)THEN
               RHOBIN(I)=0.D0
           ELSE
               RHOBIN(I)=(XM-XMM)/DVOL
*              RHOBIN(I) = XMRBIN(I)/Pi43/R(I)**3
           END IF
       END IF
*
 340   CONTINUE
*
       DO 350 IX=1,IBIN,2
*       Use radially sorted binaries in pairs - exclude binaries to be dissolved
       IB1=INAME(IX)
       IF(IEV(IB1).LT.-1)GOTO 350
       IF(IX+1.GT.IBIN)GOTO 350
       IB2=INAME(IX+1)
       IF(IEV(IB2).LT.-1)GOTO 350
*        Determine harder binary
       IF(EB(IB2).GT.EB(IB1))THEN
       IB1=INAME(IX+1)
       IB2=INAME(IX)
       END IF
*        No relaxation if one of the binaries has minimum rmin
*      IF(IB1.EQ.IMCRIT.OR.IB2.EQ.IMCRIT)THEN
*      PRINT*,' Binaries ',IB1,IB2,' excl relax IMCRIT ',IMCRIT
*      GOTO 350
*      END IF
*
       I=ISH(IB1)
       I1=I-1
       K=ISH(IB2)
       K1=K-1
*
       LRELAX=INAMEX(IB1)+INAMEX(IB2).GT.0
*
*        binary binary Relaxation due?
       IF(LRELAX)THEN
*        Compute local binary density
       XNPART = RHOBIN(I)/XMBIN
*
       CALL RELAX(IB1,IB2,IPT1,IPT2,VR2,VT2,SM2,XNPART)
*
*       Reset IEV after printout in relax
*
       I=ISH(IB1)
       I1=I-1
       K=ISH(IB2)
       K1=K-1
*       Reset relaxation interval for both binaries
       DTBACC(IB1) = C12*(TBIN(IB1)-TBIN(IB2))
       DTBACC(IB2) = C12*(TBIN(IB2)-TBIN(IB1))
*
       U1=PHI(I1)-PHTID+PHIBIN(I1)+
     *(PHI(I)+PHIBIN(I)-PHI(I1)-PHIBIN(I1))/(R(I)-R(I1))*(RB(IB1)-R(I1))
       U2=PHI(K1)-PHTID+PHIBIN(K1)+
     *(PHI(K)+PHIBIN(K)-PHI(K1)-PHIBIN(K1))/(R(K)-R(K1))*(RB(IB2)-R(K1))
       ECMOL1 = XMBIN1*(C12*(VR(IB1)**2+VT(IB1)**2)+U1)
       ECMOL2 = XMBIN2*(C12*(VR(IB2)**2+VT(IB2)**2)+U2)
*
       QQ1 = RBMAX(IB1)/RBMIN(IB1)
       QQ2 = RBMAX(IB2)/RBMIN(IB2)
*     WRITE(6,1123)NAMEB(IB1),QQ1,RBMIN(IB1),RB(IB1),RBMAX(IB1)
*1123 FORMAT(' Bin-Rx N(IB1)=',I6,' QQ1=',1PD9.2,' rmin,rb,rmax=',3D9.2)
*     WRITE(6,1124)NAMEB(IB2),QQ2,RBMIN(IB2),RB(IB2),RBMAX(IB2)
*1124 FORMAT(' Bin-Rx N(IB2)=',I6,' QQ2=',1PD9.2,' rmin,rb,rmax=',3D9.2)
*
       IF(IPT1.GT.0)THEN
       IF(IPT1.EQ.1)PRINT*,' Binary 1 newpos1 4rel escape ',IB1,' N=',
     *        NAMEB(IB1),' EB=',EB(IB1),' NESC=',NESC
       IF(IPT1.EQ.2)PRINT*,' Binary 1 newpos2 4rel escape ',IB1,' N=',
     *        NAMEB(IB1),' EB=',EB(IB1),' NESC=',NESC
      XNRBES = XNRBES + 1.D0
      XMRBES = XMRBES + BODY1(IB1) + BODY2(IB1)
      XERBES = XERBES + ECMOL1
      XBRBES = XBRBES + EB(IB1)
       NESC = NESC + 1
       IESC(NESC) = IB1
       NAMEB(IB1) = - ABS(NAMEB(IB1))
       IPT1 = 0
       GOTO 350
       END IF
*
       IF(IPT2.GT.0)THEN
       IF(IPT1.EQ.1)PRINT*,' Binary 2 newpos1 4rel escape ',IB2,' N=',
     *        NAMEB(IB2),' EB=',EB(IB2),' NESC=',NESC
       IF(IPT1.EQ.2)PRINT*,' Binary 2 newpos2 4rel escape ',IB2,' N=',
     *        NAMEB(IB2),' EB=',EB(IB2),' NESC=',NESC
      XNRBES = XNRBES + 1.D0
      XMRBES = XMRBES + BODY1(IB2) + BODY2(IB2)
      XERBES = XERBES + ECMOL2
      XBRBES = XBRBES + EB(IB2)
       NESC = NESC + 1
       IESC(NESC) = IB2
       NAMEB(IB2) = - ABS(NAMEB(IB2))
       IPT2 = 0
       GOTO 350
       END IF
*
       END IF
*        End bin-bin relaxation
 350   CONTINUE
*
*        Now check for strong encounters
       IF(LS(24))THEN
*
       DO 355 IX=1,IBIN,2
*
       IF(IX+1.GT.IBIN)GOTO 355
*
*        Use radially sorted binaries in pairs
*        Exclude dissolved or delayed binaries to be dissolved
*        Take safe double condition on IEV and NAMEB
       IB1=INAME(IX)
       IF(IEV(IB1).LT.-1)GOTO 355
       IF(NAMEB(IB1).LT.0)GOTO 355
       IB2=INAME(IX+1)
       IF(IEV(IB2).LT.-1)GOTO 355
*        Determine harder binary
       IF(EB(IB2).GT.EB(IB1))THEN
       IF(NAMEB(IB2).LT.0)GOTO 355
       IB1=INAME(IX+1)
       IB2=INAME(IX)
       END IF
*
       I=ISH(IB1)
       I1=I-1
       K=ISH(IB2)
       K1=K-1
*        In case of multi-mass system pick component
       JB1=ICO(IB1)
       JB2=ICO(IB2)
       J=JB2
*        Determine sigmas for kT
      SIGR21 = DEXP(X(I20+JB1,I)-X(I00+JB1,I))
      SIGT21 = DEXP(X(I02+JB1,I)-X(I00+JB1,I))        
      SIGR22 = DEXP(X(I20+JB2,I)-X(I00+JB2,I))
      SIGT22 = DEXP(X(I02+JB2,I)-X(I00+JB2,I))    
      XM1 = XMIND(JB1)
      XM2 = XMIND(JB2)
      EKT1 = C13*(SIGR21 + 2.D0*SIGT21)*XM1
      EKT2 = C13*(SIGR22 + 2.D0*SIGT22)*XM2
*        Compute encounter probability
       XRAN = RANF()
       FI = 2.D0*PI*XRAN
       VTX2 = VT(IB2)*COS(FI)
       VTY2 = VT(IB2)*SIN(FI)
       VREL=DSQRT((VR(IB1)-VR(IB2))**2+VTY2**2+(VTX2-VT(IB1))**2)
       XMBIN1=BODY1(IB1)+BODY2(IB1)
       XMBIN2=BODY1(IB2)+BODY2(IB2)
*      RCRIT4=2.5D0*SEMIA(IB2)
*      RCRIT4=5.0D0*SEMIA(IB2)
*
*     compute  RCRIT from PMAX.
*     calculate PMAX using formula given by Hut and Bahcall ApJ 268, 1983         
*     and Bacon and Sigurdsson astro-ph96/03036
*
       XXMT = XMBIN1 + XMBIN2
       XXMIU = (XMBIN1*XMBIN2)/XXMT
       VC2 = 2.D0*GRAV/XXMIU*(EB(IB1) + EB(IB2))
       VC2 = DSQRT(VREL/VC2)
*
       IF(LPLAN) THEN
*    Use 100*a*(1+e) for planetary system run
          PMAX = 100.D0*SEMIA(IB2)*(1.0D0 + ECC(IB2))
       ELSE
          PMAX = (5.D0/VC2 + 0.6D0*(1.D0 + ECC(IB2)))*SEMIA(IB2)
       END IF
*
*     New to speed up the program
*       IF(PMAX.GT.10.D0*SEMIA(IB2)) PMAX = 10.D0*SEMIA(IB2) 
*
       XXXY = GRAV*XXMT/(VREL*VREL)
       RCRIT4 = XXXY*(DSQRT(1.D0+(PMAX/XXXY)**2) - 1.D0)
*
       AREA=PI*RCRIT4**2*(1.D0+2.D0*GRAV*(XMBIN1+XMBIN2)/RCRIT4/VREL**2)
       PBB=RHOBIN(I)/XMBIN1*VREL*AREA/BREM
*
       IMIN=IBSTA(IB2)+1
       IMAX=IBSTA(IB2)+ILEN(IB2)
       RMIN=RBMIN(IB2)
       RMAX=RBMAX(IB2)
*       TTRB=(RMAX-RMIN)/DABS(VR(IB2))
       TTRB = TORB(IB1)
       IF(RB(IB1).GT.RMIN.AND.RB(IB1).LT.RMAX)THEN
       PBSTOD=GRAV*(XMBIN1+XMBIN2)/VREL**2*RCRIT4/
     *   2.D0/RB(IB1)**2/BREM/TTRB*VREL/ABS(VR(IB1))
       ELSE
       PBSTOD=0.D0
       END IF
*       
      XRAN=RANF()
*     test
*       PBB =0.D0
*       PBB1 = PBB
*       PBB = PBSTOD
*
*        Compare with random number for encounter probability
       IF(XRAN.LT.PBB)THEN
       ITOT4 = ITOT4 + 1
*
       WRITE(6,557)TIME,ITOT4,IDISS4,PBB,PBSTOD
 557   FORMAT(1X,10('-'),'START QUAD  ',5('-'),'T=',1P,D10.2,5('-'),
     *  ' ITOT4=',I8,' IDISS4=',I8,10('-'),/,' PBB,PBSTOD=',1P,2D12.4)
*
*      PRINT*,' quad PBB,PBSTOD,RHOB=',PBB1,PBSTOD,RHOBIN(I),' GFOC=',
*    *    AREA/PI/RCRIT4**2  
*
       IEV(IB1) = 4
       IEV(IB2) = 4
       IEVENT = 4
       U1=PHI(I1)-PHTID+PHIBIN(I1)+
     *(PHI(I)+PHIBIN(I)-PHI(I1)-PHIBIN(I1))/(R(I)-R(I1))*(RB(IB1)-R(I1))
       U2=PHI(K1)-PHTID+PHIBIN(K1)+
     *(PHI(K)+PHIBIN(K)-PHI(K1)-PHIBIN(K1))/(R(K)-R(K1))*(RB(IB2)-R(K1))
       ECMOL1 = XMBIN1*(C12*(VR(IB1)**2+VT(IB1)**2)+U1)
       ECMOL2 = XMBIN2*(C12*(VR(IB2)**2+VT(IB2)**2)+U2)
*
*       IMIN=IBSTA(IB2)+1
*       IMAX=IBSTA(IB2)+ILEN(IB2)
*       RMIN=RBMIN(IB2)
*       RMAX=RBMAX(IB2)
*       TTRB=(RMAX-RMIN)/DABS(VR(IB2))
*       IF(RB(IB1).GT.RMIN.AND.RB(IB1).LT.RMAX)THEN
*       PBSTOD=GRAV*(XMBIN1+XMBIN2)/VREL**2*RCRIT4/
*     *   2.D0/RB(IB1)**2/BREM/TTRB
*       ELSE
*       PBSTOD=0.D0
*       END IF
*
*       PRINT*,' quad PBB,PBSTOD,RHOB=',PBB,PBSTOD,RHOBIN(I),' GFOC=',
*     *    AREA/PI/RCRIT4**2
*
       EBOLD1 = EB(IB1)
       EBOLD2 = EB(IB2)
       EBSUM = EB(IB1) + EB(IB2)
*        Start determination of close encounters
       IF(.NOT.LS(28))THEN
*        Take approximate cross sections see Gao et al. 1991, ApJ 370, 567
 2002  CONTINUE
       IREJ=0
       IACC=0
*        Calculate f(z)
       XRW = RANF()
       FXW = XRW**C34
*
       XRAN = RANF()
*        Accept interaction
       IF(XRAN.LT.FXW)THEN
*
       IACC=1
       DELTA = DSQRT(2.D0/7.D0*(1.D0-XRW)/XRW)
*       One quarter of the energy to the binary
       ERCOIL = DELTA*EBSUM
       VB = ERCOIL*C12*XMBIN1/(XMBIN1+XMBIN2)
       EB(IB1) = EBSUM*(1.D0+DELTA)
       IEV(IB1) = 4
       IEVENT = 4
*       Three quarters of the energy to the singles
       XKICK(J)=XKICK(J)+C34*ERCOIL
       XHEAT(J)=XHEAT(J)+C34*ERCOIL
*      PRINT*,' Quad IB=',IB1,
*    *   ' N=',NAMEB(IB1),' Heat added=',C34*ERCOIL
*
*       Calculate Kick on Binary
*
       CALL KICK(IB1,VB)
*       Force Relaxation at next timestep
       DTBACC(IB1) = 1.D30
*       Force Mass Update at next timestep
       IBKICK = 1
       ECMB1 = C12*XMBIN1*(VR(IB1)**2+VT(IB1)**2)
       ECMB2 = C12*XMBIN2*(VR(IB2)**2+VT(IB2)**2)
       ETBIN1 =  ECMB1 + XMBIN1*U1
*       Put second binary energy to zero because here it is always removed.
       ETBIN2 = 0.D0
       BMRAT = 1.D0
*
       WRITE(6,1011)IB1,NAMEB(IB1),TIME,EB(IB1),EBOLD1,DELTA,ERCOIL
 1011  FORMAT(1X,' Bin ',I6,' N=',I6,' T,EB aft/bef quad kick=',1P,
     *  3D13.5,' DELTA,ERCOIL=',2D13.5) 
       IDISS4 = IDISS4 + 1
       WRITE(6,999)TIME,IB2,NAMEB(IB2),IDISS4,EBOLD2,
     &  EBOLD2/BMRAT/EKT2
 999   FORMAT(1X,' Dissolution by Quad T,IB2,N,IDISS4,EBold,',
     *   'EBold/BMRAT/kT=',1P,D13.5,3I6,2D13.5)
       CALL FLUSH(6)
*
*       Dissolve second binary
      XN4BIS = XN4BIS + 1.D0
      XM4BIS = XM4BIS + BODY1(IB2) + BODY2(IB2)
      XE4BIS = XE4BIS + ECMOL2
      XB4BIS = XB4BIS + EB(IB2)
*
       IEV(IB2) = -4
*
*       Two single stars get disrupted binaries cm-energy
*       plus fraction of 3/4*DeltaE
       XRAN=RANF()
       ESNGL1=XRAN*(ECMOL2+C34*ERCOIL)
       ESNGL2=(1.D0-XRAN)*(ECMOL2+C34*ERCOIL)
*       Do not use single star energies from 4b integration
*       (compare escaper check below)
       ESB1 = 0.D0
       ESB2 = 0.D0
*       Set energy of dissolved binary to -infty for single escaper condition
       EB(IB2) = -1.D30
*                                                     
       ELSE
*       Monte-Carlo rejection for close encounter
       IREJ=IREJ+1
       END IF
*
       IF(IACC.EQ.0)GOTO 2002
*
*       Four body integration for hard binary encounters
*       instead of approximate cross sections
       ELSE
*
      VROLD1 = VR(IB1)
      VTOLD1 = VT(IB1)    
      VROLD2 = VR(IB2)
      VTOLD2 = VT(IB2)    
*
       TCOMP1 = ETIME(TARRAY)
       CALL QUMAIN(IB1,IB2,EKT1,EKT2,TIME)
       TCOMP2 = ETIME(TARRAY)
       ttquad = ttquad + TCOMP2-TCOMP1
       WRITE(6,888)TCOMP2-TCOMP0,TCOMP2-TCOMP1,tttrip,ttquad
 888   FORMAT(1X,' CPU tot 4bnow trip quad=',1P,4D13.5)
*
*       End kick determination by kick or qumain
      XMBIN1 = BODY1(IB1) + BODY2(IB1)
      XMBIN2 = BODY1(IB2) + BODY2(IB2)
      XMQUAD = XMBIN1 + XMBIN2
      ECMB1 = C12*XMBIN1*(VR(IB1)**2+VT(IB1)**2)
      ECMB2 = C12*XMBIN2*(VR(IB2)**2+VT(IB2)**2) 
      ETBIN1 = ECMB1 + XMBIN1*U1
      ETBIN2 = ECMB2 + XMBIN2*U2
*       Dissolve second binary
*       In case of planets (mass ratio < 1.e-3) use v^2 criterion
      BMRAT = MIN(BODY1(IB2),BODY2(IB2))/
     &              MAX(BODY1(IB2),BODY2(IB2))
      IF(BMRAT.LT.1.D-3)THEN
      BMKT = 0.1D0*BMRAT
      ELSE
      BMKT = 1.D0
      BMRAT = 1.D0
      END IF
*
       IF(EB(IB2).LT.0.1D0*BMKT*EKT2)THEN
*
       ECMOL2 = C12*XMBIN2*(VROLD2**2 + VTOLD2**2)
       XN4BIS = XN4BIS + 1.D0
       XM4BIS = XM4BIS + BODY1(IB2) + BODY2(IB2)
       XE4BIS = XE4BIS + ECMOL2
       XB4BIS = XB4BIS + EB(IB2)          
       IEV(IB2) = -4
*
       IDISS4 = IDISS4 + 1
       WRITE(6,999)TIME,IB2,NAMEB(IB2),IDISS4,EBOLD2,EBOLD2/BMRAT/EKT2
           XRAN = RANF()
*       If binary is dissolved distribute energy randomly between singles.
           ESB1 = XRAN*ECMB2
           ESB2 = (1.D0-XRAN)*ECMB2
           ESNGL1 = ESB1 + BODY1(IB2)*U1
           ESNGL2 = ESB2 + BODY2(IB2)*U2
           ERCOIL = 0.D0
       ELSE
           ESNGL1 = 0.D0
           ESNGL2 = 0.D0
*       Force relaxation of second binary next time
           DTBACC(IB2) = 1.D30
       END IF
*       End selection of 4b integration or approximate cross sections
      END IF
*       Force relaxation of remaining binary next time
       DTBACC(IB1) = 1.D30
*
       WRITE(6,998)TIME,IB1,IB2,ETBIN1,ETBIN2,ECMB1,ECMB2
 998   FORMAT(' T,IB1,2=',1P,D13.5,2I6,' C.M.: Etot bin1,2=',
     *    2D13.5,' Ekin bin1,2=',2D13.5)
*
*       Add mass to single stars if they do not escape
       IF(ESNGL1.LT.0.D0)THEN
       XMASS(J) = XMASS(J) + BODY1(IB2)
       XQMASS(J) = XQMASS(J) + BODY1(IB2)
*       Add heating for case of 4b integration and single star creation
*       (note that for Gao case ESB's are zero!)
       XKICK(J) = XKICK(J) + ESB1
       XHEAT(J) = XHEAT(J) + ESB1
*
       IQUAD = 1
       WRITE(6,997)IB2,TIME,ESNGL1,ESB1
 997   FORMAT(1X,' Single Star 1 from 4b dissolution of ',I6,
     *  ' T,Etot,Ekin=',1P,3D13.5)
      XN4SBI = XN4SBI + 1.D0
      XM4SBI = XM4SBI + BODY1(IB2)
      XE4SBI = XE4SBI + ESNGL1
*       Delayed removal of binary only if both singles remain in system.
*
       ELSE
*       If second binary yet exists there cannot be any single star escapers
        IF(EB(IB2).LE.0.D0)THEN
*       Take out kick heating in case of single star escaper
*       (this is only relevant for Gao case)
           XKICK(J)=XKICK(J) - XRAN*C34*ERCOIL
           XHEAT(J)=XHEAT(J) - XRAN*C34*ERCOIL      
       WRITE(6,996)IB2,TIME,ESNGL1,ESB1
 996   FORMAT(1X,' Single Star Escaper 1 from 4b dissolution of ',I6,
     *  ' T,Etot,Ekin=',1P,3D13.5)
*
         IEV(IB2)=IEV(IB2)+1
*
          XN4SES = XN4SES + 1.D0
          XM4SES = XM4SES + BODY1(IB2)
          XE4SES = XE4SES + ESNGL1
        END IF
          END IF
*
       IF(ESNGL2.LT.0.D0)THEN
       XMASS(J) = XMASS(J) + BODY2(IB2)
       XQMASS(J) = XQMASS(J) + BODY2(IB2)
*       Add heating for case of 4b integration and single star creation
*       (note that for Gao case ESB's are zero!) 
       XKICK(J) = XKICK(J) + ESB2
       XHEAT(J) = XHEAT(J) + ESB2        
       IQUAD = 1
       WRITE(6,995)IB2,TIME,ESNGL2,ESB2
 995   FORMAT(1X,' Single Star 2 from 4b dissolution of ',I6,
     *  ' T,Etot,Ekin=',1P,3D13.5)
      XN4SBI = XN4SBI + 1.D0
      XM4SBI = XM4SBI + BODY2(IB2)
      XE4SBI = XE4SBI + ESNGL2
*
       ELSE
*       If second binary yet exists there cannot be any single star escapers
        IF(EB(IB2).LE.0.D0)THEN 
*       Take out kick heating in case of single star escaper
*       (this is only relevant for Gao case)
          XKICK(J)=XKICK(J) - (1.D0-XRAN)*C34*ERCOIL
          XHEAT(J)=XHEAT(J) - (1.D0-XRAN)*C34*ERCOIL
       WRITE(6,994)IB2,TIME,ESNGL2,ESB2
 994   FORMAT(1X,' Single Star Escaper 2 from 4b dissolution of ',I6,
     *  ' T,Etot,Ekin=',1P,3D13.5)
*       Immediate removal of binary if one single star escapes.
          IEV(IB2)=IEV(IB2)+1
*
         XN4SES = XN4SES + 1.D0
         XM4SES = XM4SES + BODY2(IB2)
         XE4SES = XE4SES + ESNGL2
        END IF
       END IF
*       Delayed removal of binary only if both singles remain in system.
       IF(IEV(IB2).EQ.-4)THEN
       NESC2 = NESC2 + 1
       IESC2(NESC2) = IB2
*
       ELSE IF(IEV(IB2).GT.-4.AND.IEV(IB2).LT.-1) THEN
*       Immediate removal of binary if one single star escapes.
       NESC = NESC + 1
       IESC(NESC) = IB2
       NAMEB(IB2) = - ABS(NAMEB(IB2))
*
       END IF
*
*         Check first binary for dissolution or escape
*
        IF(EB(IB1).LE.0.D0)THEN
       WRITE(6,993)IB1,NAMEB(IB1),TIME,EB(IB1)/BMRAT/EKT1,EB(IB1)
 993   FORMAT(' Bin 1',I6,' N=',I6,' dissolved by quad T,EB/BMRAT/EKT,',
     *     'EB=',1P,3D13.5)
        ELSE
       WRITE(6,992)IB1,NAMEB(IB1),TIME,EB(IB1)/BMRAT/EKT1,EB(IB1)
 992   FORMAT(' Bin 1',I6,' N=',I6,' survived quad T,EB/BMRAT/EKT,EB=',
     *     1P,3D13.5)
*       Check quad kick escape for first binary (relevant also for Gao case)
*       only if it survived.
*
       IF(ETBIN1.GT.0.D0)THEN
       NESC = NESC + 1
       WRITE(6,991)IB1,NAMEB(IB1),TIME,EB(IB1),NESC
 991   FORMAT(' Bin 1',I6,' N=',I6,' quad kick escape T,EB,NESC=',
     *    1P,2D13.5,I6)
      XN4BES = XN4BES + 1.D0
      XM4BES = XM4BES + BODY1(IB1) + BODY2(IB1)
      XE4BES = XE4BES + ETBIN1
      XB4BES = XB4BES + EB(IB1)
       IESC(NESC) = IB1
       NAMEB(IB1) = - ABS(NAMEB(IB1))
       GOTO 357
       END IF
*
       END IF
*
*         Check second binary for dissolution or escape
*      If second binary yet exists there cannot be any single star escapers
        IF(EB(IB2).LE.0.D0)THEN
       WRITE(6,990)IB2,NAMEB(IB2),TIME,EB(IB2)/BMRAT/EKT2,EB(IB2)
 990   FORMAT(' Bin 2',I6,' N=',I6,' dissolved by quad T,EB/BMRAT/EKT,',
     *    'EB=',1P,3D13.5)
        ELSE
       WRITE(6,989)IB2,NAMEB(IB2),TIME,EB(IB2)/BMRAT/EKT2,EB(IB2)
 989   FORMAT(' Bin 2',I6,' N=',I6,' survived quad T,EB/BMRAT/EKT,EB=',
     *     1P,3D13.5)
*
*       Check quad kick escape for second binary (only in QUMAIN case relevent)
*       only if it survived.
*
       IF(ETBIN2.GT.0.D0)THEN
       NESC = NESC + 1
       WRITE(6,988)IB2,NAMEB(IB2),TIME,EB(IB2),NESC
 988   FORMAT(' Bin 2',I6,' N=',I6,' quad kick escape T,EB,NESC=',
     *    1P,2D13.5,I6)
      XN4BES = XN4BES + 1.D0
      XM4BES = XM4BES + BODY1(IB2) + BODY2(IB2)
      XE4BES = XE4BES + ETBIN2
      XB4BES = XB4BES + EB(IB2)
       IESC(NESC) = IB2
       NAMEB(IB2) = - ABS(NAMEB(IB2))
       GOTO 357
       END IF                
*
        END IF
*
*       Call newpos to determine new binary mass distribution via orbit
      DEFANG = 0.D0
      IPT1 = 0
      CALL NEWPOS(IB1,IPT1,DEFANG)
*       Newpos could detect binary close to escape at r=rnj - should escape
       IF(IPT1.GT.0)THEN
       NESC = NESC + 1
       IF(IPT1.EQ.1)PRINT*,' Binary ',IB1,' N=',NAMEB(IB1),' EB=',
     *  EB(IB1),'  escape after kick NESC=',NESC
       IF(IPT1.EQ.2)PRINT*,' Binary ',IB1,' N=',NAMEB(IB1),' EB=',
     *  EB(IB1),'  newpos escape after kick NESC=',NESC
      XNRBES = XNRBES + 1.D0
      XMRBES = XMRBES + BODY1(IB1) + BODY2(IB1)
      XERBES = XERBES + ETBIN1
      XBRBES = XBRBES + EB(IB1)
       IESC(NESC) = IB1
       NAMEB(IB1) = - ABS(NAMEB(IB1))
       IPT1 = 0
       END IF                                       
*       End escape checking for first binary
*
*       Start escape checking for second binary
*       Call newpos for second binary if this was not dissolved
      IF(EB(IB2).GT.0.D0)THEN
          DEFANG = 0.D0
          IPT2 = 0
          CALL NEWPOS(IB2,IPT2,DEFANG)   
*       Newpos could detect binary close to escape at r=rnj - should escape
          IF(IPT2.GT.0)THEN
            NESC = NESC + 1
       IF(IPT2.EQ.1)PRINT*,' Binary ',IB2,' N=',NAMEB(IB2),' EB=',
     *  EB(IB2),'  escape after kick NESC=',NESC
       IF(IPT2.EQ.2)PRINT*,' Binary ',IB2,' N=',NAMEB(IB2),' EB=',
     *  EB(IB2),'  newpos escape after kick NESC=',NESC
            XNRBES = XNRBES + 1.D0
            XMRBES = XMRBES + BODY1(IB2) + BODY2(IB2)
            XERBES = XERBES + ETBIN2
            XBRBES = XBRBES + EB(IB2)
            IESC(NESC) = IB2
            NAMEB(IB2) = - ABS(NAMEB(IB2))
            IPT2 = 0
          END IF                           
*        End treatment of second binary
      END IF
*
 357  CONTINUE
       PRINT*,'----------END QUAD---------------------------'
*       End close 4b encounter selected for this binary pair
      END IF
*
 355   CONTINUE
*       End close 4b encounters
       END IF
*
       IF(NESC.GT.0)THEN
       IBESC = 1
       PRINT*,' Removal of Binaries after QUAD, ',
     *  ' NESC=',NESC,' list:',(IESC(III),III=1,NESC)
       CALL ESCAPE(NESC,IESC)
       END IF           
*
*       Note on yet to be dissolved binaries (delayed removal) from both
*       TRIPLE and QUAD, because they produced two singles in the system
*
       IF(NESC2.GT.0)THEN
       PRINT*,' Delayed removal in close 3b/4b encounters NESC2=',
     *   ' names:',(NAMEB(IESC2(III)),III=1,NESC2)
       END IF            
*
*         Compute Binary Potential to save changes for energy balance
*         and for initialization in case of initial model
       IPOT = 0
       CALL PHIMAS(IPOT)
*
*          Accumulate binary potential changes by relaxation, kicks,
*          and binary escapers in the change of external
*          energy (potential) XEERR only for the energy balance
*
         DO 1100 J=1,NCOMP
*
         IF(NMOD.EQ.0)THEN
         DO 650 I=1,NJ
 650     PHIOLD(I) = PHIBIN(I)
         ELSE
*      DPBIN1 = 1*m*DPhibin for standard processes
         DPBIN1=0.D0
*
         DO 660 I=2,NJ
         IM = I-1
         DVOL = PI43*(R(I)**3-R(I-1)**3)
         DPHI = PHIBIN(I) - PHIOLD(I)
         DPBIN1 = DPBIN1 + DPHI*DEXP(X(I00+J,I))*DVOL
 660     CONTINUE
*
         DPTOT = DPTOT + DPBIN1
         XEERR = XEERR + DPBIN1
*
         END IF
*
*         Update old potential used to compute changes in next call
       DO 630 I=1,NJ
       PHIOLD(I) = PHIBIN(I)
 630   CONTINUE
*
 1100  CONTINUE
*
*         Now late removal dissolved binaries from the binary system
       DO 450 IB = 1,IBIN
       IF(IEV(IB).EQ.-4)THEN
       NAMEB(IB) = -ABS(NAMEB(IB))
*      PRINT*,' Binary IB,NAME=',IB,NAMEB(IB),' late removed '
       END IF
 450   CONTINUE
*
       IF(NESC2.GT.0)THEN
       IBESC = 1
       CALL ESCAPE(NESC2,IESC2)
       END IF
*         
*         Now look for binary formation
*
       RPREV = 0.D0
       XN3B=0.D0
*
       DO 400 J=1,NCOMP
*
        SIGR2=DEXP(X(I20+J,1)-X(I00+J,1))
        SIGT2=DEXP(X(I02+J,1)-X(I00+J,1))
        SIG2=C13*(SIGR2+2.D0*SIGT2)
      ECORE=0.D0
      XMCORE=0.D0
      E2CORE=0.D0
*
       DO 399 I=2,NJ
       DVOL=PI43*(R(I)**3-R(I-1)**3)
       IF(R(I).LE.RCORE)THEN
       ECORE=ECORE+C12*(DEXP(X(I20+J,I))+2.D0*DEXP(X(I02+J,I)))*DVOL
       XMCORE=XMCORE+DEXP(X(I00+J,I))*DVOL
       END IF
*
       IF(R(I).LE.R(I2C))THEN
       E2CORE=E2CORE+C12*(DEXP(X(I20+J,I))+2.D0*DEXP(X(I02+J,I)))*DVOL
       END IF
*
 399   CONTINUE
       XNCORE=XMCORE/XMIND(J)
       XN3OLD=0.D0
       XMCUM = 0.D0
       XN3B = 0.D0
       LPRINT=.TRUE.
*
       DO 401 I=2,NJ-1
*
*       Do not form binaries outside half-mass radius
       IF(DEXP(X(IMR,I) - X(IMR,NJ)).GT.0.5D0)GOTO 402
*
        SIGR2=DEXP(X(I20+J,I)-X(I00+J,I))
        SIGT2=DEXP(X(I02+J,I)-X(I00+J,I))
       SIG2=C13*(SIGR2+2.D0*SIGT2)
       XE1=GRAV*XMIND(J)/SIG2
       XE2=DEXP(X(I00+J,I))/XMIND(J)
*
      DVOL=PI43*(R(I)**3-R(I-1)**3)
      XMCUM = XMCUM + DEXP(X(IMR,I))-DEXP(X(IMR,I-1))
      XN3FR = 0.9D0*XE1**5*XE2**3*DSQRT(SIG2)*DVOL/BREM
      XN3B=XN3B+XN3FR
      RPREV=R(I)
*
       XRAN = RANF()
*
       LCHECK = I.GT.10.AND.XN3FR.LT.XN3OLD
       IF(((LCHECK.AND.LPRINT).OR.I.EQ.130).AND.MOD(IMOD,10).EQ.0)THEN
       LPRINT=.FALSE.
       END IF
*      IF(XFEED(J).NE.0.D0.AND.MOD(IMOD,10).EQ.0)THEN
*      PRINT*,' I,J,ESHELL=',I,J,C32*SIG2*DEXP(X(I00+J,I))*DVOL,
*    *   ' EFEED=',XFEED(J),' XN3FR,XN3B=',XN3FR,XN3B
*      END IF
       XN3OLD = XN3FR
*
*  IF BINARIES ARE PLANETARY SYSTEMS HERE DO NOT CREATE BINARIES
*
       IF(LPLAN)XN3B = 0.D0
*
*    Only check binary formation if more than 4m have accumulated
*    Take 2 times the binary, which has to remain in single star component
*    Take care in multi-mass system
       IF(4.D0*XMIND(J).LT.XMCUM)THEN
*    Following line for test purposes only
*       IF(IBINT.LT.100.AND.IMOD.GT.50.AND.MOD(IMOD,4).EQ.0)XN3B=1.D0
       IF(XRAN.LT.XN3B)THEN
*
       IBIN=IBIN+1
       IBINT=IBINT+1
       IEVENT=1
       ISHX = I
       IFORM = 1
*
 199   CONTINUE
*
*     LPR=NMOD.GT.960
*     IF(LPR)PRINT*,' Before Check binary formation NMOD=',NMOD
       RB(IBIN)=R(ISHX)*0.99D0
       ISH(IBIN)=ISHX
       ICO(IBIN)=J
       IEV(IBIN)=2
       BODY1(IBIN)=XMIND(J)
       BODY2(IBIN)=XMIND(J)
       XMBIN = BODY1(IBIN) + BODY2(IBIN)
       SIZE1(IBIN)=RIND(J)
       SIZE2(IBIN)=RIND(J)
       NAMEB(IBIN)=IBINT
*        First check binary mass and location
       RHONEW = (BODY1(IBIN)+BODY2(IBIN))/PI43/RB(IBIN)**3
       RHOLOC = DEXP(X(I00+J,ISHX))
*
       IF(RHONEW.GT.RHOLOC)THEN
       ISHX = ISHX + 1
       GOTO 199
       END IF
*
       PRINT*,' Binary Number ',IBIN,' N=',NAMEB(IBIN),
     *        ' created at R=',R(ISHX)*0.99D0,
     *        ' Component ',J,' Shell ',ISHX,' IMOD=',IMOD
       PRINT*,' Binary XN3B=',XN3B,' XRAN=',XRAN
       PRINT*,' Binary Amount of Heating =',3.D0*XMIND(J)*SIG2
       PRINT*,' Binary Energy EB=',-3.D0*XMIND(J)*SIG2
*       Choose random (thermalized) eccentricity.
       ECC2 = RANF()
*      ECC(IBIN) = SQRT(ECC2)
       ECC(IBIN) = 0.D0
       EB(IBIN)=C103*3.D0*XMIND(J)*SIG2
*       Apply 3kt to single stars (1.5 kt to c.m. binary see below)
*       Apply 10kt to single stars (10 kt to c.m. binary see below)
       XKICK(J)=XKICK(J)+C103*3.D0*XMIND(J)*SIG2
       XHEAT(J)=XHEAT(J)+C103*3.D0*XMIND(J)*SIG2
*
       SEMIA(IBIN) = C12*GRAV*BODY1(IBIN)*BODY2(IBIN)/EB(IBIN)
*        Movement of Binary in Equipartition
      VESC = DSQRT(-2.D0*(PHI(ISHX)-PHTID+PHIBIN(ISHX)))
*
*  Test for check of dynamical friction
*     IF(IBIN.EQ.1)THEN
*     VR(IBIN)=DSQRT(-1.0D0*(PHI(ISHX)-PHTID+PHIBIN(ISHX)))
*     VT(IBIN)=1.D-1
*     VR2=VR(IBIN)
*     VT2=1.D-1
*     ELSE
      VR(IBIN)=DSQRT(C12*SIGR2)
      VT(IBIN)=DSQRT(SIGT2)
*     END IF
*       Force Relaxation next time
      IF(RB(IBIN).LT.R(ITRANS))THEN
      TBIN(IBIN) = TRXLOC(I2C)/1.D1
      ELSE IF(RB(IBIN).LT.R(IHALF))THEN
      TBIN(IBIN) = TRXLOC(I2C)/1.D1
      ELSE
      TBIN(IBIN) = TRXLOC(I2C)/1.D1
      END IF
*
      DTBACC(IBIN) = 1.D30
      PRINT*,' IB=',IBIN,' torb=',torb(ibin),' trxloc=',trxloc(i),
     *  ' tbin=',tbin(ibin),' mg tcr=',
     *   rb(ibin)/sqrt(vr(ibin)**2+vt(ibin)**2)
      PRINT*,' IB=',IBIN,' i,rb,rcore,15*rcore=',
     *   i,rb(ibin),rcore,15.*rcore
      PRINT*,' creation of ',IBIN,' at RHO=',DEXP(X(I00+J,1))
*
*       Total amount of energy to distribute 4.5 kt (3kt from
*       single stars plus 1.5 kt binary c.m. in equipartition)
*      2*Esingle + Ebind = Ecm-bin-equip + Ekick-bin + Eheat-single
*      2*3/2kT     3kT     3/2kT           |----DELTAE ----|
*                                              4.5 kT
*      Eheat-single = 2/3 Ekick-bin (lin. momentum) = 3 kT (see above)
*
       DELTAE = 4.5D0*XMIND(J)*SIG2
       VB = DELTAE*XMIND(J)/(XMIND(J)+BODY1(IBIN)+BODY2(IBIN))
*
*       Calculate isotropic Kick on Binary for test purpose
*      VR(IBIN) = DSQRT(VR(IBIN)**2+C23*VB)
*      IF(RANF().LT.0.5D0)VR(IBIN) = -VR(IBIN)
*      VT(IBIN) = DSQRT(VT(IBIN)**2+2.D0*C23*VB) 
*
       CALL KICK(IBIN,VB)
*
       ECMNEW = XMBIN*(C12*(VR2**2+VT2**2)+PHI(ISHX)-PHTID+PHIBIN(ISHX))
       PRINT*,' Binary ',IBIN,' N=',NAMEB(IBIN),
     *        ' Start with VR,VT=',VR(IBIN),VT(IBIN)
       PRINT*,' Binary ',IBIN,' N=',NAMEB(IBIN),
     *        ' CM Energy=',ECMNEW
*
*       Call newpos to determine smooth binary mass distribution via orbit
      DEFANG = 0.D0
      IPT1 = 0
      CALL NEWPOS(IBIN,IPT1,DEFANG)
*
*        Update gas model quantities for newly formed binary
      XMASS(J) = XMASS(J) - BODY1(IBIN) - BODY2(IBIN)
      XFMASS(J) = XFMASS(J) - BODY1(IBIN) - BODY2(IBIN)
      XN3SBI = XN3SBI + 2.D0
      XM3SBI = XM3SBI + 2.D0*XMIND(J)
      XE3SBI = XE3SBI + ECMNEW
      XB3SBI = xB3SBI + EB(IBIN)
*
       END IF
*
*       Reset accumulated mass and probability after 4m are reached
       XMCUM = 0.D0
       XN3B = 0.D0
*
       END IF
*
 401   CONTINUE
*
 402   CONTINUE
*
 400   CONTINUE
*         Update Potential and Mass of Binaries once more.
       IPOT=0
       CALL PHIMAS(IPOT)
*         Do not reduce timestep for planetary runs
       IF(.NOT.LPLAN) THEN
*
*         Reduce timestep if event occurred with empty reservoir
       ERTOT=0.D0
       DO 200 K=1,NCOMP
 200   ERTOT=ERTOT+EREM(K)
       IF(IEVENT.EQ.1.AND.ERTOT.EQ.0.D0)THEN
       BREM=1.D1*BREM
       PRINT*,' IEVENT=',IEVENT,' DT reduction DTOLD=',
     *  1.D1/BREM,' DTNEW=',1.D0/BREM
       END IF
       IF(IEVENT.EQ.2)THEN
       BREM=1.D1*BREM
       PRINT*,' IEVENT=',IEVENT,' DT reduction DTOLD=',
     *  1.D1/BREM,' DTNEW=',1.D0/BREM
       CALL FLUSH(6)
       END IF
*
       END IF
*
*        Compute the energy change due to binary formation and dissolution.
         DO 1200 J=1,NCOMP
*
         DPBIN2=0.D0
*
         DO 860 I=2,NJ
         IM = I-1
         DVOL = PI43*(R(I)**3-R(I-1)**3)
         DPHI = PHIBIN(I) - PHIOLD(I)
         DPBIN2 = DPBIN2+C23*DPHI*DEXP(X(I00+J,I))*DVOL
 860     CONTINUE
         DPTOT = DPTOT + DPBIN2
         XEERR = XEERR + DPBIN2
*
*         Update old potential used to compute changes in next call
       DO 830 I=1,NJ
       PHIOLD(I) = PHIBIN(I)
 830   CONTINUE
*
 1200  CONTINUE
*         Sum standard energy generation for normalization
*
         GF5=GRAV**5*27.D0*DSQRT(3.D0)
        DYA = 0.D0
        RHOT = 0.D0
*
        DO 80 K=1,NCOMP
        DYA = DYA + DEXP(VX(I20+K,1))+2.D0*DEXP(VX(I02+K,1))
  80    RHOT = RHOT + DEXP(VX(I00+K,1))
        DYA = DYA/RHOT
*
        DO 1000 J=1,NCOMP
*
         EBCHCK(J)=0.D0
         XMCHCK(J)=0.D0
         EFCHCK(J)=0.D0
*
         DO 600 I=2,NJ
*
         DYB=0.D0
         RHOT=0.D0
*
         DO 617 K=1,NCOMP
         AFCC=DEXP(X(I02+K,I)-X(I20+K,I))
         CBTERM = DLOG(XMIND(K))
      DBTERM = DEXP(CBTERM+C52*X(I00+K,I)-C32*X(I20+K,I))/
     *            (1.D0+2.D0*AFCC)**C32
         DYB=DYB+DBTERM
         RHOT=RHOT+DEXP(X(I00+K,I))
 617     CONTINUE
*     Mass Distribution Cutoff
      IF(X(IMR,I)-X(IMR,NJ).GT.-1.609D0)THEN
      CFACM=2.D0/(1.D0+(R(I)/RIR20)**4)
      CFACE=2.D0/(1.D0+(R(I)/RIR20)**4)
      CFACF=2.D0/(1.D0+(R(I)/RIR20)**4)
      ELSE
      CFACM=1.D0
      CFACE=1.D0
      CFACF=1.D0
      END IF
*      Cutoff for Feedback Energy
*      IF(X(IMR,I)-X(IMR,NJ).GT.-0.693D0)THEN
*      CFACE=2.D0/(1.D0+(R(I)/RHALF)**2)
*      ELSE
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
         DER3=PI43*(R(I)**3-R(I-1)**3)
         DETERM=XBIN*DEXP(X(I00+J,I))*GF5*
     *      DYA*DYB**3/RHOT*DER3/BREM*CFACE
         DMTERM=XBIN*DEXP(X(I00+J,I))*GF5*
     *      DYA*DYB**3/RHOT**3*DER3/BREM*CFACM
         DFTERM=XBIN*DEXP(X(I00+J,I))*GF5*
     *      DYA*DYB**3/RHOT*DER3/BREM*CFACF
         EBCHCK(J)=EBCHCK(J)+DETERM
         XMCHCK(J)=XMCHCK(J)+DMTERM
         EFCHCK(J)=EFCHCK(J)+DFTERM
*
 600     CONTINUE
*          Mass and Energy balance for binary formation case
*          and triple or quad binary disruption with singles remaining
         IF(IFORM.EQ.1.OR.ISESC.EQ.1.OR.IQUAD.EQ.1.OR.ITRIP.EQ.1)THEN
         IF(LS(6))THEN
         PRINT*,' Program stopped because extrapolation is not allowed'
         STOP
         END IF
*
         XMASSY(J)=XBIN*XMASS(J)/XMCHCK(J)
         XKMASSY(J)=XBIN*XKMASS(J)/XMCHCK(J)
         XQMASSY(J)=XBIN*XQMASS(J)/XMCHCK(J)
         XFMASSY(J)=XBIN*XFMASS(J)/XMCHCK(J)
*
         DMTERM2=0.D0
         DMTERM3=0.D0
         DMTERM4=0.D0
         DMTERM=0.D0
         DETERM2=0.D0
         DETERM3=0.D0
         DETERM4=0.D0
*        XXX1 = 0.D0
*        XXX2 = 0.D0
*        XXX3 = 0.D0
*        XXXM = 0.D0
*
         DO 750 I=2,NJ
*
         IM=I-1
         DYB=0.D0
         RHOT=0.D0
*
         DO 717 K=1,NCOMP
         AFCC=DEXP(X(I02+K,I)-X(I20+K,I))
         CBTERM = DLOG(XMIND(K))
      DBTERM = DEXP(CBTERM+C52*X(I00+K,I)-C32*X(I20+K,I))/
     *            (1.D0+2.D0*AFCC)**C32
         DYB=DYB+DBTERM
         RHOT=RHOT+DEXP(X(I00+K,I))
*
 717     CONTINUE
      IF(X(IMR,I)-X(IMR,NJ).GT.-1.609D0)THEN
      CFACM=2.D0/(1.D0+(R(I)/RHALF)**4)
      ELSE
      RHALF=R(I)
      CFACM=1.D0
      END IF
*
*
*     LPR=NMOD.GT.960
*     IF(LPR)PRINT*,' Before Balances NMOD=',NMOD
*      Energy generation term used for normalization
         DER3=PI43*(R(I)**3-R(IM)**3)
         DMTERM=XBIN*DEXP(X(I00+J,I))*GF5*
     *      DYA*DYB**3/RHOT**3*DER3/BREM*CFACM
*      Now add up mass loss due to single star kick escapers
         DMTEMP3 = XKMASSY(J)*DMTERM/XBIN
*      Now add up mass gain due to disruption of binaries by quad
         DMTEMP4 = XQMASSY(J)*DMTERM/XBIN
*
         DMTEMP2 = XFMASSY(J)*DMTERM/XBIN
*        DMTEMP2=0.D0
*        DO 720 IB=IBIN,1,-1
*        IF(IEV(IB).EQ.2.AND.DMBIN(IB,I).GT.0.D0)THEN
*        DMTEMP2=DMTEMP2-DMBIN(IB,I)
*        PRINT*,' I,IB,DRHO=',I,IB,DMTEMP/DER3,' RHO=',DEXP(X(I00+J,I))
*        END IF
*720     CONTINUE
*
         DMTERM2=DMTERM2+DMTEMP2
         DMTERM3=DMTERM3+DMTEMP3
         DMTERM4=DMTERM4+DMTEMP4
         SIGR2=DEXP(X(I20+J,I)-X(I00+J,I))
         SIGT2=DEXP(X(I02+J,I)-X(I00+J,I))
         DETEMP2=C12*(SIGR2+2.D0*SIGT2)*DMTEMP2
         DETEMP3=C12*(SIGR2+2.D0*SIGT2)*DMTEMP3
         DETEMP4=C12*(SIGR2+2.D0*SIGT2)*DMTEMP4
         PHIAV=PHI(I)-PHTID+2.D0*PHIBIN(I)
         PHIAV2=PHI(I)-PHTID+PHIBIN(I)
*        XXX1 = XXX1 + C12*(SIGR2+2.D0*SIGT2)*DMTEMP2
*        XXX2 = XXX2 + C12*(PHI(I)-PHTID)*DMTEMP2
*        XXX3 = XXX3 + PHIBIN(I)*DMTEMP2
*        XXXM = XXXM + DMTEMP2
         DETERM2=DETERM2+DETEMP2+C12*PHIAV*DMTEMP2
         DETERM3=DETERM3+DETEMP3+PHIAV2*DMTEMP3
         DETERM4=DETERM4+DETEMP4+C12*PHIAV*DMTEMP4
*
 750     CONTINUE
*
         END IF
*
 1000    CONTINUE
*
*          Total energy correction for binary formation
         IF(IFORM.EQ.1)THEN
         XMERR = XMERR + DMTERM2
         XEERR = XEERR + DETERM2
         XMERR2 = XMERR2 + DMTERM2
         XEERR2 = XEERR2 + DETERM2
         END IF
*          Total energy correction for single escapers
         IF(ISESC.EQ.1)THEN
         XMERR = XMERR + DMTERM3
         XEERR = XEERR + DETERM3
         XMERR3 = XMERR3 + DMTERM3
         XEERR3 = XEERR3 + DETERM3
         XINDHS(1)=XINDHS(1)+DETERM3
         END IF
*          Total energy correction for disruption of binaries by quad
         IF(IQUAD.EQ.1)THEN
         XMERR = XMERR + DMTERM4
         XEERR = XEERR + DETERM4
         XMERR4 = XMERR4 + DMTERM4
         XEERR4 = XEERR4 + DETERM4
         XINDHS(1)=XINDHS(1)+DETERM4
         END IF
*
*       PRINT*,' DPBIN1,DPBIN2,DETERM2,Sum=',DPBIN1,DPBIN2,DETERM2,
*    *     DETERM2+DPBIN1+DPBIN2,' XMASSY=',XMASSY(1)
      XXMASS = XXMASS + DEMASS
      XXENER = XXENER + DEENER
*
*        IF(IBIN.GT.0)THEN
*        PRINT*,' BINSTO-Test XXMASS=',XXMASS,' XXENER=',XXENER
*        PRINT*,' Cumulative  XMERR=',XMERR,' XEERR=',XEERR
*        END IF
*
*     IF(IFORM.EQ.1.OR.ISESC.EQ.1.OR.IQUAD.EQ.1)THEN
*        PRINT*,' Mass/Energy loss by bin formation=',DMTERM2,DETERM2
*        PRINT*,' Mass/Energy loss by single escape=',DMTERM3,DETERM3
*        PRINT*,' Mass/Energy gain by quad disruption=',DMTERM4,DETERM4
*        END IF
*
*    Limit energy generation rate by central relaxation rate
*
         DO 1001 J=1,NCOMP
*
         I=2
         SIGR=DEXP((X(I20+J,I)-X(I00+J,I))/2.D0)
         SIGT=DEXP((X(I02+J,I)-X(I00+J,I))/2.D0)
         SIG=DSQRT(C13*(SIGR**2+2.D0*SIGT**2))
         TRXCEN=SIG**3/(DEXP(X(I00+J,I))*XMIND(J)*XCTOT*CS(60))
*
         DTIME=1.D0/BREM
*        APPLY = C13*ECORE/TRXCEN*DTIME
         XAPPLY = 5.D0
         IF(IBIN.LT.1500)XAPPLY = 2.D0
         IF(IBIN.LT.750)XAPPLY = 1.D0
         IF(IBIN.lt.600)XAPPLY = 0.75D0
         APPLY = XAPPLY*C13*E2CORE/TRXLOC(I2C)*DTIME
*    Do not cool but wait for heating, remember to apply rate in eqs.
*        IF(XHEAT(J).LT.0.D0)APPLY=0.D0
*
         IF(XHEAT(J).GT.APPLY)THEN
         EREM(J) = XHEAT(J) - APPLY
         XHEAT(J) = APPLY
         ELSE
*    Take care of relaxation which can produce cooling
             IF(XHEAT(J).LT.-APPLY)THEN
             EREM(J) = XHEAT(J) + APPLY
             XHEAT(J) = -APPLY
             ELSE
             EREM(J)=0.D0
             END IF
         END IF
*
*        APPLYF = C13*ECORE/TRXCEN*DTIME
         XAPPLY = 5.D0
         IF(IBIN.LT.1500)XAPPLY = 2.D0
         IF(IBIN.LT.750)XAPPLY = 1.D0
         IF(IBIN.lt.600)XAPPLY = 0.75D0          
         APPLYF = XAPPLY*C13*E2CORE/TRXLOC(I2C)*DTIME
*
         IF(XFEED(J).GT.APPLYF)THEN
         EFREM(J) = XFEED(J) - APPLYF
         XFEED(J) = APPLYF
         ELSE
*    Take care of relaxation which can produce cooling
             IF(XFEED(J).LT.-APPLYF)THEN
             EFREM(J) = XFEED(J) + APPLYF
             XFEED(J) = -APPLYF
             ELSE
             EFREM(J)=0.D0
             END IF
         END IF
*    Normalization of prim. bin. heating relative to stat. 3b heating
*   
         XSBINY(J)=XBIN*XHEAT(J)/EBCHCK(J)
         XFBINY(J)=XBIN*XFEED(J)/EFCHCK(J)
*
*        IF(IBIN.GT.0.AND.MOD(NMOD,NPR).EQ.0)THEN
*        PRINT*,' T,XHEAT,EREM,FACTOR=',TIME,XHEAT(J),EREM(J),
*    * XSBINY(J)/XBIN
*        WRITE(6,2222)TIME,XFEED(J),EFREM(J),XFBINY(J)/XBIN
*2222  FORMAT(1X,' T=',1P,D13.5,' XF=',D13.5,' EFR=',2D13.5)
*        PRINT*,' XSBINY,XMASSY=',XSBINY(J),XMASSY(J)
*        PRINT*,' NCORE,ECORE,ECORE/TRXCEN=',XNCORE,ECORE,ECORE/TRXCEN
*        CALL FLUSH(6)
*        END IF
*
 1001    CONTINUE
*
*   Update total internal and external binary energies
*   after info.f has been called first
      IF(IFRFIN.GT.0)THEN
*       external c.m. energy of binaries
      XBBEXT = 0.D0
*       internal binding energy of binaries
      XBBINT = 0.D0
*       Compute correction terms to get total system energy from
*       total energies of singles and binaries
      XCORR = 0.D0
      TEPOT = 0.D0
      TETERM = 0.D0
*
         DO 5555 J=1,NCOMP
         EPOT(J) = 0.D0
         ETERM(J) = 0.D0
*
         DO 555 I=2,NJ
         IM=I-1
         DER3=PI43*(R(I)**3-R(IM)**3)
       PHIAV=C12*(PHI(I)+2.D0*PHIBIN(I)+PHI(IM)+2.D0*PHIBIN(IM))-PHTID
         EPOT(J)=EPOT(J)+C12*DEXP(X(I00+J,I))*PHIAV*DER3
         ETERM(J)=ETERM(J)+
     *       (2.D0*DEXP(X(I02+J,I))+DEXP(X(I20+J,I)))/2.D0*DER3
 555  XCORR = XCORR - C12*DEXP(X(I00+J,I))*PHIBIN(I)*DER3
*
         TEPOT = TEPOT + EPOT(J)
         TETERM = TETERM + ETERM(J)
*
 5555    CONTINUE
*
         DO 500 IB = 1,IBIN
*
      XBBINT = XBBINT + EB(IB)
      I = ISH(IB)
      XMBIN = BODY1(IB) + BODY2(IB)
      XBBEXT = XBBEXT + C12*XMBIN*(VR(IB)**2+VT(IB)**2 +
     *         2.D0*(PHI(I)-PHTID) + PHIBIN(I))
      XCORR = XCORR - C12*XMBIN*(PHI(I)-PHTID)
*
 500     CONTINUE
*
*     WRITE(6,777)TIME,EPOTJ,ETERMJ,XBBEXT,XCORR,
*    *   EPOTJ+ETERMJ+XBBEXT+XCORR
*777  FORMAT(' T=',1P,D13.5,' TOT BAL ',5D13.5)
*
      IF(NMOD.EQ.0)THEN
      WRITE(85,*)' Numbers'
      WRITE(85,*)' t RBES 3BES 4BES RSES 3SES 4SES 4BIS 4SBI 3SBI 3BIS'
      WRITE(86,*)' Masses'
      WRITE(86,*)' t RBES 3BES 4BES RSES 3SES 4SES 4BIS 4SBI 3SBI 3BIS'
      WRITE(87,*)' External Energies'
      WRITE(87,*)' t RBES 3BES 4BES RSES 3SES 4SES 4BIS 4SBI 3SBI 3BIS'
      WRITE(88,*)' Internal Energies'
      WRITE(88,*)' t RBES 3BES 4BES 3SBI BBINT BBEXT 3BIS'
      END IF
*
       IF(IRIT.EQ.0)THEN
       WRITE(85,185)TIME,XNRBES,XN3BES,XN4BES,XNRSES,XN3SES,XN4SES,
     *              XN4BIS,XN4SBI,XN3SBI,XN3BIS
       WRITE(86,185)TIME,XMRBES,XM3BES,XM4BES,XMRSES,XM3SES,XM4SES,
     *              XM4BIS,XM4SBI,XM3SBI,XM3BIS
       WRITE(87,185)TIME,XERBES,XE3BES,XE4BES,XERSES,XE3SES,XE4SES,
     *              XE4BIS,XE4SBI,XE3SBI,XE3BIS
       WRITE(88,185)TIME,XBRBES,XB3BES,XB4BES,XB3SBI,XBBEXT,XBBINT,
     *              XCORR,XB3BIS
 185   FORMAT(1X,1P,11D13.5)
       END IF
*
      IF(IFRFIN+45.GT.IVDIM)THEN
      PRINT*,' STOP Routine binsto.f uses too large storage of VIT '
      STOP
      END IF                          
*
      IFR = IFRFIN
*
      VIT(IFR+1) = XNRBES
      VIT(IFR+2) = XN3BES
      VIT(IFR+3) = XN4BES
      VIT(IFR+4) = XNRSES
      VIT(IFR+5) = XN3SES
      VIT(IFR+6) = XN4SES
      VIT(IFR+7) = XN4BIS
      VIT(IFR+8) = XN4SBI
      VIT(IFR+9) = XN3SBI
      VIT(IFR+10) = XN3BIS
      ZS(IFR+1) = ' XNRBES '
      ZS(IFR+2) = ' XN3BES '
      ZS(IFR+3) = ' XN4BES '
      ZS(IFR+4) = ' XNRSES '
      ZS(IFR+5) = ' XN3SES '
      ZS(IFR+6) = ' XN4SES '
      ZS(IFR+7) = ' XN4BIS '
      ZS(IFR+8) = ' XN4SBI '
      ZS(IFR+9) = ' XN3SBI '     
      ZS(IFR+10) = ' XN3BIS '
*
      IFR = IFR + 10
*
      VIT(IFR+1) = XMRBES
      VIT(IFR+2) = XM3BES
      VIT(IFR+3) = XM4BES
      VIT(IFR+4) = XMRSES
      VIT(IFR+5) = XM3SES
      VIT(IFR+6) = XM4SES
      VIT(IFR+7) = XM4BIS
      VIT(IFR+8) = XM4SBI
      VIT(IFR+9) = XM3SBI
      VIT(IFR+10) = XM3BIS
      ZS(IFR+1) = ' XMRBES '
      ZS(IFR+2) = ' XM3BES '
      ZS(IFR+3) = ' XM4BES '
      ZS(IFR+4) = ' XMRSES '
      ZS(IFR+5) = ' XM3SES '
      ZS(IFR+6) = ' XM4SES '
      ZS(IFR+7) = ' XM4BIS '
      ZS(IFR+8) = ' XM4SBI '
      ZS(IFR+9) = ' XM3SBI '         
      ZS(IFR+10) = ' XM3BIS '
*
      IFR = IFR + 10
*
      VIT(IFR+1) = XERBES
      VIT(IFR+2) = XE3BES
      VIT(IFR+3) = XE4BES
      VIT(IFR+4) = XERSES
      VIT(IFR+5) = XE3SES
      VIT(IFR+6) = XE4SES
      VIT(IFR+7) = XE4BIS
      VIT(IFR+8) = XE4SBI
      VIT(IFR+9) = XE3SBI
      VIT(IFR+10) = XE3BIS
      ZS(IFR+1) = ' XERBES '
      ZS(IFR+2) = ' XE3BES '
      ZS(IFR+3) = ' XE4BES '
      ZS(IFR+4) = ' XERSES '
      ZS(IFR+5) = ' XE3SES '
      ZS(IFR+6) = ' XE4SES '
      ZS(IFR+7) = ' XE4BIS '
      ZS(IFR+8) = ' XE4SBI '
      ZS(IFR+9) = ' XE3SBI '         
      ZS(IFR+10) = ' XE3BIS '
*
      IFR = IFR + 10
                            
      VIT(IFR+1) = XBRBES
      VIT(IFR+2) = XB3BES
      VIT(IFR+3) = XB4BES
      VIT(IFR+4) = XB3SBI
      VIT(IFR+5) = XB3BIS
*
      VIT(IFR+6) = XBBEXT
      VIT(IFR+7) = XBBINT
      VIT(IFR+8) = XCORR
      VIT(IFR+9) = RHOBAV
      ZS(IFR+1) = ' XBRBES '
      ZS(IFR+2) = ' XB3BES '
      ZS(IFR+3) = ' XB4BES '
      ZS(IFR+4) = ' XB3SBI '
      ZS(IFR+5) = ' XB3BIS '
*
      ZS(IFR+6) = ' XBBEXT '
      ZS(IFR+7) = ' XBBINT '
      ZS(IFR+8) = ' XCORR  '
      ZS(IFR+9) = ' RHOBAV ' 
*
      IFR = IFR + 9
*
* Total internal energy of escaping binaries
      VIT(IFR+1) = XBRBES+XB3BES+XB4BES
      ZS(IFR+1) = ' EE-INT '
* Total energy balance for binary and single system (former bal.f)
      VIT(IFR+2) = TEPOT + TETERM + XBBEXT + XCORR
      ZS(IFR+2) = ' TOTBAL '
* Total external energy of escaping binaries and singles (former fort.87)
      VIT(IFR+3) = XERBES+XE3BES+XE4BES+XERSES+XE3SES+XE4SES
      ZS(IFR+3) = ' EE-EXT '
* Total internal energy of binaries remaining in the system
      VIT(IFR+4) = XBBINT
      ZS(IFR+4) = ' EB-INT '
* Total number of escaping binaries
      VIT(IFR+5) = XNRBES+XN3BES+XN4BES
      ZS(IFR+5) = ' NE-EXTB'
* Total number of escaping singles
      VIT(IFR+6) = XNRSES+XN3SES+XN4SES
      ZS(IFR+6) = ' NE-EXTS' 
* Total mass of escaping binaries
      VIT(IFR+7) = XMRBES+XM3BES+XM4BES
      ZS(IFR+7) = ' ME-EXTB'
* Total mass of escaping singles
      VIT(IFR+8) = XMRSES+XM3SES+XM4SES
      ZS(IFR+8) = ' ME-EXTS'
* Total external energy of escaping binaries
      VIT(IFR+9) = XERBES+XE3BES+XE4BES
      ZS(IFR+9) = ' EE-EXTB'
* Total external energy of escaping singles
      VIT(IFR+10) = XERSES+XE3SES+XE4SES
      ZS(IFR+10) = ' EE-EXTS'
*
      IFR = IFR + 10
*
      IF(IZ00.EQ.1)THEN
          PRINT*,' Binsto: Last used IFR=',IFR,' IVDIM=',IVDIM
          PRINT*,' Binsto used APPLY-Fac=',XAPPLY
          PRINT*,' Use 1,2 for IBIN<500,1500'
          IZ00=2
          IF(IFR.GT.IVDIM)STOP
      END IF

* Refresh memorize content of VIT-vector in Unit 90
          IF(IY00.EQ.0)THEN
          IY00=1
          DO 290 IJK=0,IVDIM,6
*
          IF(IJK+6.LE.IVDIM)THEN
          IJLMAX=6
          ELSE
          IJLMAX=MOD(IVDIM,6)
          END IF
*
          WRITE(90,190)(IJK+IJL,ZS(IJK+IJL),IJL=1,IJLMAX)
 190      FORMAT(1X,6('|',I4,A9))
 290      CONTINUE
          END IF
*
          CLOSE(90)
*                                             
      END IF
*
*       Save Index of binary with minimum rmin for use in next model
      IF(IBIN.GT.0)THEN
      RMMIN = 1.D30
      DO 700 IB = 1,IBIN
      IF(RBMIN(IB).LT.RMMIN)THEN
      RMMIN = RBMIN(IB)
      IMCRIT = NAMEB(IB)
      END IF
 700  CONTINUE
      END IF
* 
*       Apply Mass Changes with relaxation rate of transition radius
*       where relaxation time is larger than dynamical time (or core radius)
      DO 800 I = 1,NJ
      VMRBIN(I) = XMRBIN(I)
*     IRELAX = I2C
*     IF(I2C.LT.ITRANS)IRELAX = ITRANS
*     VMRBIN(I) = VMRBIN(I)+(XMRBIN(I)-VMRBIN(I))/TRXLOC(IRELAX)/BREM
*     IF(DABS(VMRBIN(I)).LT.XMMIN)VMRBIN(I)=0.D0
*     IF(DABS(VMRBIN(I)-XMRBIN(I)).LT.XMMIN)VMRBIN(I)=XMRBIN(I)
 800  CONTINUE
*       Update absolute clock in case time step was changed
      TIME = TOLD + 1.D0/BREM
*
       RETURN
*
       END
