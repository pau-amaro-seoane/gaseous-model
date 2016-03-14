      SUBROUTINE INFO
C------------------------------------------------------------
C ROUTINE TO SAVE TIME-DEPENDENT INFORMATION OF EVERY
C COMPUTED MODEL --------------------------------------------
C THE INFORMATION IS STORED IN THE ARRAY VIT AND
C WRITTEN UNFORMATTEDLY ON FILE UNIT 3
C------------------------------------------------------------
C INDEX OF VIT DENOTES THE STORED QUANTITY
C CORRESPONDING TO THE FOLLOWING LIST:
C 1: TIME
C TO BE DEFINED
C
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
         INCLUDE 'compar.f'
          DIMENSION XLAGR(10),RLAGR(10,NCOMPO+2),ILAGR(11,NCOMPO+2)
          DIMENSION AVMASS(10)
          DIMENSION EPOT(NCOMPO),EKIN(NCOMPO),ETERM(NCOMPO)
          DIMENSION ETOT(NCOMPO)
      CHARACTER*9 ZS
      CHARACTER*5 ZVS(19)
      COMMON/INFOHELP/ZS(IVDIM)
      DATA XLAGR/0.01D0,0.02D0,0.05D0,0.1D0,0.2D0,0.3D0,0.4D0,
     *           0.5D0,0.75D0,0.9D0/
      DATA NLAGR/10/
C
      IF(IZ00.EQ.0)THEN
      IZ00=1
      OPEN(3,FILE='fort.3',FORM='UNFORMATTED',
     *    STATUS='UNKNOWN',IOSTAT=IERR)
      IF(IERR.NE.0)PRINT*,' ERROR IN OPEN 3 IERR=',IERR
      IF(LS(10))THEN
      DO 111 J=1,IREC
      READ(3,ERR=999,IOSTAT=IREERR)(VIT(K),K=1,IVDIM)
 111  CONTINUE
      END IF
*
      PRINT*,IREC,' Records read from Info-File'
      CALL FLUSH(6)
*
      DO 112 K=1,IVDIM
 112  ZS(K)='         '
C
      END IF
C
 888      CONTINUE
C
          VIT(1)=TIME
          ZS(1)='    TIME '
C      With black hole: "central" quantities just outside RGRAV
          ICEN=2
          IF(LS(19))THEN
 102      ICEN=ICEN+1
          IF(R(ICEN).LT.RGRAV.AND.ICEN.LT.NJ)GOTO 102
          END IF
*      Take radial index just inside RGRAV
          ICEN=ICEN-1
*
          RHO=0.D0
          DO 100 I=1,NCOMP
*      Interpolate core values at cusp boundary if central object present
          IF(LS(19))THEN
          SIGRIN=DEXP((X(I20+I,ICEN)-X(I00+I,ICEN))/2.D0)
          SIGTIN=DEXP((X(I02+I,ICEN)-X(I00+I,ICEN))/2.D0)
          SIGRUP=DEXP((X(I20+I,ICEN+1)-X(I00+I,ICEN+1))/2.D0)
          SIGTUP=DEXP((X(I02+I,ICEN+1)-X(I00+I,ICEN+1))/2.D0)
          DRHO=(DEXP(X(I00+I,ICEN+1))-DEXP(X(I00+I,ICEN)))/
     *            (R(ICEN+1)-R(ICEN))*(RGRAV-R(ICEN))
          DSIGR=(SIGRUP-SIGRIN)/(R(ICEN+1)-R(ICEN))*(RGRAV-R(ICEN))
          DSIGT=(SIGTUP-SIGTIN)/(R(ICEN+1)-R(ICEN))*(RGRAV-R(ICEN))
          ELSE
          DRHO=0.D0
          DSIGR=0.D0
          DSIGT=0.D0
          END IF
*
          RHO=RHO+DEXP(X(I00+I,ICEN))+DRHO
          VIT(I+1)=DEXP(X(I00+I,ICEN))+DRHO
 100      ZS(I+1)='    RHOC '
          RHO=RHO/DBLE(NCOMP)
C
          SIG=0.D0
          DO 110 I=1,NCOMP
          SIG=SIG+DEXP((X(I20+I,ICEN)-X(I00+I,ICEN))/2.D0)+DSIGR
          VIT(NCOMP+1+I)=DEXP((X(I02+I,ICEN)-X(I00+I,ICEN))/2.D0)+DSIGR
          ZS(NCOMP+1+I)='   SIGRC '
        VIT(2*NCOMP+1+I)=DEXP((X(I20+I,ICEN)-X(I00+I,ICEN))/2.D0)+DSIGT
 110      ZS(2*NCOMP+1+I)='   SIGTC '
          SIG=SIG/DBLE(NCOMP)
C
C INFORMATION ABOUT GLOBAL QUANTITIES CALCULATED HERE
C E.G.    THE VARIOUS FORMS OF ENERGIES
C
C-------------------------------------------------------------------------
C    Nomenclature:
C     EPOT(J),EKIN(J),ETERM(J),ETOT(J)   actual energy values for comp. J
C     TEPOT, TEKIN, TETERM, TETOT       total values
C     T0POT, T0KIN, etc. same as TEPOT, TEKIN etc. but values
C      at the beginning of the run (stored for file in vector AEI)
C-------------------------------------------------------------------------
C
         UNEN=UNM*UNR*UNR/UNT/UNT
C
         IFR=3*NCOMP+2
C
         TEKIN=0.D0
         TETERM=0.D0
         TEPOT=0.D0
         TETOT=0.D0
C
         DO 20 J=1,NCOMP
         EPOT(J)=0.D0
         EKIN(J)=0.D0
         ETERM(J)=0.D0
C
         DO 101 I=2,NJ
         IM=I-1
C
         DER3=PI43*(R(I)**3-R(IM)**3)
C
       PHIAV=C12*(PHI(I)+2.D0*PHIBIN(I)+PHI(IM)+2.D0*PHIBIN(IM))
         EPOT(J)=EPOT(J)+C12*DEXP(X(I00+J,I))*PHIAV*DER3
         EKIN(J)=EKIN(J)+(X(I10+J,I)**2/2.D0)*DEXP(X(I00+J,I))*DER3
         ETERM(J)=ETERM(J)+
     *       (2.D0*DEXP(X(I02+J,I))+DEXP(X(I20+J,I)))/2.D0*DER3
C
 101    CONTINUE
C
         ETOT(J)=ETERM(J)+EKIN(J)+EPOT(J)
C---------Sum components-----------------------
       TEKIN=EKIN(J)+TEKIN
       TETERM=ETERM(J)+TETERM
       TEPOT=EPOT(J)+TEPOT
       TETOT=ETOT(J)+TETOT
C-------------------------------------
C
         VIT(IFR+(J-1)*3)=EPOT(J)
         VIT(IFR+(J-1)*3+1)=EKIN(J)
         VIT(IFR+(J-1)*3+2)=ETERM(J)
         ZS(IFR+(J-1)*3)='    EPOT '
         ZS(IFR+(J-1)*3+1)='    EKIN '
         ZS(IFR+(J-1)*3+2)='   ETERM '
C
 20      CONTINUE
C
         IFR=IFR+3*NCOMP
         T0POT=AEI(5)
         T0KIN=AEI(6)
         T0TERM=AEI(7)
         VIT(IFR)=TEPOT-T0POT
         VIT(IFR+1)=TEKIN-T0KIN
         VIT(IFR+2)=TETERM-T0TERM
         ZS(IFR)='   DEPOT '
         ZS(IFR+1)='   DEKIN '
         ZS(IFR+2)='   DEKIN '
C-----------------Relaxation/Core Collapse Rates-------------------
         IFR=IFR+3
         I=3
         CSUM=0.D0
         DO 567 J=1,NCOMP
         SIG2=DEXP(X(I20+J,I)-X(I00+J,I))*
     *               (1.D0+2.D0*DEXP(X(I02+J,I)-X(I20+J,I)))
         TRX(J,J)=(SIG2/3.D0)**C32/
     *    CS(60)/XMIND(J)/XCOUL(J)/DEXP(X(I00+J,I))
         CCR=BREM*(X(I00+J,I)-VX(I00+J,I))*TRX(J,J)
         CSUM=CSUM+CCR
      VIT(IFR+J-1)=CCR
      ZS(IFR+J-1)=' CC-RATE '
 567     CONTINUE
*
*           Return without writing if important quantities
*           RHO,SIG,CC-Rate have not changed more than CORR,
*           but not for NMOD<10 or if 5*NRIT times it was not written
         LRET=DABS(RHO-VRHO)/RHO.LT.CORR .AND.
     *        DABS(SIG-VSIG)/SIG.LT.CORR .AND.
     *        DABS((CSUM-VCSUM)/(CSUM+EPS)).LT.CORR .AND.
     *        NMOD.GT.20. AND. IMOD-IMOLD.LT.5*NRIT
         IF(LRET)RETURN
         IMOLD=IMOD
*
C------------------Equipartition Values-----------------------------
          IFR=IFR+NCOMP
*
          VIT(IFR)=TTRX
          ZS(IFR)='    TTRX '
*
          IF(NCOMP.GT.1)THEN
          I=3 
          DO 566 J=2,NCOMP
          VIT(IFR+J-1)=XMIND(J)*DEXP(X(I20+J,I)-X(I00+J,I))/
     *     XMIND(1)/DEXP(X(I20+1,I)-X(I00+1,I)) - 1.D0
          ZS(IFR+J-1)='   EQUIP '
 566      CONTINUE 
*
          IFR = IFR + NCOMP - 1
*
          END IF
*
          VIT(IFR+1)=XMHOLE
          VIT(IFR+2)=DMTOT
          ZS(IFR+1)='   MHOLE '
          ZS(IFR+2)='   DMTOT '
*
          IFR=IFR+2
*
          DO 569 J=1,NCOMP
          DO 570 K=1,5
          VIT(IFR+(J-1)*5+K)=DMHOLE(K,J)
 570      ZS(IFR+(J-1)*5+K)='  DMHOLE '
          IFR=IFR+5
 569      CONTINUE
*
C      HOMOLOGY INVARIANTS U,V OF THE ISOTHERMAL MODEL
      I=3
      IP=I+1
      DR=DLOG(R(IP)/R(I))
      DO 55 J=1,NCOMP
      U3=(X(IMR,IP)-X(IMR,I))/DR
 49   V3=-(X(I20+J,IP)-X(I20+J,I))/DR
      W3=-(X(I02+J,IP)-X(I02+J,I))/DR
      VIT(IFR+1+3*J-3)=U3
      VIT(IFR+1+3*J-2)=V3
      VIT(IFR+1+3*J-1)=W3
      ZS(IFR+1+3*J-3)='       U '
      ZS(IFR+1+3*J-2)='       V '
      ZS(IFR+1+3*J-1)='       W '
 55   CONTINUE
*
      IFR = IFR + 3*NCOMP
      NC1 = NCOMP + 1
*        Determine Lagrangian radii for whole system
*        Here possible binaries are included
          DO 700 J = 1,NLAGR
          IP = 1
 705      I = IP
          IP = IP+1
*        If mass has been lost set rlagr to large value
          IF(IP.GT.NJ)THEN
          RLAGR(J,NC1) = 1.D30
          ELSE
          XMII = DEXP(X(IMR,I))+XMRBIN(I)
          XMIP = DEXP(X(IMR,IP))+XMRBIN(IP)
          XMICR = XLAGR(J)*(DEXP(X(IMR,NJ))+XMRBIN(NJ))
          IF(XMIP.LE.XMICR) GOTO 705
          RIP3 = R(IP)**3
          RI3 = R(I)**3
          ALPH = (RIP3 - RI3)/(XMIP-XMII)
          RLAGR(J,NC1) = RI3 + ALPH*(XMICR-XMII)
          RLAGR(J,NC1) = RLAGR(J,NC1)**C13
          END IF
*
        VIT(IFR+J) = RLAGR(J,NC1)
        WRITE(ZS(IFR+J),666)INT(100.*XLAGR(J))
 666    FORMAT(' TOR(',I3,')')
*
 700      CONTINUE
*
          IFR = IFR + NLAGR
*
          DO 800 J = 1,NLAGR
*
          IP = 1
 805      IP = IP+1
          IF(IP.GT.NJ)THEN
          ILAGR(J,NC1) = NJ
          ELSE
*        Determine mesh number ILAGR of Lagrangian Radii
          XMI = DEXP(X(IMR,IP))+XMRBIN(IP)
*        Decide whether to take initial or actual mass for Lagr. radii
*         XMICR = XLAGR(J)*(XMTOT0+XMBTOT)
          XMICR = XLAGR(J)*(DEXP(X(IMR,NJ))+XMRBIN(NJ))
          IF(XMI.LE.XMICR)GOTO 805
          ILAGR(J,NC1) = IP
          END IF
*
 800    CONTINUE
*
*        In case of primordial binaries determine their Lagr radii
      IF(LS(5))THEN
      NC2 = NCOMP + 2
*
          DO 770 J = 1,NLAGR
          IP = 1
 775      I = IP
          IP = IP+1
*        If mass has been lost set rlagr to large value
          IF(IP.GT.NJ)THEN
          RLAGR(J,NC2) = 1.D30
          ELSE
*
          XMII = XMRBIN(I)
          XMIP = XMRBIN(IP)
          XMICR = XLAGR(J)*XMRBIN(NJ)
          IF(XMIP.LE.XMICR) GOTO 775
          RIP3 = R(IP)**3
          RI3 = R(I)**3
          ALPH = (RIP3 - RI3)/(XMIP-XMII)
          RLAGR(J,NC2) = RI3 + ALPH*(XMICR-XMII)
          RLAGR(J,NC2) = RLAGR(J,NC2)**C13
          END IF
*
        VIT(IFR+J) = RLAGR(J,NC2)
        WRITE(ZS(IFR+J),677)INT(100.*XLAGR(J))
 677    FORMAT(' TBR(',I3,')')
*
 770      CONTINUE
        IFR = IFR + NLAGR
*
          DO 880 J = 1,NLAGR
*
          IP = 1
 885      IP = IP+1
          IF(IP.GT.NJ)THEN
          ILAGR(J,NC2) = NJ
          ELSE
*        Determine mesh number ILAGR of Lagrangian Radii
          XMIP = XMRBIN(IP)
          XMICR = XLAGR(J)*XMRBIN(NJ)
          IF(XMIP.LE.XMICR)GOTO 885
          ILAGR(J,NC2) = IP
          END IF
*
 880      CONTINUE
*
          ELSE
*        Advance locations in VIT vector to avoid change of position of items
          IFR = IFR + NLAGR
*
          END IF
*
          RFAC=R(3)/R(2)
*        Determine Lagrangian radii for components
          DO 680 K = 1,NCOMP
*
          DO 710 J = 1,NLAGR
          IP = 1
          XMI = 0.D0
*
 715      IP = IP+1
          I = IP-1
*        If mass has been lost set rlagr to large value
          IF(IP.GT.NJ)THEN
          RLAGR(J,K) = 1.D30
          ELSE
          RIP3 = R(IP)**3
          RI3 = R(I)**3
          XMIM = XMI
          XMI = XMI + PI43*DEXP(X(I00+K,IP))*(RIP3-RI3)
*        Decide whether to take initial or actual mass for Lagr. radii
*         IF(XMI.LE.XLAGR(J)*XMTOTI(K)) GOTO 815
          IF(XMI.LE.XLAGR(J)*XMTOT(K)) GOTO 715
*
          ALPH = (RIP3 - RI3)/(XMI-XMIM)
          RLAGR(J,K) = RI3 + ALPH*(XLAGR(J)*XMTOT(K)-XMIM)
          RLAGR(J,K) = RLAGR(J,K)**C13
          END IF
*
        VIT(IFR+J) = RLAGR(J,K)
        WRITE(ZS(IFR+J),676)K,INT(100.*XLAGR(J))
 676    FORMAT(' ',I2,'R(',I3,')')
*
 710      CONTINUE
*
          IFR = IFR + NLAGR
*
          DO 810 J = 1,NLAGR
*
          IP = 1
          XMI = 0.D0
*
 815      IP = IP+1
*        Determine mesh number ILAGR of Lagrangian Radii
          I = IP-1
*
          IF(IP.GT.NJ)THEN
          ILAGR(J,K) = NJ
          ELSE
          RIP3 = R(IP)**3
          RI3 = R(I)**3
          XMIM = XMI
          XMI = XMI + PI43*DEXP(X(I00+K,IP))*(RIP3-RI3)
*
          IF(XMI.LE.XLAGR(J)*XMTOT(K)) GOTO 815
          ILAGR(J,K) = IP 
          END IF
*
 810      CONTINUE
*
 680      CONTINUE
*
*        It follows storage of mass-weighted average
*        quantities between R(ILAGR(J))-radii
*
          TINY=1.D-30
*
          KKMAX = NCOMP + 1
          IF(NCOMP.EQ.1)KKMAX = 1
*
          DO 910 KK = 1,KKMAX
*
          DO 900 J = 1,NLAGR
*
          IF(J.EQ.1)THEN
          IGR = 2
          RDOWN=0.D0
          ELSE
          IGR = ILAGR(J-1,KK) - 1
          RDOWN=RLAGR(J-1,KK)
          END IF
*
          RUP=RLAGR(J,KK)
*
          IF(IFR+19.GT.IVDIM)THEN
          PRINT*,' INFO: VECTOR VIT HAS TOO SMALL A DIMENSION '
          PRINT*,'       CALCULATION STOPPED '
          PRINT*,' J,KK=',J,KK
          STOP
          END IF
*
          DO 850 KKK = 1,19
 850      VIT(IFR+KKK) = 0.D0
*
          SMASS = 0.D0
          XNMASS = 0.D0
*
  950     CONTINUE
*
          IGR = IGR + 1
*
          IF(R(IGR-1).LT.RDOWN)THEN
          RLOW=RDOWN
          ELSE
          RLOW=R(IGR-1)
          END IF
*
          IF(R(IGR).GE.RUP)THEN
          RHIGH=RUP
          ELSE
          RHIGH=R(IGR)
          END IF
*
          K = 0
 990      CONTINUE
* For K=NCOMP+1 take averages over all components (for whole system)
          IF(KK.EQ.NC1)THEN
          K = K + 1
          ELSE
          K = KK
          END IF
*
          HX(I00+K)=DEXP(X(I00+K,IGR))
          HX(I20+K)=DEXP(X(I20+K,IGR))
          HX(I02+K)=DEXP(X(I02+K,IGR))
          HXM(I20+K)=DEXP(VX(I20+K,IGR))
          HXM(I02+K)=DEXP(VX(I02+K,IGR)) 
          HX(I30+K)=3.D0*HX(I20+K)*X(I30+K,IGR)
          HX(I12+K)=2.D0*HX(I02+K)*X(I12+K,IGR)
          HX(I10+K)=X(I10+K,IGR)
          HXM(I30+K)=3.D0*DEXP(X(I20+K,IGR-1))*X(I30+K,IGR-1)
          HXM(I12+K)=2.D0*DEXP(X(I02+K,IGR-1))*X(I12+K,IGR-1)
          HXM(I10+K)=X(I10+K,IGR-1)
*
          DEM = PI43 * HX(I00+K)*(RHIGH**3 - RLOW**3)
          DEN = DEM/XMIND(K)
*
          VIT(IFR+1) = VIT(IFR+1) + DEM*HX(I20+K)/HX(I00+K)
          VIT(IFR+2) = VIT(IFR+2) + DEM*HX(I02+K)/HX(I00+K)
          ANI = 2.D0 - 2.D0*HX(I02+K)/HX(I20+K)
          VANI= 2.D0 - 2.D0*HXM(I02+K)/HXM(I20+K)
          VIT(IFR+3) = VIT(IFR+3) + DEM*ANI
          IF(DABS(ANI).LT.1.D-30)ANI=1.D-30
*
          PANI = 2.D0 * ANI * HX(I20+K)
          VPANI= 2.D0 *VANI * HXM(I20+K)
          AAV = DSQRT(DABS(ANI*VANI))
          PAV = DSQRT(DABS(PANI*VPANI))
          TANI=PAV/BREM/(PANI-VPANI+TINY)
          TANI2=AAV/BREM/(ANI-VANI+TINY)
          SUM = C13*(HX(I20+K)+HX(I02+K))/HX(I00+K)
          TRXX = SUM**C32/HX(I00+K)/CS(60)/XMIND(K)/XCOUL(K)
*
          VIT(IFR+4) = VIT(IFR+4) + DEM*TANI
          VIT(IFR+5) = VIT(IFR+5) + DEM*TANI2
          VIT(IFR+6) = VIT(IFR+6) + DEM*TRXX/(TANI+TINY)
          VIT(IFR+7) = VIT(IFR+7) + DEM*TRXX/(TANI2+TINY)
*
          DER = RHIGH - RLOW
          DER3 = RHIGH**3 - RLOW**3
          RAV = DSQRT(RHIGH*RLOW)
*
          QU = 2.D0*RAV*HX(I20+K)*(HX(I10+K)/RHIGH -
     *               HXM(I10+K)/R(IGR-1))/DER
*
          QFR = 3.D0/DER3*(RHIGH**2*HX(I30+K) -
     *                  RLOW**2*HXM(I30+K))
          QFT = -C32/DER3*(RHIGH**2*HX(I12+K) - 
     *                  RLOW**2*HXM(I12+K)) + 
     *       2.D0*(HX(I12+K)+HXM(I12+K))/RAV
          QF = QFR + QFT
*
          VIT(IFR+8) = VIT(IFR+8) + DEM*QU/(PANI+TINY)
          VIT(IFR+9) = VIT(IFR+9) + DEM*QFR/(PANI+TINY)
          VIT(IFR+10)= VIT(IFR+10) + DEM*QFT/(PANI+TINY)
          VIT(IFR+11)= VIT(IFR+11) + DEM*QF/(PANI+TINY)
          VIT(IFR+12)= VIT(IFR+12) + DEM*9.D0/10.D0/TRXX
          VIT(IFR+13)= VIT(IFR+13) + DEM/(TANI+TINY)
          VIT(IFR+14)= VIT(IFR+14) + DEM*QU
          VIT(IFR+15)= VIT(IFR+15) + DEM*QFR
          VIT(IFR+16)= VIT(IFR+16) + DEM*QFT
          VIT(IFR+17)= VIT(IFR+17) + DEM*HX(I30+K)
          VIT(IFR+18)= VIT(IFR+18) + DEM*HX(I12+K)
*        
          SMASS = SMASS + DEM
          XNMASS = XNMASS + DEN
* Next component
          IF(KK.EQ.NC1.AND.K.LT.NCOMP)GOTO 990
* Next radial mesh
          IF(IGR.LT.ILAGR(J,KK))GOTO 950
* Anisotropy at Lagrangian radii (no averaging)
*
          AI = 2.D0 - 2.D0*HX(I02+K)/HX(I20+K)
          AIM = 2.D0 - 2.D0*HXM(I02+K)/HXM(I20+K)
          ALPH = (AI - AIM)/(R(IGR)**3-RUP**3)
*
          VIT(IFR+19)= AIM + ALPH*(RUP**3-R(IGR-1)**3)
*
          DO 960 KKK = 1,18
 960      VIT(IFR+KKK) = VIT(IFR+KKK)/SMASS
*
          ZVS(1)='SGR2 '
          ZVS(2)='SGT2 '
          ZVS(3)=' ANI '
          ZVS(4)=' TAN '
          ZVS(5)=' TAN2'
          ZVS(6)='TX/TA'
          ZVS(7)='TX/T2'
          ZVS(8)=' QUPA'
          ZVS(9)='QFRPA'
          ZVS(10)='QFTPA'
          ZVS(11)='QF/PA'
          ZVS(12)='L/TAN'
          ZVS(13)='1/TAN'
          ZVS(14)='  QU '
          ZVS(15)=' QFR '
          ZVS(16)=' QFT '
          ZVS(17)='  VR '
          ZVS(18)='  VT '
          ZVS(19)=' ANIL'
*
          DO 665 ILK=1,19
          WRITE(ZS(IFR+ILK),667)KK,J,ZVS(ILK)
 667      FORMAT(2I2,A5)
 665      CONTINUE
*
	  AVMASS(J) = SMASS/XNMASS
          IFR = IFR+19
 900      CONTINUE
*
 910      CONTINUE
*
*           Store cumulative escape masses and energies
          DO 970 J=1,NCOMP
          VIT(IFR+1)=TMESC(J)
          VIT(IFR+2)=TEESC(J)
          VIT(IFR+3)=DMESC(J)
          VIT(IFR+4)=DEESC(J)
          ZS(IFR+1)=' TOTMESC '
          ZS(IFR+2)=' TOTEESC '
          ZS(IFR+3)='   DMESC '
          ZS(IFR+4)='   DEESC '
          IFR=IFR+4
 970      CONTINUE
*           Store cumulative binary/black hole energies
          ISTA = 100 + NCOMP
          DO 980 J=1,NCOMP
*           Only kick energies
          VIT(IFR+1)=AEI(ISTA+J)
          ZS(IFR+1)=' EKICK   '
*           Binary radial and tangential energies
          VIT(IFR+2)=AEI(ISTA+NCOMP+J)+AEI(ISTA+2*NCOMP+J)
          ZS(IFR+2)=' EBINR+T '
*           Energy diffusion by black hole
          VIT(IFR+3)=AEI(ISTA+3*NCOMP+J)
          ZS(IFR+3)=' BH:EDIFF'
*           Loss cone energy balance by black hole
          VIT(IFR+4)=AEI(ISTA+4*NCOMP+J)
          ZS(IFR+4)=' BH:LOSSC'
*           Escaper energy balance
          VIT(IFR+5)=AEI(ISTA+5*NCOMP+J)
          ZS(IFR+5)=' ENERESC '
*           Gravitational radiation from the cusp around black hole
          VIT(IFR+6)=AEI(ISTA+6*NCOMP+J)
          ZS(IFR+6)=' BH:GRRAD'
*           Indirect Heating due to single star escapers
          VIT(IFR+7)=AEI(ISTA+7*NCOMP+J)
          ZS(IFR+7)=' BI:INDHS'
*           Indirect Heating due to binary star escapers
          VIT(IFR+8)=AEI(ISTA+8*NCOMP+J)
          ZS(IFR+8)=' BI:INDHB'
*
          IFR=IFR+8
 980      CONTINUE
*
          VIT(IFR+1)=PHI(1)+PHIBIN(1)
          ZS(IFR+1)='  PHICEN '
          VIT(IFR+2)=DABS(PHI(1)+PHIBIN(1))/SIG**2
          ZS(IFR+2)='       X '
          IFR=IFR+2
*
	  DO 985 J=1,NLAGR
          VIT(IFR+J)=AVMASS(J)
          WRITE(ZS(IFR+J),984)J
 984      FORMAT(I2,' MAV   ')
 985      CONTINUE
          IFR=IFR+NLAGR
*
*           Calculation of core quantities
*
	  DO 1200 K=1,NC1
*
	  IF(K.LT.NC1)THEN
          IMIN=K
          IMAX=K
          ELSE
          IMIN=1
          IMAX=NCOMP
          END IF
*
	  SIGAV=0.D0
          RHOAV=0.D0
          XMCORE = 0.D0
          XNCORE = 0.D0
          SGCORE = 0.D0
          SMASS = 0.D0
*
      DO 1201 IK=IMIN,IMAX
      DEM = PI43*DEXP(X(I00+IK,1))*R(2)**3
      SMASS = SMASS + DEM
*          Take 3d sigma for core radius
      SIGAV = SIGAV + DEM*DSQRT(DEXP(X(I20+IK,1)-X(I00+IK,1))+
     *                          2.D0*DEXP(X(I02+IK,1)-X(I00+IK,1)))
      RHOAV = RHOAV + DEM*DEXP(X(I00+IK,1))
 1201     CONTINUE
      SIGAV=SIGAV/SMASS
      RHOAV=RHOAV/SMASS
*
	  RCC=SIGAV/DSQRT(4.D0*PI*GRAV*RHOAV)
*
	  DO 1202 IK=IMIN,IMAX
*
	  DO 1210 I=2,NJ
*
          RUP = R(I)
          IF(R(I-1).GT.RCC)GOTO 1202
          IF(R(I).GT.RCC)RUP=RCC
*
	  DEM = PI43*DEXP(X(I00+IK,I))*(RUP**3-R(I-1)**3)
         XMCORE = XMCORE + DEM
         XNCORE = XNCORE + DEM/XMIND(IK)
         SGCORE = SGCORE + DEM*DSQRT(C13*(DEXP(X(I20+IK,I)-X(I00+IK,I))+
     *                            2.D0*DEXP(X(I02+IK,I)-X(I00+IK,I))))
 1210     CONTINUE
 1202     CONTINUE
*
	  SGCORE = SGCORE/XMCORE
*
	  VIT(IFR+1)=RCC
          VIT(IFR+2)=SGCORE
          VIT(IFR+3)=XNCORE
          VIT(IFR+4)=XMCORE
          WRITE(ZS(IFR+1),1251)K
          WRITE(ZS(IFR+2),1252)K
          WRITE(ZS(IFR+3),1253)K
          WRITE(ZS(IFR+4),1254)K
 1251     FORMAT(I2,' RC    ')
 1252     FORMAT(I2,' SGC   ')
 1253     FORMAT(I2,' NC    ')
 1254     FORMAT(I2,' MC    ')
*
          IFR=IFR+4
 1200     CONTINUE
*
*           This is -DPTOT see BINSTO
          VIT(IFR+1)=-DPTOT
          ZS(IFR+1)=' -DPTOT  '
*           Cumulative mass loss due to formation of binaries,
*           3b and 4b escapers
          VIT(IFR+2)=XMERR
          ZS(IFR+2)='  XMERR  '
*           Cumulative thermal energy loss due to formation of binaries,
*           3b and 4b escapers
          VIT(IFR+3)=XEERR
          ZS(IFR+3)='  XEERR  '
*           Cum mass loss due to binary formation
          VIT(IFR+4)=XMERR2
          ZS(IFR+4)='  XMERR2 '
*           Cum energy loss due to binary formation
          VIT(IFR+5)=XEERR2
          ZS(IFR+5)='  XEERR2 '
*           Cum mass loss due to 3b kick escapers
          VIT(IFR+6)=XMERR3
          ZS(IFR+6)='  XMERR3 '
*           Cum energy loss due to 3b kick escapers
          VIT(IFR+7)=XEERR3
          ZS(IFR+7)='  XEERR3 '
*           Cum mass loss due to 4b kick escapers
          VIT(IFR+8)=XMERR4
          ZS(IFR+8)='  XMERR4 '
*           Cum energy loss due to 3b kick escapers
          VIT(IFR+9)=XEERR4
          ZS(IFR+9)='  XEERR4 '
*      Save total masses for case of tidal mass loss
          VIT(IFR+10)=DEXP(X(IMR,NJ))
          ZS(IFR+10)='   MTOT  '
          VIT(IFR+11)=XMRBIN(NJ)
          ZS(IFR+11)=' MBINTOT '
*
          IFR=IFR+11
*
          IF(NCOMP.GT.1)THEN
          I=3
          DO 568 J=2,NCOMP
          HELP=XMIND(J)*DEXP(X(I20+J,I)-X(I00+J,I))/
     *     XMIND(1)/DEXP(X(I20+1,I)-X(I00+1,I))
          VIT(IFR+J-1)=HELP
          ZS(IFR+J-1)='   EQUIP '
 568      CONTINUE
*
          IFR = IFR + NCOMP - 1
*
          END IF
*
*           Special data for collaborative experiment
          KK = NCOMP + 1
*
          XBR = 0.D0
          XBRH = 0.D0
          XBR1 = 0.D0
          XBT = 0.D0
          XBTH = 0.D0
          XBT1 = 0.D0
          SMASS = 0.D0
          XNMASS = 0.D0
*
*           Take total Lagrangian Radii
          DO 9900 J = 1,NLAGR
*
          IF(J.EQ.1)THEN
          IGR = 2
          RDOWN=0.D0
          ELSE
          IGR = ILAGR(J-1,KK) - 1
          RDOWN=RLAGR(J-1,KK)
          END IF
*
          RUP=RLAGR(J,KK)
*
 9950     CONTINUE
*
          IGR = IGR + 1
*
          IF(R(IGR-1).LT.RDOWN)THEN
          RLOW=RDOWN
          ELSE
          RLOW=R(IGR-1)
          END IF
*
         IF(R(IGR).GE.RUP)THEN
          RHIGH=RUP
          ELSE
          RHIGH=R(IGR)
          END IF
*
          K = 0
 9990     CONTINUE
          K = K + 1
* accumulate particle number times <v**2>
          DER3 = PI43*(RHIGH**3 - RLOW**3)
*
          HX(I00+K)=DEXP(X(I00+K,IGR))
          HX(I20+K)=DEXP(X(I20+K,IGR)-X(I00+K,IGR))
          HX(I02+K)=DEXP(X(I02+K,IGR)-X(I00+K,IGR))
          DERPR = DER3*HX(I20+K)/XMIND(K)
          DERPT = DER3*HX(I02+K)/XMIND(K)
*
          IF(J.LE.4)THEN
          XBR1 = XBR1 + DERPR
          XBT1 = XBT1 + DERPT
          END IF
          IF(J.LE.8)THEN
          XBRH = XBRH + DERPR
          XBTH = XBTH + DERPT
          END IF
          IF(J.LE.10)THEN
          XBR = XBR + DERPR
          XBT = XBT + DERPT
          END IF
          SMASS = SMASS + DER3*HX(I00+K)
          XNMASS = XNMASS + DER3*HX(I00+K)/XMIND(K)
*
          IF(K.LT.NCOMP)GOTO 9990
* Next radial mesh
          IF(IGR.LT.ILAGR(J,KK))GOTO 9950
* Save particle numbers and average masses
          IF(J.LE.4)THEN
          VIT(IFR+4)=SMASS/XNMASS
          XNM1=XNMASS
          END IF
          IF(J.LE.8)THEN
          VIT(IFR+5)=SMASS/XNMASS
          XNMH=XNMASS
          END IF
          IF(J.LE.10)THEN
          VIT(IFR+6)=SMASS/XNMASS
          XNM=XNMASS
          END IF
 9900      CONTINUE
* Divide by particle number to get correct average of <v**2>
          XBR1 = XBR1/XNM1
          XBT1 = XBT1/XNM1
          XBRH = XBRH/XNMH
          XBTH = XBTH/XNMH
          XBR = XBR/XNM
          XBT = XBT/XNM
          VIT(IFR+1) = 1.D0 - XBT1/XBR1
          VIT(IFR+2) = 1.D0 - XBTH/XBRH
          VIT(IFR+3) = 1.D0 - XBT/XBR
*
          ZS(IFR+1)=' BETA_10 '
          ZS(IFR+2)=' BETA_H  '
          ZS(IFR+3)=' BETA    '
          ZS(IFR+4)=' MBAR_10 '
          ZS(IFR+5)=' MBAR_H  '
          ZS(IFR+6)=' MBAR    '
*
          IFR = IFR + 6
*
* Save tidal quantities
          VIT(IFR+1) = RTIDAL
          VIT(IFR+2) = PHTID    
          ZS(IFR+1)=' RTIDAL '
          ZS(IFR+2)=' PHTIDAL'        
*
* Save last used IFR for possible use in binsto.f
	  IF(IZ00.EQ.1)THEN
          PRINT*,' Info: Last used IFR=',IFR,' IVDIM=',IVDIM
          IFRFIN = IFR
          IZ00=2
          IF(IFR.GT.IVDIM)STOP
          END IF
*
          WRITE(3,ERR=998,IOSTAT=IWRERR)(VIT(K),K=1,IVDIM)
          IREC=IREC+1
*           Store previous RHO,SIG,CC-Rate to avoid too much output
          VRHO=RHO
          VSIG=SIG
          VCSUM=CSUM
* Memorize content of VIT-vector in Unit 90
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
*
          CLOSE(90)
*
          END IF
*
          RETURN
C
 998      CONTINUE
          PRINT*,' Error occurred during Write on Unit 3 '
          PRINT*,' IOSTAT=',IWRERR
          PRINT*,' This was record number IREC=',IREC
          PRINT*,' FORTRAN STOP'
          STOP
C
 999      CONTINUE
          PRINT*,' Error occurred during Read from Unit 3 '
          PRINT*,' IOSTAT=',IREERR
          PRINT*,' Last Record found on Unit 3=',J
          PRINT*,' Program Execution continues...'
          GOTO 888
C
          END
