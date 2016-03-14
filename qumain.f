         SUBROUTINE QUMAIN(IB1,IB2,EKT1,EKT2,QTIME)
*
*        Interface for triple integration from gaseous model program
*        (from S. Aarseth NBODY6 Program)
*
*        IB1 Number of First Binary more strongly bound
*        IB2 Number of Second Binary less strongly bound
*        EKT kt in the vicinity of interacting binaries
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  M,MIJ
      REAL*8  XORB(2),VORB(2),XREL(3),VREL(3),PX(3),QX(3)
      REAL*8  XREL1(3),VREL1(3),XREL2(3),VRELB(3),EI1(3),EI2(3),
     &        XLI1(3),XLI2(3)
      REAL*8  XC1(3),XC1DOT(3),XC2(3),XC2DOT(3),XC3(3),XC3DOT(3)
      COMMON/CREG/M(4),X4(3,4),XDOT4(3,4),P(12),Q(12),TIME4,ENERGY,
     &      EPSR2,XR(9),W(9),R(6),TA(6),MIJ(6),CM(10),RMAX4,TMAX,
     &              DS,TSTEP,EPS,NSTEP4,NAME4(4),KZ15,KZ27,NREG,NFN
      COMMON/TPR/   SWITCH,GTYPE,GTYPE0
      COMMON/CONFIG/  R2(4,4),I1,I2,I3,I4
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3,IP(4)
      COMMON/CCOLL/  QK(12),PK(12),ICALL,ICOLL,NDISS4
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      COMMON/KSAVE/  K1,K2
      COMMON/EBSAVE/  EBS
      COMMON/IND6/  IND(6)
      COMMON/QQMIN/  QRMIN
      DATA  IND  /1,2,3,4,5,6/
      DATA  ISUBS /0/
      SAVE
*
      PARAMETER(NTAB=32)
      INTEGER IDUM2,K,IV,IY
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
      INTEGER ITRM4
      SAVE DEB4,ITRM4
      DATA ITRM4/0/
*
      REAL RANF
      EXTERNAL RANF
*
      INCLUDE 'params.f'
*
* 
      COMMON/BINARY/BODY1(NBINO),BODY2(NBINO),SIZE1(NBINO),SIZE2(NBINO),
     *ECC(NBINO),SEMIA(NBINO),RB(NBINO),TORB(NBINO),VR(NBINO),VT(NBINO),
     *EB(NBINO),RBMIN(NBINO),RBMAX(NBINO),XIN(NBINO),OMEG(NBINO),
     &XMRBIN(NJO),PHIBIN(NJO),
     *PHIOLD(NJO),EREM(NCOMPO),XHEAT(NCOMPO),XKICK(NCOMPO),
     *XMASS(NCOMPO),TBIN(NBINO),RHOBIN(NJO),VMRBIN(NJO)
      COMMON/BININT/INAME(NBINO),INAMEX(NBINO),
     *ISH(NBINO),ICO(NBINO),IEV(NBINO),NAMEB(NBINO),
     *IBSTA(NBINO),ILEN(NBINO)
*
      DATA ONEPI,C13,C12/3.141592654D0,0.333333333D0,0.5D0/
*       This data statement already in trmain.
*     DATA  EP  /0.04D0,0.0016D0,0.64D-4,0.256D-5/
*
*       Binary 1(stars 1 and 2) meets Binary 2 (stars 3 and 4)
*       Binary 2 has binding energy smaller than Binary 1
*       Transformation to c.m. of the quod system necessary
*
      BMRAT = MIN(BODY1(IB1),BODY2(IB1),BODY1(IB2),BODY2(IB2))/
     &           MAX(BODY1(IB1),BODY2(IB1),BODY1(IB2),BODY2(IB2))
      IF(BMRAT.LT.1.D-3)THEN
      BMKT = 0.1D0*BMRAT
      ELSE
      BMKT = 1.D0
      BMRAT = 1.D0
      END IF
*
      PRINT*,' Binary ',IB1,' N=',NAMEB(IB1),' invoked QUMAIN ',
     *  ' with Binary ',IB2,' N=',NAMEB(IB2)
      CALL FLUSH(6)
*       Here ECMOL1,2 in cluster frame
      XMBIN1=BODY1(IB1)+BODY2(IB1)
      XMBIN2=BODY1(IB2)+BODY2(IB2)
      XMQUAD = XMBIN1 + XMBIN2
      ECMOL1=XMBIN1*(VR(IB1)**2+VT(IB1)**2)/2.D0
      ECMOL2=XMBIN2*(VR(IB2)**2+VT(IB2)**2)/2.D0
      EBOLD1=EB(IB1)
      EBOLD2=EB(IB2)
      ECCOL1=ECC(IB1)
      ECCOL2=ECC(IB2)
      XINOL1=XIN(IB1)
      XINOL2=XIN(IB2)
      OMGOL1=OMEG(IB1)
      OMGOL2=OMEG(IB2)
*
      WRITE(6,887)QTIME,ECMOL1,ECMOL2,EBOLD1,EBOLD2
      WRITE(6,888)QTIME,EKT1,EBOLD1,EBOLD1/BMRAT/EKT1,
     &    SEMIA(IB1),ECC(IB1)
      WRITE(6,889)QTIME,EKT2,EBOLD2,EBOLD2/BMRAT/EKT2,
     &    SEMIA(IB2),ECC(IB2)
 887  FORMAT(1X,' Start with T,ECM1,2,EB1,2=',1P,5D13.5)
 888  FORMAT(1X,' Start T,kT1,EB1,EB1/BMRAT/kT,a,e=',1P,6D13.5)
 889  FORMAT(1X,' Start T,kT2,EB2,EB2/BMRAT/kT,a,e=',1P,6D13.5)
*
      GRAV = 1.D0
      TWOPI = 2.D0*ONEPI
*
 7891 CONTINUE
      FI = TWOPI*RANF()
      SINFI = SIN(FI)
      COSFI = COS(FI)
*
*     RGRAVT defined in the same way as RCRIT in computation of
*     interaction probability in binsto (quad interaction)
*
*      RGRAV0 = 2.5D0*SEMIA(IB2)
*      RGRAVT = RGRAV0
*
*      RGRAV1 = (EB(IB2)/EKT2)**C13*SEMIA(IB2)
*      IF(RGRAV1.GT.RGRAV0)RGRAVT=RGRAV1
*
*            NEW NEW NEW (the same procedure as in trmain.f)                                                                                    
*    calculate Pmax using formula given by Hut and Bahcal ApJ 268, 1983
*     and Bacon and Sigurdsson astro-ph96/03036
*
      VRELX = VR(IB2)-VR(IB1)
      VRELY = VT(IB2)*SINFI
      VRELZ = VT(IB2)*COSFI-VT(IB1)
      VRELT2 = VRELX**2+VRELY**2+VRELZ**2
      VRELT = DSQRT(VRELT2)
*
      VC2RAW = 2.D0*XMQUAD*(EBOLD1 + EBOLD2)/(XMBIN1*XMBIN2)
      VC2 = VC2RAW/BMRAT
*
      VVC2 = DSQRT(VRELT2/VC2)                   
*
      IF(BMRAT.LT.1.D-3)THEN 
*    Use 100*a*(1+e) for planetary system run
          PMAX = 100.D0*SEMIA(IB2)*(1.0D0 + ECC(IB2))
      ELSE
          PMAX = SEMIA(IB2)*(0.6D0*(1.0D0 + ECCOL2) + 5.0D0/VVC2)
      END IF
*
*    New - to speed up the program
*      IF(PMAX.GT.10.D0*SEMIA(IB2)) PMAX = 10.D0*SEMIA(IB2)
*
      PMAX2 = PMAX*PMAX
*
*     Calculate minimum distance RGRAVT for given PMAX                          
*                                                                               
      RXXX = GRAV*XMQUAD/VRELT2                                                 
      RGRAVT = RXXX*(DSQRT(1.D0 + (PMAX/RXXX)**2) - 1.D0)                         
      RGRAV1 = RGRAVT                                                           
      RGRAV0 = RGRAVT                                                           
*
      WRITE(6,3459) XMBIN1,XMBIN2,XMQUAD,VRELT2,VC2,VVC2,RXXX,
     &    BMRAT,VC2RAW
 3459 FORMAT(1X,'XMBIN1,2,XMQUAD,VRELT2,VC2,VVC2,RXXX,BMRAT,VC2RAW= ',
     &    1P,9E12.4)
*
      VCMTX = (XMBIN2*VR(IB2) + XMBIN1*VR(IB1))/XMQUAD
      VCMTY = XMBIN2*VT(IB2)*SINFI/XMQUAD
      VCMTZ = (XMBIN2*VT(IB2)*COSFI + XMBIN1*VT(IB1))/XMQUAD
      ECMTOT = XMQUAD*(VCMTX**2+VCMTY**2+VCMTZ**2)/2.D0
*
      VCM1X = -XMBIN2/XMQUAD*VRELX
      VCM1Y = -XMBIN2/XMQUAD*VRELY
      VCM1Z = -XMBIN2/XMQUAD*VRELZ      
      VCM2X = XMBIN1/XMQUAD*VRELX
      VCM2Y = XMBIN1/XMQUAD*VRELY
      VCM2Z = XMBIN1/XMQUAD*VRELZ  
*
*       Save kinetic energy of binaries in 4b c.m. frame before encounter
*       at infinity.
      ECMOL1 = C12*XMBIN1*(VCM1X**2+VCM1Y**2+VCM1Z**2)
      ECMOL2 = C12*XMBIN2*(VCM2X**2+VCM2Y**2+VCM2Z**2)  
      ECMOLD = ECMOL1 + ECMOL2
*     PRINT*,' After Transformation ECMTOT=',ECMTOT,' ECM=',ECMOLD
*     PRINT*,' After Transformation Sum=',ECMTOT+ECMOLD
*      
*      Determine max impact parameter for rmin = 2 rgrav
*      PMAX2 = RGRAVT*(RGRAVT + 2.D0*GRAV*XMQUAD/
*     *   VRELT/VRELT)
*
      KZ15 = 2
      KZ27 = 2
      NAME4(1) = 1
      NAME4(2) = 2
      NAME4(3) = 3
      NAME4(4) = 4
*
*       first binary
*
      IF (BODY1(IB1).LT.BODY2(IB1)) THEN
          M(1) = BODY1(IB1)
          M(2) = BODY2(IB1)
          SIZE(1) = 0.D0
          SIZE(2) = 0.D0
      ELSE
          M(1) = BODY2(IB1)
          M(2) = BODY1(IB1)
          SIZE(1) = 0.D0
          SIZE(2) = 0.D0
      END IF
*
*       second binary
*     
      IF (BODY1(IB2).LT.BODY2(IB2)) THEN
          M(3) = BODY1(IB2)
          M(4) = BODY2(IB2)
          SIZE(3) = 0.D0
          SIZE(4) = 0.D0
      ELSE
          M(3) = BODY2(IB2)
          M(4) = BODY1(IB2)
          SIZE(3) = 0.D0
          SIZE(4) = 0.D0
      END IF
*
      DO 2 K = 1,7
          CM(K) = 0.0D0
    2 CONTINUE
*
      PP2 = RANF()*PMAX2
      PP = DSQRT(PP2)
*
*      Select start distance for quadruple for approximate perturbation
*      of 1.e-3 and assume worst case of apocentres aligned.
*      WARNING! Check this carefully for multi-mass case!! Definitely wrong!!
*      XMUE = M(3)*M(4)/XMBIN2
*      RX = 3.D0*SEMIA(IB2)*(1.D0+ECC(IB2))*
*     *                          (1.D0+1.D1*(XMBIN1/XMUE)**C13)
*
*     Select distance to start quad for approximate force perturbation
*     gama_crit=1.e-3 and assume the worst case of apocentre directed to single 
*     star. Determination of Rx based on ksinit.f from Nbody. Modifications:
*     1) for m3=0 Rx=A*a(1+e)
*     2) A is choosen to 3
*     3) perturber: binary1
*
      XMUE = 8.D0*XMBIN1/XMBIN2
      RX = 3.D0*SEMIA(IB2)*(1.D0+ECC(IB2))*(1.D0+1.D1*(XMUE)**C13)
*
*     determine the minimum distance for actual impact parameter pp
*
      RGRAVX = RXXX*(DSQRT(1.D0 + (PP/RXXX)**2) - 1.D0)  
*      Take care if RX based on perturbation gives too small value
      IF(RX.LT.2.D0*RGRAVX) RX = 2.D0*RGRAVX
*
      WRITE(6,2003)QTIME,RX,RGRAV0,RGRAVX,PP,DSQRT(PMAX2)
 2003 FORMAT(1X,' Start with T,RX=',1P,2D13.5,' RGRAV0,X=',
     * 2D13.5,' PP,PMAX=',2D13.5)
*
      VRELT2 = VRELT*VRELT
      ECC2 = 1.D0 + (PP*VRELT2/(GRAV*XMQUAD))**2
      ECCT = DSQRT(ECC2)
      XJ = PP*VRELT
*     PRINT*,' ECCT=',ECCT,' PP=',PP,' PMAX=',DSQRT(PMAX2),
*    * ' RX=',RX,' RGRAVT=',RGRAVT
*     PRINT*,' Arg Theta=',(XJ**2/(RX*GRAV*XMQUAD)-1.D0),
*    * (XJ**2/(RX*GRAV*XMQUAD))/ECCT-1.D0/ECCT
      THETAX = DACOS((XJ**2/(RX*GRAV*XMQUAD)-1.D0)/ECCT)
      PSI = DATAN(PP*VRELT2/GRAV/XMQUAD)
      ALPHA = ONEPI - THETAX - PSI
      DEGR1 = ONEPI/180.
*     PRINT*,' ALPHA, THETAX,PSI,PI=',ALPHA/DEGR1,THETAX/DEGR1,
*    * PSI/DEGR1,ONEPI/DEGR1,VRELT
*
      VTX = PP*VRELT/RX
      VRX = -DSQRT(VRELT2 - VTX**2 + 2.D0*GRAV*XMQUAD/RX)
*
*       Velocity and Position Components of binaries 1 and 2
*
      XRELP = RX*DCOS(ALPHA)
      YRELP = RX*DSIN(ALPHA)
*
      XVRELP = (VRX*DCOS(ALPHA) - VTX*DSIN(ALPHA))
      YVRELP = (VRX*DSIN(ALPHA) + VTX*DCOS(ALPHA))
*
*     PRINT*,' XRELP,YRELP=',XRELP,YRELP,' XVRELP,YVRELP=',
*    *    XVRELP,YVRELP,' VRX,VTX=',VRX,VTX
*
      EPOT2 = -GRAV*XMBIN2*XMBIN1/RX
*     PRINT*,' ***XXXXXX E CHECK****'
      VVREL2=XVRELP**2 + YVRELP**2
      ERELOL=XMBIN2*XMBIN1/XMQUAD*VVREL2/2.D0
*     PRINT*,' ECMTOT, ERELOL, Sum+EPOT=',ECMTOT,ERELOL,
*    *   ECMTOT+ERELOL+EPOT2
*     CALL FLUSH(6)
*
      VRELXY = DSQRT(VRELX**2+VRELY**2)
      XXI = DACOS(VRELZ/VRELT)
      SINXI = SIN(XXI)
      COSXI = COS(XXI)
*
      PPSI = DACOS(VRELX/VRELXY)
      COSPSI = COS(PPSI)
      SINPSI = SIN(PPSI)
*
*     PRINT*,' XI, SIN,COS=',XXI,SINXI,COSXI
*     PRINT*,' PSI, SIN,COS=',PPSI,SINPSI,COSPSI
*
      XXREL = XRELP*SINXI*COSPSI - YRELP*SINPSI
      YREL = XRELP*SINXI*SINPSI + YRELP*COSPSI
      ZREL = XRELP*COSXI 
*
      XVREL = XVRELP*SINXI*COSPSI - YVRELP*SINPSI
      YVREL = XVRELP*SINXI*SINPSI + YVRELP*COSPSI
      ZVREL = XVRELP*COSXI 
*
*     cm of binary 2 stored in 3
*
      X4(1,3) = XMBIN1/XMQUAD*XXREL
      X4(2,3) = XMBIN1/XMQUAD*YREL
      X4(3,3) = XMBIN1/XMQUAD*ZREL
*
      XDOT4(1,3) = XMBIN1/XMQUAD*XVREL
      XDOT4(2,3) = XMBIN1/XMQUAD*YVREL
      XDOT4(3,3) = XMBIN1/XMQUAD*ZVREL
*
*     PRINT*,' Before QUAD CALL    ',-VRELT,0.D0
*     PRINT*,' cm Velocity Bin 2=',XDOT4(1,3),XDOT4(2,3),XDOT4(3,3)
*     PRINT*,' cm Position Bin 2=',X4(1,3),X4(2,3),X4(3,3)
*
*     cm of binary 1 in 1
*
      X4(1,1) = -XMBIN2/XMQUAD*XXREL
      X4(2,1) = -XMBIN2/XMQUAD*YREL
      X4(3,1) = -XMBIN2/XMQUAD*ZREL
*
      XDOT4(1,1) = -XMBIN2/XMQUAD*XVREL
      XDOT4(2,1) = -XMBIN2/XMQUAD*YVREL
      XDOT4(3,1) = -XMBIN2/XMQUAD*ZVREL
*
*     PRINT*,' cm Velocity Bin 1=',XDOT4(1,3),XDOT4(2,3),XDOT4(3,3)
*     PRINT*,' cm Position Bin 1=',X4(1,3),X4(2,3),X4(3,3)
*
*     PRINT*,' ***SECOND E CHECK****'
      VVREL2=(XDOT4(1,1)-XDOT4(1,3))**2+
     *       (XDOT4(2,1)-XDOT4(2,3))**2+
     *       (XDOT4(3,1)-XDOT4(3,3))**2
      ERELOL=XMBIN1*XMBIN2/XMQUAD*VVREL2/2.D0
      EBINA2=XMBIN2*(XDOT4(1,3)**2+XDOT4(2,3)**2+XDOT4(3,3)**2)/2.D0
      EBINA1=XMBIN1*(XDOT4(1,1)**2+XDOT4(2,1)**2+XDOT4(3,1)**2)/2.D0
      ENERG2=EBINA1+EBINA2-EB(IB1)-EB(IB2)
*     PRINT*,' ECMTOT, ERELOL, Sum=',ECMTOT,ERELOL,
*    *   ECMTOT+ERELOL+EPOT2
*     PRINT*,' EBINA2,EBINA1,Sum+ECMTOT+EPOT=',EBINA2,EBINA1,
*    *   ECMTOT+EBINA2+EBINA1+EPOT2
*
*       Randomize perihelion, node & inclination first binary.
          PEI = TWOPI*RANF()
          OMEGA = TWOPI*RANF()
          ZI = 0.25*TWOPI*RANF()
*
*       Set transformation elements (Brouwer & Clemence p. 35).
          PX(1) = COS(PEI)*COS(OMEGA) - SIN(PEI)*SIN(OMEGA)*COS(ZI)
          QX(1) =-SIN(PEI)*COS(OMEGA) - COS(PEI)*SIN(OMEGA)*COS(ZI)
          PX(2) = COS(PEI)*SIN(OMEGA) + SIN(PEI)*COS(OMEGA)*COS(ZI)
          QX(2) =-SIN(PEI)*SIN(OMEGA) + COS(PEI)*COS(OMEGA)*COS(ZI)
          PX(3) = SIN(PEI)*SIN(ZI)
          QX(3) = COS(PEI)*SIN(ZI)
*
*       Specify relative motion at apocentre and sum binding energy.
          XORB(1) = SEMIA(IB1)*(1.0 + ECC(IB1))
          XORB(2) = 0.0
          VORB(1) = 0.0
          VORB(2) = SQRT(XMBIN1*(1.0D0 - ECC(IB1))/
     *                         (SEMIA(IB1)*(1.0D0 + ECC(IB1))))
*
*       Transform to relative variables.
          DO 401 K = 1,3
              XREL(K) = PX(K)*XORB(1) + QX(K)*XORB(2)
              VREL(K) = PX(K)*VORB(1) + QX(K)*VORB(2)
 401    CONTINUE
*
        DIST1 = DSQRT(XREL(1)**2+XREL(2)**2+XREL(3)**2)
        VREL2 = VREL(1)**2+VREL(2)**2+VREL(3)**2
        SEMIX1 = 2.D0/DIST1-VREL2/XMBIN1
        SEMIX1 = 1.D0/SEMIX1
        EBX1 = C12*(XMBIN1/2.D0)**2/SEMIX1
*       PRINT*,' EBX1 bef Quad=',EBX1
*
          I1 = 1
          I2 = 2
*
*       Set global variables for each component.
          DO 500 K = 1,3
              X4(K,I1) = X4(K,I1) + M(I2)*XREL(K)/XMBIN1
              X4(K,I2) = X4(K,I1) - XREL(K)
              XDOT4(K,I1) = XDOT4(K,I1) + M(I2)*VREL(K)/XMBIN1
              XDOT4(K,I2) = XDOT4(K,I1) - VREL(K)
 500    CONTINUE
*
*     PRINT*,' Masses =',M(1),M(2),ECMOLD
*     PRINT*,' Velocity 1=',XDOT4(1,1),XDOT4(2,1),XDOT4(3,1)
*     PRINT*,' Velocity 2=',XDOT4(1,2),XDOT4(2,2),XDOT4(3,2)
*     PRINT*,' Position 1=',X4(1,1),X4(2,1),X4(3,1)
*     PRINT*,' Position 2=',X4(1,2),X4(2,2),X4(3,2)
*
*       Randomize perihelion, node & inclination second binary.
          PEI = TWOPI*RANF()
          OMEGA = TWOPI*RANF()
          ZI = 0.25*TWOPI*RANF()
*
*       Set transformation elements (Brouwer & Clemence p. 35).
          PX(1) = COS(PEI)*COS(OMEGA) - SIN(PEI)*SIN(OMEGA)*COS(ZI)
          QX(1) =-SIN(PEI)*COS(OMEGA) - COS(PEI)*SIN(OMEGA)*COS(ZI)
          PX(2) = COS(PEI)*SIN(OMEGA) + SIN(PEI)*COS(OMEGA)*COS(ZI)
          QX(2) =-SIN(PEI)*SIN(OMEGA) + COS(PEI)*COS(OMEGA)*COS(ZI)
          PX(3) = SIN(PEI)*SIN(ZI)
          QX(3) = COS(PEI)*SIN(ZI)
*
*       Specify relative motion at apocentre and sum binding energy.
          XORB(1) = SEMIA(IB2)*(1.0 + ECC(IB2))
          XORB(2) = 0.0
          VORB(1) = 0.0
          VORB(2) = SQRT(XMBIN2*(1.0D0 - ECC(IB2))/
     *                         (SEMIA(IB2)*(1.0D0 + ECC(IB2))))
*
*       Transform to relative variables.
          DO 400 K = 1,3
              XREL(K) = PX(K)*XORB(1) + QX(K)*XORB(2)
              VREL(K) = PX(K)*VORB(1) + QX(K)*VORB(2)
 400    CONTINUE
*
        DIST1 = DSQRT(XREL(1)**2+XREL(2)**2+XREL(3)**2)
        VREL2 = VREL(1)**2+VREL(2)**2+VREL(3)**2
        SEMIX1 = 2.D0/DIST1-VREL2/XMBIN1
        SEMIX1 = 1.D0/SEMIX1
        EBX1 = C12*(XMBIN1/2.D0)**2/SEMIX1
*       PRINT*,' EBX2 bef Quad=',EBX1   
*
          I3 = 3
          I4 = 4
*
*       Set global variables for each component.
          DO 501 K = 1,3
              X4(K,I3) = X4(K,I3) + M(I4)*XREL(K)/XMBIN2
              X4(K,I4) = X4(K,I3) - XREL(K)
              XDOT4(K,I3) = XDOT4(K,I3) + M(I4)*VREL(K)/XMBIN2
              XDOT4(K,I4) = XDOT4(K,I3) - VREL(K)
 501    CONTINUE
*
*     PRINT*,' Masses =',M(3),M(4),ECMOLD
*     PRINT*,' Velocity 3=',XDOT4(1,3),XDOT4(2,3),XDOT4(3,3)
*     PRINT*,' Velocity 4=',XDOT4(1,4),XDOT4(2,4),XDOT4(3,4)
*     PRINT*,' Position 3=',X4(1,3),X4(2,3),X4(3,3)
*     PRINT*,' Position 4=',X4(1,4),X4(2,4),X4(3,4)
*
*       Transform to the local c.m. reference frame.
      DO 4 L = 1,4
          CM(7) = CM(7) + M(L)
          DO 3 K = 1,3
              CM(K) = CM(K) + M(L)*X4(K,L)
              CM(K+3) = CM(K+3) + M(L)*XDOT4(K,L)
    3     CONTINUE
    4 CONTINUE
*
*       Set c.m. coordinates & velocities for triple system.
      DO 5 K = 1,6
          CM(K) = CM(K)/CM(7)
    5 CONTINUE
*
*     PRINT*,' CM: ',(CM(K),K=1,7)
*     PRINT*,' Now Call QUAD with Coordinates:'
*     PRINT*,' Binary 1-1=',X4(1,1),X4(2,1),X4(3,1)
*     PRINT*,' Binary 1-2=',X4(1,2),X4(2,2),X4(3,2)
*     PRINT*,' Binary 2-1=',X4(1,3),X4(2,3),X4(3,3)
*     PRINT*,' Binary 2-2=',X4(1,4),X4(2,4),X4(3,4)
*     PRINT*,' Velocities:'
*     PRINT*,' Binary 1-1 DOT=',XDOT4(1,1),XDOT4(2,1),XDOT4(3,1)
*     PRINT*,' Binary 1-2 DOT=',XDOT4(1,2),XDOT4(2,2),XDOT4(3,2)
*     PRINT*,' Binary 2-1 DOT=',XDOT4(1,3),XDOT4(2,3),XDOT4(3,3)
*     PRINT*,' Binary 2-2 DOT=',XDOT4(1,4),XDOT4(2,4),XDOT4(3,4)
*
*        calculate total energy (like as in NEWSYS in NODY6)
*
      T = 0.0D0
      V = 0.0D0
      DO 10 I = 1,4
          T = T + C12*M(I)*(XDOT4(1,I)**2 + XDOT4(2,I)**2 +
     *               XDOT4(3,I)**2)
          IF (I.EQ.4) GO TO 10
          J1 = I + 1
          DO 9 J = J1,4
                  RR = (X4(1,J) - X4(1,I))**2 +
     *            (X4(2,J) - X4(2,I))**2 + (X4(3,J) - X4(3,I))**2
              V = V - M(I)*M(J)/DSQRT(RR)
    9     CONTINUE
   10 CONTINUE
*
      ENERGY = T + V
      GAM = T - V
*        put into CM(8) for later use in quad
      CM(8) = ENERGY
*     PRINT*,' ENERGY,ENERG2=',ENERGY,ENERG2
*     CALL FLUSH(6)
*
*
*     DIST = 0.D0
*     VDOT = 0.D0
*     DO 450 K=1,3
*     XC1(K)=(M(1)*X4(K,1)+M(2)*X4(K,2))/XMBIN1
*     XC1DOT(K)=(M(1)*XDOT4(K,1)+M(2)*XDOT4(K,2))/XMBIN1
*     XC2(K)=(M(3)*X4(K,3)+M(4)*X4(K,4))/XMBIN2
*     XC2DOT(K)=(M(3)*XDOT4(K,3)+M(4)*XDOT4(K,4))/XMBIN2
*     XREL(K) = XC1(K)-XC2(K)
*     VREL(K) = XC1DOT(K)-XC2DOT(K)
*     DIST = DIST + XREL(K)*XREL(K)
*     VDOT = VDOT + VREL(K)*VREL(K)
*450  CONTINUE
*     DIST=SQRT(DIST)
*     VDOT=SQRT(VDOT)
*
*     PRINT*,' Binary 1 CM=',(XC1(K),K=1,3)
*     PRINT*,' Binary 1 CMDOT=',(XC1DOT(K),K=1,3)
*     PRINT*,' Binary 2 CM=',(XC2(K),K=1,3)
*     PRINT*,' Binary 2 CMDOT=',(XC2DOT(K),K=1,3)
*
      I = 0
*
*     PRINT*,' Before CALL QUAD CM7=',CM(7)
*     PRINT*,' RX/VRELT=',RX/VRELT,' TORBBIN1=',
*    *       2.D0*ONEPI*SEMIA(IB1)**1.5D0/DSQRT(XMBIN1)
*     PRINT*,' quma: cm7=',M(1),X4(1,1),XDot4(1,1),P(1),Q(1),CM(7)
*     PRINT*,' quma: ',xr(1),w(1),r(1),ta(1),mij(1),name4(1)
*     CALL FLUSH(6)
*
*       Evaluate orbital elements.
      RM1 = 0.D0
      RM2 = 0.D0
      VVX1 = 0.D0
      VVX2 = 0.D0
      VVX3 = 0.D0
      VVX4 = 0.D0
      VCMBIN = 0.D0
      RM13 = 0.D0
      RM14 = 0.D0
      RM23 = 0.D0
      RM24 = 0.D0
*
      DO 775 K = 1,3
          VVX1 = VVX1 + XDOT4(K,I1)**2
          VVX2 = VVX2 + XDOT4(K,I2)**2
          VVX3 = VVX3 + XDOT4(K,I3)**2
          VVX4 = VVX4 + XDOT4(K,I4)**2
          RM1 = RM1 + (X4(K,I1)-X4(K,I2))**2
          RM2 = RM2 + (X4(K,I3)-X4(K,I4))**2
          RM13 = RM13 + (X4(K,I1)-X4(K,I3))**2
          RM14 = RM14 + (X4(K,I1)-X4(K,I4))**2  
          RM23 = RM23 + (X4(K,I2)-X4(K,I3))**2  
          RM24 = RM24 + (X4(K,I2)-X4(K,I4))**2  
*
  775 CONTINUE
      VCMBIN2 = C12*VCMBIN2/(M(I1)+M(I2))
      VVX1 = 0.5D0*M(I1)*VVX1
      VVX2 = 0.5D0*M(I2)*VVX2
      VVX3 = 0.5D0*M(I3)*VVX3
      VVX4 = 0.5D0*M(I4)*VVX4
*
      RM1 = DSQRT(RM1)
      RM13 = DSQRT(RM13)
      RM14 = DSQRT(RM14)
*
      RM2 = DSQRT(RM2)
      RM23 = DSQRT(RM23)
      RM24 = DSQRT(RM24)
*
      EGX12 = -M(I1)*M(I2)/RM1
      EGX13 = -M(I1)*M(I3)/RM13
      EGX14 = -M(I1)*M(I4)/RM14
      EGX23 = -M(I2)*M(I3)/RM23
      EGX24 = -M(I2)*M(I4)/RM24
      EGX34 = -M(I3)*M(I4)/RM2
*
      FX12 = -EGX12/RM1
      FX13 = -EGX13/RM13
      FX14 = -EGX14/RM14
      FX23 = -EGX23/RM23
      FX24 = -EGX24/RM24
      FX34 = -EGX34/RM2
*
      ETTX = VVX1 + VVX2 + VVX3 + VVX4 +
     *       EGX12 + EGX13 + EGX14 + EGX23 + EGX24 + EGX34
      ETTX1 = ETTX - EB(IB1)
      ETTX2 = ETTX1 - EB(IB2)
*
*     WRITE(6,1000)QTIME,FX12,FX13,FX14,FX23,FX24,FX34,ETTX,ETTX1,ETTX2
*1000 FORMAT(' T,Forces before QUAD=',1P,7(D12.5,2X),/,
*    *       ' ETOT1-3 before QUAD=',3(D12.5,1X))
*
*     PRINT*,' EGX12,13,14,23,24,34=',EGX12,EGX13,EGX14,EGX23,
*    *   EGX24,EGX34
*
*     PRINT*,' ETTX=',ETTX
*           
      CALL FLUSH(6)
*
      CALL QUAD(I,RX,QTIME,EKT1)
*       Return from Quad with I1,I2 indices of closest pair and RX=RGRAV.
*
      RGRAV = RX
      XMBIN1 = M(I1) + M(I2)
      XMBIN2 = M(I3) + M(I4)
      XMQUAD = XMBIN1 + XMBIN2    
*
*       Save binary c.m. data for later use (in system where triple rests)
      DO 74 K=1,3
      XC1(K)=(M(I1)*X4(K,I1)+M(I2)*X4(K,I2))/XMBIN1
      XC1DOT(K)=(M(I1)*XDOT4(K,I1)+M(I2)*XDOT4(K,I2))/XMBIN1
      XC2(K)=(M(I3)*X4(K,I3)+M(I4)*X4(K,I4))/XMBIN2
      XC2DOT(K)=(M(I3)*XDOT4(K,I3)+M(I4)*XDOT4(K,I4))/XMBIN2
 74   CONTINUE              
*
*       Evaluate orbital elements.
      VREL21 = 0.0D0
      RDOT1 = 0.0D0
      VREL22 = 0.0D0
      RDOT2 = 0.0D0
      RI21 = 0.0D0
      RI22 = 0.0D0
*
      DO 75 K = 1,3
         XREL1(K) = X4(K,I1) - X4(K,I2)
         VREL1(K) = XDOT4(K,I1) - XDOT4(K,I2)
         XREL2(K) = X4(K,I3) - X4(K,I4)
         VRELB(K) = XDOT4(K,I3) - XDOT4(K,I4)
*
          RDOT1 = RDOT1 + XREL1(K)*VREL1(K)
          VREL21 = VREL21 + VREL1(K)**2
          RDOT2 = RDOT2 + XREL2(K)*VRELB(K)
          VREL22 = VREL22 + VRELB(K)**2
          RI21 = RI21 + XREL1(K)**2
          RI22 = RI22 + XREL2(K)**2
*
   75 CONTINUE
*
*       Construct Runge-Lenz vector (Heggie & Rasio, 1995, IAU174, Eq.(5)).
          EI21 = 0.0D0
          XLI21 = 0.0D0
          CHECK1 = 0.0D0
          EI22 = 0.0D0
          XLI22 = 0.0D0
          CHECK2 = 0.0D0
*
          DO 76 K = 1,3
              KP1 = 1 + MOD(K,3)
              KP2 = 1 + MOD(K+1,3)
              EI1(K) = (VREL21*XREL1(K) - RDOT1*VREL1(K))/XMBIN1 -
     &                                        XREL1(K)/DSQRT(RI21)
              EI2(K) = (VREL22*XREL2(K) - RDOT2*VRELB(K))/XMBIN2 -
     &                                        XREL2(K)/DSQRT(RI22)
              XLI1(K) = XREL1(KP1)*VREL1(KP2)-XREL1(KP2)*VREL1(KP1)
              XLI2(K) = XREL2(KP1)*VRELB(KP2)-XREL2(KP2)*VRELB(KP1)
              EI21 = EI21 + EI1(K)**2
              EI22 = EI22 + EI2(K)**2
              XLI21 = XLI21 + XLI1(K)**2
              XLI22 = XLI22 + XLI2(K)**2
              CHECK1 = CHECK1 + EI1(K)*XLI1(K)
              CHECK2 = CHECK2 + EI2(K)*XLI2(K)
   76     CONTINUE
*
          XIARG1 = DSQRT((XLI1(1)**2+XLI1(2)**2)/XLI21)
          IF (XLI1(3).GT.0.D0) THEN
              XINCL1 = ASIN(XIARG1)
              ELSE
              XINCL1 = TWOPI/2.D0 - ASIN(XIARG1)
          END IF
*
          XOARG1 = EI1(2)/DSQRT(EI1(1)**2+EI1(2)**2)
          IF (EI1(1).GT.0.D0) THEN
          IF (EI1(2).GT.0.D0) XOMEG1 = ASIN(XOARG1)
          IF (EI1(2).LT.0.D0) XOMEG1 = TWOPI - ASIN(-XOARG1)
          ELSE
          IF (EI1(2).LT.0.D0) XOMEG1 = TWOPI/2.D0 + ASIN(-XOARG1)
          IF (EI1(2).GT.0.D0) XOMEG1 = TWOPI/2.D0 - ASIN(XOARG1)
          END IF
*
          XIARG2 = DSQRT((XLI2(1)**2+XLI2(2)**2)/XLI22)
          IF (XLI2(3).GT.0.D0) THEN
              XINCL2 = ASIN(XIARG2)
              ELSE
              XINCL2 = TWOPI/2.D0 - ASIN(XIARG2)
          END IF
*
          XOARG2 = EI2(2)/DSQRT(EI2(1)**2+EI2(2)**2)
          IF (EI2(1).GT.0.D0) THEN
          IF (EI2(2).GT.0.D0) XOMEG2 = ASIN(XOARG2)
          IF (EI2(2).LT.0.D0) XOMEG2 = TWOPI - ASIN(-XOARG2)
          ELSE
          IF (EI2(2).LT.0.D0) XOMEG2 = TWOPI/2.D0 + ASIN(-XOARG2)
          IF (EI2(2).GT.0.D0) XOMEG2 = TWOPI/2.D0 - ASIN(XOARG2)
          END IF
*
      XIN(IB1) = XINCL1
      OMEG(IB1) = XOMEG1
      XIN(IB2) = XINCL2
      OMEG(IB2) = XOMEG2
*       Determine semi-major axis & eccentricity of inner binary.
*
      RM1 = DSQRT(R2(I1,I2))
      RM2 = DSQRT(R2(I3,I4))
*
      SEMI1 = 2.0D0/RM1 - VREL21/XMBIN1
      SEMI1= 1.0/SEMI1
      SEMI2 = 2.0D0/RM2 - VREL22/XMBIN2
      SEMI2= 1.0/SEMI2
*       Save new binary data for gaseous model.
      EB(IB1) = C12*M(I1)*M(I2)/SEMI1
      EB(IB2) = C12*M(I3)*M(I4)/SEMI2            
*
      ECC(IB1) = SQRT((1.0D0 - RM1/SEMI1)**2 + RDOT1**2/(SEMI1*XMBIN1))
      ECC(IB2) = SQRT((1.0D0 - RM2/SEMI2)**2 + RDOT2**2/(SEMI2*XMBIN2))
*
      SEMIA(IB1) = SEMI1
      SEMIA(IB2) = SEMI2
      BODY1(IB1) = M(I1)
      BODY2(IB1) = M(I2)
      BODY1(IB2) = M(I3)
      BODY2(IB2) = M(I4)
*
*     SIZE1(IB1) = SIZE(I1)
*     SIZE2(IB1) = SIZE(I2)
*     SIZE1(IB2) = SIZE(I3)
*     SIZE2(IB2) = SIZE(I4)
*               
      WRITE(6,3001)QTIME,I1,I2,RM1,SEMI1,ECC(IB1),EB(IB1),EKT1
 3001 FORMAT(1X,' After Quad Bin 1 I1,I2,RM,SEMI,E,EB,EKT=',
     &   1P,D13.5,2I5,5D13.5)
      WRITE(6,3002)QTIME,I3,I4,RM2,SEMI2,ECC(IB2),EB(IB2),EKT2
 3002 FORMAT(1X,' After Quad Bin 2 I3,I4,RM,SEMI,E,EB,EKT=',
     &   1P,D13.5,2I5,5D13.5)
      CALL FLUSH(6)
*
      VVX1 = 0.D0
      VVX2 = 0.D0
      VVX3 = 0.D0
      VVX4 = 0.D0
      VCMBIN2 = 0.D0
*
      DO 875 K = 1,3
          VCMBIN2 = VCMBIN2 + (M(I1)*XDOT4(K,I1)+M(I2)*XDOT4(K,I2))**2
          VVX1 = VVX1 + XDOT4(K,I1)**2
          VVX2 = VVX2 + XDOT4(K,I2)**2
          VVX3 = VVX3 + XDOT4(K,I3)**2
          VVX4 = VVX4 + XDOT4(K,I4)**2
 875  CONTINUE
*
      VCMBIN2 = C12*VCMBIN2/(M(I1)+M(I2))
      VVX1 = C12*M(I1)*VVX1
      VVX2 = C12*M(I2)*VVX2
      VVX3 = C12*M(I3)*VVX3
      VVX4 = C12*M(I4)*VVX4
*
      RM13 = DSQRT(R2(I1,I3))
      RM14 = DSQRT(R2(I1,I4))
      RM32 = DSQRT(R2(I3,I2))
      RM24 = DSQRT(R2(I2,I4))
*
      EGX12 = -M(I1)*M(I2)/RM1
      EGX13 = -M(I1)*M(I3)/RM13
      EGX14 = -M(I1)*M(I4)/RM14
      EGX23 = -M(I2)*M(I3)/RM32
      EGX24 = -M(I2)*M(I4)/RM24
      EGX34 = -M(I3)*M(I4)/RM2
*
      FX12 = -EGX12/RM1
      FX13 = -EGX13/RM13
      FX14 = -EGX14/RM14
      FX23 = -EGX23/RM32
      FX24 = -EGX24/RM24
      FX34 = -EGX34/RM2
*
      ETTX = VVX1 + VVX2 + VVX3 + VVX4 +
     *       EGX12 + EGX13 + EGX14 + EGX23 + EGX24 + EGX34
      ETTX1 = ETTX - EB(IB1)
      ETTX2 = ETTX1 - EB(IB2)
*
*     PRINT*,' EKIN1-4=',VVX1,VVX2,VVX3,VVX4
*     PRINT*,' EGX12,13,14,23,24,34=',EGX12,EGX13,EGX14,EGX23,
*    *   EGX24,EGX34
*
*       Determine new c.m. energy by subtracting Delta Eb (note Eb>0).
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
         ECMNEW = ECMOL1 + ECMOL2 + EB(IB1) - EBOLD1 - EBOLD2
      ELSE
         ECMNEW = ECMOL1 + ECMOL2 + EB(IB1) - EBOLD1 + EB(IB2) - EBOLD2
      END IF
*
      IDEST = 0
      IF(EB(IB2).LT.0.D0)IDEST = 1
      IF(ECMNEW.LT.0.D0)IDEST = -1
*
      WRITE(6,3003)NAMEB(IB1),NAMEB(IB2),QTIME,EB(IB1),EB(IB2),EBOLD1,
     &   EBOLD2,ECC(IB1),ECC(IB2),ECCOL1,ECCOL2,ECMOL1,ECMOL2,
     &   BODY1(IB1),BODY2(IB1),BODY1(IB2),BODY2(IB2),
     &   VRELT,RGRAV,RGRAV0,RGRAV1,PP,IDEST,QRMIN,XINOL1,XINOL2,
     &   OMGOL1,OMGOL2,XINCL1,XINCL2,XOMEG1,XOMEG2
 3003 FORMAT(' Cross 4b N;T;Eb(1,2)n,o;ecc(1,2)n,o;Ecmo=',2I6,1P,
     & 11D13.5,' M1-4,Vrel,Rgr1-3,p,icase,2angles(old,new)=',9D13.5,
     & I3,9D13.5)
      CALL FLUSH(6)
*     WRITE(6,1001)QTIME,FX12,FX13,FX14,FX23,FX24,FX34,IDEST,
*    *    ETTX,ETTX1,ETTX2,IDEST
*1001 FORMAT(' T,Forces after QUAD=',1P,7(D12.5,2X),I2,/,
*    *       ' ETOT1-3 after QUAD=',3(D12.5,1X),I2)
*       
      IF(ECMNEW.LT.0.D0)THEN
      PRINT*,' Warning 4b possibly bound ECMNEW/kT1,ECMNEW=',
     *       ECMNEW/EKT1,ECMNEW
      PRINT*,' Warning 4b possibly bound ECMNEW/kT2,ECMNEW=',
     *       ECMNEW/EKT2,ECMNEW
      PRINT*,' Warning 4b assumed zero c.m. energy in 4b frame '
*
      IF(ECMNEW.LT.-0.1D0*EKT1)THEN
            ISUBS = ISUBS + 1
            PRINT*,' 4b Subsystem more than 0.1kT, ISUBS=',ISUBS
            WRITE(6,2000)((X4(K,J),K=1,3),(XDOT4(K,J),K=1,3),J=1,4)
 2000       FORMAT(1P,4(1X,6(D12.5,2X),/))
*       Evaluate orbital elements of possible outer binary
      VREL2 = 0.0D0
      RDOT = 0.0D0
      RMCM3 = 0.0D0
      RMCM4 = 0.0D0 
*       Check for separations to 3rd and 4th body
      DO 84 K = 1,3
         RMCM3 = RMCM3 + (XC1(K) - X4(K,I3))**2
         RMCM4 = RMCM4 + (XC1(K) - X4(K,I4))**2
 84   CONTINUE
*
      RMCM3 = DSQRT(RMCM3)
      RMCM4 = DSQRT(RMCM4)   
*
      EGCM3 = -XMBIN1*M(I3)/RMCM3
      EGCM4 = -XMBIN1*M(I4)/RMCM4
*
      FXCM3 = -EGCM3/RMCM3
      FXCM4 = -EGCM4/RMCM4
*       Select third body with larger force on inner binary
      IF(FXCM3.GT.FXCM4)THEN
         IMAX = I3
         RMCM = RMCM3
         XMTRIP = XMBIN1 + M(I3)
      ELSE
         IMAX = I4
         RMCM = RMCM4
         XMTRIP = XMBIN1 + M(I4)
      END IF
*
      DO 85 K = 1,3
         RDOT = RDOT +
     &               (XC1(K) - X4(K,IMAX))*(XC1DOT(K) - XDOT4(K,IMAX))
         VREL2 = VREL2 + (XC1DOT(K) - XDOT4(K,IMAX))**2
   85 CONTINUE
*
*       Determine semi-major axis & eccentricity of outer binary.
      SEMITR = 2.0D0/RMCM - VREL2/XMTRIP
      SEMITR = 1.0/SEMITR
*       Save new binary data for gaseous model.
      EBTRIP = C12*M(IMAX)*XMBIN1/SEMITR
*
      E2 = (1.0D0 - RMCM/SEMITR)**2 + RDOT**2/(SEMITR*XMTRIP)
      IF(E2.GT.0.D0)THEN
      ECCTRP = SQRT(E2)
      ELSE
      ECCTRP = 1.D30
      END IF
*
      PRINT*,' 4b H-Triple? T,IMAX,EB/kT,e,a=',QTIME,IMAX,EBTRIP/EKT1,
     *   ECCTRP,SEMITR
*
*       Check stability criterion of Mardling & Aarseth 2000
*
      IF(ECCTRP.LT.1.D0)THEN
      XCRIT = 2.8D0*SEMIA(IB1)*((1.D0+M(IMAX)/XMBIN1)*(1.D0+ECCTRP)/
     *    DSQRT(1.D0-ECCTRP))**4.D-1
*       Pericentre Distance of outer binary
      YCRIT = SEMITR*(1.D0-ECCTRP)
      PRINT*,' 4b H-Triple? T,Rpcrit,Xcrit=',QTIME,YCRIT,XCRIT
      END IF
*       Determine c.m. of new triple
      DO 86 K=1,3
      XC3(K)=(XMBIN1*XC1(K)+M(IMAX)*X4(K,IMAX))/XMTRIP
      XC3DOT(K)=(XMBIN1*XC1DOT(K)+M(IMAX)*XDOT4(K,IMAX))/XMTRIP
 86   CONTINUE
*
      ECMTRP = 0.D0
      ESNGL4 = 0.D0
      IF(I3.EQ.IMAX)THEN
        ISNGL = I4
      ELSE
        ISNGL = I3
      END IF
*
      RTR4 = 0.D0
*
      DO 87 K=1,3                        
      RTR4 = RTR4 + (XC3(K) - X4(K,ISNGL))**2
      ECMTRP = ECMTRP + XC3DOT(K)**2
      ESNGL4 = ESNGL4 + XDOT4(K,ISNGL)**2
 87   CONTINUE
*
      RTR4 = DSQRT(RTR4)
      EGRTR4 = -XMTRIP*M(ISNGL)/RTR4
*
      ECMTRP = C12*XMTRIP*ECMTRP
      ESNGL4 = C12*M(ISNGL)*ESNGL4
      ECMNEWX = ECMTRP + ESNGL4
      VBIN1X = DSQRT(2.D0*XMTRIP/XMQUAD/M(ISNGL)*ECMNEWX)
      VBIN2X = DSQRT(2.D0*M(ISNGL)/XMQUAD/XMTRIP*ECMNEWX)   
*
      PRINT*,' 4b H-Triple T,EGRTR4/kT,ECMTRP/kT,ESNGL4/kT=',
     *    QTIME,EGRTR4/EKT1,ECMTRP/EKT1,ESNGL4/EKT1
      PRINT*,' 4b H-Triple T,VBIN1,2X=',QTIME,VBIN1X,-VBIN2X
         END IF
*        For correct energy balance add binding energy of triple to that
*        of stronger binary, and assume single star reaches just infinity
      EB(IB1) = EB(IB1) - ECMNEW  
      SEMIA(IB1) = C12*M(I1)*M(I2)/EB(IB1)
      ECMNEW = 0.D0
      END IF
*                              
*       Determine modulus of sing/bin velocity at infinity from energy balance.
      VBIN1 = DSQRT(2.D0*XMBIN2/XMQUAD/XMBIN1*ECMNEW)
      VBIN2 = DSQRT(2.D0*XMBIN1/XMQUAD/XMBIN2*ECMNEW)
*
*     WRITE(6,4001)TIME4,ECMOLD,ECMNEW,VBIN1,-VBIN2
*4001 FORMAT(1X,' End with TIME4,ECMOLD,NEW,vels=',1P,5(D12.5,2X))
*
      q1 =ranf() - 0.5
      costeta = -2.d0*q1
      sinteta = sqrt(1.d0 - costeta**2)
*       Determine vector velocity by two random angles, add 3b c.m. motion.
      dv = vbin1
      q2 = 2.d0*onepi*ranf()
      sinffi = sin(q2)
      cosffi = cos(q2)
      vr(ib1) = vcmtx + dv*costeta
      vtsy = vcmty + dv*sinteta*cosffi
      vtsz = vcmtz + dv*sinteta*sinffi
      vt(ib1) = dsqrt(vtsy**2 + vtsz**2)
*       single and binary move in opposite direction here.
      costeta = -costeta
      sinffi = -sinffi
      cosffi = -cosffi
      dv = vbin2
      vr(ib2) = vcmtx + dv*costeta
      vtsy = vcmty + dv*sinteta*cosffi
      vtsz = vcmtz + dv*sinteta*sinffi
      vt(ib2) = dsqrt(vtsy**2 + vtsz**2)
*
*     print*,' vbin1,vbin2,vcm=',dsqrt(vr(ib1)**2+vt(ib1)**2),
*    *          dsqrt(vr(ib2)**2+vt(ib2)**2),
*    *          dsqrt(vcmtx**2+vcmty**2+vcmtz**2)
*
      DIST = 0.D0
      VDOT = 0.D0
      DO 402 K=1,3
      XREL(K) = XC1(K)-XC2(K)
      VREL(K) = XC1DOT(K)-XC2DOT(K)
      DIST = DIST + XREL(K)*XREL(K)
      VDOT = VDOT + VREL(K)*VREL(K)
 402  CONTINUE
      DIST=SQRT(DIST)
      VDOT=SQRT(VDOT)
*
      EPOTNE = -GRAV*XMBIN1*XMBIN2/DIST
      ERELNE=XMBIN1*XMBIN2/XMQUAD*VDOT**2/2.D0
      EBINA1=XMBIN1*(XC1DOT(1)**2+XC1DOT(2)**2+XC1DOT(3)**2)/2.D0
      EBINA2=XMBIN2*(XC2DOT(1)**2+XC2DOT(2)**2+XC2DOT(3)**2)/2.D0
*     PRINT*,' ECMTOT, EREOLOL, Sum+EPOT=',ECMTOT,ERELNE,
*    *   ECMTOT+ERELNE+EPOTNE
*     PRINT*,' ESNGL3,EBINA3,Sum+ECMTOT+EPOT=',ESNGL3,EBINA3,
*    *   ECMTOT+ESNGL3+EBINA3+EPOTNE
*
      EBNEW1 = EB(IB1)
      EBNEW2 = EB(IB2)
*
*      Cumulative DED/EB analysis if binary not dissolved
*
      IF(EBNEW1.GT.0.D0)THEN            
      ITRM4 = ITRM4 + 1
      DEB4 = DBLE(ITRM4-1)*DEB4/DBLE(ITRM4) + 
     *       (EBNEW1-EBOLD1)/EBOLD1/DBLE(ITRM4)
      END IF
*
      IF(EBNEW2.GT.0.D0)THEN
      ITRM4 = ITRM4 + 1
      DEB4 = DBLE(ITRM4-1)*DEB4/DBLE(ITRM4) +
     *       (EBNEW2-EBOLD2)/EBOLD2/DBLE(ITRM4)
      END IF        
*
      WRITE(6,800)QTIME,DEB4,(EBNEW1-EBOLD1)/EBOLD1,
     *    (EBNEW2-EBOLD2)/EBOLD2,itrm4
 800  FORMAT(1X,1P,' T= ',D12.5,' 4b cum= ',D12.5,' DEB1/EB1= ',D12.5,
     *  ' DEB2/EB2= ',D12.5,' itrm4=',I8)    
*
      RETURN
*
      END
