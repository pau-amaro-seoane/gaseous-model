         SUBROUTINE TRMAIN(IBIN,IC,VR1,VT1,XM,XR,EKT,TTIME)
*
*        Interface for triple integration from gaseous model program
*        (from S. Aarseth NBODY6 Program)
*
*        IBIN Number of Binary
*        IC   Field star mass group index
*        VR1, VT1 Radial and tangential velocity of field star
*        XM,XR Field star mass and radius
*
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  M
      REAL*8  XORB(2),VORB(2),XREL(3),VREL(3),PX(3),QX(3),EI(3),XLI(3)
      REAL*8  XCM(3),XCMDOT(3)
      COMMON/AZREG/  TIME3,TMAX,Q(8),P(8),R1,R2,R3,ENERGY,M(3),X3(3,3),
     &               XDOT3(3,3),CM(10),C11,C12A,C19,C20,C24,C25,
     &               NSTEP3,NAME3(3),KZ15,KZ27
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3
      COMMON/AZCOLL/  RK(3),QK(8),PK(8),ICALL,ICOLL,NDISS3
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      PARAMETER(NTAB=32)
      INTEGER IDUM2,K,IV,IY
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
      INTEGER ITRM3
      SAVE DEB3,ITRM3
      DATA ITRM3/0/
      DATA  ISUBS /0/ 
      SAVE
*
      REAL RANF
      EXTERNAL RANF
C     Common block BINARY
C     contains binary data for stochastic 3-B binaries
      INCLUDE 'params.f'
C     Common block BINARY
C     contains binary data for stochastic 3-B binaries
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
      DATA  EP  /0.04D0,0.0016D0,0.64D-4,0.256D-5/
*       Field star 1 meets Binary 2,3
*       Transformation to c.m. of the triple system necessary
*
      BMRAT = MIN(BODY1(IBIN),BODY2(IBIN))/
     &              MAX(BODY1(IBIN),BODY2(IBIN))
      IF(BMRAT.LT.1.D-3)THEN
      BMKT = 0.1D0*BMRAT 
      ELSE
      BMKT = 1.D0
      BMRAT = 1.D0
      END IF
*
      PRINT*,' Binary ',IBIN,' N=',NAMEB(IBIN),
     *    ' invoked TRMAIN with field star '
      CALL FLUSH(6)
*     PRINT*,' VR,VT=',VR1,VT1
      ESNGOL=XM*(VR1**2+VT1**2)/2.D0
      XMBIN=BODY1(IBIN)+BODY2(IBIN)
      XMTRIP = XMBIN + XM
      EBINOL=XMBIN*(VR(IBIN)**2+VT(IBIN)**2)/2.D0
      EBOLD=EB(IBIN)
      ECCOLD=ECC(IBIN)
      XINOLD=XIN(IBIN)
      OMEGOL=OMEG(IBIN)
*
      WRITE(6,888)TTIME,ESNGOL,EBINOL,EBOLD
      WRITE(6,889)TTIME,EKT,EBOLD,EBOLD/BMRAT/EKT,SEMIA(IBIN),ECC(IBIN)
 888  FORMAT(1X,' Start with T,ECM1,2,EB=',1P,4D13.5)
 889  FORMAT(1X,' Start T,kT,EB,EB/BMRAT/kT,a,e=',1P,6D13.5)
*
      GRAV = 1.D0
      TWOPI = 2.D0*ONEPI
*
 7891 CONTINUE
      FI = TWOPI*RANF()
      SINFI = SIN(FI)
      COSFI = COS(FI)
*
*    calculate Pmax using formula given by Hut and Bahcal ApJ 268, 1983
*
      VRELX = VR1-VR(IBIN)
      VRELY = VT1*SINFI
      VRELZ = VT1*COSFI-VT(IBIN)
      VRELT2 = VRELX**2+VRELY**2+VRELZ**2
      VRELT = DSQRT(VRELT2)
*
      VC2RAW = BODY1(IBIN)*BODY2(IBIN)*XMTRIP/(XM*XMBIN*SEMIA(IBIN))
      VC2 = VC2RAW/BMRAT
*   take care for rare case we come here with an unbound ''binary'' (planet)
      IF(VC2.LT.0.D0)THEN
      VVC2 = 1.D0
      ELSE
      VVC2 = DSQRT(VRELT2/VC2)
      END IF
*
      IF(BMRAT.LT.1.D-3)THEN 
*    Use 100*a*(1+e) for planetary system run
         PMAX = 100.D0*SEMIA(IBIN)*(1.0D0 + ECC(IBIN))
      ELSE
         PMAX = SEMIA(IBIN)*(0.6D0*(1.0D0 + ECC(IBIN)) + 4.0D0/VVC2)
      END IF
*
*    New to speed up the program
*      IF(PMAX.GT.10.D0*SEMIA(IBIN)) PMAX = 10.D0*SEMIA(IBIN)
*
      PMAX2 = PMAX*PMAX
*
*     Calculate minimum distance RGRAVT for given PMAX
*
      RXXX = GRAV*XMTRIP/VRELT2
      RGRAVT = RXXX*(DSQRT(1.D0 + (PMAX/RXXX)**2) -1.D0)
      RGRAV1 = RGRAVT
      RGRAV0 = RGRAVT
*
      WRITE(6,3459) XMBIN,XM,XMTRIP,VRELT2,VC2,VVC2,RXXX,BMRAT,VC2RAW
 3459 FORMAT(1X,'XMBIN,XM,XMTRIP,VRELT2,VC2,VVC2,RXXX,BMRAT,VC2RAW= ',
     &   1P9E12.4)
*
*      RGRAV0 = 2.5D0*SEMIA(IBIN)
*      RGRAVT = RGRAV0
*       Take out because confusing
*      RGRAV1 = (EB(IBIN)/EKT)**C13*SEMIA(IBIN)
*      IF(RGRAV1.GT.RGRAV0)RGRAVT=RGRAV1
*
      VCMTX = (XM*VR1 + XMBIN*VR(IBIN))/XMTRIP
      VCMTY = XM*VT1*SINFI/XMTRIP
      VCMTZ = (XM*VT1*COSFI + XMBIN*VT(IBIN))/XMTRIP
      ECMTOT = XMTRIP*(VCMTX**2+VCMTY**2+VCMTZ**2)/2.D0
*
      VCM1X = -XM/XMTRIP*VRELX
      VCM1Y = -XM/XMTRIP*VRELY
      VCM1Z = -XM/XMTRIP*VRELZ
      VCM2X = XMBIN/XMTRIP*VRELX
      VCM2Y = XMBIN/XMTRIP*VRELY
      VCM2Z = XMBIN/XMTRIP*VRELZ
*
*       Save kinetic energy of binary and single in 3b c.m. frame 
*       before encounter (at infinity).
      ECMOL1 = C12*XMBIN*(VCM1X**2+VCM1Y**2+VCM1Z**2)
      ECMOL2 = C12*XM*(VCM2X**2+VCM2Y**2+VCM2Z**2)
      ECMOLD = ECMOL1 + ECMOL2    
*      Determine max impact parameter for rmin = 2 rgrav
*      PMAX2 = RGRAVT*(RGRAVT + 2.D0*GRAV*(XMBIN+XM)/
*     *   VRELT/VRELT)
*
      KZ15 = 2
      KZ27 = 2
      NAME3(1) = 1
      NAME3(2) = 2
      NAME3(3) = 3
*
      M(1) = XM
      SIZE(1) = 0.D0
      IF (BODY1(IBIN).LT.BODY2(IBIN)) THEN
          M(2) = BODY1(IBIN)
          M(3) = BODY2(IBIN)
          SIZE(2) = 0.D0
          SIZE(3) = 0.D0
      ELSE
          M(2) = BODY2(IBIN)
          M(3) = BODY1(IBIN)
          SIZE(2) = 0.D0
          SIZE(3) = 0.D0
      END IF
*
      DO 2 K = 1,7
          CM(K) = 0.0D0
    2 CONTINUE
*
      PP2 = RANF()*PMAX2
      PP = DSQRT(PP2)
*
*      Select start distance for triple for approximate force perturbation
*      of 1.e-3 and assume worst case of apocentre directed to single star.
*      WARNING! Check this carefully for multi-mass case!! Definitely wrong!!
*      XMUE = M(2)*M(3)/XMBIN
*      RX = 1.D1*SEMIA(IBIN)*(1.D0+ECC(IBIN))*(1.D0+1.D1*(XM/XMUE)**C13)
*
*     Select distance to starr triple for approximate force perturbation
*     gama_crit=1.e-3 and assume the worst case of apocentre directed to single
*     star. Determination of Rx based on ksinit.f from Nbody. Modifications:
*     1) for m3=0 Rx=A*a(1+e)
*     2) A is choosen to 10
*
      XMUE = 8.D0*XM/XMBIN      
      RX = 1.D1*SEMIA(IBIN)*(1.D0+ECC(IBIN))*(1.D0+1.D1*(XMUE)**C13)
*
*     determine the minimum distance for actual impact parameter pp
*
      RGRAVX = RXXX*(DSQRT(1.D0 + (PP/RXXX)**2) - 1.D0)  
*      Take care if RX based on perturbation gives too small value
      IF(RX.LT.2.D0*RGRAVX) RX = 2.D0*RGRAVX
* 
      WRITE(6,2003)TTIME,RX,RGRAV0,RGRAV1,PP,DSQRT(PMAX2)
 2003 FORMAT(1X,' Start with T,RX=',1P,2D13.5,' RGRAV0,X=',
     * 2D13.5,' PP,PMAX=',2D13.5)
*
      VRELT2 = VRELT*VRELT
      ECC2 = 1.D0 + (PP*VRELT2/(GRAV*XMTRIP))**2
      ECCT = DSQRT(ECC2)
      XJ = PP*VRELT
*     PRINT*,' ECCT=',ECCT,' PP=',PP,' PMAX=',DSQRT(PMAX2),
*    * ' RX=',RX,' RGRAVT=',RGRAVT
*     PRINT*,' Arg Theta=',(XJ**2/(RX*GRAV*XMTRIP)-1.D0),
*    * (XJ**2/(RX*GRAV*XMTRIP))/ECCT-1.D0/ECCT
      THETAX = DACOS((XJ**2/(RX*GRAV*XMTRIP)-1.D0)/ECCT)
      PSI = DATAN(PP*VRELT2/GRAV/XMTRIP)
      ALPHA = ONEPI - THETAX - PSI
      DEGR1 = ONEPI/180.
*     PRINT*,' ALPHA, THETAX,PSI,PI=',ALPHA/DEGR1,THETAX/DEGR1,
*    * PSI/DEGR1,ONEPI/DEGR1,VRELT
*
      VTX = PP*VRELT/RX
      VRX = -DSQRT(VRELT2 - VTX**2 + 2.D0*GRAV*XMTRIP/RX)
*
*       Velocity Components of field star

      XRELP = RX*DCOS(ALPHA)
      YRELP = RX*DSIN(ALPHA)
*
      XVRELP = (VRX*DCOS(ALPHA) - VTX*DSIN(ALPHA))
      YVRELP = (VRX*DSIN(ALPHA) + VTX*DCOS(ALPHA))
*
*     PRINT*,' XRELP,YRELP=',XRELP,YRELP,' XVRELP,YVRELP=',
*    *    XVRELP,YVRELP,' VRX,VTX=',VRX,VTX
*
      EPOT2 = -GRAV*XM*XMBIN/RX
*     PRINT*,' ***XXXXXX E CHECK****'
      VVREL2=XVRELP**2 + YVRELP**2
      ERELOL=XM*XMBIN/XMTRIP*VVREL2/2.D0
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
      X3(1,1) = XMBIN/XMTRIP*XXREL
      X3(2,1) = XMBIN/XMTRIP*YREL
      X3(3,1) = XMBIN/XMTRIP*ZREL
*
      XDOT3(1,1) = XMBIN/XMTRIP*XVREL
      XDOT3(2,1) = XMBIN/XMTRIP*YVREL
      XDOT3(3,1) = XMBIN/XMTRIP*ZVREL
*
*     PRINT*,' Before TRIPLE CALL    ',-VRELT,0.D0
*     PRINT*,' Velocity Sin=',XDOT3(1,1),XDOT3(2,1),XDOT3(3,1)
*     PRINT*,' Position Sin=',X3(1,1),X3(2,1),X3(3,1)
*       Velocity Components of binary c.m.
*
      X3(1,2) = -XM/XMTRIP*XXREL
      X3(2,2) = -XM/XMTRIP*YREL
      X3(3,2) = -XM/XMTRIP*ZREL
*
      XDOT3(1,2) = -XM/XMTRIP*XVREL
      XDOT3(2,2) = -XM/XMTRIP*YVREL
      XDOT3(3,2) = -XM/XMTRIP*ZVREL
*
*     PRINT*,' Velocity Bin=',XDOT3(1,2),XDOT3(2,2),XDOT3(3,2)
*     PRINT*,' Position Bin=',X3(1,2),X3(2,2),X3(3,2)
*
*     PRINT*,' ***SECOND E CHECK****'
      VVREL2=(XDOT3(1,1)-XDOT3(1,2))**2+
     *       (XDOT3(2,1)-XDOT3(2,2))**2+
     *       (XDOT3(3,1)-XDOT3(3,2))**2
      ERELOL=XM*XMBIN/XMTRIP*VVREL2/2.D0
      ESNGL2=XM*(XDOT3(1,1)**2+XDOT3(2,1)**2+XDOT3(3,1)**2)/2.D0
      EBINA2=XMBIN*(XDOT3(1,2)**2+XDOT3(2,2)**2+XDOT3(3,2)**2)/2.D0
*     PRINT*,' ECMTOT, ERELOL, Sum=',ECMTOT,ERELOL,
*    *   ECMTOT+ERELOL+EPOT2
*     PRINT*,' ESNGL2,EBINA2,Sum+ECMTOT+EPOT=',ESNGL2,EBINA2,
*    *   ECMTOT+ESNGL2+EBINA2+EPOT2
*
*       Randomize perihelion, node & inclination. (from binpop, NBODY6)
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
          XORB(1) = SEMIA(IBIN)*(1.0 + ECC(IBIN))
          XORB(2) = 0.0
          VORB(1) = 0.0
          VORB(2) = SQRT(XMBIN*(1.0D0 - ECC(IBIN))/
     *                         (SEMIA(IBIN)*(1.0D0 + ECC(IBIN))))
*
*       Transform to relative variables.
          DO 40 K = 1,3
              XREL(K) = PX(K)*XORB(1) + QX(K)*XORB(2)
              VREL(K) = PX(K)*VORB(1) + QX(K)*VORB(2)
   40     CONTINUE
*
          I1 = 2
          I2 = 3
*
*       Set global variables for each component.
          DO 50 K = 1,3
              X3(K,I1) = X3(K,I1) + M(I2)*XREL(K)/XMBIN
              X3(K,I2) = X3(K,I1) - XREL(K)
              XDOT3(K,I1) = XDOT3(K,I1) + M(I2)*VREL(K)/XMBIN
              XDOT3(K,I2) = XDOT3(K,I1) - VREL(K)
   50     CONTINUE
*
*     PRINT*,' Masses =',M(1),M(2),M(3)
*     PRINT*,' Velocity 1=',XDOT3(1,1),XDOT3(2,1),XDOT3(3,1)
*     PRINT*,' Velocity 2=',XDOT3(1,2),XDOT3(2,2),XDOT3(3,2)
*     PRINT*,' Velocity 3=',XDOT3(1,3),XDOT3(2,3),XDOT3(3,3)
*     PRINT*,' Position 1=',X3(1,1),X3(2,1),X3(3,1)
*     PRINT*,' Position 2=',X3(1,2),X3(2,2),X3(3,2)
*     PRINT*,' Position 3=',X3(1,3),X3(2,3),X3(3,3)
*
*
*       Transform to the local c.m. reference frame.
      DO 4 L = 1,3
          CM(7) = CM(7) + M(L)
          DO 3 K = 1,3
              CM(K) = CM(K) + M(L)*X3(K,L)
              CM(K+3) = CM(K+3) + M(L)*XDOT3(K,L)
    3     CONTINUE
    4 CONTINUE
*
*       Set c.m. coordinates & velocities for triple system.
      DO 5 K = 1,6
          CM(K) = CM(K)/CM(7)
    5 CONTINUE
*
*     PRINT*,' CM: ',(CM(K),K=1,6)
*     PRINT*,' Now Call TRIPLE with Coordinates:'
*     PRINT*,' Particle X=',X3(1,1),X3(2,1),X3(3,1)
*     PRINT*,' Binary 1X=',X3(1,2),X3(2,2),X3(3,2)
*     PRINT*,' Binary 2X=',X3(1,3),X3(2,3),X3(3,3)
*     PRINT*,' Velocities:'
*     PRINT*,' Particle XDOT=',XDOT3(1,1),XDOT3(2,1),XDOT3(3,1)
*     PRINT*,' Binary 1XDOT=',XDOT3(1,2),XDOT3(2,2),XDOT3(3,2)
*     PRINT*,' Binary 2XDOT=',XDOT3(1,3),XDOT3(2,3),XDOT3(3,3)
*
*     DIST = 0.D0
*     VDOT = 0.D0
*     DO 400 K=1,3
*     XCM(K)=(M(2)*X3(K,2)+M(3)*X3(K,3))/(M(2)+M(3))
*     XCMDOT(K)=(M(2)*XDOT3(K,2)+M(3)*XDOT3(K,3))/(M(2)+M(3))
*     XREL(K) = X3(K,1)-XCM(K)
*     VREL(K) = XDOT3(K,1)-XCMDOT(K)
*     DIST = DIST + XREL(K)*XREL(K)
*     VDOT = VDOT + VREL(K)*VREL(K)
*400  CONTINUE
*     DIST=SQRT(DIST)
*     VDOT=SQRT(VDOT)
*
      VVX1 = XDOT3(1,1)**2 + XDOT3(2,1)**2 + XDOT3(3,1)**2
      VVX2 = XDOT3(1,2)**2 + XDOT3(2,2)**2 + XDOT3(3,2)**2 
      VVX3 = XDOT3(1,3)**2 + XDOT3(2,3)**2 + XDOT3(3,3)**2 
*
      DBIN2 = (X3(1,2)-X3(1,3))**2 + (X3(2,2)-X3(2,3))**2 +
     *        (X3(3,2)-X3(3,3))**2
      RKJ2 = (X3(1,2)-X3(1,1))**2 + (X3(2,2)-X3(2,1))**2 +
     *       (X3(3,2)-X3(3,1))**2 
      RLJ2 = (X3(1,3)-X3(1,1))**2 + (X3(2,3)-X3(2,1))**2 +
     *       (X3(3,3)-X3(3,1))**2
      FK1 = -M(1)*((X3(1,2)-X3(1,1))/RKJ2**1.5D0 -
     *             (X3(1,3)-X3(1,1))/RLJ2**1.5D0)
      FK2 = -M(1)*((X3(2,2)-X3(2,1))/RKJ2**1.5D0 -
     *             (X3(2,3)-X3(2,1))/RLJ2**1.5D0)
      FK3 = -M(1)*((X3(3,2)-X3(3,1))/RKJ2**1.5D0 -
     *             (X3(3,3)-X3(3,1))/RLJ2**1.5D0)
      FFK = DSQRT(FK1**2+FK2**2+FK3**2)
      GAMMA = FFK * DBIN2/XMBIN
*     PRINT*,' Binary CM=',(XCM(K),K=1,3)
*     PRINT*,' Binary CMDOT=',(XCMDOT(K),K=1,3)
*     PRINT*,' Binary Single R,VREL,GAMMA=',DIST,VDOT,GAMMA
*
      RM1 = DSQRT(DBIN2)
      RM2 = DSQRT(RKJ2)
      RM3 = DSQRT(RLJ2)
*
      EGX23 = - M(2)*M(3)/RM1
      EGX12 = - M(1)*M(2)/RM2  
      EGX13 = - M(1)*M(3)/RM3  
*
      FX23 = -EGX23/RM1
      FX12 = -EGX12/RM2      
      FX13 = -EGX13/RM3      
*
      ETTX = VVX1 + VVX2 + VVX3 + EGX23 + EGX12 + EGX12
      ETTX1 = ETTX - EB(IBIN)
*
*
*     WRITE(6,1000)TTIME,FX23,FX12,FX13,ETTX,ETTX1
*1000 FORMAT(' T,Forces before TRIPLE=',1P,4(D12.5,2X),/,
*    *       ' ETOT1-2 before TRIPLE=',2(D12.5,1X))
*
      CALL TRANS3(0)
*
      I = 0
*
      NAMB = NAMEB(IBIN)
*
      CALL FLUSH(6)
*
      CALL TRIPLE(I,RGRAVT,NAMB,TTIME,TRMIN,EKT)
*       Rgrav saved from triple
      RGRAV = RGRAVT
*
      IMIN = 1
      IF(R2.LT.R1) IMIN = 2
      IMAX = 3 - IMIN
*
      RM = MIN(R1,R2)
      XMBIN = M(IMIN) + M(3)
      XM = M(IMAX)
      XMTRIP = XM + XMBIN       
*
*       Save binary c.m. data for later use (in system where triple rests)
      DO 74 K=1,3
      XCM(K)=(M(IMIN)*X3(K,IMIN)+M(3)*X3(K,3))/(M(IMIN)+M(3))
      XCMDOT(K)=(M(IMIN)*XDOT3(K,IMIN)+M(3)*XDOT3(K,3))/(M(IMIN)+M(3))
 74   CONTINUE      
*
*       Evaluate orbital elements.
      VREL2 = 0.0D0
      RDOT = 0.0D0
      RI2 = 0.0D0
*           
      DO 75 K = 1,3
         XREL(K) = X3(K,3) - X3(K,IMIN)
         VREL(K) = XDOT3(K,3) - XDOT3(K,IMIN)
         RI2 = RI2 + XREL(K)**2
         RDOT = RDOT + XREL(K)*VREL(K)
         VREL2 = VREL2 + VREL(K)**2
   75 CONTINUE
*
*       Construct Runge-Lenz vector (Heggie & Rasio, 1995, IAU174, Eq.(5)).
          EI2 = 0.0D0
          XLI2 = 0.0D0
          CHECK = 0.0D0
          DO 76 K = 1,3
              KP1 = 1 + MOD(K,3)
              KP2 = 1 + MOD(K+1,3)
              EI(K) = (VREL2*XREL(K) - RDOT*VREL(K))/XMBIN -
     &                                        XREL(K)/DSQRT(RI2)
              XLI(K) = XREL(KP1)*VREL(KP2)-XREL(KP2)*VREL(KP1)
              EI2 = EI2 + EI(K)**2
              XLI2 = XLI2 + XLI(K)**2
              CHECK = CHECK + EI(K)*XLI(K)
   76     CONTINUE
*
          XIARG = DSQRT((XLI(1)**2+XLI(2)**2)/XLI2)
          IF (XLI(3).GT.0.D0) THEN
              XINCL = ASIN(XIARG)
              ELSE
              XINCL = TWOPI/2.D0 - ASIN(XIARG)
          END IF
*
          XOARG = EI(2)/DSQRT(EI(1)**2+EI(2)**2)
          IF (EI(1).GT.0.D0) THEN
          IF (EI(2).GT.0.D0) XOMEG = ASIN(XOARG)
          IF (EI(2).LT.0.D0) XOMEG = TWOPI - ASIN(-XOARG)
          ELSE
          IF (EI(2).LT.0.D0) XOMEG = TWOPI/2.D0 + ASIN(-XOARG)
          IF (EI(2).GT.0.D0) XOMEG = TWOPI/2.D0 - ASIN(XOARG)
          END IF
*
      XIN(IBIN) = XINCL
      OMEG(IBIN) = XOMEG
*       Determine semi-major axis & eccentricity of inner binary.
      SEMI = 2.0D0/RM - VREL2/XMBIN
      SEMI = 1.0/SEMI
*       Save new binary data for gaseous model.
      EB(IBIN) = 0.5D0*M(IMIN)*M(3)/SEMI      
*
      E2 = (1.0D0 - RM/SEMI)**2 + RDOT**2/(SEMI*XMBIN)
      IF(E2.GT.0.D0)THEN
      ECC(IBIN) = SQRT(E2)
      ELSE
      ECC(IBIN) = 1.D30
      END IF
*
      SEMIA(IBIN) = SEMI
      BODY1(IBIN) = M(IMIN)
      BODY2(IBIN) = M(3)
*     SIZE1(IBIN) = SIZE(IMIN)
*     SIZE2(IBIN) = SIZE(3)
*               
      WRITE(6,3001)TTIME,IMIN,RM,SEMI,ECC(IBIN),EB(IBIN),EKT
 3001 FORMAT(1X,' After Triple Bin T,IMIN,RM,SEMI,E,EB,EKT=',
     &   1P,D13.5,I5,5D13.5)
      CALL FLUSH(6)
*
      VVX1 = 0.D0
      VVX2 = 0.D0
      VVX3 = 0.D0
      RM1 = 0.D0
      RM2 = 0.D0
      RM3 = 0.D0
      VCMBIN2 = 0.D0
*
      DO 875 K = 1,3
         VCMBIN2 = VCMBIN2 + (M(IMIN)*XDOT3(K,IMIN)+M(3)*XDOT3(K,3))**2
         VVX1 = VVX1 + XDOT3(K,1)**2
         VVX2 = VVX2 + XDOT3(K,2)**2
         VVX3 = VVX3 + XDOT3(K,3)**2
         RM1 = RM1 + (X3(K,IMIN)-X3(K,3))**2
         RM2 = RM2 + (X3(K,IMAX)-X3(K,3))**2
         RM3 = RM3 + (X3(K,IMIN)-X3(K,IMAX))**2
 875  CONTINUE
*
      VCMBIN2 = C12*VCMBIN2/(M(IMIN)+M(3))
      VVX1 = C12*M(1)*VVX1
      VVX2 = C12*M(2)*VVX2
      VVX3 = C12*M(3)*VVX3
      RM1 = DSQRT(RM1)
      RM2 = DSQRT(RM2)
      RM3 = DSQRT(RM3)
*
      EGX1 = -M(IMIN)*M(3)/RM1
      EGX2 = -M(IMAX)*M(3)/RM2
      EGX3 = -M(IMIN)*M(IMAX)/RM3
*
      FX1 = -EGX1/RM1
      FX2 = -EGX2/RM2
      FX3 = -EGX3/RM3
*
      ETTX = VVX1 + VVX2 + VVX3 + EGX1 + EGX2 + EGX3
      ETTX1 = ETTX - EB(IBIN)
*
*       Determine new c.m. energy by subtracting Delta Eb (note Eb>0).
*       In case of planets (mass ratio < 1.e-3) use v^2 criterion
      BMRAT = MIN(BODY1(IBIN),BODY2(IBIN))/
     &              MAX(BODY1(IBIN),BODY2(IBIN))
      IF(BMRAT.LT.1.D-3)THEN
      BMKT = 0.1D0*BMRAT
      ELSE
      BMKT = 1.D0 
      BMRAT = 1.D0
      END IF
*
      IF(EB(IBIN).LT.0.1D0*BMKT*EKT)THEN
         ECMNEW = ECMOLD + EB(IBIN)
         IDEST = 1
      ELSE
         ECMNEW = ECMOLD + EB(IBIN) - EBOLD
         IDEST = 0
      END IF  
*       Possibly bound triple
      IF(ECMNEW.LT.0)IDEST = -1
*
      WRITE(6,3002)NAMEB(IBIN),TTIME,EB(IBIN),EBOLD,ECC(IBIN),ECCOLD,
     &   ECMOL1,ECMOL2,
     &   BODY1(IBIN),BODY2(IBIN),XM,VRELT,RGRAV,RGRAV0,RGRAV1,PP,IDEST,
     &   TRMIN,XINOLD,OMEGOL,XINCL,XOMEG
 3002 FORMAT(' Cross 3b N;T;Ebn,o;eccn,o;Ecmo=',I6,1P,7D13.5,
     & ' M1-3,Vrel,Rgr1-3,p,icase,trmin,2angles(old,new)=',8D13.5,
     &   I3,5D13.5)
      CALL FLUSH(6)
*     WRITE(6,1001)TTIME,FX1,FX2,FX3,IDEST,ETTX,ETTX1,IDEST
*1001 FORMAT(' T,Forces after TRIPLE=',1P,4(D12.5,2X),I2,/,
*    *       ' ETOT1-2 after TRIPLE=',2(D12.5,1X),I2)
*
      IF(ECMNEW.LT.0.D0)THEN
      PRINT*,' Warning 3b possibly bound subsystem ECMNEW/kT=',
     *  ECMNEW/EKT,ECMNEW
      PRINT*,' Warning 3b assumed zero c.m. energy in 3b frame ' 
         IF(ECMNEW.LT.-0.1D0*EKT)THEN
            ISUBS = ISUBS + 1
            PRINT*,' 3b Subsystem more than 0.1kT, T,ISUBS=',
     *      TTIME,ISUBS
            WRITE(6,2000)((X3(K,J),K=1,3),(XDOT3(K,J),K=1,3),J=1,3)
 2000       FORMAT(1P,3(1X,6(D12.5,2X),/))
*       Evaluate orbital elements of possible outer binary.
      VREL2 = 0.0D0
      RDOT = 0.0D0
      RMCM = 0.0D0
*
      DO 85 K = 1,3
         RMCM = RMCM + (XCM(K) - X3(K,IMAX))**2
         RDOT = RDOT +
     &               (XCM(K) - X3(K,IMAX))*(XCMDOT(K) - XDOT3(K,IMAX))
         VREL2 = VREL2 + (XCMDOT(K) - XDOT3(K,IMAX))**2
   85 CONTINUE
*
      RMCM = DSQRT(RMCM)
*       Determine semi-major axis & eccentricity of outer binary.
      SEMITR = 2.0D0/RMCM - VREL2/XMTRIP
      SEMITR = 1.0/SEMITR
*       Save new binary data for gaseous model.
      EBTRIP = C12*M(IMAX)*XMBIN/SEMITR
*
      ECCTRP = SQRT((1.0D0 - RMCM/SEMITR)**2 + RDOT**2/(SEMITR*XMTRIP))
*
      PRINT*,'3b H-Triple? T,EB/kT,e,a=',TTIME,EBTRIP/EKT,ECCTRP,SEMITR
*
         END IF        
*        For correct energy balance add binding energy of triple to that
*        of remaining binary, and assume single star reaches just infinity
      EB(IBIN) = EB(IBIN) - ECMNEW
      SEMIA(IBIN) = C12*M(IMIN)*M(3)/EB(IBIN)
      ECMNEW = 0.D0
      END IF
*
*       Determine modulus of sing/bin velocity at infinity from energy balance.
      VSNGL = DSQRT(2.D0*XMBIN/XMTRIP/XM*ECMNEW)
      VBINA = DSQRT(2.D0*XM/XMTRIP/XMBIN*ECMNEW)
*
*     PRINT*,' End 3b with TIME3,ECMOLD,NEW,vels=',TIME3,
*    *      ECMOLD,ECMNEW,VSNGL,-VBINA
      q1 =ranf() - 0.5
      costeta = -2.d0*q1
      sinteta = sqrt(1.d0 - costeta**2)
*       Determine vector velocity by two random angles, add 3b c.m. motion.
      dv = vsngl
      q2 = 2.d0*onepi*ranf()
      sinffi = sin(q2)
      cosffi = cos(q2)
      vr1 = vcmtx + dv*costeta
      vtsy = vcmty + dv*sinteta*cosffi
      vtsz = vcmtz + dv*sinteta*sinffi
      vt1 = dsqrt(vtsy**2 + vtsz**2)
*       single and binary move in opposite direction here.
      costeta = -costeta
      sinffi = -sinffi
      cosffi = -cosffi
      dv = vbina
      vr(ibin) = vcmtx + dv*costeta
      vtsy = vcmty + dv*sinteta*cosffi
      vtsz = vcmtz + dv*sinteta*sinffi
      vt(ibin) = dsqrt(vtsy**2 + vtsz**2)
*
*     print*,' vsngl,vbina,vcm=',dsqrt(vr1**2+vt1**2),
*    *          dsqrt(vr(ibin)**2+vt(ibin)**2),
*    *          dsqrt(vcmtx**2+vcmty**2+vcmtz**2)
*
*     PRINT*,' After TRIPLE with Coordinates (Binary 3 ',IMIN,':'
*     PRINT*,' Particle X=',X3(1,IMAX),X3(2,IMAX),X3(3,IMAX)
*     PRINT*,' Binary 1X=',X3(1,IMIN),X3(2,IMIN),X3(3,IMIN)
*     PRINT*,' Binary 2X=',X3(1,3),X3(2,3),X3(3,3)
*     PRINT*,' Velocities:'
*     PRINT*,' Particle XDOT=',XDOT3(1,IMAX),XDOT3(2,IMAX),XDOT3(3,IMAX)
*     PRINT*,' Binary 1XDOT=',XDOT3(1,IMIN),XDOT3(2,IMIN),XDOT3(3,IMIN)
*     PRINT*,' Binary 2XDOT=',XDOT3(1,3),XDOT3(2,3),XDOT3(3,3)
*
      DIST = 0.D0
      VDOT = 0.D0
      DO 402 K=1,3
      XREL(K) = X3(K,IMAX)-XCM(K)
      VREL(K) = XDOT3(K,IMAX)-XCMDOT(K)
      DIST = DIST + XREL(K)*XREL(K)
      VDOT = VDOT + VREL(K)*VREL(K)
 402  CONTINUE
      DIST=SQRT(DIST)
      VDOT=SQRT(VDOT)
*
*     PRINT*,' ***THIRD E CHECK****'
      VVREL2=(XDOT3(1,1)-XCMDOT(1))**2+
     *       (XDOT3(2,1)-XCMDOT(2))**2+
     *       (XDOT3(3,1)-XCMDOT(3))**2
      XREL2 = (X3(1,1)-XCM(1))**2+
     *       (X3(2,1)-XCM(2))**2+
     *       (X3(3,1)-XCM(3))**2
      EPOTNE = -GRAV*XM*XMBIN/DSQRT(XREL2)
      ERELNE=XM*XMBIN/XMTRIP*VVREL2/2.D0
      ESNGL3=XM*(XDOT3(1,1)**2+XDOT3(2,1)**2+XDOT3(3,1)**2)/2.D0
      EBINA3=XMBIN*(XCMDOT(1)**2+XCMDOT(2)**2+XCMDOT(3)**2)/2.D0
*     PRINT*,' ECMTOT, EREOLOL, Sum+EPOT=',ECMTOT,ERELNE,
*    *   ECMTOT+ERELNE+EPOTNE
*     PRINT*,' ESNGL3,EBINA3,Sum+ECMTOT+EPOT=',ESNGL3,EBINA3,
*    *   ECMTOT+ESNGL3+EBINA3+EPOTNE
*
      EBNEW = EB(IBIN)
*
*      Cumulative DED/EB analysis if binary not dissolved
*
      IF(EBNEW.GT.0.D0)THEN
*
      ITRM3 = ITRM3 + 1
      DEB3OL = DEB3
      DEB3 = DBLE(ITRM3-1)*DEB3/DBLE(ITRM3) + 
     *       (EBNEW-EBOLD)/EBOLD/DBLE(ITRM3)
*      END IF
      WRITE(6,800)TTIME,DEB3,(EBNEW-EBOLD)/EBOLD,itrm3
 800  FORMAT(1X,1P,' T= ',D12.5,' 3b cum= ',D12.5,' DEB/EB= ',D12.5,
     *  ' itrm3=',I8)
*
      END IF
*
*     VR(IBIN) = VCMTX + XCMDOT(1)
*     VTYY = VCMTY + XCMDOT(2)
*     VTZZ = VCMTZ + XCMDOT(3)
*     VT(IBIN) = DSQRT(VTYY**2+VTZZ**2)
*     VR1 = VCMTX + XDOT3(1,IMAX)
*     VTYY = VCMTY + XDOT3(2,IMAX)
*     VTZZ = VCMTZ + XDOT3(3,IMAX)
*     VT1 = DSQRT(VTYY**2+VTZZ**2)
*
*     PRINT*,' 1 bin:',- XM/XMTRIP*VREL(1),XCMDOT(1)
*     PRINT*,' 2 bin:',- XM/XMTRIP*VREL(2),XCMDOT(2)
*     PRINT*,' 3 bin:',- XM/XMTRIP*VREL(3),XCMDOT(3)
*     PRINT*,' 1 sng:',+ XMBIN/XMTRIP*VREL(1),XDOT3(1,IMAX)
*     PRINT*,' 2 sng:',+ XMBIN/XMTRIP*VREL(2),XDOT3(2,IMAX)
*     PRINT*,' 3 sng:',+ XMBIN/XMTRIP*VREL(3),XDOT3(3,IMAX)
*
*     PRINT*,' ******FOURTH CHECK*****'
*     ESNGNE = XM*(VR1**2+VT1**2)/2.D0
*     EBINNE = XMBIN*(VR(IBIN)**2+VT(IBIN)**2)/2.D0
*     PRINT*,' ESNGNE,EBINNE,Sum=',ESNGNE,EBINNE,ESNGNE+EBINNE+EPOTNE
*     PRINT*,' Binary Single VR,VT=',VR(IBIN),VT(IBIN),VR1,VT1
*     PRINT*,' Binary Single R,VREL=',DIST,VDOT
*     PRINT*,' Binary Energy Bind =',EB(IBIN)
*
      RETURN
*
      END
