      SUBROUTINE TRIPLE(ISUB,RX,NAMB,TTIME,TRMIN,EKT)
*
*
*       Three-body regularization.
*       --------------------------
*
*       Method of Aarseth & Zare, Celestial Mechanics 10, 185.
*       ......................................................
*
*
      IMPLICIT REAL*8  (A-H,M,O-Z)
      PARAMETER  (NCMAX=10)
      COMMON/AZREG/  TIME3,TMAX,Q(8),P(8),R1,R2,R3,ENERGY,M(3),X3(3,3),
     &               XDOT3(3,3),CM(10),C11,C12A,C19,C20,C24,C25,
     &               NSTEP3,NAME3(3),KZ15,KZ27
      COMMON/CLOSE/  RIJ(4,4),RCOLL,QPERI,SIZE(4),ECOLL3
      COMMON/AZCOLL/  RK(3),QK(8),PK(8),ICALL,ICOLL,NDISS3
      COMMON/BSSAVE/  EP(4),DSC,FACM,TFAC,ITFAC,JC
      COMMON/AZCONF/  ICONF(3)
      DATA ONEPI,C13,C12/3.141592654D0,0.333333333D0,0.5D0/
      REAL*8  Y(17)
      SAVE
*
*
*       Main variables for AZ regularization
*       ************************************
*
*       ---------------------------------------------------------------------
*       CM(1-7) Coordinates, velocities & total mass of triple system.
*       CM(9)   Binary energy change (copied to SBCOLL in START3).
*       CM(10)  Relative energy error (not used).
*       C11     Inverse mass factor for DERQP3 (also C12A,C19,C20,C24,C25).
*       ENERGY  Twice the initial total energy.
*       ICALL   Pericentre indicator (activated if MIN(R1,R2) < RSTAR).
*       ICOLL   Collision or capture indicator (activated in routine DERQP3).
*       KZ15    Triple option copied in START3 (full output if > 1).
*       KZ27    Tidal capture & collision option (copied in START3).
*       M       Particle mass (CM(7) = M(1) + M(2) + M(3)).
*       NAME3   Particle identity (initialized to global name).
*       NREG    Number of regularization switches.
*       NSTEP3  Number of DIFSY3 calls.
*       ICONF   Three-body configuration array (initialized to 1, 2, 3).
*       P       Regularized momenta.
*       Q       Regularized coordinates.
*       QPERI   Pericentre distance of closest two-body motion.
*       R1      Distance between M(1) and M(3).
*       R2      Distance between M(2) and M(3).
*       R3      Distance between M(1) and M(2) (not regularized).
*       RCOLL   Minimum two-body separation.
*       RCRIT   Escape test criterion (RGRAV or harmonic mean with RMAXS).
*       RGRAV   Gravitational radius (sum (M(I)*M(J))/ABS(0.5*ENERGY)).
*       RIJ     Minimum pairwise separations.
*       RSTAR   Pericentre check distance (= 0.2*RGRAV).
*       SIZE    Individual stellar radius (collisions or tidal capture).
*       TIME3   Local physical time in scaled units.
*       TCR     Local crossing time ((sum M(I))**2.5/ABS(ENERGY)**1.5).
*       TOL     Tolerance for DIFSY3 (1.0E-10 is recommended).
*       X3      Particle coordinates (X3(3,I) is Z-component).
*       XDOT3   Velocity components (XDOT3(3,I) is Z-component).
*       ---------------------------------------------------------------------
*
*
*       Save termination indicator and check for restart.
      ITERM = ISUB
*
*       Specify tolerance & DSC for DIFSY3 (relative energy error < TOL).
      TOL = 1.0D-10
      DSC = 1.0
*
*       Initialize diagnostic & COMMON variables.
      R12MIN = 100.0
      R3MIN = 100.0
      RCOLL = 100.0
      ZMI1 = 100.0
      ZMI2 = 100.0
      ZMI3 = 100.0
      DO 10 J = 1,3
          ICONF(J) = J
          DO 5 K = 1,3
              RIJ(J,K) = 100.0
    5     CONTINUE
   10 CONTINUE
      ICALL = 0
      ICOLL = 0
      NEXT = 0
      NDISS3 = 0
      ECOLL3 = 0.0D0
      JC = 0
      ITFAC = 0
      NSTEP3 = 0
      ITRY = 0
      IQ = 0
      NREG = 0
      RMAX3 = RX
*
*       Initialize local time & regularized time.
      TIME3 = 0.0D0
      TAU3 = 0.0D0
      Y(17) = 0.0D0
*       Specify the number of first-order equations for the integrator.
      NEQ = 17
*
*       Transform to regularized variables.
      CALL TRANS3(1)
*
*       Define gravitational radius and pericentre check distance.
      RGRAV = (M(1)*M(2) + M(1)*M(3) + M(2)*M(3))/(0.5D0*ABS(ENERGY))
      RSTAR = 0.2*RGRAV
*       Introduce average mass factor for Bulirsch-Stoer integrator.
      FACM = (M(1)*M(2) + M(1)*M(3) + M(2)*M(3))/3.0
*
*       Ensure that unperturbed boundary exceeds system size.
      R1 = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      R2 = Q(5)**2 + Q(6)**2 + Q(7)**2 + Q(8)**2
*
*       Modify escape distance criterion in case of hierarchical system.
      RCRIT = 5.D0*RGRAV
*
*       Specify local crossing time (also meaningful if ENERGY > 0).
      TCR = (M(1) + M(2) + M(3))**2.5/ABS(ENERGY)**1.5
*
*       Define a nominal crossing time in near-parabolic cases.
      RMAX = MAX(R1,R2)
      TSTAR = RMAX*SQRT(RMAX/CM(7))
*
*       Determine the smallest two-body time-scale from parabolic orbit.
      IM = 1
      RM = R1
      IF (R2.LT.R1) THEN
          IM = 2
          RM = R2
      END IF
      VP2 = 2.0*(M(IM) + M(3))/RM
      TP = RM/SQRT(VP2)
      TSTAR = MIN(TP,TSTAR)
*
*       Set step for time transformation DT/DTAU = R1*R2/(R1 + R2)**0.5.
      TPR = R1*R2/SQRT(R1 + R2)
      DTAU3 = MIN(TCR,TSTAR)*TOL**0.1/TPR      
*     IF(DTAU3.GT.1.D2)DTAU3 = 1.D2
*
*       Form initial binding energy of the binary.
      VREL2 = (XDOT3(1,2) - XDOT3(1,3))**2 +
     &        (XDOT3(2,2) - XDOT3(2,3))**2 +
     &        (XDOT3(3,2) - XDOT3(3,3))**2
*
      EB0 = (0.5D0*VREL2/(M(2) + M(3)) - 1.0/R2)*M(2)*M(3)
*
*     PRINT*,' ***HERE IS TRIPLE****'
*     WRITE(6,1111)NAMB,TIME3,R1,R2
*     PRINT*,' TIME3,DTAU3,ENERGY,EB0=',TIME3,DTAU3,ENERGY,EB0
*     CALL FLUSH(6)
*
*       Initialize input array for the integrator.
   15 DO 20 K = 1,8
          Y(K) = Q(K)
          Y(K+8) = P(K)
   20 CONTINUE
      Y(17) = TIME3
*
*       Advance the equations of motion by Bulirsch-Stoer integrator.
   30 CONTINUE
*       Avoid any external termination criterion
      ITERM = 0
*
*     IF(MOD(NSTEP3,100).EQ.0)PRINT*,' TIME3,DTAU3=',TIME3,DTAU3
      CALL DIFSY3(NEQ,TOL,DTAU3,TAU3,Y)
*
*       Copy regularized coordinates, momenta & time to COMMON variables.
      DO 40 K = 1,8
          Q(K) = Y(K)
          P(K) = Y(K+8)
   40 CONTINUE
      TIME0 = TIME3
      TIME3 = Y(17)
*
*       Update relative distances (NB! Not quite latest value).
      R1 = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      R2 = Q(5)**2 + Q(6)**2 + Q(7)**2 + Q(8)**2
*
*       Check minimum two-body separations and increase step counter.
      R3MIN = MIN(R3MIN,R3)
      RM = MIN(R1,R2)
      R12MIN = MIN(R12MIN,RM)
      NSTEP3 = NSTEP3 + 1
*
*       Check minimum pairwise separations (also in DERQP3).
      RK(1) = R1
      RK(2) = R2
      RK(3) = R3
*
*     IF(MOD(NSTEP3,10).EQ.0)THEN
*     WRITE(6,1111)NAMB,TIME3,R1,R2,R3
 1111 FORMAT(' NAMB TIME3 R1-3=',I5,1P,4E12.4)
*     PRINT*,' TIME3,DTAU3,ENERGY=',TIME3,DTAU3,ENERGY
*     END IF
*
*       Consider pairs 1-2, 1-3 & 2-3 (initial names = 1, 2, 3).
      DO 44 K = 1,3
          DO 42 L = K+1,3
              I = ICONF(K)
              J = ICONF(L)
*       Use cyclic loop index (3,1,2) for distances R3, R1 & R2.
              KK = K - 1
              IF (KK.EQ.0) KK = 3
              RIJ(I,J) = MIN(RIJ(I,J),RK(KK))
              RIJ(J,I) = MIN(RIJ(J,I),RK(KK))
   42     CONTINUE
   44 CONTINUE
*
*       Update smallest moment of inertia during normal forward integration.
      IF (TIME3.GT.TIME0.AND.JC.EQ.0) THEN
          ZMI = RM**2
          ZMI1 = ZMI2
          ZMI2 = ZMI3
          ZMI3 = ZMI
          DTAU0 = DTAU3
      END IF
*
*       Switch on search indicator during iteration or just after pericentre.
      IF (ICOLL.LT.0) ICALL = 1
      IF (RM.LT.RSTAR.AND.NSTEP3.GT.NEXT) THEN
          IF (ZMI3.GT.ZMI2.AND.ZMI2.LT.ZMI1) THEN
              ICALL = 1
          END IF
      END IF
*
*       Check for tidal dissipation or collision during last step.
      IF (ICOLL.LE.0) GO TO 48
*
*       Restore the minimum configuration.
      DO 45 K = 1,8
          Q(K) = QK(K)
          P(K) = PK(K)
   45 CONTINUE
*
*       Set two-body separations.
      R1 = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      R2 = Q(5)**2 + Q(6)**2 + Q(7)**2 + Q(8)**2
*
*       Delay next search a few steps to avoid the same pericentre.
      NEXT = NSTEP3 + 2
*
*       Distinguish between tidal energy loss and collision (ICOLL holds IM).
      IM = ICOLL
      ICOLL = 0
      IF (QPERI.LT.2.7*MAX(SIZE(IM),SIZE(3))) THEN
          IF (QPERI.LT.0.75*(SIZE(IM) + SIZE(3))) THEN
*
*       Combine internal energy and external c.m. kinetic energy.
              H3 = 0.5D0*ENERGY + 0.5D0*CM(7)*(CM(4)**2 + CM(5)**2 +
     &                                                    CM(6)**2)
*
        PRINT*,' SIZES=',(SIZE(K),K=1,3)
        PRINT*,' Collision! '
*       Form composite body and begin KS regularization of new pair.
*       Create composite body for gaseous model here
              GO TO 100
          ELSE
*
*       Modify variables due to tidal effects and check stability parameter.
              CALL QPMOD3(IM,ITERM)
*       Ensure positive step from pericentre and switch off search indicator.
              DTAU3 = ABS(DTAU0)
              ICALL = 0
*       Update relative distances (NB! Not quite latest value).
      R1 = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      R2 = Q(5)**2 + Q(6)**2 + Q(7)**2 + Q(8)**2
              IF (ITERM.LT.0) GO TO 90
*       Initialize input array and continue integration.
              GO TO 15
          END IF
      END IF
*
*       See whether switching of reference body is desirable.
   48 IF (R3.GT.RM) GO TO 70
*
*       Use a simple distance test to determine new reference body IMIN.
      IMIN = 1
      IF (R2.LT.R1) IMIN = 2
*
*       Transform to physical variables and rename the exchanged particles.
      CALL TRANS3(2)
*
      DO 50 K = 1,3
          TEMP1 = X3(K,3)
          TEMP2 = XDOT3(K,3)
          X3(K,3) = X3(K,IMIN)
          XDOT3(K,3) = XDOT3(K,IMIN)
          X3(K,IMIN) = TEMP1
          XDOT3(K,IMIN) = TEMP2
   50 CONTINUE
*
      TEMP1 = M(3)
      M(3) = M(IMIN)
      M(IMIN) = TEMP1
      TEMP2 = SIZE(3)
      SIZE(3) = SIZE(IMIN)
      SIZE(IMIN) = TEMP2
      NAME33 = NAME3(3)
      NAME3(3) = NAME3(IMIN)
      NAME3(IMIN) = NAME33
      I3 = ICONF(3)
      ICONF(3) = ICONF(IMIN)
      ICONF(IMIN) = I3
*
*       Transform back to regularized variables and initialize input array.
      CALL TRANS3(4)
      DO 60 K = 1,8
          Y(K) = Q(K)
          Y(K+8) = P(K)
   60 CONTINUE
*
*       Update regularization counter at the end of switching procedure.
      NREG = NREG + 1
*
*       Set consistent relative distances and minimum separation.
      R1 = Q(1)**2 + Q(2)**2 + Q(3)**2 + Q(4)**2
      R2 = Q(5)**2 + Q(6)**2 + Q(7)**2 + Q(8)**2
      RM = MIN(R1,R2)
*
*       Terminate triple integration if R3 > specified perturber limit.
   70 CONTINUE
*     CALL TRANS3(3)
*     PRINT*,' Now Call TRIPLE with Coordinates:'
*     PRINT*,' Particle X=',X3(1,1),X3(2,1),X3(3,1)
*     PRINT*,' Binary 1X=',X3(1,2),X3(2,2),X3(3,2)
*     PRINT*,' Binary 2X=',X3(1,3),X3(2,3),X3(3,3)
*     PRINT*,' Velocities:'
*     PRINT*,' Particle XDOT=',XDOT3(1,1),XDOT3(2,1),XDOT3(3,1)
*     PRINT*,' Binary 1XDOT=',XDOT3(1,2),XDOT3(2,2),XDOT3(3,2)
*     PRINT*,' Binary 2XDOT=',XDOT3(1,3),XDOT3(2,3),XDOT3(3,3)
*
*     IF (R3.GT.5.D1*RGRAV) GO TO 90
*       Use RGRAVT inherited from trmain.
*     IF (R3.GT.1.D1*RGRAVT) GO TO 90
      IF (ITERM.LT.0) GO TO 74
*
*       Check for hierarchical stability in case of tidal capture.
      IF (NDISS3.GT.0) THEN
          ISKIP = 10
          AN = INT(NSTEP3/FLOAT(ISKIP)) - FLOAT(NSTEP3)/FLOAT(ISKIP)
*       Evaluate stability parameter every ISKIP step.
          IF (ABS(AN).LT.0.01) THEN
              CALL STABL3(ITERM)
              IF (ITERM.LT.0) GO TO 90
          END IF
      END IF
*
*       See if the configuration permits testing of escape condition.
*     GO TO 30
*     IF (R1 + R2 + R3.LT.3.0*RCRIT) GO TO 30
*
*       Obtain current physical variables.
      CALL TRANS3(2)
*
*       Identify index of second binary component & escaper.
      IMIN = 1
      IF (R2.LT.R1) IMIN = 2
      I = 3 - IMIN
*
*       Skip escape test if distant body is approaching the system c.m.
      RIDOT = X3(1,I)*XDOT3(1,I) + X3(2,I)*XDOT3(2,I) +
     &                             X3(3,I)*XDOT3(3,I)
*     GO TO 30
*     IF (RIDOT.LT.0.0) GO TO 30
*
*       Set distance & radial velocity of body #I with respect to binary.
      MB = M(IMIN) + M(3)
      FAC = CM(7)/MB
      RI = SQRT(X3(1,I)**2 + X3(2,I)**2 + X3(3,I)**2)
      RIDOT = FAC*RIDOT/RI
      RI = FAC*RI
*
*       Check the escape criterion due to Standish (Celes. Mech. 4, 44).
      RATIO = RGRAV/(MB*RI)
      VCRIT2 = 2.0*CM(7)*(1.0/RI + M(3)*M(IMIN)*RATIO**2/(RI - RGRAV))
*     GO TO 30
*     IF (RIDOT**2.LT.VCRIT2.OR.RI.LT.RGRAV) GO TO 30
*
*       Define basic variables for termination of tidal capture event.
   74 IF (ITERM.LT.0.AND.NDISS3.GT.0) THEN
          CALL TRANS3(2)
          IMIN = 1
          IF (R2.LT.R1) IMIN = 2
          I = 3 - IMIN
          RI = R3
          RM = MIN(R1,R2)
          CALL EREL3(IM,EBS,SEMI)
          MB = M(IMIN) + M(3)
          FAC = CM(7)/MB
      END IF
*
*       Evaluate orbital elements.
      VREL2 = 0.0D0
      RDOT = 0.0D0
      VI2 = 0.0D0
      VRELCM = 0.D0
      RDOTCM = 0.D0
      RCM = 0.0D0
*            
      DO 75 K = 1,3
          XCM = (M(IMIN)*X3(K,IMIN)+M(3)*X3(K,3))/MB
          XCMDOT = (M(IMIN)*XDOT3(K,IMIN)+M(3)*XDOT3(K,3))/MB
          RCM = RCM + (XCM - X3(K,I))**2
          VRELCM = VRELCM + (XCMDOT - XDOT3(K,I))**2
          RDOTCM = RDOTCM + (XCM - X3(K,I))*(XCMDOT - XDOT3(K,I))    
          RDOT = RDOT +
     &               (X3(K,3) - X3(K,IMIN))*(XDOT3(K,3) - XDOT3(K,IMIN))
          VREL2 = VREL2 + (XDOT3(K,3) - XDOT3(K,IMIN))**2
          VI2 = VI2 + XDOT3(K,I)**2
   75 CONTINUE
      RCM = DSQRT(RCM)
*
      SEMICM = 2.0D0/RCM - VRELCM/(MB+M(I))
      SEMICM = 1.0/SEMICM
      EBCM = -0.5D0*MB*M(I)/SEMICM
      E2 = (1.0D0-RCM/SEMICM)**2+RDOTCM**2/(SEMICM*(MB+M(I)))
      IF(E2.GT.0.D0)THEN
      ECCCM = SQRT(E2)
      ELSE
      ECCCM = 1.D30
      END IF
*              
*
*       Form outer semi-major axis (not used).
      VI2 = VI2*FAC**2
      SEMI1 = 2.0D0/RI - VI2/CM(7)
      SEMI1 = 1.0/SEMI1
*
*       Determine semi-major axis & eccentricity of inner binary.
      RB = RM
      SEMI = 2.0D0/RB - VREL2/MB
      SEMI = 1.0/SEMI
      E2 = (1.0D0 - RB/SEMI)**2 + RDOT**2/(SEMI*MB)
      IF(E2.GT.0.D0)THEN  
      ECC = SQRT(E2)
      ELSE
      ECC = 1.D30
      END IF
*
      IF(SEMI.LT.0.D0.AND.RM.GT.1.2D0*RMAX3)THEN
      PRINT*,' triple: term dissol SEMI,RM,RX=',SEMI,RM,RX
      GO TO 90
      END IF
*       Emergency termination if runaway dissolution.
      IF(SEMI.GT.1.D0)THEN
*        PRINT*,' triple: Term 1 due NSTEP3,TIME3,R3,RMAX3=',
*    *         NSTEP3,TIME3,R3,RMAX3
         WRITE(6,888)TTIME,SEMI,ECC,RM,RX
 888     FORMAT(1X,' triple: T,SEMI,ECC,RM,RX=',1P,5D13.5)
      GO TO 90
      END IF
*
      IF(SEMI.GT.0.D0)THEN
*
*       Check triple left over
      IF(ECCCM.LT.1.D0)THEN
      XCRIT = 2.8D0*SEMI*((1.D0+M(I)/MB)*(1.D0+ECCCM)/
     *    DSQRT(1.D0-ECCCM))**4.D-1
*       Pericentre Distance of outer binary
      YCRIT = SEMICM*(1.D0-ECCCM)
*       Terminate if stable triple acc. to MA is left over
      IF(YCRIT.GT.XCRIT)
     *PRINT*,' triple: Find MA-crit Rpcrit,Xcrit=',YCRIT,XCRIT
      END IF          
*
*     define RMAX3 in the same way as RX in trmain.f
*      
*      XMUE = M(IMIN)*M(3)/MB
*      RMAX3 = 1.2D1*SEMI*(1.D0+ECC)*
*     *                          (1.D0+1.D1*(M(I)/XMUE)**C13)   
      XMUE = 8.D0*M(3)/MB
      RMAX3 = 1.2D1*SEMI*(1.D0+ECC)*(1.D0+1.D1*(XMUE)**C13)

      END IF
*
      IF(R3.GT.RMAX3)THEN
*
         IF(SEMI.LT.0.D0.AND.EBCM.GT.-0.1D0*EKT) THEN
         WRITE(6,887)TTIME,SEMI,ECC,RM,RMAX3,EBCM
 887     FORMAT(1X,' triple: Term due to EBCM>0,TTIME,SEMI,ECC,',
     *     ' RM,RMAX3,EBCM=',1P,6D13.5)
         GO TO 90
         END IF
*
         IF(SEMI.LT.0.D0.AND.RM.LE.1.2D0*RMAX3) THEN
*        PRINT*,' triple: Cont 1 due NSTEP3,TIME3,R3,RMAX3=',
*    *         NSTEP3,TIME3,R3,RMAX3
         WRITE(6,889)TTIME,SEMI,ECC,RM,1.2D0*RMAX3
 889     FORMAT(1X,' triple: T,SEMI,ECC,RM,1.2RMAX3=',1P,5D13.5)
         GO TO 30
         END IF
*     PRINT*,' triple: Term 2 due to NSTEP3,TIME3,R3,RMAX3=',
*    *         NSTEP3,TIME3,R3,RMAX3
      WRITE(6,886)TTIME,SEMI,ECC,SEMICM,ECCCM,RM
 886  FORMAT(1X,' triple: T,SEMI,ECC,SEMICM,ECCCM,RM=',1P,6D13.5)
      GO TO 90
      ELSE
      IF(RMAX3.GT.5.D2)THEN
*     PRINT*,' triple: Cont 2 due to NSTEP3,TIME3,R3,RMAX3=',
*    *         NSTEP3,TIME3,R3,RMAX3
      WRITE(6,886)TTIME,SEMI,ECC,SEMICM,ECCCM,RM
      END IF
      GO TO 30
      END IF
*
*       Delay new regularization if binary is near small pericentre.
      ITRY = ITRY + 1
      IF(NSTEP3.LT.100000)GO TO 30
*     IF (RB.LT.SEMI.AND.ITRY.LE.10) GO TO 30
*
*       Set binary energy ratio for triple system.
      EB = -M(IMIN)*M(3)/(SEMI*ENERGY)
*       Form relative perturbation on inner binary due to body #I.
      GB = 2.0*M(I)/MB*(RB/RI)**3
*
*       Print final configuration and terminate integration of escaper.
      WRITE (6,80)  NAME3(IMIN), NAME3(3), MB, SEMI, ECC, EB, GB, RB,
     &              M(I), RI
   80 FORMAT (/,' BINARY ',2I5,'  MB =',F7.4,'  A =',1PE8.1,'  E =',
     &               0PF5.2,'  EB =',F5.2,'  GB =',1PE8.1,'  RB =',E8.1,
     &                  '  MI =',0PF7.4,'  RI =',1PE8.1)
*
*       Transform to physical variables for termination.
   90 CALL TRANS3(3)
*
*       Identify second component and obtain osculating elements.
      IMIN = 1
      IF (R2.LT.R1) IMIN = 2
      RB = MIN(R1,R2)
      VREL2 = (XDOT3(1,IMIN) - XDOT3(1,3))**2 +
     &        (XDOT3(2,IMIN) - XDOT3(2,3))**2 +
     &        (XDOT3(3,IMIN) - XDOT3(3,3))**2
*
*       Form the semi-major axis & energy of closest pair.
      SEMI = 2.0D0/RB - VREL2/(M(IMIN) + M(3))
      SEMI = 1.0/SEMI
      EB = -0.5D0*M(IMIN)*M(3)/SEMI
*
*       Avoid new KS regularization inside half semi-major axis.
*     ITRY = ITRY + 1
*     IF (RB.LT.ABS(SEMI).AND.ITRY.LE.10) THEN
*         NEXT = NSTEP3 + 10
*         GO TO 15
*     END IF
*
*       Check optional print diagnostics of triple integration.
      IF (KZ15.GT.1) THEN
          TC = TIME3/TCR
*       Print relative energy change of binary (including exchange).
          DB = (EB - EB0)/EB0
          WRITE (6,95)  NREG, R12MIN, R3MIN, R3, RGRAV, TC, NSTEP3, DB,
     &                  NAME3(3-IMIN)
   95     FORMAT (' END TRIPLE    NREG =',I3,'  MIN(RB&R3) =',1PE8.1,
     &                 E9.1,'  R3 =',E8.1,'  RG =',E8.1,'  TC =',E12.5,
     &                 '  STEPS =',I8,'  DB =',1P,E13.5,'  I =',I5)
          WRITE (6,91) TTIME,MIN(R12MIN,R3MIN),R1,R2,R3
   91     FORMAT (' END TRIPLE T,DMIN=',1P,2E13.5,' R1-3=',3E13.5)
          TRMIN = MIN(R12MIN,R3MIN)
      END IF
*
*       Form binary energy change (set in SBCOLL).
      CM(9) = EB - EB0
*
  100 CONTINUE
      ISUB = ITERM
*       Return RGRAV for diagnostics in TRMAIN
      RX = RGRAV
      RETURN
*
      END
