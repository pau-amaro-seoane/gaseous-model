          SUBROUTINE DENS(JJ,AMS)
C*******************************************************************
C  
C         Compute Initial Density Profile for Component JJ
C
C        Parameter NOL(JJ) chooses different models
C        EXP(JJ) power law exponent (only for NOL(JJ)=2
C        RSC(JJ) scale radius
C
C        Note: Parameters AMS scales with the total mass;
C              DENS is called by INITST two times:
C              first with arbitrary AMS;
C              second with AMS scaled in order to reach XMTOT(J)
C                           as desired total mass
C******************************************************************
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
         INCLUDE 'compar.f'
         CHARACTER*7 NAME,NAMED
C
           RR=RSC(JJ)
           EXPO=EXP(JJ)
           N=IABS(NOL(JJ))
*      Code for King model: NOL > 100
          IF(N.GT.100)THEN
          W0 = DBLE(N-100)
          N=4
          END IF
*
           IF(N.GT.20)N=N-20
C
          GOTO(1,2,3,4,5,6) N
 4        CONTINUE
*
          IUNIT = 70 + NOL(JJ) - 100
          WRITE(NAME,100)IUNIT
 100      FORMAT('fort.',I2)
          PRINT*,' King Model Data to be read from file ',NAME
          OPEN(IUNIT,FILE=NAME)
*
          READ(IUNIT,*)W0,CONC,XIRAT,RKING,RTIDAL
          READ(IUNIT,'(A7)')NAMED
*
          DO 10 J=2,NJ
          IF(J.GT.NJ)GOTO 18
          READ(IUNIT,*,ERR=19,END=14)R(J),XD,RRJ,XD2,VJ
          IF(J.EQ.2)THEN
          PHI0 = -XD
          VJ0 = VJ
*       Determine central sigma from King data
          SIGC2 = -PHI0/W0*VJ0
          PRINT*,' King model: W0,CONC=',W0,CONC,' RK,RT=',RKING,RTIDAL
          PRINT*,' King model: PHI0=',PHI0,' SIGC2=',SIGC2
          END IF
*
          RHOM = RRJ*AMS/XMTOT(JJ)
          SIG2 = SIGC2 * VJ/VJ0
*         PRINT*,' J,R,RHO,SIG2=',J,R(J),RHOM,SIG2
*         IF(IZ00.EQ.0.AND.RHOM.LT.1.D-4*DEXP(X(I00+JJ,2)))GOTO 14
          X(I00+JJ,J)=DLOG(RRJ*AMS/XMTOT(JJ))
          X(I20+JJ,J)=X(I00+JJ,J)+DLOG(SIG2)
          X(I02+JJ,J)=X(I20+JJ,J)
 10       CONTINUE
*      Pressure determination redundant here - more detailed in initst.f
 14       CONTINUE
*      Assume total mass is normalized to one here
          PHTID = -GRAV/RTIDAL
          PRINT*,'PHTID,J,XMTT=',PHTID,J,1.0
*           Cut the model or extend it outside tidal radius
          IF(LS(27))THEN
          IF(IZ00.EQ.0)NJ = J - 1
          IZ00 = 1
          END IF
*
          JM = J-1
          PRINT*,JM,' Data Points read from King Model data'
          PRINT*,' First Rho=',DEXP(X(I00+JJ,2)),' Last Rho=',
     *      RHOM,' Last R=',R(JM),' AMS=',AMS,' XMT=',XMTOT(JJ),
     *      ' NJ=',NJ,' PHTID=',PHTID
          REWIND IUNIT
*
          IF(J.LE.NJ)THEN
          DXDR = 3.D0
          DO 17 IJ=J,NJ
          R(IJ) = R(IJ-1)*STEU3
          RRJ = R(JM)/R(IJ)
          X(I00+JJ,IJ)=DLOG(RHOM) + DXDR*DLOG(RRJ)
          X(I20+JJ,IJ)=DLOG(SIG2*RHOM) + (DXDR+1.D0)*DLOG(RRJ)
          X(I02+JJ,IJ)=X(I20+JJ,IJ)
 17       CONTINUE
          END IF
*
 18       CONTINUE
          R(1) = 0.D0
          X(I00+JJ,1) = X(I00+JJ,2)
          X(I20+JJ,1) = X(I20+JJ,2)
          X(I02+JJ,1) = X(I02+JJ,2)
          RRJ = R(JM)/R(NJ)/STEU3
          XOUT(JJ)=2.D0*X(I00+JJ,NJ)-X(I00+JJ,NJ-1)
          PROUT(JJ)=2.D0*X(I20+JJ,NJ)-X(I20+JJ,NJ-1)
          PTOUT(JJ)=PROUT(JJ)
*
          RETURN
C********************************************************************
C      NOL=5:
C      n=5 Polytrope rho propto 1/(1+(r/rsc)**2)**5/2
C      Logarithmic variables
C--------------------------------------------------------------------
 5        FACD=1.D0/PI43/RR*AMS/RR/RR
C
          DO 20 J=2,NJ
          JM=J-1
          RRJ=((R(JM)+R(J))/RR/2.D0)**2+1.D0
          IF(J.EQ.2)RRJ=((R(J)/STEU3+R(J))/RR/2.D0)**2+1.D0
          X(I00+JJ,J)=DLOG(FACD)-25.D-1*DLOG(RRJ)
 20       CONTINUE
C           Outer Value
          ROUT=R(NJ)*STEU1
          RRJ=((R(NJ)+ROUT)/RR/2.D0)**2+1.D0
          XOUT(JJ)=DLOG(FACD)-25.D-1*DLOG(RRJ)
C
          X(I00+JJ,1)=X(I00+JJ,2)
C           For tidal cut of Plummer's model use W0 input parameter
C           as potential at tidal radius
          RETURN
C*******************************************************************
C       NOL=2:
C
C      constant density from centre to RSC
C      outwards power law with exponent EXP
C      logarithmic variables
C--------------------------------------------------------------------
  2       CONTINUE
          FACC=R(NJ)/R(2)-2.D0/3.D0
          FACD=1.D0/PI43/RR*AMS/RR/RR
C
          X(I00+JJ,1)=DLOG(FACD)
          DO 50 J=2,NJ
          RI=RR
          IF(R(J).LE.RR) RI=R(J)
          IF(DABS(EXPO).GT.0.D0)THEN
          X(I00+JJ,J)=DLOG(FACD)+EXPO*(DLOG(RI)-DLOG(R(J)))
          ELSE
          X(I00+JJ,J)=DLOG(FACD)
          END IF
  50      CONTINUE
C
          RETURN
C********************************************************************
C       NOL=1:
C       so called Hubble-profile
C       rho propto 1/(1+(r/rsc)**2)**3/2
C----------------------------------------------------------------------
  1       CONTINUE
          FACD=1.D0/PI43/RR*AMS/RR/RR
          X(I00+JJ,1)=DLOG(FACD)
          X(I00+JJ,2)=DLOG(FACD)-15.D-1*DLOG((R(2)/RR)**2+1.D0)
C
          DO 15 J=3,NJ
          JM=J-1
          RRJ=((R(JM)+R(J))/RR/2.D0)**2+1.D0
          X(I00+JJ,J)=DLOG(FACD)-15.D-1*DLOG(RRJ)
 15       CONTINUE
C
          RETURN
C
  3       CONTINUE
C
          RETURN
C********************************************************************
C       NOL=6:
C       so called Hernquist-profile
C       rho propto 1/r/(rsc+r)**3
C----------------------------------------------------------------------
  6       CONTINUE
          FACD=AMS/2.D0/PI*RR
          RFAC=R(3)/R(2)
          RLOW=R(2)/RFAC
          X(I00+JJ,1)=DLOG(FACD/RLOW/(RLOW+RR)**3)
          RLOW=(RLOW+R(2))/2.D0
          X(I00+JJ,2)=DLOG(FACD/RLOW/(RLOW+RR)**3)
C
          DO 16 J=3,NJ
          JM=J-1
          RRJ=(R(JM)+R(J))/2.D0
          X(I00+JJ,J)=DLOG(FACD/RRJ/(RRJ+RR)**3)
 16       CONTINUE
C
          RETURN
 19       PRINT*,' READ ERROR from ',NAME
          STOP
          END
