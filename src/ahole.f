      SUBROUTINE AHOLE
C     This subprogram contains all operations concerning central fixed
C     star-accreting black hole. It is called by main program star and
C     (some entries) by other subroutines.
C
C     Variables and their meanings:
C
C     XMHOLE: Mass of (fixed) central black hole
C     EFFI: efficiency of mass energy conversion at tidal radius
C     RTIDE(NCOMP): Tidal disruption radius for each stellar component
C     RHOSTA(NCOMP): Average stellar density "    "     "    "(used for RTIDE)
C     XBIND(NCOMP): Average stellar binding energy "    "    "(used for RTIDE)
C     SMHOLE: mass of (moving) central black hole (not yet supported)
C     ROSZ: oscillation radius of (moving) central black hole (n.y.supp.)
C     DMHOLE(5,NCOMP): field containing various contributions to
C     star accretion rate of black hole
C     DMHOLE(1) : energy diffusion of orbits
C     DMHOLE(3) : loss cone accretion rate
C     DMHOLE(4) : gas accretion rate (not yet supported)
C     DMHOLE(5) : loss cone accretion rate for enhanced loss cone
C
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
      INCLUDE 'compar.f'
      LOGICAL LOAD(NCOMPO)
      DIMENSION YLC(NCOMPO,NJO),YLC2(NCOMPO,NJO)
C
      DO 1000 J=1,NCOMP
C
         IF(IMOD.GT.0)THEN
            XMHOLE=XMHOLE+DMHOLE(1,J)/BREM
            DELTM(1,J)=DELTM(1,J)+DMHOLE(1,J)/BREM
            DMTOT=DMTOT+DMHOLE(1,J)/BREM
         END IF
C
         RTIDE(J)=(2.D0*XMHOLE/(PI43*RHOSTA(J)*XBIND(J)))**C13
C
         I=2
C
         SIG2=DEXP(X(I20+J,I)-X(I00+J,I))*
     *        (1.D0+2.D0*DEXP(X(I02+J,I)-X(I20+J,I)))
         TRX(J,J)=(C13*SIG2)**C32/
     *        CS(60)/XMIND(J)/XCOUL(J)/DEXP(X(I00+J,I))
         DMHOLE(1,J)=PI43*DEXP(X(I00+J,I))*R(I)**3/TRX(J,J)
C
 1000 CONTINUE
C
      K=0
C
      DO 50 I=1,NJ
         IF(DEXP(X(IMR,I)).LT.XMHOLE) K = K+1
 50   CONTINUE
C
      DM=DEXP(X(IMR,K+1)) - DEXP(X(IMR,K))
      DR=R(K+1) - R(K)
      RGRAV=R(K)+(2.D0*XMHOLE-DEXP(X(IMR,K)))*DR/DM
C
      RETURN
C
      ENTRY LHOLE
C
C     Computation of loss cone swallowing rate and critical radius estimate
C
      DO 1001 J=1,NCOMP
         LCRIT=.TRUE. ! Modif M.Freitag 11.10.03
C
         IF(IMOD.GT.0)THEN
            XMHOLE=XMHOLE+DMHOLE(3,J)/BREM
            DELTM(3,J)=DELTM(3,J)+DMHOLE(3,J)/BREM
            DMTOT=DMTOT+DMHOLE(3,J)/BREM
         END IF
C
         DMHOLE(3,J)=0.D0
         IF(IMOD.EQ.0)XBETA=1.D0
C
         DO 99 I=2,NJ
            SIG2=DEXP(X(I20+J,I)-X(I00+J,I))*
     *           (1.D0+2.D0*DEXP(X(I02+J,I)-X(I20+J,I)))
            TRX(J,J)=(C13*SIG2)**C32/
     *           CS(60)/XMIND(J)/XCOUL(J)/DEXP(X(I00+J,I))
            TCR=R(I)/DEXP((X(I20+J,I)-X(I00+J,I))/2.D0)
C
            TOUT=TCR
*     Prepare HX for use in routine LOSSCO
            HX(I00+J)=X(I00+J,I)
            HX(I20+J)=X(I20+J,I)
            HX(I02+J)=X(I02+J,I)
            HX(I10+J)=X(I10+J,I)
*
            RLC=RTIDE(J)
            IF(R(I).LT.1.D1*R(2))RLC=RLC*R(I)/1.D1/R(2)
*
            CALL LOSSCO(XMHOLE,RLC,J)
*
            ERARG=DSQRT(4.D0*TRX(J,J)*DOMEGA/TOUT/XBETA)

            IF(ERARG.LT.1.D1)THEN
               EARG=ERARG
               PF=DERF(EARG)
            ELSE
               PF=1.D0
            END IF


*            IF((ERARG.GT.1.D0).AND.LCRIT)THEN ! Modif M. Freitag 11.10.03
             IF((ERARG.LT.2.D0).AND.LCRIT)THEN ! Modif P. Amaro-Seoane
                                               ! and R. Spurzem
                                               ! Beijing, 12.08.2011
                                               ! 1- Equation (23) Amaro-Seoane et al 2004
                                               ! (sqrt and factor 4 were missing in Marc's correction)
                                               ! 2- ERARG is _decreasing_
               KL=I
               PRINT*,' CRITICAL RADIUS AT LHOLE I=',I,' RCRIT=',R(I)
               LCRIT=.FALSE.
            END IF

            XK=1.D0+XBETA*TOUT/DOMEGA/PF/TRX(J,J)
            IF(IMOD.EQ.0)Y(J,I)=TOUT/TRX(J,J)/DOMEGA/PF/XK
*     Start with empty loss cone
*     IF(IMOD.EQ.0)Y(J,I)=0.D0
            IF(Y(J,I).GT.1.D0)Y(J,I)=1.D0
            DMHOLE(3,J)=DMHOLE(3,J)+PI43*(R(I)**3-R(I-1)**3)*
     *           Y(J,I)*DEXP(X(I00+J,I))*DOMEGA/TOUT*PF
c$$$            IF((PF.GE.1.D0).AND.LCRIT)THEN
c$$$               KL=I
c$$$               PRINT*,' CRITICAL RADIUS AT LHOLE I=',I,' RCRIT=',R(I)
c$$$               LCRIT=.FALSE.
c$$$            END IF
 99      CONTINUE
         RCRIT(J)=R(KL)
C
 1001 CONTINUE
C
      IF(IMOD.EQ.0)THEN
         DO 1010 J=1,NCOMP
            WRITE(6,215)J,RTIDE(J),RCRIT(J),DMHOLE(1,J),DMHOLE(3,J),
     *           -GRAV*XMHOLE/RTIDE(J)
 1010    CONTINUE
 215     FORMAT(' START HOLE COMP. ',I3,' RTIDE=',1P,D11.4,
     *        ' RCRIT=',D11.4,' DMHOLE(1)=',D11.4,' DMHOLE(3)=',D11.4,
     *        ' PHIRT=',D11.4)
      END IF
C
C
C**** Loss cone update between the timesteps
C**** Solve diffusion equation for filling degree
C**** and advect loss-cone stars
C
C     number of orbits until trapping
*     do not use XALPHA here because now it is used for tidal field
      YALPHA=1.D0
      XBETA=1.D0
C**   CALCULATE AND STORE DOMEGA FOR J,J-1,J+1=VDOI,VDOM,DOMEGA********
      VDOM=0.D0
C
      DO 1002 JC=1,NCOMP
C     diminish RLC for inner radii for stability (small error in MDOT)
C
         IF(IMOD.GT.0)THEN
C
            I=2
*     Prepare HX for use in routine LOSSCO
            HX(I00+JC)=X(I00+JC,I)
            HX(I20+JC)=X(I20+JC,I)
            HX(I02+JC)=X(I02+JC,I)
            HX(I10+JC)=X(I10+JC,I)
*
            RLC=RTIDE(JC)
            RLC=RTIDE(JC)/1.D1
*
            CALL LOSSCO(XMHOLE,RLC,JC)
C
            VDOI=DOMEGA
C
         END IF
C
         IF(MOD(NMOD,NCOL).EQ.0)WRITE(6,106)
C
         DO 159 J=2,NJ
C
            DV=PI43*(R(J)**3-R(J-1)**3)
C
            IF(J.GT.2)THEN
               VDOM=VDOI
               VDOI=DOMEGA
            END IF
C
            I=MIN(J+1,NJ-1)
*     Prepare HX for use in routine LOSSCO
            HX(I00+JC)=X(I00+JC,I)
            HX(I20+JC)=X(I20+JC,I)
            HX(I02+JC)=X(I02+JC,I)
            HX(I10+JC)=X(I10+JC,I)
*
*     diminish RLC for inner radii for stability (small error in MDOT)
            RLC=RTIDE(JC)
            IF(R(I).LT.1.D1*R(2))RLC=RLC*R(I)/1.D1/R(2)
*
            CALL LOSSCO(XMHOLE,RLC,JC)
*
            IF(IMOD.GT.0)THEN
C
C**** END OF PREPARATION OF DOMEGA FOR THREE NEIGHBOURING GRIDS*****
C**** LOSS CONE TRANSPORT******************************************
               FIM=1.D0
               IF(X(I10+JC,J-1).GT.0)FIM=0.D0
               FIP=1.D0-FIM
               GIM=1.D0
               IF(X(I10+JC,J).GT.0)GIM=0.D0
               GIP=1.D0-GIM
*
               HX(I00+JC)=DEXP(X(I00+JC,J))
               HXM(I00+JC)=DEXP(X(I00+JC,J-1))
               HXP(I00+JC)=DEXP(X(I00+JC,J+1))
*
          DIVPK=R(J-1)**2*X(I10+JC,J-1)*VDOI*(FIP*HXM(I00+JC)*Y(JC,J-1)+
     *              FIM*HX(I00+JC)*Y(JC,J))-R(J)**2*X(I10+JC,J)*DOMEGA*
     *              (GIP*HX(I00+JC)*Y(JC,J)+GIM*HXP(I00+JC)*Y(JC,J+1))
*     No lc transport for test purposes
               YLC(JC,J)=Y(JC,J)
*     YLC(JC,J)=DEXP(VX(I00+JC,J))/HX(I00+JC)*Y(JC,J)+
*     *                    DIVPK/HX(I00+JC)/VDOI/DV
            END IF
C
C******END LOSS-CONE TRANSPORT*************************************
C
            SIGR2=DEXP(X(I20+JC,J)-X(I00+JC,J))
            SIGR=DSQRT(SIGR2)
            SIGT2=DEXP(X(I02+JC,J)-X(I00+JC,J))
            SIG2=SIGR2*(1.D0+2.D0*DEXP(X(I02+JC,J)-X(I20+JC,J)))
            TRX(JC,JC)=(C13*SIG2)**C32/
     *           DEXP(X(I00+JC,J))/CS(60)/XMIND(JC)/XCOUL(JC)
            TCR=R(J)/SIGR
            TOUT=TCR
C
C***  XCF=SQRT OF RATIO OF TOUT/TIN**XCF1=SAME FOR NTR=1*************
            XCF1=DSQRT(TOUT*XBETA/TRX(JC,JC)/DOMEGA)
            XCF=XCF1/DSQRT(YALPHA)
            IF(XCF.LT.1.D0)ICRIT=J
            IF(XCF1.LT.1.D0)ICRIT1=J
            ERARG=2.D0/XCF
            IF(ERARG.LT.1.D1)THEN
               EARG=ERARG
               PF=DERF(EARG)
            ELSE
               PF=1.D0
            END IF
C
            XK=1.D0+TOUT*XBETA/YALPHA/TRX(JC,JC)/DOMEGA/PF
            ERARG=2.D0/XCF1
            IF(ERARG.LT.1.D1)THEN
               EARG=ERARG
               PF1=DERF(EARG)
            ELSE
               PF1=1.D0
            END IF
C***  Update of lc-filling-factor in YLC; storage on Y(I) later***
            IF(IMOD.GT.0)THEN
               XEXP=DEXP(-YALPHA*PF*XK/TOUT/BREM)
               YLC(JC,J)=YLC(JC,J)*XEXP+XCF*XCF/PF/XK*(1.D0-XEXP)
               IF(YLC(JC,J).GT.1.D0)YLC(JC,J)=1.D0
               IF(YLC(JC,J).LE.0.D0)YLC(JC,J)=1.D-30
            END IF
C
C
C     IF(LSTE(19))THEN
C     W(200+J)=PF
C     IF(W(200+J).GT.1.D0)W(200+J)=1.D0
C     DRHO=Y(JC,J)*DEXP(X(I00+JC,J))/TRX(JC,JC)
C     DRHO1=DRHO
C     END IF
C
C
            DRHO=Y(JC,J)*DEXP(X(I00+JC,J))*DOMEGA/TOUT
            DRHO1=PF1*DRHO
            DRHO=YALPHA*PF*DRHO
C**** DRHO1: Strike rate for one crossing, DRHO: trapping rate for ntr**
            DMDT=DRHO*DV
            THLC=DATAN(VLC*SQRT(2.D0*SIGT2)/SIGR)
            THLCO=DSQRT(PI*DOMEGA)
C***  THLC=Frank/Rees-Model, THLCO: Velocity Integration with SB-distr.***
            THD=DSQRT(TOUT/TRX(JC,JC)/YALPHA)
            TIN=TRX(JC,JC)*DOMEGA
            XER=EPSR/DOMEGA
            XET=EPST/DOMEGA
C*****Begin information printing block**********************************
 106        FORMAT(1X,' **(JC)**#1:I= ','#2:K(T)=**** ','#3:MDOTLC=** ',
     *           '#4:XER=***** ','#5:XET=***** ','#6:PF(1TCR)= ',
     *           '#7:TH-D=**** ','#8:TH-LC=*** ','#9:THLC-OM=* ',
     *           '#10:ULC/SR=** ','#11:VLC/ST=** ',
     *           ' #12:TIN=**** ',' #13:TOUT=*** ',' #14:DOMEGA=* ',
     *           '#15:ERARG=** ','#16:R(J)=** ')
            IF((MOD(NMOD,NCOL)+MOD(J,NWRITE)).EQ.0)WRITE(6,107)
     *           JC,J,Y(JC,J),DMDT,XER,XET,PF,
     *           THD,THLC,THLCO,ULC,VLC,TIN,TOUT,DOMEGA,
     *           ERARG,R(J)
 107        FORMAT(1X,I4,I4,1P,16(' ',D11.4))
C*******End information printing block**********************************
C
 159     CONTINUE
C
         IF(IMOD.EQ.0)RETURN
C
C*****Store new filling of loss-cone************************************
         DO 1559 I=2,NJ
            Y(JC,I)=YLC(JC,I)
 1559    continue
C
 1002 CONTINUE
C
      RETURN
C*********************************************************************
      ENTRY LHOLEN
C
C     Computation of loss cone swallowing rate and critical radius estimate
C     for loss cone enhanced by small non-sphericity of the system
C
      LCRIT=.TRUE.
      DO 2001 J=1,NCOMP
C
         DMHOLE(5,J)=0.D0
         IF(IMOD.EQ.0)XBETA=1.D0
C
         DO 299 I=2,NJ
            SIG2=DEXP(X(I20+J,I)-X(I00+J,I))*
     *           (1.D0+2.D0*DEXP(X(I02+J,I)-X(I20+J,I)))
            TRX(J,J)=(C13*SIG2)**C32/
     *           CS(60)/XMIND(J)/XCOUL(J)/DEXP(X(I00+J,I))
            TCR=R(I)/DEXP((X(I20+J,I)-X(I00+J,I))/2.D0)
*     Prepare HX for use in routine LOSSCO
            HX(I00+J)=X(I00+J,I)
            HX(I20+J)=X(I20+J,I)
            HX(I02+J)=X(I02+J,I)
            HX(I10+J)=X(I10+J,I)
*
            RLC=RTIDE(J)
            IF(R(I).LT.1.D1*R(2))RLC=RLC*R(I)/1.D1/R(2)
*
            CALL LOSSCO(XMHOLE,RLC,J)
C
C     Change TOUT due to small non-sphericity of the system
C
C            IF(DEXP(X(IMR,I)).LT.XMHOLE)THEN
C               TCR=DSQRT(R(I)**3/GRAV/XMHOLE)
C               VMEAN=VLCTH/2.D0*(VLCTH-VLC)+VLCPH/2.D0*(VLCPH-VLC)
C               DEN=DSQRT(1-R(I)/GRAV/XMHOLE*(VLCPH**2+VLCTH**2)/4.D0)
C               DEELL=R(I)/GRAV/XMHOLE*VMEAN/DEN
C               TOUT=TCR*XMHOLE/DEXP(X(IMR,I))*5.D0/6.D0/ELL*DEELL
C            ELSE
C               TOUT=TCR*5.D0/6.D0/ELL
C            END IF
            IF(DEXP(X(IMR,I)).LT.XMHOLE)THEN
               TCR=DSQRT(R(I)**3/GRAV/XMHOLE)
               VMEAN=VLCTH/2.D0*(VLCTH-VLC)+VLCPH/2.D0*(VLCPH-VLC)
*     Pau: sqrt(...) invalid sometimes, check first
               IF (R(I)/GRAV/XMHOLE*(VLCPH**2+VLCTH**2)/4.D0.GT.1.0)
     &         THEN
                  TOUT=TCR*5.D0/6.D0/ELL
               ELSE
                  DEN=DSQRT(1-R(I)/GRAV/XMHOLE*(VLCPH**2+VLCTH**2)/4.D0)
                  DEELL=R(I)/GRAV/XMHOLE*VMEAN/DEN
                  TOUT=TCR*XMHOLE/DEXP(X(IMR,I))*5.D0/6.D0/ELL*DEELL
               END IF
            ELSE
               TOUT=TCR*5.D0/6.D0/ELL
            END IF



C
            ERARG=DSQRT(4.D0*TRX(J,J)*DOMEG1/TOUT/XBETA)

            IF(ERARG.LT.1.D1)THEN
               EARG=ERARG
               PF=DERF(EARG)
            ELSE
               PF=1.D0
            END IF

*            IF((ERARG.GT.1.D0).AND.LCRIT)THEN ! Modif M. Freitag 11.10.03
             IF((ERARG.LT.2.D0).AND.LCRIT)THEN ! Modif P. Amaro-Seoane
                                               ! and R. Spurzem
                                               ! Beijing, 12.08.2011
                                               ! 1- Equation (23) Amaro-Seoane et al 2004
                                               ! (sqrt and factor 4 were missing in Marc's correction)
                                               ! 2- ERARG is _decreasing_
               KL=I
               PRINT*,' CRITICAL RADIUS AT LHOLE I=',I,' RCRIT=',R(I)
               LCRIT=.FALSE.
            END IF

            XK=1.D0+XBETA*TOUT/DOMEG1/PF/TRX(J,J)
            IF(IMOD.EQ.0)Y2(J,I)=TOUT/TRX(J,J)/DOMEG1/PF/XK
            IF(Y2(J,I).GT.1.D0)Y2(J,I)=1.D0
            DMHOLE(5,J)=DMHOLE(5,J)+PI43*(R(I)**3-R(I-1)**3)*
     *           Y2(J,I)*DEXP(X(I00+J,I))*DOMEG1/TOUT*PF
C             IF((PF.GT.1.D0).AND.LCRIT)THEN
C                KL=I
C                PRINT*,' CRITICAL RADIUS AT LHOLE I=',I,' RCRIT=',R(I)
C                LCRIT=.FALSE.
C             END IF
 299     CONTINUE
C
         XMHOLE=XMHOLE+DMHOLE(5,J)/BREM
         DELTM(5,J)=DELTM(5,J)+DMHOLE(5,J)/BREM
         DMTOT=DMTOT+DMHOLE(5,J)/BREM
         RCRIT(J)=R(KL)
C
 2001 CONTINUE
C
      IF(IMOD.EQ.0)THEN
         DO 2010 J=1,NCOMP
            WRITE(6,2215)J,RTIDE(J),RCRIT(J),DMHOLE(1,J),DMHOLE(5,J)
 2010    CONTINUE
 2215    FORMAT(' START HOLE COMP. ',I3,' RTIDE=',1P,D11.4,
     *        ' RCRIT=',D11.4,' DMHOLE(1)=',D11.4,' DMHOLE(5)=',D11.4)
      END IF
C
C
C**** Loss cone update between the timesteps
C**** Solve diffusion equation for filling degree
C**** and advect loss-cone stars
C
C     number of orbits until trapping
      YALPHA=1.D0
      XBETA=1.D0
C**   CALCULATE AND STORE DOMEGA FOR J,J-1,J+1=VDOI,VDOM,DOMEGA********
      VDOM=0.D0
C
      DO 2002 JC=1,NCOMP
C
         IF(IMOD.GT.0)THEN
C
            I=2
*     Prepare HX for use in routine LOSSCO
            HX(I00+JC)=X(I00+JC,I)
            HX(I20+JC)=X(I20+JC,I)
            HX(I02+JC)=X(I02+JC,I)
            HX(I10+JC)=X(I10+JC,I)
*
C     diminish RLC for inner radii for stability (small error in MDOT)
            RLC=RTIDE(JC)
            IF(R(I).LT.1.D1*R(2))RLC=RLC*R(I)/1.D1/R(2)
*
            CALL LOSSCO(XMHOLE,RLC,JC)
C
            VDOI=DOMEG1
C
         END IF
C
         IF(MOD(NMOD,NCOL).EQ.0)WRITE(6,206)
C
         DO 259 J=2,NJ
C
            DV=PI43*(R(J)**3-R(J-1)**3)
C
            IF(J.GT.2)THEN
               VDOM=VDOI
               VDOI=DOMEG1
            END IF
C
            I=MIN(J+1,NJ-1)
*     Prepare HX for use in routine LOSSCO
            HX(I00+JC)=X(I00+JC,I)
            HX(I20+JC)=X(I20+JC,I)
            HX(I02+JC)=X(I02+JC,I)
            HX(I10+JC)=X(I10+JC,I)
*
*     diminish RLC for inner radii for stability (small error in MDOT)
            RLC=RTIDE(JC)
            IF(R(I).LT.1.D1*R(2))RLC=RLC*R(I)/1.D1/R(2)
*
            CALL LOSSCO(XMHOLE,RLC,JC)
*
            IF(IMOD.GT.0)THEN
C
C**** END OF PREPARATION OF DOMEGA FOR THREE NEIGHBOURING GRIDS*****
C**** LOSS CONE TRANSPORT******************************************
               FIM=1.D0
               IF(X(I10+JC,J-1).GT.0)FIM=0.D0
               FIP=1.D0-FIM
               GIM=1.D0
               IF(X(I10+JC,J).GT.0)GIM=0.D0
               GIP=1.D0-GIM
*
               HX(I00+JC)=DEXP(X(I00+JC,J))
               HXM(I00+JC)=DEXP(X(I00+JC,J-1))
               HXP(I00+JC)=DEXP(X(I00+JC,J+1))
*
         DIVPK=R(J-1)**2*X(I10+JC,J-1)*VDOI*(FIP*HXM(I00+JC)*Y2(JC,J-1)+
     *              FIM*HX(I00+JC)*Y2(JC,J))-R(J)**2*X(I10+JC,J)*DOMEG1*
     *              (GIP*HX(I00+JC)*Y2(JC,J)+GIM*HXP(I00+JC)*Y2(JC,J+1))
*     No lc transport for test purposes
               YLC2(JC,J)=Y(JC,J)
*     YLC2(JC,J)=DEXP(VX(I00+JC,J))/HX(I00+JC)*Y2(JC,J)+
*     *                    DIVPK/HX(I00+JC)/VDOI/DV
            END IF
C
C******END LOSS-CONE TRANSPORT*************************************
C
            SIGR2=DEXP(X(I20+JC,J)-X(I00+JC,J))
            SIGR=DSQRT(SIGR2)
            SIGT2=DEXP(X(I02+JC,J)-X(I00+JC,J))
            SIG2=SIGR2*(1.D0+2.D0*DEXP(X(I02+JC,J)-X(I20+JC,J)))
            TRX(JC,JC)=(C13*SIG2)**C32/
     *           DEXP(X(I00+JC,J))/CS(60)/XMIND(JC)/XCOUL(JC)
            TCR=R(J)/SIGR
C     Correction of TOUT due to small non-sphericity of the system
C
            IF(DEXP(X(IMR,I)).LT.XMHOLE)THEN
               TCR=DSQRT(R(I)**3/GRAV/XMHOLE)
               VMEAN=VLCTH/2.D0*(VLCTH-VLC)+VLCPH/2.D0*(VLCPH-VLC)
*     Pau: sqrt(...) invalid sometimes, check first
               IF (R(I)/GRAV/XMHOLE*(VLCPH**2+VLCTH**2)/4.D0.GT.1.0)
     &         THEN
                  TOUT=TCR*5.D0/6.D0/ELL
               ELSE
                  DEN=DSQRT(1-R(I)/GRAV/XMHOLE*(VLCPH**2+VLCTH**2)/4.D0)
                  DEELL=R(I)/GRAV/XMHOLE*VMEAN/DEN
                  TOUT=TCR*XMHOLE/DEXP(X(IMR,I))*5.D0/6.D0/ELL*DEELL
               END IF
            ELSE
               TOUT=TCR*5.D0/6.D0/ELL
            END IF

c            IF(DEXP(X(IMR,I)).LT.XMHOLE)THEN
c               TCR=DSQRT(R(I)**3/GRAV/XMHOLE)
c               VMEAN=VLCTH/2.D0*(VLCTH-VLC)+VLCPH/2.D0*(VLCPH-VLC)
c               DEN=DSQRT(1-R(I)/GRAV/XMHOLE*(VLCPH**2+VLCTH**2)/4.D0)
c               DEELL=R(I)/GRAV/XMHOLE*VMEAN/DEN
c               TOUT=TCR*XMHOLE/DEXP(X(IMR,I))*5.D0/6.D0/ELL*DEELL
c            ELSE
c               TOUT=TCR*5.D0/6.D0/ELL
c            END IF
C
C***  XCF=SQRT OF RATIO OF TOUT/TIN**XCF1=SAME FOR NTR=1*************
            XCF1=DSQRT(TOUT*XBETA/TRX(JC,JC)/DOMEG1)
            XCF=XCF1/DSQRT(YALPHA)
            IF(XCF.LT.1.D0)ICRIT=J
            IF(XCF1.LT.1.D0)ICRIT1=J
            ERARG=2.D0/XCF
            IF(ERARG.LT.1.D1)THEN
               EARG=ERARG
               PF=DERF(EARG)
            ELSE
               PF=1.D0
            END IF
C
            XK=1.D0+TOUT*XBETA/YALPHA/TRX(JC,JC)/DOMEG1/PF
            ERARG=2.D0/XCF1
            IF(ERARG.LT.1.D1)THEN
               EARG=ERARG
               PF1=DERF(EARG)
            ELSE
               PF1=1.D0
            END IF
C***  Update of lc-filling-factor in YLC; storage on Y2(I) later***
            IF(IMOD.GT.0)THEN
               XEXP=DEXP(-YALPHA*PF*XK/TOUT/BREM)
               YLC2(JC,J)=YLC2(JC,J)*XEXP+XCF*XCF/PF/XK*(1.D0-XEXP)
               IF(YLC2(JC,J).GT.1.D0)YLC2(JC,J)=1.D0
               IF(YLC2(JC,J).LE.0.D0)YLC2(JC,J)=1.D-30
            END IF
C
            DRHO=Y2(JC,J)*DEXP(X(I00+JC,J))*(DOMEG1-DOMEGA)/TOUT
            DRHO1=PF1*DRHO
            DRHO=YALPHA*PF*DRHO
C     END IF
C**** DRHO1: Strike rate for one crossing, DRHO: trapping rate for ntr**
            DMDT=DRHO*DV
            THLCTH=DATAN(VLCTH*SQRT(2.D0*SIGT2)/SIGR)
            THLCPH=DATAN(VLCPH*SQRT(2.D0*SIGT2)/SIGR)
C***  THLC=Frank/Rees-Model, THLCO: Velocity Integration with SB-distr.***
            THD=DSQRT(TOUT/TRX(JC,JC)/YALPHA)
            TIN=DOMEG1*TRX(JC,JC)
            XER1=EPSR1/DOMEG1
            XET1=EPST1/DOMEG1
C*****Begin information printing block**********************************
 206        FORMAT(1X,' **(JC)**#1:I= ','#2:K(T)=**** ','#3:MDOTLC=** ',
     *           '#4:XER1=**** ','#5:XET1=**** ','#6:PF(TOUT)= ',
     *           '#7:TH-D=**** ','#8:THLCTH=** ',
     *           '#9:THLCPH=** ','#10:VLCTH/ST= ','#11:VLCPH/ST= ',
     *           ' #12:TIN=**** ',' #13:TOUT=*** ','#14:DOMEG1=** ',
     *           '#15:ERARG=** ', '#16:R(J)=** ')
            IF((MOD(NMOD,NCOL)+MOD(J,NWRITE)).EQ.0)WRITE(6,207)
     *           JC,J,Y2(JC,J),DMDT,XER1,XET1,PF,
     *           THD,THLCTH,THLCPH,VLCTH,VLCPH,TIN,TOUT,DOMEG1,
     *           ERARG,R(J)
 207        FORMAT(1X,I4,I4,1P,16(' ',D11.4))
C*******End information printing block**********************************
C
 259     CONTINUE
C
         IF(IMOD.EQ.0)RETURN
C
C*****Store new filling of loss-cone************************************
         DO 2559 I=2,NJ
            Y2(JC,I)=YLC2(JC,I)
 2559    continue
C
 2002 CONTINUE
C
      RETURN
C*******************************************************************
      ENTRY NHOLE
C
C     CONSTRUCT LOADED POLYTROPE OF INDEX NOL-20 FOR START MODEL
C
      LTOT=.FALSE.
*
      DO 1003 J=1,NCOMP
*     Construct loaded polytrope only if NOL(J)>20 selected
         LOAD(J)=.FALSE.
         IF(NOL(J).GT.20)THEN
            LOAD(J)=.TRUE.
            NOL(J)=NOL(J)-20
         END IF
         LTOT=LTOT.OR.LOAD(J)
*
         IF(.NOT.LOAD(J))GOTO 1003
*
         KL2=0
 400     KL2=KL2+1
         IF(DEXP(X(IMR,KL2)).LT.XMHOLE)GOTO 400
         PRINT*,' KL2=',KL2
         XSC=X(I00+J,KL2)
         SIGSC=DEXP(X(I20+J,KL2)-XSC)
         GAM=1.D0+1.D0/DBLE(NOL(J))
         FACG=GAM/(GAM+1.D0)
         GEXP=1.D0/GAM
*
         RFAC=R(3)/R(2)
         VMKL2=PI43*(R(2)/RFAC)**3*DEXP(X(I00+J,1))
         DO 399 I=2,KL2
            RIM=R(I-1)
            IF(I.EQ.2)RIM=R(2)/RFAC
            VMKL2=VMKL2+PI43*(R(I)**3-RIM**3)*DEXP(X(I00+J,I))
 399     CONTINUE
C
         X(I00+J,1)=XSC+GEXP*DLOG(1.D0+FACG*R(KL2)/R(2)*RFAC)
         DLTM=PI43*(R(2)/RFAC)**3*DEXP(X(I00+J,1))
         DO 401 I=2,KL2
            RIM=R(I-1)
            IF(I.EQ.2)RIM=R(2)/RFAC
            X(I00+J,I)=XSC+GEXP*DLOG(1.D0+FACG*R(KL2)/R(I))
            DLTM=DLTM+PI43*(R(I)**3-RIM**3)*DEXP(X(I00+J,I))
 401     CONTINUE
C
         DEM=(DLTM-VMKL2+XMTOT(J))/XMTOT(J)
         PRINT*,' J=',J,' DLTM=',DLTM,'VMKL2=',VMKL2,' DEM=',DEM
C
         DO 402 I=1,NJ
            X(I00+J,I)=X(I00+J,I)-DLOG(DEM)
 402     CONTINUE
C
 1003 CONTINUE
C
      IF(.NOT.LTOT)RETURN
*
      RHOSUM=0.D0
      DO 251 J=1,NCOMP
         RHOSUM=RHOSUM+DEXP(X(I00+J,1))
 251  CONTINUE
*
      X(IMR,1)=DLOG(PI43)+DLOG(RHOSUM)+3.D0*DLOG(R(2)/RFAC)
*
      IMMAX=0
      DO 450 I=2,NJ
         RHOSUM=0.D0
         DO 250 J=1,NCOMP
            RHOSUM=RHOSUM+DEXP(X(I00+J,I))
 250     CONTINUE
*
         IF(I.EQ.2)THEN
            DER=DLOG(RFAC)
            RAV=C12*DLOG(R(I)*R(I)/RFAC)
            RM=R(I)/RFAC
         ELSE
            DER=DLOG(R(I)/R(I-1))
            RAV=C12*DLOG(R(I)*R(I-1))
            RM=R(IM)
         END IF
*
         IMASS=0
         AMAX=0.D0
 135     IMASS=IMASS+1
         IF(IMASS.EQ.1)THEN
            EFAC=DEXP(3.D0*DLOG(RM)-X(IMR,I-1))
            VMASS=0.D0
         ELSE
            XMAV=C12*(X(IMR,I)+X(IMR,I-1))
            EFAC=DEXP(3.D0*RAV-XMAV)
            VMASS=X(IMR,I)
         END IF

         X(IMR,I)=X(IMR,I-1)+PI4*RHOSUM*DER*EFAC
         XMCORR=X(IMR,I)-VMASS
         IF(IMASS.GT.100)THEN
            PRINT*,' Failure of initial total mass iteration at:'
            PRINT*,' I=',I
            STOP
         END IF
*
         IF(DABS(XMCORR).GT.EPS)GOTO 135
         IF(IMASS.GT.IMMAX)THEN
            IMMAX=IMASS
            KMAX=I
         END IF
*
 450  CONTINUE
*
      PRINT*,' NHOLE: Comp. ',J,' Fin Max IMASS=',IMMAX,
     *     ' at mesh ',KMAX
*
      RETURN
C
      END
