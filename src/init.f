          SUBROUTINE INIT
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C        INITIALIZE THE STARTING MODELS FOR EACH COMPONENT
C
C       FIRST THE DENSITY PROFILES ARE CALCULATED (DENS)
C       SECOND THE TOTAL MASS DISTRIBUTION -
C       THIRD VELOCITIES AND VEL.DISPERSIONS ARE CALC.:
C
C      AFTER CONSTRUCTING THE INITIAL MODELS SOME
C      INFORMATION IS DISPLAYED.
C      Logarithmic variables
C
         IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
         INCLUDE 'compar.f'
C
C     INITIAL VALUES OF ALL VARIABLES

      LCUT=.FALSE.
*
      DO 5 JJJ=1,NJ
      DO 55 KKK=1,NG
 55     X(KKK,JJJ)=-30.
      DO 56 MMM=1,NCOMP
 56     Y(MMM,JJJ)=-30.
 5    CONTINUE
C
         DO 100 J=1,NCOMP
C
         IAIM=0
         XMFAC=1.D0
 120     IAIM=IAIM+1
         AMS=XMTOT(J)*XMFAC
         CALL DENS(J,AMS)
C  INTEGRATE MASS AND CORRECT AMS TO GET WISHED MASS XMTOT(I)
*
         FP(J,1)=DLOG(PI43)+X(I00+J,1)+3.D0*DLOG(R(2)/STEU3)
*
         IMMAX=0
         DO 150 I=2,NJ
*
         IF(I.EQ.2)THEN
         DER=DLOG(STEU3)
         RAV=C12*DLOG(R(I)*R(I)/STEU3)
         RM=R(I)/STEU3
         ELSE
         DER=DLOG(R(I)/R(I-1))
         RAV=C12*DLOG(R(I)*R(I-1))
         RM=R(IM)
         END IF
*
         IMASS=0
         AMAX=0.D0
 130     IMASS=IMASS+1
         IF(IMASS.EQ.1)THEN
         EFAC=DEXP(3.D0*DLOG(RM)-FP(J,I-1))
         VMASS=0.D0
         ELSE
         XMAV=C12*(FP(J,I)+FP(J,I-1))
         EFAC=DEXP(3.D0*RAV-XMAV)
         VMASS=FP(J,I)
         END IF
C
       FP(J,I)=FP(J,I-1)+PI4*DEXP(X(I00+J,I))*DER*EFAC
         XMCORR=FP(J,I)-VMASS
         IF(IMASS.GT.100)THEN
         PRINT*,' Failure of initial mass iteration at:'
         PRINT*,' J=',J,' I=',I
         STOP
         END IF
         IF(DABS(XMCORR).GT.EPS)GOTO 130
*
         IF(IMASS.GT.IMMAX)THEN
         IMMAX=IMASS
         KMAX=I
         END IF
*
 150     CONTINUE
         PRINT*,' Comp. ',J,' Now Max IMASS=',IMMAX,' at mesh ',KMAX
*
         XMFAC=XMTOT(J)/DEXP(FP(J,NJ))
C
         IF(DABS(XMFAC-1.D0).GT.EPS)GOTO 120
C
 100     CONTINUE
C
         RHOSUM=0.D0
         DO 251 J=1,NCOMP
 251     RHOSUM=RHOSUM+DEXP(X(I00+J,1))
*
         X(IMR,1)=DLOG(PI43)+DLOG(RHOSUM)+3.D0*DLOG(R(2)/STEU3)
*
         IMMAX=0
         DO 200 I=2,NJ
         RHOSUM=0.D0
         DO 250 J=1,NCOMP
 250     RHOSUM=RHOSUM+DEXP(X(I00+J,I))
*
         IF(I.EQ.2)THEN
         DER=DLOG(STEU3)
         RAV=C12*DLOG(R(I)*R(I)/STEU3)
         RM=R(I)/STEU3
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
*
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
 200     CONTINUE
*
         PRINT*,' Comp. ',J,' Fin Max IMASS=',IMMAX,' at mesh ',KMAX
C
C-----------END OF MASS INTEGRATION-----------------------
         IF(LS(16))CALL NHOLE
*
 330     CONTINUE
*
         IF(LS(5).OR.LS(27))THEN
*          Calculation of gravitational potential
*          FP(), FM() used as temporary storage
*          necessary here in case of primordial binaries
*          Start at outer boundary
         PHI(NJ)=-GRAV*DEXP(X(IMR,NJ))/R(NJ)
*
         DO 340 I=1,NJ-1
 340     PHI(I)=0.D0
*
         DO 350 J=1,NCOMP
         FM(J,1)=0.D0
         DO 360 I=2,NJ
         DER3=R(I)**3-R(I-1)**3
         FM(J,I)=FM(J,I-1)+DEXP(X(I00+J,I))*DER3
 360     CONTINUE
*
         DER3=R(NJ)**3-R(NJ-1)**3
         RAV=(R(NJ)+R(NJ-1))/2.D0
         FP(J,NJ)=DEXP(X(I00+J,NJ))*DER3/RAV
         DO 370 I=NJ-1,2,-1
         DER3=R(I)**3-R(I-1)**3
         RAV=(R(I)+R(I-1))/2.D0
         FP(J,I)=FP(J,I+1)+DEXP(X(I00+J,I))*DER3/RAV
 370     CONTINUE
*
         PHI(1)=PHI(1)-PI43*GRAV*FP(J,2)
         RFAC=R(3)/R(2)
         IF(LS(19))PHI(1)=PHI(1)-GRAV*XMHOLE/R(2)*RFAC
         DO 380 I=2,NJ-1
         RAV=(R(I)+R(I-1))/2.D0
         PHI(I)=PHI(I)-PI43*GRAV*(FM(J,I)/R(I)+FP(J,I+1))
         IF(LS(19))PHI(I)=PHI(I)-GRAV*XMHOLE/R(I)
 380     CONTINUE
*
 350     CONTINUE
*
*       Tidal cutoff
*
*       IF(.NOT.LCUT)THEN
*
*       DO 390 I = NJ,2,-1
*       PHII=-GRAV*(DEXP(X(IMR,I))+XMRBIN(I))/R(I)
*       PRINT*,' I=',I,' PHII,PHTID=',PHII,PHTID
*       IF(PHII.LT.PHTID)GOTO 401
*390    CONTINUE
*
*401    NJSAV = I
*       LCUT = .TRUE.
*       PRINT*,' Tidal Cutoff at NJ=',NJSAV
*
*       GOTO 330
*       END IF
*
        END IF
*
*       Locate primordial binaries
      IF(LS(5))THEN
      DO 450 IB=1,IBIN
      DO 455 I=2,NJ
      IF(RB(IB).GT.R(I))ISH(IB)=I
 455  CONTINUE
      DMBIN(IB,ISH(IB)) = DMBIN(IB,ISH(IB)) + BODY1(IB) + BODY2(IB)
 450  CONTINUE
*
      XMRBIN(1) = 0.D0
      DO 460 I=2,NJ
      DMBT = 0.D0
      DO 465 IB=1,IBIN
 465  DMBT = DMBT + DMBIN(IB,I)
 460  XMRBIN(I) = XMRBIN(I-1) + DMBT
*
      IPOT=0
      CALL PHIMAS(IPOT)
*       Call newpos to determine smooth binary mass distribution via orbit
*
      DO 470 IB=1,IBIN
      DEFANG = 0.D0
      IPT1 = 0
*
      CALL NEWPOS(IB,IPT1,DEFANG)
      CALL PHIMAS(IPOT)
*
         I=IBSTA(IB)
         J=ICO(IB)
         IF(J.NE.1)THEN
         PRINT*,' J is in error IB=',IB,'J=',J
         STOP 'J is in error'
         END IF
         IF(I.LT.1.OR.I.GT.NJ)THEN
         PRINT*,' I is in error IB=',IB,'I=',I,' RBMAX/MIN=',
     &    RBMAX(IB),RBMIN(IB),' v=',VR(IB),VT(IB),' body=',
     &    body1(ib),body2(ib), ' RB= ',RB(IB)
*        STOP 'I is in error'
         END IF
*
*        SIG2=DEXP(X(I20+J,I)-X(I00+J,I))*
*    *             (1.D0+2.D0*DEXP(X(I02+J,I)-X(I20+J,I)))
*        TORB(IB)=(SIG2/3.D0)**C32/
*    *    CS(60)/XMIND(J)/XCOUL(J)/DEXP(X(I00+J,I))
*       Use torb to store trxmin of binary
*     TBIN(IB) = 1.D0/BREM
      TBIN(IB) = TORB(IB)/3.D0
*
      if(mod(ib,300).eq.0.or.ibin.lt.500)then
      print*,' binary ',ib,' initialized, masses=',body1(ib),body2(ib)
      call flush(6)
      end if
*
 470  CONTINUE
*
      END IF
*
      CALL INITST
C
      DO 400 I=1,NG
      DO 410 J=1,NJ
 410  VX(I,J)=X(I,J)
 400  CONTINUE
C
C TO FOLLOW: INFORMATION DISPLAY ABOUT START MODELS
C
*     IF(LS(27))NJ = NJSAV
*
      RETURN
      END
