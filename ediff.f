      SUBROUTINE EDIFF
*     Star accretion due to energy diffusion near tidal radius
*     Distribution of close encounter energy with black hole in core
*     Loss of mass and energy due to gravitational radiation
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*     Only called for second grid point except if
*     gravitational radiation term is present or core energy distribution
      INCLUDE 'compar.f'
C     
*     Distribution of energy diffusion energy partly in core
      IF(R(I).GT.RGRAV.AND.XBIN.NE.0.D0)THEN
         DO 1001 J=1,NCOMP
            CTERM=-C23*XBINY(J)*DEXP(HX(I00+J)*AEXP)
            G(I20+J)=G(I20+J)+CTERM/DEXP(HX(I20+J))
            G(I02+J)=G(I02+J)+CTERM/DEXP(HX(I02+J))
            D(I20+J,I20+J)=D(I20+J,I20+J)-CTERM/DEXP(HX(I20+J))
            D(I02+J,I02+J)=D(I02+J,I02+J)-CTERM/DEXP(HX(I02+J))
            D(I20+J,I00+J)=D(I20+J,I00+J)+AEXP*CTERM/DEXP(HX(I20+J))
            D(I02+J,I00+J)=D(I02+J,I00+J)+AEXP*CTERM/DEXP(HX(I02+J))
            IF(LS(25))PRINT*,' EDIFF distributed J=',J,' G20=',
     *           CTERM/DEXP(HX(I20+J)),
     *           ' G02=',CTERM/DEXP(HX(I02+J))
 1001    CONTINUE
      END IF
*     
      IF(R(I).LT.RGRAV.AND.LS(20))THEN
*     
         DO 1002 J=1,NCOMP
            TGRAV=CS(25)*(R(I)/XMIND(J))**2*
     *           R(I)**2/XMHOLE/(1.D0+XMHOLE/XMIND(J))
            G(I00+J)=G(I00+J)+1.D0/TGRAV
            G(I20+J)=G(I20+J)+1.D0/TGRAV
            G(I02+J)=G(I02+J)+1.D0/TGRAV
 1002    CONTINUE
*     
      END IF
      IF(.NOT.L2)RETURN
*     
      DO 1000 J=1,NCOMP
C     
         DM=1.D0/TRX(J,J)
C     
         G(I00+J) = G(I00+J) + DM
C     
         D(I00+J,I00+J)=D(I00+J,I00+J)+2.5D0*DM
         IF(LS(13))THEN
            D(I00+J,I20+J)=D(I00+J,I20+J)-C32*DM
         ELSE
            D(I00+J,I20+J)=D(I00+J,I20+J)
     *           +DM*(-C32+3.D0*AFAC(J)/(1.D0+2.D0*AFAC(J)))
            D(I00+J,I02+J)=D(I00+J,I02+J)
     *           -DM*3.D0*AFAC(J)/(1.D0+2.D0*AFAC(J))
         END IF
C     
C-----------Energy equations
C     
         SIGR2=DEXP(HX(I20+J)-HX(I00+J))
         XLPT=C23*GRAV*XMHOLE/RTIDE(J)/SIGR2
*     if energy is distributed in core it must be taken from here
         IF(AEXP.GT.0.D0)XLPT=(1.D0-XBIN)*XLPT
         DENR=DM*(1.D0-XLPT)
*     
         G(I20+J)=G(I20+J)+DENR
C     
         D(I20+J,I00+J)=D(I20+J,I00+J)+2.5D0*DENR-DM*XLPT
         D(I20+J,I20+J)=D(I20+J,I20+J)-C32*DENR+DM*XLPT
*     
         IF(.NOT.LS(13))THEN
            D(I20+J,I20+J)=D(I20+J,I20+J)
     *           +3.D0*DENR*AFAC(J)/(1.D0+2.D0*AFAC(J))
            D(I20+J,I02+J)=D(I20+J,I02+J)
     *           -3.D0*DENR*AFAC(J)/(1.D0+2.D0*AFAC(J))
*     
            SIGT2=DEXP(HX(I02+J)-HX(I00+J))
            XLPT=C23*GRAV*XMHOLE/RTIDE(J)/SIGT2
            DENT=DM*(1.D0-XLPT)
*     
            G(I02+J)=G(I02+J)+DENT
C     
            D(I02+J,I00+J)=D(I02+J,I00+J)+2.5D0*DENT-DM*XLPT
            D(I02+J,I02+J)=D(I02+J,I02+J)+DM*XLPT
     *           -3.D0*DENT*AFAC(J)/(1.D0+2.D0*AFAC(J))
            D(I02+J,I20+J)=D(I02+J,I20+J)-C32*DENR
     *           +3.D0*DENT*AFAC(J)/(1.D0+2.D0*AFAC(J))
*     
         END IF
*     
*     Diagnostic printout
         IF(LS(25))PRINT*,'  EDIFF: COMP-',J,' G00=',DM,
     *        ' G20=',DENR,' G02=',DENT
*     
 1000 CONTINUE
C     
      RETURN
C     
      END
