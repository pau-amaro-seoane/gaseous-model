      SUBROUTINE FORMB3
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
      INCLUDE 'compar.f'
      COMMON/EBTEST/EBALA,EBEQ,EBEQ2,DEMASS,DEENER,DEMASS2
      DIMENSION DYB00(NCOMPO),DYB20(NCOMPO),DYB02(NCOMPO)
*
*       Energy generation by three body binary formation and hardening
*       Take Lee et al. formula
*       Possible radial/tangential difference in energy generation
*         by using factor AEXP
*
      C52=2.5D0
      DYB=0.0D0
      RHOT=0.D0
*
      DO 901 J = 1,NCOMP
*
      RHOT = RHOT + DEXP(HX(I00+J))
* 
      CBTERM = DLOG(XMIND(J))
      DBTERM = DEXP(CBTERM+C52*HX(I00+J)-C32*HX(I20+J))/
     *            (1.D0+2.D0*AFAC(J))**C32
      DYB = DYB + DBTERM
      DYB00(J)=C52*DBTERM
      DYB20(J)=-C32*DBTERM
*
      IF(.NOT.LS(13))THEN 
      DYB20(J)=DYB20(J)+3.D0*DBTERM*AFAC(J)/(1.D0+2.D0*AFAC(J))
      DYB02(J)=-3.D0*DBTERM*AFAC(J)/(1.D0+2.D0*AFAC(J))
      ELSE
      DYB02(J)=0.D0
      END IF
*
 901  CONTINUE
*
      DO 1010 J = 1,NCOMP
*      Take corrected XBIN-value
      CBIN=DLOG(18.D0*DSQRT(3.D0)*XBINY(J))+5.D0*DLOG(GRAV)
      CBR=CBIN+AEXP*VX(I00+J,1)
      CBT=CBIN-AEXP*VX(I00+J,1)
*
* Radial energy generation (DYA determined in GINIT central sigma)
*
      C3 = DEXP(CBR+(1.D0-AEXP)*HX(I00+J)-HX(I20+J))*DYA*DYB**3/RHOT
*
      G(I20+J) = G(I20+J) - C3
*
      DER3=PI43*(R(I)**3-R(I-1)**3)
      EBXX = C3*DEXP(HX(I20+J))*DER3/BREM
      EBEQ2 = EBEQ2 + C3*C32*DEXP(HX(I20+J))*DER3/BREM
*     
      D(I20+J,I00+J) = D(I20+J,I00+J)-(1.D0-AEXP)*C3
      D(I20+J,I20+J) = D(I20+J,I20+J)+C3
*
      DO 1011 K = 1,NCOMP
*
      D(I20+J,I00+K) = D(I20+J,I00+K) + 
     *  C3*DEXP(HX(I00+K))/RHOT  - 3.D0*C3/DYB*DYB00(K)
      D(I20+J,I20+K) = D(I20+J,I20+K) - 3.D0*C3/DYB*DYB20(K) 
      D(I20+J,I02+K) = D(I20+J,I02+K) - 3.D0*C3/DYB*DYB02(K)
*
 1011 CONTINUE
*
*        Tangential energy generation
*
      IF(.NOT.LS(13))THEN
*
      CT3 = DEXP(CBT+(1.D0+AEXP)*HX(I00+J)-HX(I02+J))*DYA*DYB**3/RHOT
*
      EBXX = C12*(EBXX + 2.D0*CT3*DEXP(HX(I02+J))*DER3/BREM)
      EBEQ = EBEQ + EBXX
*
      G(I02+J) = G(I02+J) - CT3
*
      D(I02+J,I00+J) = D(I02+J,I00+J)-(1.D0+AEXP)*CT3
      D(I02+J,I02+J) = D(I02+J,I02+J)+CT3
*
      DO 1012 K = 1,NCOMP
*
      D(I02+J,I00+K) = D(I02+J,I00+K) +
     *  CT3*DEXP(HX(I00+K))/RHOT - 3.D0*CT3/DYB*DYB00(K)
      D(I02+J,I20+K) = D(I02+J,I20+K) - 3.D0*CT3/DYB*DYB20(K)
      D(I02+J,I02+K) = D(I02+J,I02+K) - 3.D0*CT3/DYB*DYB02(K)
*
 1012 CONTINUE
*
      END IF
*
*       Diagnostic output
*
      IF(LS(25))THEN
      PRINT*,' FORMB3 LEE: J=',J,' C3 (RAD) = ',C3,
     *' CT3 (TANG) = ',CT3
      END IF
*
 1010 CONTINUE
*
      RETURN
*
      END
