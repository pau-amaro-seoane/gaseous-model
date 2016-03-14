      SUBROUTINE FORMB2
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
      INCLUDE 'compar.f'
      DIMENSION DYB00(NCOMPO),DYB20(NCOMPO),DYB02(NCOMPO)
      DIMENSION PFAC(NCOMPO,NCOMPO),PSUM(NCOMPO)
      DIMENSION XFAC(NCOMPO,NCOMPO),SRFAC(NCOMPO,NCOMPO)
      DIMENSION DYBA(NCOMPO)
*
*
*       Take care to define DBTERM, DYB if using this routine 
*
*
*       Energy generation by tidal two body binary formation
*       Lee & Ostriker 1986,1993
*
      C52=2.5D0
      C14=0.25D0
      XNTOT=0.0D0
      RHOT=0.D0
*
      DO 800 J=1,NCOMP
      XNTOT = XNTOT + DEXP(VX(I00+J,I))/XMIND(J)
      RHOT = RHOT + DEXP(VX(I00+J,I))
 800  CONTINUE
      XMQ=RHOT/XNTOT
*
      DO 900 J=1,NCOMP
      PSUM(J)=1.D0+2.D0*AFAC(J)
      DO 950 K=1,NCOMP
      PFAC(J,K)=DEXP(HX(I20+K)-HX(I20+J))
*      Note inverse order in PFAC and XFAC
      XFAC(J,K)=DEXP(HX(I00+J)-HX(I00+K))
      SRFAC(J,K)=PFAC(J,K)*XFAC(J,K)
 950  CONTINUE
 900  CONTINUE
*
      DO 901 J = 1,NCOMP
      DYBA(J) = 0.0D0
*
      DO 902 K = 1,NCOMP
*
      VESC2 = 2.D0*GRAV*XMIND(J)/RIND(J)
      DYBA(J) = DYBA(J) + 
     *    C12*DEXP(HX(I00+J)+HX(I20+K))*
     *    CR2B(J,K)*RIND(J)**2/XMIND(J)/XMIND(K)*VESC2**(-ALP2B(J,K))*
     *    (PSUM(J)+SRFAC(J,K)*PSUM(K))**(C12-ALP2B(J,K))
 902  CONTINUE
*       Take care to define DBTERM if using this routine
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
*
* Radial energy generation (DYA determined in GINIT central sigma)
*
      C3 = DEXP(CBIN+HX(I00+J)-HX(I20+J))*DYA*DYB**3/RHOT
*
      G(I20+J) = G(I20+J) - C3
*
      D(I20+J,I00+J) = D(I20+J,I00+J)-C3
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
      CT3 = DEXP(CBIN+HX(I00+J)-HX(I02+J))*DYA*DYB**3/RHOT
*
      G(I02+J) = G(I02+J) - CT3
*
      D(I02+J,I00+J) = D(I02+J,I00+J)-CT3
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
