          SUBROUTINE I10EQ
       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
*        Equation for (rho*r/pr)*du/dt
*
         INCLUDE 'compar.f'
*
      IF(IZ00.EQ.0)THEN
      IZ00=1
      PRINT*,' I10EQ: Upstream/downstream differencing '
      END IF
*
      DO 1000 J=1,NCOMP
*
      CFAC=1.D0
      IF(LS(17))CFAC=0.D0
*
      UFAC=CFAC*DEXP(XAVR(J)+XLR-PAVR(J))
      XMFAC=GRAV*DEXP(HX(IMR)+XAVR(J)-XLR-PAVR(J))
      XMHFAC=0.D0
      IF(LS(19))XMHFAC=GRAV*XMHOLE*DEXP(XAVR(J)-XLR-PAVR(J))
*         Use slowly varying mass
      IF(LS(3))XMHFAC=XMHFAC+GRAV*VMRBIN(I)*DEXP(XAVR(J)-XLR-PAVR(J))
      DPDR=(HXP(I20+J)-HX(I20+J))/DERAV
      DUDR=FP(J,I)*(HX(I10+J)-HXM(I10+J))/DER+
     *   FM(J,I)*(HXP(I10+J)-HX(I10+J))/DERP
      DUDT=BREM*(X(I10+J,I)-VX(I10+J,I))+HX(I10+J)/R(I)*DUDR
*
      G(I10+J)=UFAC*DUDT + XMFAC + XMHFAC + DPDR
*
      D(I10+J,I10+J)=UFAC*(BREM/XIMPL+DUDR/R(I)+HX(I10+J)/R(I)*
     *  (FP(J,I)/DER-FM(J,I)/DERP))
      C(I10+J,I10+J)=-UFAC*HX(I10+J)/R(I)*FP(J,I)/DER
      D(I10+J,I20+J)=-1.D0/DERAV-C12*(UFAC*DUDT+XMFAC+XMHFAC)
      D(I10+J,I00+J)=C12*(UFAC*DUDT+XMFAC+XMHFAC)
      D(I10+J,IMR)=XMFAC
*
      IF(.NOT.LNJ)THEN
      E(I10+J,I10+J)=UFAC*HX(I10+J)/R(I)*FM(J,I)/DERP
      E(I10+J,I20+J)=1.D0/DERAV-C12*(UFAC*DUDT+XMFAC+XMHFAC)
      E(I10+J,I00+J)=C12*(UFAC*DUDT+XMFAC+XMHFAC)
      END IF
*
*        Anisotropy dependent terms
      IF(.NOT.LS(13))THEN
*
      AFACAV=DEXP(PAVT(J)-PAVR(J))
*
      G(I10+J)=G(I10+J)+2.D0*(1.D0-AFACAV)
*
      D(I10+J,I20+J)=D(I10+J,I20+J)+AFACAV
      D(I10+J,I02+J)=-AFACAV
*
        IF(.NOT.LNJ)THEN
        E(I10+J,I20+J)=E(I10+J,I20+J)+AFACAV
        E(I10+J,I02+J)=-AFACAV
        END IF
*
      END IF
*
*        Diagnostic output
*
      IF(LS(25))THEN
      PRINT*,' STARS: DUDT=',BREM*(X(I10+J,I)-VX(I10+J,I)),
     *  ' ADV.TERM=',HX(I10+J)/R(I)*DUDR,
     *  ' UFAC*DUDT=',UFAC*DUDT,' XMFAC=',XMFAC,' DPDR=',
     *  DPDR,' XMHFAC=',XMHFAC
      IF(.NOT.LS(13))PRINT*,' ANISO=',2.D0*(1.D0-AFACAV)
      END IF
*
*        Optional choice of trivial equation
*
      G(I10+J)=FX(I10+J)*G(I10+J)+
     * (1.D0-FX(I10+J))*BREM*(X(I10+J,I)-VX(I10+J,I))
*
      DO 1001 K=1,NG
      C(I10+J,K)=FX(I10+J)*C(I10+J,K)
      D(I10+J,K)=FX(I10+J)*D(I10+J,K)
      IF(.NOT.LNJ)E(I10+J,K)=FX(I10+J)*E(I10+J,K)
 1001 CONTINUE
*
      D(I10+J,I10+J)=FX(I10+J)*D(I10+J,I10+J)+
     *   (1.D0-FX(I10+J))*BREM
*
 1000 CONTINUE
*
          RETURN
*
          END
