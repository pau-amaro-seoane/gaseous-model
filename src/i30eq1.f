          SUBROUTINE I30EQ1
       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
*        Equation for dlog((pr+pt)/rho)/dlogr
*
         INCLUDE 'compar.f'
*
      C14=0.25D0
      CLAM=1.8D0*PI4*GRAV/CS(21)
*
      C1=XMHOLE/DEXP(HX(IMR))
      DMRC1=-C1
*
      DO 1000 J=1,NCOMP
*       Take care for isotropic calculation
      IF(LS(13))THEN
      XLAINV=CLAM*(1.D0+C1)/(3.D0)**C12*5.D0*
     *    DEXP(XLR+C12*PAVR(J)-C12*XAVR(J)+CRX-
     *    DLOG(XMIND(J))-DLOG(XCOUL(J)))
      DSIG2=(HXP(I20+J)-HX(I20+J)-HXP(I00+J)+HX(I00+J))/DERAV
      ELSE
*
      AFACP=DEXP(HXP(I02+J)-HXP(I20+J))
      AFACAV=DEXP(PAVT(J)-PAVR(J))
      XLAINV=CLAM*(1.D0+C1)/(1.D0+2.D0*AFACAV)**C12*
     *                      (3.D0+2.D0*AFACAV)*
     *    DEXP(XLR+C12*PAVR(J)-C12*XAVR(J)+CRX-
     *    DLOG(XMIND(J))-DLOG(XCOUL(J)))
      DSIG2=(HXP(I20+J)-HX(I20+J)-HXP(I00+J)+HX(I00+J)+
     *   DLOG(1.D0+2.D0*AFACP)-DLOG(1.D0+2.D0*AFAC(J)))/DERAV
      END IF
*
      FTERM=XLAINV*HX(I30+J)
*
      G(I30+J)=DSIG2 + FTERM
*
      D(I30+J,I30+J)=XLAINV
      D(I30+J,IMR)=FTERM*DMRC1/(C1+1.D0)
*
      D(I30+J,I20+J)=-1.D0/DERAV+C14*FTERM
      D(I30+J,I00+J)=1.D0/DERAV-C14*FTERM
      IF(.NOT.LNJ)THEN
      E(I30+J,I20+J)=1.D0/DERAV+C14*FTERM
      E(I30+J,I00+J)=-1.D0/DERAV-C14*FTERM
      END IF
*
      IF(.NOT.LS(13))THEN
      D(I30+J,I20+J)=D(I30+J,I20+J)+
     *    2.D0/DERAV*AFAC(J)/(1.D0+2.D0*AFAC(J))+
     *    C12*FTERM*AFACAV/(1.D0+2.D0*AFACAV)-
     *    FTERM*AFACAV/(3.D0+2.D0*AFACAV)
      D(I30+J,I02+J)=-2.D0/DERAV*AFAC(J)/(1.D0+2.D0*AFAC(J))-
     *    C12*FTERM*AFACAV/(1.D0+2.D0*AFACAV)+
     *    FTERM*AFACAV/(3.D0+2.D0*AFACAV)
*
      IF(.NOT.LNJ)THEN
      E(I30+J,I20+J)=E(I30+J,I20+J)-2.D0/DERAV*AFACP/(1.D0+2.D0*AFACP)+
     *    C12*FTERM*AFACAV/(1.D0+2.D0*AFACAV)-
     *    FTERM*AFACAV/(3.D0+2.D0*AFACAV)
      E(I30+J,I02+J)=2.D0/DERAV*AFACP/(1.D0+2.D0*AFACP)-
     *    C12*FTERM*AFACAV/(1.D0+2.D0*AFACAV)+
     *    FTERM*AFACAV/(3.D0+2.D0*AFACAV)
      END IF
      END IF
*
*        Transport velocities equal
*
      G(I12+J)=BREM*(HX(I12+J)-HX(I30+J))
      D(I12+J,I12+J)=BREM
      D(I12+J,I30+J)=-BREM
*
*        Diagnostic output
*
      IF(LS(25))THEN
      PRINT*,' I30EQ1 Comp. ',J,' DSIG2=',DSIG2,
     *  ' FTERM=',FTERM,' XLAINV=',XLAINV
      END IF
*
 1000    CONTINUE
*
          RETURN
*
          END
