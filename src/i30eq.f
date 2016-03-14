          SUBROUTINE I30EQ
       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
*        Equation for dlog(pr/rho)/dlogr
*
         INCLUDE 'compar.f'
*
      C14=0.25D0
*      Different Lambda for rad./tang. energy transport if XLFAC.NE.0.D0
      CLAM=PI4*GRAV/CS(21)/(1.D0+XLFAC)
      C1=XMHOLE/DEXP(HX(IMR))
      DMRC1=-C1
*
      DO 1000 J=1,NCOMP
*
      AFACAV=DEXP(PAVT(J)-PAVR(J))
      XLAINV=CLAM*(1.D0+C1)*(1.D0+2.D0*AFACAV)**C32*
     *    DEXP(XLR+C12*PAVR(J)-C12*XAVR(J)+CRX-
     *    DLOG(XMIND(J))-DLOG(XCOUL(J)))
      DSIG2=(HXP(I20+J)-HX(I20+J)-HXP(I00+J)+HX(I00+J))/DERAV
      FTERM=XLAINV*HX(I30+J)
*
      G(I30+J)=DSIG2 + FTERM
*
      D(I30+J,I30+J)=XLAINV
      D(I30+J,IMR)=FTERM*DMRC1/(C1+1.D0)
*
      D(I30+J,I20+J)=-1.D0/DERAV-C32*FTERM*AFACAV/(1.D0+2.D0*AFACAV)+
     *    C14*FTERM
      D(I30+J,I02+J)=C32*FTERM*AFACAV/(1.D0+2.D0*AFACAV)
      D(I30+J,I00+J)=1.D0/DERAV-C14*FTERM
*
      IF(.NOT.LNJ)THEN
      E(I30+J,I20+J)=1.D0/DERAV-C32*FTERM*AFACAV/(1.D0+2.D0*AFACAV)+
     *    C14*FTERM
      IF(.NOT.LS(13))E(I30+J,I02+J)=C32*FTERM*AFACAV/(1.D0+2.D0*AFACAV)
      E(I30+J,I00+J)=-1.D0/DERAV-C14*FTERM
      END IF
*
*        Diagnostic output
*
      IF(LS(25))THEN
      PRINT*,' I30EQ Comp. ',J,' DSIG2=',DSIG2,
     *  ' FTERM=',FTERM,' XLAINV=',XLAINV
      END IF
*
 1000    CONTINUE
*
          RETURN
*
          END
