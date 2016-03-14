          SUBROUTINE I12EQ
       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
*        Equation for dlog(pt/rho)/dlogr
*
         INCLUDE 'compar.f'
*
      C14=0.25D0
      C34=0.75D0
*      Different Lambda for rad./tang. energy transport if XLFAC.NE.0.D0
      CLAM=PI4*GRAV/CS(21)*(1.D0+XLFAC)
      C1=XMHOLE/DEXP(HX(IMR))
      DMRC1=-C1
*
      DO 1000 J=1,NCOMP
*
      AFACAV=DEXP(PAVR(J)-PAVT(J))
      XLAINV=CLAM*(1.D0+C1)*(2.D0+AFACAV)**C32*
     *    DEXP(XLR+C12*PAVT(J)-C12*XAVR(J)+CRX-
     *    DLOG(XMIND(J))-DLOG(XCOUL(J)))
      DSIGT2=(HXP(I02+J)-HX(I02+J)-HXP(I00+J)+HX(I00+J))/DERAV
      FTERM=XLAINV*HX(I12+J)
*
      G(I12+J)=DSIGT2 + FTERM
*
      D(I12+J,I12+J)=XLAINV
      D(I12+J,IMR)=FTERM*DMRC1/(C1+1.D0)
*
      D(I12+J,I02+J)=-1.D0/DERAV-C34*FTERM*AFACAV/(2.D0+AFACAV)+
     *    C14*FTERM
      D(I12+J,I20+J)=C34*FTERM*AFACAV/(2.D0+AFACAV)
      D(I12+J,I00+J)=1.D0/DERAV-C14*FTERM
*
      IF(.NOT.LNJ)THEN
      E(I12+J,I02+J)=1.D0/DERAV-C34*FTERM*AFACAV/(2.D0+AFACAV)+
     *    C14*FTERM
      E(I12+J,I20+J)=C34*FTERM*AFACAV/(2.D0+AFACAV)
      E(I12+J,I00+J)=-1.D0/DERAV-C14*FTERM
      END IF
*
*        Diagnostic output
*
      IF(LS(25))THEN
      PRINT*,' I12EQ Comp. ',J,' DSIGT2=',DSIGT2,
     *  ' FTERM=',FTERM,' XLAINV=',XLAINV
      END IF
*
 1000    CONTINUE
*
          RETURN
*
          END
