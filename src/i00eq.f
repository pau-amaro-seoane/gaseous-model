       SUBROUTINE I00EQ
       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
*        Equation for dlog rho /dt
*
         INCLUDE 'compar.f'
*
      DO 1000 J=1,NCOMP
*
      DXDR=(XAVR(J)-XAVRM(J))/DER
      DUDR=(HX(I10+J)-HXM(I10+J))/DER
*
      G(I00+J)=BREM*(X(I00+J,I)-VX(I00+J,I))+
     * (UAV(J)*DXDR+DUDR+2.D0*UAV(J))/RAVI
*
      D(I00+J,I00+J)=BREM/XIMPL
      IF(.NOT.LNJ)E(I00+J,I00+J)=C12*UAV(J)/DER/RAVI
      C(I00+J,I00+J)=-C12*UAV(J)/DER/RAVI
*
      D(I00+J,I10+J)=1.D0/DER/RAVI+C12*DXDR/RAVI+1.D0/RAVI
      C(I00+J,I10+J)=-1.D0/DER/RAVI+C12*DXDR/RAVI+1.D0/RAVI
*
*        Stellar escapers
*
      IF(LS(11))THEN
      G(I00+J)=G(I00+J)+DRHESC(J)
      END IF
*
*       Optional choice of trivial equation
*
      G(I00+J)=FX(I00+J)*G(I00+J)+
     *  (1.D0-FX(I00+J))*BREM*(X(I00+J,I)-VX(I00+J,I))
      DO 1001 K=1,NG
      C(I00+J,K)=FX(I00+J)*C(I00+J,K)
      D(I00+J,K)=FX(I00+J)*D(I00+J,K)
      E(I00+J,K)=FX(I00+J)*E(I00+J,K)
 1001 CONTINUE
      D(I00+J,I00+J)=D(I00+J,I00+J)*FX(I00+J)+
     *  (1.D0-FX(I00+J))*BREM
*
*       Diagnostic output
*
      IF(LS(25))THEN
      PRINT*,' COMP-',J,' DLGRHODT=',
     *  BREM*(X(I00+J,I)-VX(I00+J,I)),
     *  ' DXDR-T.=',UAV(J)*DXDR/RAVI,
     *  ' DUDR-T.=',DUDR/RAVI,
     *  ' 2*U/R=',2.D0*UAV(J)/RAVI,
     *  ' DRHESC=',DRHESC(J)
      END IF
*
 1000 CONTINUE
*
          RETURN
*
          END
