       SUBROUTINE I2IEQ
       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
*        Equation for dlog pr /dt
*
         INCLUDE 'compar.f'
*
      DO 1000 J=1,NCOMP
      DPDR=(PAVR(J)-PAVRM(J))/DER
      VELO=UAV(J)+3.D0*VRAV(J)
*
      G(I20+J)=BREM*(X(I20+J,I)-VX(I20+J,I))+VELO/RAVI*DPDR+
     *     3.D0/RAVI*(HX(I10+J)-HXM(I10+J)+HX(I30+J)-HXM(I30+J))/DER+
     *     2.D0/RAVI*(UAV(J)+VRAV(J))
*
      D(I20+J,I20+J)=BREM/XIMPL
      C(I20+J,I20+J)=-C12*VELO/RAVI/DER
*
      D(I20+J,I10+J)=C12/RAVI*DPDR+3.D0/RAVI/DER+1.D0/RAVI
      C(I20+J,I10+J)=C12/RAVI*DPDR-3.D0/RAVI/DER+1.D0/RAVI
*
*        Switch off fluxes if LS(9)=.TRUE.
      IF(.NOT.LS(9))THEN
      D(I20+J,I30+J)=C32/RAVI*DPDR+3.D0/RAVI/DER+1.D0/RAVI
      C(I20+J,I30+J)=C32/RAVI*DPDR-3.D0/RAVI/DER+1.D0/RAVI
      END IF
*
      IF(.NOT.LNJ)E(I20+J,I20+J)=C12*VELO/RAVI/DER
*
*         Isotropy
*
      G(I02+J)=(HX(I02+J)-HX(I20+J))
      D(I02+J,I02+J)=1.D0
      D(I02+J,I20+J)=-1.D0
*
*        Stellar Escapers
*
      IF(LS(11))THEN
      G(I20+J)=G(I20+J)+DERESC(J)
      END IF
*
*        Optional choice of trivial equation
*
      G(I20+J)=FX(I20+J)*G(I20+J)+
     * (1.D0-FX(I20+J))*BREM*(X(I20+J,I)-VX(I20+J,I))
      DO 1001 K=1,NG
      C(I20+J,K)=FX(I20+J)*C(I20+J,K)
      D(I20+J,K)=FX(I20+J)*D(I20+J,K)
      E(I20+J,K)=FX(I20+J)*E(I20+J,K)
 1001 CONTINUE
      D(I20+J,I20+J)=FX(I20+J)*D(I20+J,I20+J)+
     *  (1.D0-FX(I20+J))*BREM
*
*        Diagnostic output
*
       IF(LS(25))THEN
       PRINT*,' COMP. ',J,' DPDT=',
     *    BREM*(X(I20+J,I)-VX(I20+J,I)),
     * ' U/R*DPDR=',UAV(J)/RAVI*DPDR,' 3*VR/R*DPDR=',
     *  3.D0*VRAV(J)/RAVI*DPDR,' 3/R*DUDR=',
     *  3.D0/RAVI/DER*(HX(I10+J)-HXM(I10+J)),' 3/R*DVRDR=',
     *  3.D0/RAVI/DER*(HX(I30+J)-HXM(I30+J)),' 2*U/R=',
     *  2.D0*UAV(J)/RAVI,' 2*VR/R=',2.D0*VRAV(J)/RAVI,
     *  ' COLLTERM=',CT,
     * ' DEESC=',ESTERM
       END IF
*
 1000 CONTINUE
*
          RETURN
*
          END
