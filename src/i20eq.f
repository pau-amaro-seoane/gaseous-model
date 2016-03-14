       SUBROUTINE I20EQ
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
     *     2.D0/RAVI*(UAV(J)+3.D0*VRAV(J)-2.D0*VTAV(J)*AFAC(J))
*
      D(I20+J,I20+J)=BREM/XIMPL+4.D0*VTAV(J)/RAVI*AFAC(J)
      D(I20+J,I02+J)=-4.D0*VTAV(J)/RAVI*AFAC(J)
      C(I20+J,I20+J)=-C12*VELO/RAVI/DER
      IF(.NOT.LNJ)E(I20+J,I20+J)=C12*VELO/RAVI/DER
*
      D(I20+J,I10+J)=C12/RAVI*DPDR+3.D0/RAVI/DER+1.D0/RAVI
      C(I20+J,I10+J)=C12/RAVI*DPDR-3.D0/RAVI/DER+1.D0/RAVI
*
*        Switch off fluxes if LS(9)=.TRUE.
      IF(.NOT.LS(9))THEN
      D(I20+J,I30+J)=C32/RAVI*DPDR+3.D0/RAVI/DER+3.D0/RAVI
      C(I20+J,I30+J)=C32/RAVI*DPDR-3.D0/RAVI/DER+3.D0/RAVI
      D(I20+J,I12+J)=-2.D0/RAVI*AFAC(J)
      C(I20+J,I12+J)=-2.D0/RAVI*AFAC(J)
      END IF
*
*        Relaxation terms
      RXT=TRX(J,J)
*
      IF(.NOT.LS(15))THEN
*
*     VAFAC = DEXP(VX(I02+J,I)-VX(I20+J,I))
      XFTAN = 1.D0
*     IF(VAFAC.GT.2.D0)XFTAN=1.D0/(1.D0+(VAFAC-2.D0)**2)
*
      CT=6.D-1/XTAN/RXT*(1.D0-AFAC(J))*XFTAN
*
      G(I20+J)=G(I20+J)+CT
      D(I20+J,I20+J)=D(I20+J,I20+J)+6.D-1/XTAN/RXT*AFAC(J)+
     *  3.D0*CT*AFAC(J)/(1.D0+2.D0*AFAC(J))-C32*CT
      D(I20+J,I02+J)=D(I20+J,I02+J)-6.D-1/XTAN/RXT*AFAC(J)-
     *  3.D0*CT*AFAC(J)/(1.D0+2.D0*AFAC(J)) 
      D(I20+J,I00+J)=D(I20+J,I00+J)+25.D-1*CT
C
      END IF
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
     *  2.D0*UAV(J)/RAVI,' 6*VR/R=',6.D0*VRAV(J)/RAVI,
     *  ' -4*VT*PT/PR/R=',-4.D0*VTAV(J)/RAVI*AFAC(J),
     *  ' COLLTERM=',CT,' RXT=',RXT,
     * ' DERESC=',DERESC(J)
      PRINT*,' D20/20=',6.D-1/XTAN/RXT*AFAC(J)+
     *  3.D0*CT*AFAC(J)/(1.D0+2.D0*AFAC(J))-C32*CT
      PRINT*,' D20/02=',-6.D-1/XTAN/RXT*AFAC(J)-
     *  3.D0*CT*AFAC(J)/(1.D0+2.D0*AFAC(J))
      PRINT*,' D20/00=',+25.D-1*CT
       END IF
*
 1000 CONTINUE
*
          RETURN
*
          END
