       SUBROUTINE I02EQ
       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
*        Equation for dlog pt /dt
*
         INCLUDE 'compar.f'
*
      DO 1000 J=1,NCOMP
      DPTDR=(PAVT(J)-PAVTM(J))/DER
      VELO=UAV(J)+VTAV(J)
*
      G(I02+J)=BREM*(X(I02+J,I)-VX(I02+J,I))+VELO/RAVI*DPTDR+
     *     1.D0/RAVI*(HX(I10+J)-HXM(I10+J)+HX(I12+J)-HXM(I12+J))/DER+
     *     4.D0/RAVI*(UAV(J)+VTAV(J))
*
      D(I02+J,I02+J)=BREM/XIMPL
      C(I02+J,I02+J)=-C12*VELO/RAVI/DER
      IF(.NOT.LNJ)E(I02+J,I02+J)=C12*VELO/RAVI/DER
*
      D(I02+J,I10+J)=C12/RAVI*DPTDR+1.D0/RAVI/DER+2.D0/RAVI
      C(I02+J,I10+J)=C12/RAVI*DPTDR-1.D0/RAVI/DER+2.D0/RAVI
*
*        Switch off fluxes if LS(9)=.TRUE.
      IF(.NOT.LS(9))THEN
      D(I02+J,I12+J)=C12/RAVI*DPTDR+1.D0/RAVI/DER+2.D0/RAVI
      C(I02+J,I12+J)=C12/RAVI*DPTDR-1.D0/RAVI/DER+2.D0/RAVI
      END IF
*
*        Relaxation terms
      IF(.NOT.LS(15))THEN
*
      RXT=(2.D0+ATFAC(J))**C32*
     * DEXP(C32*HX(I02+J)-25.D-1*HX(I00+J)+CRX
     *                   -DLOG(XMIND(J))-DLOG(XCOUL(J)))
*
*     VAFAC = DEXP(VX(I02+J,I)-VX(I20+J,I))
      XFTAN = 1.D0
*     IF(VAFAC.GT.2.D0)XFTAN=1.D0/(1.D0+(VAFAC-2.D0)**2)
*
      CT=3.D-1/XTAN/RXT*(ATFAC(J)-1.D0)*XFTAN
      G(I02+J)=G(I02+J)-CT
      D(I02+J,I20+J)=D(I02+J,I20+J)-3.D-1/XTAN/RXT*ATFAC(J)+
     *    C32*CT*ATFAC(J)/(2.D0+ATFAC(J))
      D(I02+J,I02+J)=D(I02+J,I02+J)+3.D-1/XTAN/RXT*ATFAC(J)-
     *    C32*CT*ATFAC(J)/(2.D0+ATFAC(J))+C32*CT
      D(I02+J,I00+J)=D(I02+J,I00+J)-25.D-1*CT
*
      END IF
*
*        Stellar Escapers
*
      IF(LS(11))THEN
      G(I02+J)=G(I02+J)+DETESC(J)
      END IF
*
*        Optional choice of trivial equation
*
      G(I02+J)=FX(I02+J)*G(I02+J)+
     *  (1.D0-FX(I02+J))*BREM*(X(I02+J,I)-VX(I02+J,I))
      DO 1001 K=1,NG
      C(I02+J,K)=FX(I02+J)*C(I02+J,K)
      D(I02+J,K)=FX(I02+J)*D(I02+J,K)
 1001 E(I02+J,K)=FX(I02+J)*E(I02+J,K)
      D(I02+J,I02+J)=FX(I02+J)*D(I02+J,I02+J)+
     *  (1.D0-FX(I02+J))*BREM
*
*        Diagnostic output
*
       IF(LS(25))THEN
       PRINT*,' COMP. ',J,' TANG. DPDT=',
     *       BREM*(X(I02+J,I)-VX(I02+J,I)),
     *  ' U/R*DPTDR=',UAV(J)/RAVI*DPTDR,' VT/R*DPTDR=',
     *  VTAV(J)/RAVI*DPTDR,' 1/R*DUDR=',
     *   1.D0/RAVI*(HX(I10+J)-HXM(I10+J))/DER,' 1/R*DVTDR=',
     *   1.D0/RAVI*(HX(I12+J)-HXM(I12+J))/DER,' 4*U/R=',
     *   4.D0*UAV(J)/RAVI,' 4*VT/R=',4.D0*VTAV(J)/RAVI,
     *   ' COLLTERM=',-CT,' RXT=',RXT,
     *  ' DETESC=',DETESC(J)
      PRINT*,' D02/20=',-3.D-1/XTAN/RXT*ATFAC(J)+
     *    C32*CT*ATFAC(J)/(2.D0+ATFAC(J))
      PRINT*,' D02/02=',+3.D-1/XTAN/RXT*ATFAC(J)-
     *    C32*CT*ATFAC(J)/(2.D0+ATFAC(J))+C32*CT
      PRINT*,' D02/00=',-25.D-1*CT
       END IF
*
 1000 CONTINUE
*
          RETURN
*
          END
