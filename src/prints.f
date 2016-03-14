      SUBROUTINE PRINTS
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
         INCLUDE 'compar.f'
C
      DTIME=1.D0/BREM
      LWRIT=.TRUE.
      IF(LS(23).AND.(ITER.LT.ITMAX/2).AND.(MOD(IMOD,NPR).NE.0).
     * AND.(IMOD.GT.6))LWRIT=.FALSE.
        LP=.FALSE.
C
      I=2
      IP=I+1
      RHO=0.D0
      SIG=0.D0
      DO 25 K=1,NCOMP
      SIG=SIG+DSQRT(C13*(DEXP(X(I20+K,I))+2.D0*DEXP(X(I02+K,I)))/
     *   DEXP(X(I00+K,I)))  
 25   RHO=RHO+DEXP(X(I00+K,I))
      SIG=SIG/DBLE(NCOMP)/UNU
      RHO=RHO/UNRH
C
      T=TIME
      DT=1.D0/BREM
      IF(ITER.GT.99)ITER=99
C
      IF(LWRIT)THEN
      CHANGE=CS(17)
      IF(NMOD.EQ.1)CHANGE=0.D0
      WRITE(6,5)T,DT,RHO,SIG,
     * NMOD,IMOD,CHANGE,ITER,ITVEC,(JG(K),K=1,MIN(9,NPOS))
      CALL FLUSH(6)
      IF(LS(19))WRITE(6,7)XMHOLE,DMTOT,(DMHOLE(K,1),K=1,5)
*        Further output only for multi-mass systems
      IF(NCOMP.GT.1)THEN
      DO 100 J=1,NCOMP
      SIGR=DEXP((X(I20+J,2)-X(I00+J,2))/2.D0)
      SIGT=DEXP((X(I02+J,2)-X(I00+J,2))/2.D0)
      WRITE(6,6)J,T,DEXP(X(I00+J,2))/UNRH,SIGR/UNU,SIGT/UNU,
     *       XMIND(J)/UNM*(SIGR/UNU)**2,
     *       XMIND(J)/UNM*(SIGT/UNU)**2
 100  CONTINUE
      END IF
C
      END IF
C
 5        FORMAT(1X,'T=',1P,D12.4,' DT=',D9.2,' RHO=',D12.4,
     * ' SIG=',D12.4,'*',I8,'*',I8,' C=',
     * D9.2,' IT=',I2,'/',I2,' JG=',9I3)
 6        FORMAT(1X,' COMP=',I3,' T=',1P,D12.5,' RHO=',1P,D9.2,
     *              ' SIGR=',D9.2,
     *              ' SIGT=',D9.2,' M*SIGR2=',D9.2,' M*SIGT2/2=',D9.2)
 7        FORMAT(1X,' XMHOLE=',1P,D9.2,' DMTOT=',D9.2,
     *              ' DMHOLE(1-5,1)=',5D9.2)
C
*        Calculate velocity of grid movement CS(19)
*        and grid number of core radius IS(12)
*        from core radius of first component
*
      CS(19)=0.D0
      IS(12)=2
*
      IF(LS(26))THEN
      J=1
      SUM1=DEXP(X(I20+J,1)-X(I00+J,1))*
     *               (1.D0+2.D0*DEXP(X(I02+J,1)-X(I20+J,1)))
*   Take into account binaries for core radius
        IF(LS(5))RHOT = DEXP(X(I00+J,1)) + RHOBAV
      RCORE=DSQRT(SUM1/(4.D0*PI*GRAV*RHOT))
      CS(19)=CS(19)-
     *    RCORE/3.D0*BREM*(X(I00+J,1)-VX(I00+J,1))
      IS(12)=0
      DO 5041 IRAD=1,NJ
      IF(R(IRAD).GT.RCORE)GOTO 5042
 5041 CONTINUE
 5042 CONTINUE
      IS(12)=IRAD-1
      IF(LWRIT)WRITE(6,8)IS(12),CS(19)
   8  FORMAT(1X,1P,' IGR=',I5,' UGR=',D12.5)
      END IF
*
                  RETURN
                  END
