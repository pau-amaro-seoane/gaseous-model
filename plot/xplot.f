          SUBROUTINE XPLOT(ILAUF,ICOUNT)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
          INCLUDE 'compar.f'
          INCLUDE 'plotsplo.f'
C
          IF(ILAUF.GT.IMODA)GOTO 1111
          JS=6*(ICOUNT-1)
          IU=5*(ICOUNT-1)
          NP=NJ-1
          IMIN=2
C
          DO 111 K=1,5
          II=IWORK(JS+K+1)
C
          DO 112 JC=1,NCOMP
          XMASS(JC)=0.D0
C
          DO 110 IO=IMIN,NP
          IF(II.EQ.0)VEC(K,IO)=R(IO)/UNIT(IU+K)
          IF(II.EQ.IMR)VEC(K,IO)=X(IMR,IO)/UNIT(IU+K)
          IF(II.EQ.I00+JC)VEC(K,IO)=X(I00+JC,IO)/UNIT(IU+K)
          IF(II.EQ.I20+JC)
     *  VEC(K,IO)=DSQRT(X(I20+JC,IO)/X(I00+JC,IO)/UNIT(IU+K))
          IF(II.EQ.I02+JC)
     *  VEC(K,IO)=DSQRT(X(I02+JC,IO)/X(I00+JC,IO)/UNIT(IU+K))
          IF(II.EQ.I10+JC)VEC(K,IO)=DABS(X(I10+JC,IO))/UNIT(IU+K)
          IF(II.EQ.I30+JC)VEC(K,IO)=DABS(X(I30+JC,IO))/UNIT(IU+K)
          IF(II.EQ.I12+JC)VEC(K,IO)=DABS(X(I12+JC,IO))/UNIT(IU+K)
          IF(II.EQ.NG+JC)THEN
          XMASS(JC)=XMASS(JC)+PI43*X(I00+JC,IO)*(R(IO)**3-R(IO-1)**3)
          VEC(K,IO)=XMASS(JC)/SUNM/UNIT(IU+K)
          END IF
          IF(II.EQ.NG+NCOMP+1)VEC(K,IO)=AEI(2)/YEAR
 110      CONTINUE
C
      GOTO(1,2,3,4,5)IABS(II)-10
C
 1        CONTINUE
          INCLUDE 'plsplo11.f'
          GOTO 9900
 2        CONTINUE
          INCLUDE 'plsplo12.f'
          GOTO 9900
 3        CONTINUE
          INCLUDE 'plsplo13.f'
          GOTO 9900
 4        CONTINUE
          INCLUDE 'plsplo14.f'
          GOTO 9900
 5        CONTINUE
          INCLUDE 'plsplo15.f'
 9900     CONTINUE
C
 112      CONTINUE
 111      CONTINUE
C
          GOTO 9500
C
 1111     CONTINUE
C
          JS=6*(ICOUNT+NPLOT-1)
          IU=5*(ICOUNT+NPLOT-1)
C
          DO 121 K=1,5
C
          II=IWORK(JS+K+1)
          ITP=NMOMIN
          DO 120 IO=1,NTP-1
          ITP=ITP+NINT
 1221     CONTINUE
          IF(II.EQ.0)VEC(K,IO)=VST(ITP,1)/UNIT(IU+K)
          IF(II.GT.0)VEC(K,IO)=VST(ITP,II)/UNIT(IU+K)
C
 120      CONTINUE
 1222     CONTINUE
          IF(ITP.GT.NDATA)NTP=IO-1
C
          GOTO(21,22,23,24,25)IABS(II)-20
C
 21       CONTINUE
          INCLUDE 'plsplo21.f'
          GOTO 9800
 22       CONTINUE
          INCLUDE 'plsplo22.f'
          GOTO 9800
 23       CONTINUE
          INCLUDE 'plsplo23.f'
          GOTO 9800
 24       CONTINUE
          INCLUDE 'plsplo24.f'
          GOTO 9800
 25       CONTINUE
          INCLUDE 'plsplo25.f'
C
 9800     CONTINUE
 121      CONTINUE
C
 9500     CONTINUE
C
          KMIN=IMIN
          KMAX=NP
          IF(ILAUF.GT.IMODA)THEN
          KMIN=1
          KMAX=NTP-1
          END IF
          DO 9501 IO=KMIN,KMAX
          DO 9502 K=1,5
          IF(LLOG(JS+K))THEN
          IF(DABS(VEC(K,IO)).LT.1.D-50)VEC(K,IO)=1.D-50
          VEC(K,IO)=DLOG10(DABS(VEC(K,IO)))
          VEC(K,IO)=DINT(1.D5*VEC(K,IO))/1.D5
          END IF
 9502     CONTINUE
 9501     CONTINUE
          RETURN
 999      PRINT*,' Error in Read from INFO-File in XPLOT'
          STOP
C
          END
