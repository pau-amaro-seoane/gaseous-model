       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
       IREC=0
       IWRITE=0
C
       PRINT*,' Enter Log Start/Stop Time:'
       READ*,TSTART,TSTOP
C
 1999  IREC=IREC+1
C
       IF(MOD(IREC,500).EQ.0)PRINT*,IREC,' Records read from unit 50...'
       IF(IREC.GT.1)THEN
       VTIME=TIME
       ELSE
       VTIME=0.D0
       END IF
C
       READ(50,*,ERR=700,END=700)TIME,TTRX,RHO,SIG,XI
C
       IF(XI.GT.0.01.OR.XI.LT.-0.01)GOTO 1999
       IF(XI.GT.0.005.OR.XI.LT.-0.008)
     *    PRINT*,' IREC=',IREC,' Bad XI=',XI,' T=',1.D1**TIME
       IF(TIME.LT.TSTART)GOTO 1999
C
       IWRITE=IWRITE+1
*
       IF(TIME.GT.TSTOP)GOTO 700
C
       WRITE(51,5005)TIME,TTRX,RHO,SIG,XI
 5005 FORMAT(1P,5(D14.7,1X))
C
       GOTO 1999
C
 700   CONTINUE
       PRINT*,IWRITE-1,' Records written to unit 51 '
       STOP
       END
