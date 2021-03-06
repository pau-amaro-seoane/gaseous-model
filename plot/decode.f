      SUBROUTINE DECODE(IJK,ISET,ILAUF)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER*4(I-K,M-N),LOGICAL(L)
      REAL*4 VVEC(2000),XVEC(1000)
          INCLUDE 'compar.f'
          INCLUDE 'equiv.f'
          INCLUDE 'plotsplo.f'
*     LOGICAL LST(2000,IVDIM)
C
      IF(ILAUF.GT.IMODA)GOTO 2222
C
      LPR=.FALSE.
      NX=NJ
      DO 150 J=1,NG
      CALL DCODE(20,XVEC,NX,IERR,ILKIN)
      IF(IERR.NE.0)PRINT*,' Error in DCODE Source Line ',ILKIN
      DO 160 K=1,NX
 160  X(J,K)=1.D1**DBLE(XVEC(K))
C
      IF(LPR)WRITE(6,400)J,XVEC(2),(XVEC(K),K=NX/4,NX,NX/4)
      IF(LPR)WRITE(6,399)J,X(J,2),(X(J,K),K=NX/4,NX,NX/4)
 400  FORMAT(1X,' XVEC(',I3,',2,50;100;150;200)=',1P,5D9.2,' DECODED')
 399  FORMAT(1X,' X   (',I3,',2,50;100;150;200)=',1P,5D9.2)
 150  CONTINUE
C
      CALL DCODE(20,XVEC,NX,IERR,ILKIN)
      IF(IERR.NE.0)PRINT*,' Error in DCODE Source Line ',ILKIN
C
      DO 161 K=1,NX
 161  R(K)=1.D1**DBLE(XVEC(K))
      IF(LPR)WRITE(6,402)XVEC(2),(XVEC(K),K=NX/4,NX,NX/4)
      IF(LPR)WRITE(6,404)R(2),(R(K),K=NX/4,NX,NX/4)
 402  FORMAT(1X,' XVEC(',3X,',2,50;100;150;200)=',1P,5D9.2,' DECODED')
 404  FORMAT(1X,' R   (',3X,',2,50;100;150;200)=',1P,5D9.2)
C
      LERR=.FALSE.
      DO 450 IK=1,4
      ILKIN=ILKIN+1
      READ(20,'(5F11.6)')(AEI(K),K=5*(IK-1)+1,5*IK)
 450  CONTINUE
C
      DO 451 IK=1,20
 451  AEI(IK)=1.D1**AEI(IK)
C
      ILKIN=ILKIN+1
      READ(20,'(10I7)')(IEI(K),K=1,10)
      ILKIN=ILKIN+1
      READ(20,'(10I7)')(IEI(K),K=11,20)
C
      LERR=LERR.OR.(IERR.NE.0)
      IF(IERR.NE.0)PRINT*,' ERROR AT J=',J,' IN DECODING DATA'
C
      IF(LPR)WRITE(6,403)(AEI(K),K=10,100,10)
      IF(LPR)WRITE(6,405)(IEI(K),K=1,20)
 403  FORMAT(1X,' AEI=',1P,5D9.2,/,6X,5D9.2,/)
 405  FORMAT(1X,' IEI 1-20=',10(I5,1X),/,12X,10(I5,1X))
C
      IF(LERR)GOTO 888
      RETURN
C
 2222 CONTINUE
C
      IF(IABS(ISET).EQ.1)REWIND 21
      DO 101 JK=1,IVDIM
      CALL DCODE(21,VVEC,NDATA,IERR,ILKIN)
*
*     PRINT*,' JK=',JK,' VVEC decoded NDATA=',NDATA
*     READ(21,1001)(LST(K,JK),K=1,NDATA)
*     PRINT*,' JK=',JK,' LST read'
*1001 FORMAT(78L1)
*
      LMATCH=JK.EQ.IWORK(1).OR.JK.EQ.IWORK(2).
     *    OR.JK.EQ.IWORK(3).OR.JK.EQ.IWORK(4).OR.JK.EQ.IWORK(5)
      LMATCH=LMATCH.OR.JK.GT.140
      IF(LPR)THEN
      WRITE(6,*)(IWORK(K),K=1,5)
      WRITE(6,602)
     *NDATA,JK,ILKIN,(K,VVEC(K),K=1,NDATA,10)
602   FORMAT(1X,1P,I3,1X,' JK=',I3,' ILKIN=',I7,/,
     * (I3,D9.2,I3,D9.2,I3,D9.2,I3,D9.2,I3,D9.2,I3,D9.2,/))
      END IF
      IF(IERR.NE.0)THEN
      PRINT*,' Error in READ VIT by DCODE ',ILKIN,' JK,IERR=',JK,IERR
      ISET=-ISET
      RETURN
      END IF
C
      LFAUL=.FALSE.
      DO 102 IK=1,NDATA
      BARG = VVEC(IK)
      IF(BARG.GT.1.D2)THEN
      LFAUL = .TRUE.
      BARG=1.D2
      END IF
      VST(IK,JK)=1.D1**DBLE(BARG)
*     IF(.NOT.LST(IK,JK))VST(IK,JK)=-VST(JK,IK)
 102  CONTINUE
      IF(LFAUL)PRINT*,' JK=',JK,' overflow in VIT-values '
 101  CONTINUE
C
      PRINT*,' ISET=',ISET,' Successful DECODE of ',NDATA,' VIT-Blocks '
      NTP=NDATA/NINT
      NTP2=(NMOMAX-NMOMIN)/NINT
      IF(NTP2.LT.NTP)NTP=NTP2
      PRINT*,' To be compressed to ',NTP,' Plot Points '
C
      RETURN
C
 888  PRINT*,' ERROR IN DCODE ENCOUNTERED - STOP'
      STOP
      END
