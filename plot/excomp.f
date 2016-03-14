      SUBROUTINE EXCOMP(ISET,ILAUF)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
          INCLUDE 'compar.f'
          INCLUDE 'plotsplo.f'
C
          IF(ILAUF.LE.IMODA)NMAX=NPLOT
          IF(ILAUF.GT.IMODA)NMAX=NTIME
          PRINT*,' This is EXCOMP with NMAX=',NMAX,' ILA=',ILAUF
          PRINT*,'         IMODA=',IMODA
          DO 22 ICOUNT=1,NMAX
          PRINT*,' Now XPLOT called with ICOUNT=',ICOUNT
          CALL XPLOT(ILAUF,ICOUNT)
C
          ICNT=ICOUNT
          IF(ILAUF.GT.IMODA)ICNT=NPLOT+ICOUNT
          IMIN=2
          IF(ILAUF.GT.IMODA)IMIN=1
          IF((ILAUF.EQ.1).OR.(ILAUF.GT.IMODA))THEN
          IF(ISET.EQ.1)THEN
          DO 200 K=1,5
          PARY(2*K-1,ICNT)=VEC(K,IMIN)
          PARY(2*K,ICNT)=PARY(2*K-1,ICNT)
  200     CONTINUE
          END IF
          END IF
C
          IMAX=NJ-1
          IF(ILAUF.GT.IMODA)IMAX=NTP
          DO 222 J=IMIN,IMAX
          DO 250 K=1,5
          IF(VEC(K,J).GT.PARY(2*K-1,ICNT))PARY(2*K-1,ICNT)=VEC(K,J)
          IF(VEC(K,J).LT.PARY(2*K,ICNT))PARY(2*K,ICNT)=VEC(K,J)
 250      CONTINUE
 222      CONTINUE
C
          DO 105 IL=2,10,2
      IF((PARY(IL-1,ICNT).EQ.0.D0).AND.(PARY(IL,ICNT).EQ.0.D0))GOTO 105
          DPY=PARY(IL-1,ICNT)-PARY(IL,ICNT)
          XMAX=DMAX1(DABS(PARY(IL-1,ICNT)),DABS(PARY(IL,ICNT)))
          IF(DPY.LT.XMAX/1.D1)DPY=XMAX/1.D1
          PARY(IL-1,ICNT)=PARY(IL-1,ICNT)+DPY/2.D1
          PARY(IL,ICNT)=PARY(IL,ICNT)-DPY/2.D1
          DPY=PARY(IL-1,ICNT)-PARY(IL,ICNT)
      IF(LPR)WRITE(6,205)PARY(IL,ICNT),PARY(IL-1,ICNT),DPY/1.D1
 205  FORMAT(' SCALE PAR MIN=',1PD11.4,' MAX=',D11.4,' INT=',D11.4)
          XMAX=DABS(PARY(IL-1,ICNT))
          SGN=XMAX/PARY(IL-1,ICNT)
          IMX=1
          IF(DABS(PARY(IL,ICNT)).GT.XMAX)THEN
          XMAX=DABS(PARY(IL,ICNT))
          SGN=XMAX/PARY(IL,ICNT)
          IMX=-1
          END IF
          IM=0
          IN=0
           IF(XMAX.GT.1.D3)THEN
 106        IM=IM+1
            XMAX=XMAX/1.D1
            DPY=DPY/1.D1
            IF(XMAX.GT.1.D3)GOTO 106
           ELSE
 108        CONTINUE
            IF(XMAX.GE.1.D2)GOTO 107
            IN=IN+1
            XMAX=XMAX*1.D1
            DPY=DPY*1.D1
            GOTO 108
           END IF
 107      CONTINUE
           IF(IMX.EQ.1)THEN
           IF(PARY(IL-1,ICNT).GE.0.D0)XMAX=DFLOAT(INT(XMAX)+1)
           IF(PARY(IL-1,ICNT).LT.0.D0)XMAX=DFLOAT(INT(XMAX)-1)
           END IF
           IF(IMX.EQ.-1)THEN
           IF(PARY(IL,ICNT).GE.0.D0)XMAX=DFLOAT(INT(XMAX)-1)
           IF(PARY(IL,ICNT).LT.0.D0)XMAX=DFLOAT(INT(XMAX)+1)
           END IF
 
           IP=0
 111       IP=IP+10
           IF(DFLOAT(IP).LT.DPY)GOTO 111
           DPY=DFLOAT(IP)
           IF(IM.NE.0)THEN
           XMAX=XMAX*1.D1**IM
           DPY=DPY*1.D1**IM
           END IF
           IF(IN.NE.0)THEN
           XMAX=XMAX/1.D1**IN
           DPY=DPY/1.D1**IN
           END IF
        IF(IMX.EQ.1)THEN
        PARY(IL-1,ICNT)=XMAX*SGN
        PARY(IL,ICNT)=PARY(IL-1,ICNT)-DPY
        END IF
        IF(IMX.EQ.-1)THEN
        PARY(IL,ICNT)=XMAX*SGN
        PARY(IL-1,ICNT)=PARY(IL,ICNT)+DPY
        END IF
        WRITE(6,205)PARY(IL,ICNT),PARY(IL-1,ICNT),
     * (PARY(IL-1,ICNT)-PARY(IL,ICNT))/1.D1
        JS=6*(ICNT-1)
        II=IWORK(JS+IL/2)
        IJ=IWORK(JS+IL/2+1)
        DO 1111 KJI=1,NCOMP
        LEQUA=LEQUA.OR.((IJ.EQ.I02+KJI).AND.(II.EQ.I20+KJI))
 1111   CONTINUE
        IF(LEQUA)THEN
        PRINT*,' STELLAR V-DISPERSIONS PLOT NEW PARS CALCULATED:'
        PARY(IL-3,ICNT)=DMAX1(PARY(IL-3,ICNT),PARY(IL-1,ICNT))
        PARY(IL-1,ICNT)=PARY(IL-3,ICNT)
        PARY(IL-2,ICNT)=MIN(PARY(IL-2,ICNT),PARY(IL,ICNT))
        PARY(IL,ICNT)=PARY(IL-2,ICNT)
        WRITE(6,205)PARY(IL-2,ICNT),PARY(IL-3,ICNT),
     * (PARY(IL-3,ICNT)-PARY(IL-2,ICNT))/1.D1
        WRITE(6,205)PARY(IL,ICNT),PARY(IL-1,ICNT),
     * (PARY(IL-1,ICNT)-PARY(IL,ICNT))/1.D1
          END IF
 105      CONTINUE
C
 22       CONTINUE
          RETURN
          END
