      IF(II.NE.-14)GOTO 9900
      PRINT*,' PLSPLO14 called '
      DO 514 IO=IMIN+1,NP-1
          JCC=1
          TINY=1.D-50
          DEM = PI43 * X(I00+JCC,IO)*(R(IO)**3 - R(IO-1)**3)
          DER = R(IO) - R(IO-1)
          DER3 = R(IO)**3 - R(IO-1)**3
          RAV = (R(IO)+R(IO-1))/2.D0
          UAVE = (X(I10+JCC,IO)+X(I10+JCC,IO-1))/2.D0
          PANI=X(I20+JCC,IO)-X(I02+JCC,IO)/2.D0
*
      QU=2.D0*X(I20+JCC,IO)*((X(I10+JCC,IO)-X(I10+JCC,IO-1))/DER -
     *              UAVE/RAV)
*
          QFR = 3.D0/DER3*(R(IO)**2*X(I30+JCC,IO) -
     *                  R(IO-1)**2*X(I30+JCC,IO-1))
          QFT1 = -C32/DER3*(R(IO)**2*X(I12+JCC,IO) -
     *                  R(IO-1)**2*X(I12+JCC,IO-1)) 
          QFT2 = 2.D0*(X(I12+JCC,IO)+X(I12+JCC,IO-1))/RAV
          QF = QFR + QFT1 + QFT2 
          XCOUL(JCC)=DLOG(.11D4)
          XMIND(JCC)=1.D-4
          CS(60)=16.D0/9.D0*DSQRT(PI)
         TRXX=((X(I20+JCC,IO)+X(I02+JCC,IO))/X(I00+JCC,IO)/3.D0)**C32/
     *    CS(60)/XMIND(JCC)/XCOUL(JCC)/X(I00+JCC,IO)
*
          IF(K.EQ.3)VEC(K,IO) = PANI/(QF+TINY)
          IF(K.EQ.2)VEC(K,IO) = 10.D0/9.D0*TRXX
 514    CONTINUE
          VEC(K,IMIN)=VEC(K,IMIN+1)
          VEC(K,NP)=VEC(K,NP-1)
