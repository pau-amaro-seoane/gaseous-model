      IF(II.NE.-11)GOTO 9900
      DO 511 IO=IMIN+1,NP-1
          JCC=1
          VEC(2,IO) = 2.D0-2.D0*X(I02+JCC,IO)/X(I20+JCC,IO)
 511    CONTINUE
          VEC(2,IMIN)=VEC(2,IMIN+1)
          VEC(2,NP)=VEC(2,NP-1)
