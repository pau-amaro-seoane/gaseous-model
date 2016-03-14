      IF(II.NE.-12)GOTO 9900
      PRINT*,' PLSPLO12 called '
      DO 512 IO=IMIN+1,NP-1
        JCC=1
          TINY=1.D-50
        DLNR=DLOG(R(IO))-DLOG(R(IO-1))
        DLNRH=DLOG(X(I00+JCC,IO))-DLOG(X(I00+JCC,IO-1))
*
          VEC(2,IO) = DLNRH/(DLNR+TINY)
 512    CONTINUE
          VEC(2,IMIN)=VEC(2,IMIN+1)
          VEC(2,NP)=VEC(2,NP-1)
