          SUBROUTINE SLOPE(XY,NFIT,STAT)
          IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
          DIMENSION XY(30,2),STAT(3)
C
C
          SUMX=0.D0
          SUMY=0.D0
          SUMXY=0.D0
          SUMX2=0.D0
          SUMY2=0.D0
C
          DO 1 J=1,NFIT
          SUMX=SUMX+XY(J,1)
          SUMY=SUMY+XY(J,2)
          SUMX2=SUMX2+XY(J,1)*XY(J,1)
          SUMY2=SUMY2+XY(J,2)*XY(J,2)
          SUMXY=SUMXY+XY(J,1)*XY(J,2)
 1        CONTINUE
C
          STAT(1)=SUMXY - SUMX*SUMY/DBLE(NFIT)
          STAT(1)=STAT(1)/(SUMX2 - SUMX**2/DBLE(NFIT))
C
          STAT(3)=(SUMY - STAT(1)*SUMX)/DBLE(NFIT)
C
          DSUM2=0.D0
          DO 2 J=1,NFIT
          HILF=STAT(3) + STAT(1)*XY(J,1)
          SS=XY(J,2)-HILF
          DSUM2=DSUM2 + SS*SS
 2        CONTINUE
C
          STAT(2)=DSQRT(DSUM2/DBLE(NFIT-1))
C
          RETURN
          END
