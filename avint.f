      SUBROUTINE AVINT(XAV,X,Y,N,XLO,XUP,IND)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  THIS PROGR. COMPUTES AN APPROX. TO THE INTEGRAL Y(X) OVER THE
C  INTERVAL XLO-XUP. X,Y ARE ONE-DIMENSIONAL ARRAAYS OF N ELEMENTS
C  WHERE N>=3. THE X(I) ARE USUALLY UNEQUALLY SPACED OUT AND MUST BE
C  DISTINCT AND IN ASCENDING ORDER. IF THIS IS NOT THE CASE OR IF
C  N<3, THE PROGR. RETURNS WITH IND=0, OTHERWISE IND=1.
C  THERE ARE NO RESTRICTIONS ON XLO AND XUP EXCEPT THAT OF ACCURACY
C  WHICH REQUIRES THEM TO BE CLOSE TO THE INTERVAL (X(I),X(N)).
C  REFERENCE: METHODS OF NUMERICAL INTEGRATION, BY P.J.DAVIS AND
C  P.RABINOWITZ, ACADEMIC PRESS NEW YORK, LONDON 1975
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IMPLICIT REAL*8 (A-H,O-Z),INTEGER(I-N)
      DIMENSION X(1),Y(1)
      IND=0
          XAV=0.E0
      IF(N .LT. 3) RETURN
      DO 10 I=2,N
      IF(X(I) .LE. X(I-1)) RETURN
   10 CONTINUE
      SUM =0.
      IF(XLO .LE. XUP) GOTO 5
      SYL = XUP
      XUP = XLO
      XLO = SYL
      IND = -1
      GOTO 6
    5 IND = 1
      SYL = XLO
    6 IB = 1
      J = N
      DO 1 I=1,N
      IF(X(I) .GE. XLO) GOTO 7
    1 IB=IB+1
    7 IB = MAX0(2,IB)
      IB = MIN0(IB,N-1)
      DO 2 I=1,N
      IF(XUP .GE. X(J)) GOTO 8
    2 J=J-1
    8 J = MIN0(J,N-1)
      J = MAX0(IB,J-1)
      DO 3 JM=IB,J
      X1 = X(JM-1)
      X2 = X(JM)
      X3 = X(JM+1)
      TERM1 = Y(JM-1)/((X1-X2)*(X1-X3))
      TERM2 = Y(JM)/((X2-X1)*(X2-X3))
      TERM3 = Y(JM+1)/((X3-X1)*(X3-X2))
      A = TERM1+TERM2+TERM3
      B = -(X2+X3)*TERM1-(X1+X3)*TERM2-(X1+X2)*TERM3
      C = X2*X3*TERM1+X1*X3*TERM2+X1*X2*TERM3
      IF(JM .GT. IB) GOTO 14
      CA = A
      CB = B
      CC = C
      GOTO 15
   14 CA = .5*(A+CA)
      CB = .5*(B+CB)
      CC = .5*(C+CC)
   15 SUM = SUM + CA*(X2**3-SYL**3)/3.+CB*.5*(X2**2-SYL**2)+CC*(X2-SYL)
      CA = A
      CB = B
      CC = C
    3 SYL = X2
      XAV = SUM+CA*(XUP**3-SYL**3)/3.+CB*.5*(XUP**2-SYL**2)
     F  +CC*(XUP-SYL)
      IF(IND .EQ. 1) RETURN
      IND = 1
      SYL = XUP
      XUP = XLO
      XLO = SYL
      XAV = -XAV
      RETURN
      END
