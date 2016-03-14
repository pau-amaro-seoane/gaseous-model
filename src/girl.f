      SUBROUTINE GIRL(A,N,M)
C******************************************************************
C
C    CALLED BY: HENYEY
C
C    Matrix Inversion by Gauss elimination procedure
C    Some Matrix Operations necessary for Henyey
C    See Manual for Use of HENYEY
C------------------------------------------------------------------
C    Note: this program is more stable against 
C          badly conditioned matrices as some
C          library routines from IMSL/ESSL etc.
C  
C
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
      DIMENSION A(*)
      NPM = N+M
      DO  25  J = 1,N
         NJ = (J-1)*N
         JJ = NJ + J
         J1 = J + 1
         AMAX =DABS (A(JJ))
         JM = J
      IF( J1 - N )   11,11,16
   11 DO  13  I = J1,N
         IJ = NJ + I
      IF(DABS (A(IJ)) - AMAX )   13,13,12
   12    AMAX =DABS (A(IJ))
         JM = I
   13 CONTINUE
      IF( JM - J )   40,16,14
   14    I1 = JM + NJ
         I2 = JJ
      DO  15  I = J,NPM
         ZWI = A(I1)
         A(I1) = A(I2)
         A(I2) = ZWI
         I1 = I1 + N
   15    I2 = I2 + N
   16 IF( A(JJ) )   20,40,20
   20 DO  23  I = 1,N
      IF( I - J )   21,23,21
   21    IJ = NJ + I
         IK = NJ + I
         JK = JJ
         FAKTOR = - A(IJ)/A(JJ)
      DO  22  K = J1,NPM
         JK = JK + N
         IK = IK + N
   22    A(IK) = A(IK) + FAKTOR * A(JK)
   23 CONTINUE
         JK = JJ
         FAKTOR = 1./A(JJ)
      DO  24  K = J1,NPM
         JK = JK + N
   24    A(JK) = A(JK) * FAKTOR
   25 CONTINUE
      RETURN
   40 CONTINUE
      WRITE(6,100)
  100 FORMAT(' ERROR EXIT IN SUBROUTINE GIRL')
      M = -M
      CALL FLUSH(6)
      RETURN
       END
