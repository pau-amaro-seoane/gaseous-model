      SUBROUTINE OUTER
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
      INCLUDE 'compar.f'
*
*       Extrapolation of velocities at outer boundary
*
      DO 1000 J=1,NCOMP
*
      IF(ITER.LT.5)THEN
      IF(HXM(I10+J).GT.0.D0)THEN
      IF(LS(30))UFAC=DEXP(HXM(I00+J)-HX(I00+J))*(R(IM)/R(I))**2
      IF(LS(31))UFAC=0.D0
      ELSE
      UFAC=0.D0
      END IF
      END IF
*
      G(I10+J)=HX(I10+J)-UFAC*HXM(I10+J)
      D(I10+J,I00+J)=UFAC*HXM(I10+J)
      C(I10+J,I00+J)=-UFAC*HXM(I10+J)
      D(I10+J,I10+J)=1.D0
      C(I10+J,I10+J)=-UFAC
*
      IF(ITER.LT.5)THEN
      IF(HXM(I30+J).GT.0.D0)THEN
      IF(LS(30))UFAC=DEXP(HXM(I20+J)-HX(I20+J))*(R(IM)/R(I))**2
      IF(LS(31))UFAC=0.D0
      ELSE
      UFAC=0.D0
      END IF
      END IF
*
      G(I30+J)=HX(I30+J)-UFAC*HXM(I30+J)
      D(I30+J,I30+J)=1.D0
      C(I30+J,I30+J)=-UFAC
*
      IF(ITER.LT.5)THEN
      IF(HXM(I12+J).GT.0.D0)THEN
      IF(LS(30))UFAC=DEXP(HXM(I02+J)-HX(I02+J))*(R(IM)/R(I))**2
      IF(LS(31))UFAC=0.D0
      ELSE
      UFAC=0.D0
      END IF
      END IF
*
      G(I12+J)=HX(I12+J)-UFAC*HXM(I12+J)
      D(I12+J,I12+J)=1.D0
      C(I12+J,I12+J)=-UFAC
*
 1000 CONTINUE
*
      RETURN
*
      END
