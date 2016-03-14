          SUBROUTINE INNER
       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
         INCLUDE 'compar.f'
C
          RFAC=R(I+2)/R(I+1)
          G(IMR)=HX(IMR)-HXP(IMR)+3.D0*DLOG(RFAC)
          D(IMR,IMR)=1.D0
          E(IMR,IMR)=-1.D0
C
          DO 1000 J=1,NCOMP
          IF(LEQ(I00+J))THEN
          G(I00+J)=HX(I00+J)-HXP(I00+J)
          D(I00+J,I00+J)=1.D0
          E(I00+J,I00+J)=-1.D0
          END IF
C
          IF(LEQ(I20+J))THEN
          G(I20+J)=HX(I20+J)-HXP(I20+J)
          D(I20+J,I20+J)=1.D0
          E(I20+J,I20+J)=-1.D0
          END IF
C
           IF(LEQ(I02+J))THEN
          G(I02+J)=HX(I02+J)-HXP(I02+J)
          D(I02+J,I02+J)=1.D0
          E(I02+J,I02+J)=-1.D0
           END IF
C
           IF(LEQ(I10+J))THEN
           G(I10+J)=HX(I10+J)-HXP(I10+J)/RFAC
           D(I10+J,I10+J)=1.D0
           E(I10+J,I10+J)=-1.D0/RFAC
           END IF
C
           IF(LEQ(I30+J))THEN
           G(I30+J)=HX(I30+J)-HXP(I30+J)/RFAC
           D(I30+J,I30+J)=1.D0
           E(I30+J,I30+J)=-1.D0/RFAC
           END IF
C
           IF(LEQ(I12+J))THEN
           G(I12+J)=HX(I12+J)-HXP(I12+J)/RFAC
           D(I12+J,I12+J)=1.D0
           E(I12+J,I12+J)=-1.D0/RFAC
           END IF
C
 1000      CONTINUE
C
          RETURN
          END
