          SUBROUTINE IMREQ
       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
*        Equation for dlog Mr /dlogr
*
         INCLUDE 'compar.f'
*
       RHO=0.D0
       DO 100 J=1,NCOMP
       RHO=RHO+FX(I00+J)*DEXP(HX(I00+J))+
     *       (1.D0-FX(I00+J))*DEXP(VX(I00+J,I))
 100   CONTINUE
*
      XMAV=C12*(HX(IMR)+HXM(IMR))
      ETERM = PI4*RHO*DEXP(3.D0*RAV-XMAV)
*
      G(IMR)=(HX(IMR)-HXM(IMR))/DER-ETERM
*
      D(IMR,IMR)=1.D0/DER+C12*ETERM
      C(IMR,IMR)=-1.D0/DER+C12*ETERM
      DO 110 J=1,NCOMP
      IF(LEQ(I00+J))D(IMR,I00+J)=-ETERM/RHO*DEXP(HX(I00+J))
 110  CONTINUE
*
*        Diagnostic output
*
      IF(LS(25))THEN
      PRINT*,' DLGM/DLGR=',(HX(IMR)-HXM(IMR))/DER,' -ETERM=',-ETERM
      END IF
*
          RETURN
*
          END
