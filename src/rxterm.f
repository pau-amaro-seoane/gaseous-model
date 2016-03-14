      SUBROUTINE RXTERM                                                         
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)                    
*
*      Equipartition terms due to small angle grav. encounter effects
*      Field and target particle distribution function anisotropic
*
      INCLUDE 'compar.f'                                                 
C
      DIMENSION PSUM(NCOMPO),PTSUM(NCOMPO)
        COMMON/EQUITEST/EEQUI(NCOMPO),XEEQUI(NCOMPO),XEQEX(NCOMPO)
*      Calculate mass ratios only once per run 
*
      C14=0.25D0
*
      DO 900 J=1,NCOMP
      PSUM(J)=1.D0+2.D0*AFAC(J)
      PTSUM(J)=ATFAC(J)+2.D0
 900  CONTINUE
*
      DO 1000 J=NCOMP,1,-1
*
      DO 1010 K=NCOMP,1,-1
*
      IF(K.EQ.J)GOTO 1010
*
      XMQ=XMIND(K)/XMIND(J)
      PFAC=DEXP(HX(I20+K)-HX(I20+J))
      PTFAC=DEXP(HX(I02+K)-HX(I02+J))
*      Note inverse order in PFAC and XFAC
      XFAC=DEXP(HX(I00+J)-HX(I00+K))
      XRFAC=1.D0/XFAC
      SRFAC=PFAC*XFAC
      STFAC=PTFAC*XFAC
      SFAC=PSUM(K)/PSUM(J)*SRFAC
      TFAC=PTSUM(K)/PTSUM(J)*STFAC
*
      DTR20J=TRX(J,K)*(C32-3.D0*AFAC(J)/PSUM(J)+
     *   C32*SFAC/(1.D0+SFAC)*(2.D0*AFAC(J)/PSUM(J)-1.D0))
      DTT20J=DTR20J
*     DTT02J=TRX(J,K)*C32*(1.D0-ATFAC(J)/PTSUM(J)+
*    *   TFAC/(1.D0+TFAC)*(ATFAC(J)/PTSUM(J)-1.D0))
      DTR02J=TRX(J,K)*(3.D0*AFAC(J)/PSUM(J)-
     *   3.D0*SFAC/(1.D0+SFAC)*AFAC(J)/PSUM(J))
      DTT02J=DTR02J
*     DTT20J=TRX(J,K)*C32*(ATFAC(J)/PTSUM(J)-
*    *   TFAC/(1.D0+TFAC)*ATFAC(J)/PTSUM(J))
      DTR00J=TRX(J,K)*(-2.5D0+C32*SFAC/(1.D0+SFAC))
      DTT00J=DTR00J
*     DTT00J=TRX(J,K)*(-2.5D0+C32*TFAC/(1.D0+TFAC))
      DTR20K=TRX(J,K)*C32*SFAC/(1.D0+SFAC)*
     *  (-2.D0*AFAC(K)/PSUM(K)+1.D0)
      DTT20K=DTR20K
*     DTT02K=TRX(J,K)*C32*TFAC/(1.D0+TFAC)*
*    *  (-ATFAC(K)/PTSUM(K)+1.D0)
      DTR02K=TRX(J,K)*3.D0*SFAC/(1.D0+SFAC)*AFAC(K)/PSUM(K)
      DTT02K=DTR02K
*     DTT20K=TRX(J,K)*C32*TFAC/(1.D0+TFAC)*ATFAC(K)/PTSUM(K)
      DTR00K=-TRX(J,K)*C32*SFAC/(1.D0+SFAC)
      DTT00K=DTR00K
*     DTT00K=-TRX(J,K)*C32*TFAC/(1.D0+TFAC)
*
      C1=C14*(XRFAC*PSUM(J)-XMQ*PFAC*PSUM(K))
      CT1=C14*(XRFAC*PTSUM(J)-XMQ*PTFAC*PTSUM(K))
*      Option to switch off non-classical anisotropy dependent equipartition
      IF(LS(8))THEN
      C2=0.D0
      C3=0.D0
      CT2=0.D0
      CT3=0.D0
      ELSE
      C2=(XRFAC*(1.D0-AFAC(J))*(1.D0+C32*XMQ)-
     *   PFAC*(1.D0-AFAC(K))*(C32+XMQ))/5.D0
      CT2=(XRFAC*(ATFAC(J)-1.D0)*(1.D0+C32*XMQ)-
     *   PTFAC*(ATFAC(K)-1.D0)*(C32+XMQ))/5.D0
      C3=0.3D0*(XRFAC*(1.D0-AFAC(J))+PFAC*(1.D0-AFAC(K)))*
     *    (1.D0+XMQ)*SFAC/(1.D0+SFAC)
      CT3=0.3D0*(XRFAC*(ATFAC(J)-1.D0)+PTFAC*(ATFAC(K)-1.D0))*
     *    (1.D0+XMQ)*TFAC/(1.D0+TFAC)
      END IF
*
      D20J=AFAC(J)*XRFAC*(1.D0+C32*XMQ)/5.D0+
     *    PFAC*(1.D0-AFAC(K))*(C32+XMQ)/5.D0+
     *    0.3D0*(XRFAC*AFAC(J)-PFAC*(1.D0-AFAC(K)))*
     *       (1.D0+XMQ)*
     *     SFAC/(1.D0+SFAC)+C3/(1.D0+SFAC)*
     *     (2.D0*AFAC(J)/PSUM(J)-1.D0) 
      DT02J=-ATFAC(J)*XRFAC*(1.D0+C32*XMQ)/5.D0+
     *    PTFAC*(ATFAC(K)-1.D0)*(C32+XMQ)/5.D0-
     *    0.3D0*(XRFAC*ATFAC(J)+PTFAC*(ATFAC(K)-1.D0))*
     *       (1.D0+XMQ)*
     *     TFAC/(1.D0+TFAC)-CT3/(1.D0+TFAC)*
     *     ATFAC(J)/PTSUM(J)
*
      IF(LS(8))THEN
      D20J=0.D0
      DT02J=0.D0
      END IF
*
      DRD20J=-C12*XRFAC*AFAC(J)+C14*XMQ*PFAC*PSUM(K)+D20J
      DTD02J=-C12*(C12*XRFAC*ATFAC(J)-
     *                          C12*XMQ*PTFAC*PTSUM(K)+DT02J)
*
      D20K=-PFAC*(C32+XMQ)/5.D0+
     *    0.3D0*PFAC*(1.D0+XMQ)*SFAC/(1.D0+SFAC)+
     *              C3/(1.D0+SFAC)*(-2.D0*AFAC(K)/PSUM(K)+1.D0)
      DT02K=PTFAC*(C32+XMQ)/5.D0-
     *    0.3D0*PTFAC*(1.D0+XMQ)*TFAC/(1.D0+TFAC)+
     *              CT3/(1.D0+TFAC)*ATFAC(K)/PTSUM(K)
*
      IF(LS(8))THEN
      D20K=0.D0
      DT02K=0.D0
      END IF
*
      DRD20K=C12*XMQ*PFAC*(AFAC(K)-C12*PSUM(K))+D20K
      DTD02K=C12*(XMQ*PTFAC*(C12*ATFAC(K)-C12*PTSUM(K))-DT02K)
*
      D02J=-AFAC(J)*XRFAC*(1.D0+C32*XMQ)/5.D0
     *   -0.3D0*AFAC(J)*XRFAC*(1.D0+XMQ)*
     *    SFAC/(1.D0+SFAC)-
     *                      C3/(1.D0+SFAC)*2.D0*AFAC(J)/PSUM(J)
      DT20J=ATFAC(J)*XRFAC*(1.D0+C32*XMQ)/5.D0
     *   +0.3D0*ATFAC(J)*XRFAC*(1.D0+XMQ)*
     *    TFAC/(1.D0+TFAC)
     *           +CT3/(1.D0+TFAC)*(2.D0*ATFAC(J)/PTSUM(J)-1.D0)
*
      IF(LS(8))THEN
      D02J=0.D0
      DT20J=0.D0
      END IF
*
      DRD02J=C12*AFAC(J)*XRFAC+D02J
      DTD20J=C12*(C12*ATFAC(J)*XRFAC-DT20J)
*
      D02K=PFAC*AFAC(K)*(C32+XMQ)/5.D0
     *   -0.3D0*PFAC*AFAC(K)*(1.D0+XMQ)*
     *    SFAC/(1.D0+SFAC)+
     *                 C3/(1.D0+SFAC)*2.D0*AFAC(K)/PSUM(K)
      DT20K=-PTFAC*ATFAC(K)*(C32+XMQ)/5.D0
     *   +0.3D0*PTFAC*ATFAC(K)*(1.D0+XMQ)*
     *    TFAC/(1.D0+TFAC)+
     *                 CT3/(1.D0+TFAC)*(-ATFAC(K)/PTSUM(K)+1.D0)
*
      IF(LS(8))THEN
      D02K=0.D0
      DT20K=0.D0
      END IF
*
      DRD02K=-C12*AFAC(K)*XMQ*PFAC+D02K
      DTD20K=-C12*(C12*ATFAC(K)*XMQ*PTFAC+DT20K)
*
      D00J=-XRFAC*(1.D0-AFAC(J))*(1.D0+C32*XMQ)/5.D0
     *    -0.3D0*XRFAC*(1.D0-AFAC(J))*(1.D0+XMQ)*
     *     SFAC/(1.D0+SFAC)+C3/(1.D0+SFAC)
      DT00J=-XRFAC*(ATFAC(J)-1.D0)*(1.D0+C32*XMQ)/5.D0
     *    -0.3D0*XRFAC*(ATFAC(J)-1.D0)*(1.D0+XMQ)*
     *     TFAC/(1.D0+TFAC)+CT3/(1.D0+TFAC)
*
      IF(LS(8))THEN
      D00J=0.D0
      DT00J=0.D0
      END IF
*
      DRD00J=-C14*XRFAC*PSUM(J)+D00J
      DTD00J=-C14*(XRFAC*PTSUM(J)+2.D0*DT00J)
*
      DRD00K=-DRD00J
      DTD00K=-DTD00J
*
*      Scale for the moment the timescale as anisotropy decay timescale
      CRTERM=(C1+C2+C3)/TRX(J,K)/XTEQ
      CTTERM=(CT1-(CT2+CT3)/2.D0)/TRX(J,K)/XTEQ
*
      G(I20+J)=G(I20+J)+CRTERM
*
      D(I20+J,I20+J)=D(I20+J,I20+J)+(DRD20J/XTEQ-CRTERM*DTR20J)/TRX(J,K)
      D(I20+J,I02+J)=D(I20+J,I02+J)+(DRD02J/XTEQ-CRTERM*DTR02J)/TRX(J,K)
      D(I20+J,I00+J)=D(I20+J,I00+J)+(DRD00J/XTEQ-CRTERM*DTR00J)/TRX(J,K)
*
      D(I20+J,I20+K)=(DRD20K/XTEQ-CRTERM*DTR20K)/TRX(J,K)
      D(I20+J,I02+K)=(DRD02K/XTEQ-CRTERM*DTR02K)/TRX(J,K)
      D(I20+J,I00+K)=(DRD00K/XTEQ-CRTERM*DTR00K)/TRX(J,K)
*
      G(I02+J)=G(I02+J)+CTTERM
*
      D(I02+J,I20+J)=D(I02+J,I20+J)+
     *   (DTD20J/XTEQ-CTTERM*DTT20J)/TRX(J,K)
      D(I02+J,I02+J)=D(I02+J,I02+J)+
     *   (DTD02J/XTEQ-CTTERM*DTT02J)/TRX(J,K)
      D(I02+J,I00+J)=D(I02+J,I00+J)+
     *   (DTD00J/XTEQ-CTTERM*DTT00J)/TRX(J,K)
*
      D(I02+J,I20+K)=(DTD20K/XTEQ-CTTERM*DTT20K)/TRX(J,K)
      D(I02+J,I02+K)=(DTD02K/XTEQ-CTTERM*DTT02K)/TRX(J,K)
      D(I02+J,I00+K)=(DTD00K/XTEQ-CTTERM*DTT00K)/TRX(J,K)
*
*       Diagnostic output
*
      IF(LS(25))THEN
      PRINT*,' RXTERM J,K=',J,K
      PRINT*,' Check Aniso Balanced Terms: RAD =',DEXP(HX(I20+J))*
     * (C2+C3),' TANG =',DEXP(HX(I02+J))*(CT2+CT3)
      PRINT*,' Check Energy Balanced Terms: RAD =',DEXP(HX(I20+J))*
     *  C1/TRX(J,K),' TANG=',DEXP(HX(I02+J))*CT1/TRX(J,K)
      PRINT*,' trxjk=',TRX(J,K),' trxkj*q=',TRX(K,J)*XMIND(K)/
     *  XMIND(J)*DEXP(HX(I00+K)-HX(I00+J))
      PRINT*,' CRTERM,CTTERM=',CRTERM,CTTERM
      END IF
*
 1010 CONTINUE
*
 1000 CONTINUE
*
      RETURN                                                                    
C                                                                               
      END                                                                       
