      SUBROUTINE ESCAPE(NESC,IESC)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
         INCLUDE 'compar.f'
      DIMENSION IESC(*)
*
      IBACT = 0
      XMMIN = 1.D-13
*
      IBNEW = IBIN
*
      DO 100 IB = 1,IBIN
*
      XMBIN = BODY1(IB) + BODY2(IB)
*
      IF(NAMEB(IB).LT.0)THEN
*
*     Save and correct for Escaper Data
*
      ibnew = ibnew - 1
      xmrbin(1) = 0.d0
      dmttot = 0.0d0
      ekt = xmind(1)*dexp(x(i20+1,1)-x(i00+1,1))
*
      write(50,600)time,nameb(ib),ico(ib),iev(ib),eb(ib),ecc(ib),
     *  semia(ib),rb(ib),body1(ib),body2(ib),size1(ib),size2(ib),
     *  vr(ib),vt(ib),rcore,ekt
 600  format(1x,1p,d12.4,3i5,12d12.4)
*
      depot = 0.d0
      do 500 i=2,nj
      dmttot = dmttot + dmbin(ib,i)
      depot = depot - dmbin(ib,i)*(phi(i)+phibin(i))
      xmrbin(i) = xmrbin(i) - dmttot
      if(xmrbin(i).lt.xmmin)xmrbin(i) = 0.d0
      if(dabs(xmrbin(i)-ibnew*xmbin).lt.xmmin)xmrbin(i) = ibnew*xmbin
      if(dabs(xmrbin(i)-xmrbin(i-1)).lt.xmmin)xmrbin(i)=xmrbin(i-1)
 500  continue
*       add indirect heating to energy reservoir
      j=ico(ib)    
*       no ind heat if quad event without single escapers
      if(iev(ib).ne.-4)then
      if(iev(ib).eq.-3)then
      xindhb(j) = xindhb(j) + c12*depot
      else
      xindhb(j) = xindhb(j) + depot
      end if
      end if
*
*     xheat(j) = xheat(j) + depot
*
      ELSE
*
      IBACT = IBACT + 1
*
      IF(IBACT.LT.IB)THEN
      RB(IBACT) = RB(IB)
      ISH(IBACT) = ISH(IB)
      ICO(IBACT) = ICO(IB)
      IEV(IBACT) = IEV(IB)
      BODY1(IBACT) = BODY1(IB)
      BODY2(IBACT) = BODY2(IB)
      SIZE1(IBACT) = SIZE1(IB)
      SIZE2(IBACT) = SIZE2(IB)
      NAMEB(IBACT) = NAMEB(IB)
      ECC(IBACT) = ECC(IB)
      EB(IBACT) = EB(IB)
      SEMIA(IBACT) = SEMIA(IB)
      VR(IBACT) = VR(IB)
      VT(IBACT) = VT(IB)
      TORB(IBACT) = TORB(IB)
      TBIN(IBACT) = TBIN(IB)
      INAME(IBACT) = INAME(IB)
      INAMEX(IBACT) = INAMEX(IB)
      IBSTA(IBACT) = IBSTA(IB)
      ILEN(IBACT) = ILEN(IB)
*
      DO 200 IN = 1,NJ
      DMBIN(IBACT,IN) = DMBIN(IB,IN)
 200  CONTINUE
      END IF
*
      END IF
*
 100  CONTINUE
*
      DO 110 IB = IBIN-NESC+1,IBIN
      DO 120 IN = 1,NJ
 120  DMBIN(IB,IN) = 0.D0
 110  CONTINUE
*
      IBIN = IBIN - NESC
      NESC = 0
      IPOT = 0
      CALL PHIMAS(IPOT)
*
      RETURN
*
      END
