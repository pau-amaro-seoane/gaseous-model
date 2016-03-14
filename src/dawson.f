      SUBROUTINE DAFAC(AA,BB,XFAC,XRFAC,XTFAC)
      REAL*8 AA,BB,XFAC,XRFAC,XTFAC,TINY,GFAC,GARG,PI,DAWSON
      REAL*8 GXX1,GXX2,HXX1,HXX2
      INTEGER IC1,IC2,IC3,ICTOT
      DATA PI/3.141592654D0/
      DATA IC1,IC2,IC3/0,0,0/
      EXTERNAL DAWSON,GXX1,GXX2,HXX1,HXX2
*       Determination of escaper phase space constants
      XFAC = 1.D0 - DERF(AA)
      XRFAC = 1.D0 - DERF(AA) + 2.D0*AA/DSQRT(PI)*DEXP(-AA**2)
      XTFAC = 1.D0 - DERF(AA)
*
      TINY = 1.D-8
*       Limit for isotropy and non-apocentre criterion or a=b
      IF(AA.LE.BB*(1.D0+TINY).AND.AA.GE.BB*(1.D0-TINY))THEN
*
      IC1 = IC1 + 1
      XFAC = XFAC + 2.D0*AA/DSQRT(PI)*DEXP(-AA**2)
      XRFAC = XRFAC + 4.D0/3.D0*AA**3/DSQRT(PI)*DEXP(-AA**2)
      XTFAC = XTFAC + 2.D0*AA/DSQRT(PI)*DEXP(-AA**2)
     *              + 4.D0/3.D0*AA**3/DSQRT(PI)*DEXP(-AA**2)
*
      IF(XTFAC.LT.0.D0)THEN
      PRINT*,' A,B,Xe,Xr,Xt=',AA,BB,XFAC,XRFAC,XTFAC
      PRINT*,' Terms=',- DERF(AA),2.D0*AA/DSQRT(PI)*DEXP(-AA**2),
     *    4.D0/3.D0*AA**3/DSQRT(PI)*DEXP(-AA**2)
      END IF

      END IF
*       a>b
      IF(AA.GT.BB*(1.D0+TINY))THEN
*
      IC2 = IC2 + 1
      GFAC = DSQRT(1.D0 - (BB/AA)**2)
      GARG = GFAC*AA
      XFAC = XFAC + AA/DEXP(BB**2)*GXX1(GARG)
      XERFGA = C12*DERF(GARG) - GARG/DSQRT(PI)*DEXP(-GARG**2)
      XRFAC = XRFAC + 2.D0*AA**3/DEXP(BB**2)*GXX2(GARG)
      XTFAC = XTFAC + AA/DEXP(BB**2)*(1.D0+BB**2)*GXX1(GARG)
     *              - BB**2*AA/DEXP(BB**2)*GXX2(GARG)
*
      IF(XTFAC.LT.0.D0)THEN
      PRINT*,' A,B,Xe,Xr,Xt=',AA,BB,XFAC,XRFAC,XTFAC
      PRINT*,' gfac,garg,gxx1,gxx2=',GFAC,GARG,gxx1(garg),gxx2(garg)
      PRINT*,' Terms=',- DERF(AA),
     *    + AA/DEXP(BB**2)*(1.D0+BB**2)*GXX1(GARG),
     *    - BB**2*AA/DEXP(BB**2)*GXX2(GARG)
      END IF

      END IF
*       a<b
      IF(AA.LT.BB*(1.D0-TINY))THEN
*
      IC3 = IC3 + 1
      GFAC = DSQRT((BB/AA)**2 - 1.D0)
      GARG = GFAC*AA
      XFAC = XFAC + AA/DEXP(BB**2)*HXX1(GARG)
      XRFAC = XRFAC + 2.D0*AA**3/DEXP(BB**2)*HXX2(GARG)
      XTFAC = XTFAC + AA*(1.D0+BB**2)/DEXP(BB**2)*HXX1(GARG)
     *              - AA*BB**2/DEXP(BB**2)*HXX2(GARG)
*
      IF(XTFAC.LT.0.D0)THEN
      PRINT*,' A,B,Xe,Xr,Xt=',AA,BB,XFAC,XRFAC,XTFAC
      PRINT*,' gfac,garg,hxx1,hxx2,dawson=',GFAC,GARG,HXX1(GARG),
     *   HXX2(GARG),DAWSON(GARG)
      PRINT*,' Terms=',- DERF(AA),
     *    AA*(1.D0+BB**2)/DEXP(BB**2)*HXX1(GARG)
     *    - AA*BB**2/DEXP(BB**2)*HXX2(GARG)
      END IF

      END IF
*
      ICTOT = IC1 + IC2 + IC3
*
      IF(MOD(ICTOT,50000).EQ.0)PRINT*,' dawson check ',IC1,IC2,IC3
      IF(ICTOT.GT.1000000)THEN
      IC1 = 0
      IC2 = 0
      IC3 = 0
      END IF
*
      RETURN
      END
*
      FUNCTION xierf(x)
      REAL*8 xierf,x,pi,dawson
      DATA PI/3.141592654D0/
      EXTERNAL dawson
      xierf = 2.d0/dsqrt(pi)*dexp(x**2)*dawson(x)
      return
      end
*
      FUNCTION gxx1(x)
      REAL*8 gxx1,x,pi
      DATA PI/3.141592654D0/
      parameter (a1=1.d0/3.d0,a2=0.3d0,a3=5.d0/21.d0)
*      use series expansion
      x2 = x**2
      if(dabs(x).lt.0.2)then
      gxx1 = 1.d0 - a1*x2*(1.d0 - a2*x2*(1.d0 - a3*x2))
      gxx1 = 2.d0/dsqrt(pi)*gxx1
      else
*      use explicit formula
      gxx1 = derf(x)/x
      end if
* 
      return
      end
*
      FUNCTION hxx1(x)
      REAL*8 hxx1,x,pi,xierf
      EXTERNAL xierf
      DATA PI/3.141592654D0/
      parameter (a1=1.d0/3.d0,a2=0.3d0,a3=5.d0/21.d0)
*      use series expansion
      x2 = x**2
      if(dabs(x).lt.0.2)then
      hxx1 = 1.d0 + a1*x2*(1.d0 + a2*x2*(1.d0 + a3*x2))
      hxx1 = 2.d0/dsqrt(pi)*hxx1
      else
*      use explicit formula
      hxx1 = xierf(x)/x
      end if
*
      return
      end
*
      FUNCTION gxx2(x)
      REAL*8 gxx2,x,pi
      DATA PI/3.141592654D0/
      parameter (b0=1.d0/3.d0,b1=0.2d0,b2=5.d0/14.d0,b3=7.d0/27.d0)
*      use series expansion
      x2 = x**2 
      if(dabs(x).lt.0.2)then
      gxx2 = b0 - b1*x2*(1.d0 - b2*x2*(1.d0 - b3*x2))
      gxx2 = 2.d0/dsqrt(pi)*gxx2
*      use explicit formula
      else
      gxx2 = (derf(x)/2.d0 - x/dsqrt(pi)*dexp(-x**2))/x**3
      end if
*
      return
      end
*
      FUNCTION hxx2(x)
      REAL*8 hxx2,x,pi,xierf
      EXTERNAL xierf
      DATA PI/3.141592654D0/
      parameter (b0=1.d0/3.d0,b1=0.2d0,b2=5.d0/14.d0,b3=7.d0/27.d0)
*      use series expansion
      x2 = x**2
      if(dabs(x).lt.0.2)then
      hxx2 = b0 + b1*x2*(1.d0 + b2*x2*(1.d0 + b3*x2))
      hxx2 = 2.d0/dsqrt(pi)*hxx2
*      use explicit formula
      else
      hxx2 = (-xierf(x)/2.d0 + x/dsqrt(pi)*dexp(x**2))/x**3
      end if
*
      return
      end
*
      FUNCTION dawson(x)
      INTEGER NMAX
      REAL*8 dawson,x,H,A1,A2,A3
      PARAMETER (NMAX=6,H=0.4,A1=2./3.,A2=0.4,A3=2./7.)
      INTEGER i,init,n0
      REAL*8 d1,d2,e1,e2,sum,x2,xp,xx,c(NMAX)
      SAVE init,c
      DATA init/0/
      if(init.eq.0)then
        init=1
        do 11 i=1,NMAX
          c(i)=exp(-((2.*float(i)-1.)*H)**2)
11      continue
      endif
      if(abs(x).lt.0.2)then
        x2=x**2
        dawson=x*(1.-A1*x2*(1.-A2*x2*(1.-A3*x2)))
      else
        xx=abs(x)
        n0=2*nint(0.5*xx/H)
        xp=xx-float(n0)*H
        e1=exp(2.*xp*H)
        e2=e1**2
        d1=float(n0+1)
        d2=d1-2.
        sum=0.
        do 12 i=1,NMAX
          sum=sum+c(i)*(e1/d1+1./(d2*e1))
          d1=d1+2.
          d2=d2-2.
          e1=e2*e1
12      continue
        dawson=0.5641895835*sign(exp(-xp**2),x)*sum
      endif
*       print*,' x,dawson=',x,dawson
      return
      END
