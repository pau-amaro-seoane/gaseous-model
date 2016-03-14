      subroutine orbit(im1,imin1,imax1,rmax1,rmin1,gmin1,gmax1,e1,a1)
      implicit real*8(a-h,o-z),integer(i-k,m-n),logical(l)
*
*     determine many positions for orbit integrals
*
      include 'compar.f'
      PARAMETER(NTAB=32)
      INTEGER IDUM2,IV,IY
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
      dimension pos(500),vel(500),jn(500)
      dimension dmold(njo)
      common/ttx/dmtold,dmttot
*
      real*4 ranf
*
*       if(r(ipmin).gt.rmin1)print*,' ORBIT WARNING ',ipmin,
*    *   ' rmin1,rmax1,ripmin=',rmin1,rmax1,r(ipmin)
*
        xmmin = 1.d-13
        lcheck = .true.
        isum = 0
        qmax = 1.9d0
        ivel = 0
        idiff = imax1 - imin1
        xmbin = body1(im1) + body2(im1)
*
        dmold(1) = 0.d0
        dmbin(im1,1) = 0.d0
*
        do 150 i = 2,nj
        dmold(i) = dmbin(im1,i)
        dmbin(im1,i) = 0.d0
 150    continue
*
        rmmx = rmax1
*       if(rmmx.gt.r(160))then
*       rmmx=r(160)
*       print*,' rmmx160 taken im1=',im1,time
*       print*,' rmin1,rmmx,rmax1=',rmin1,rmmx,rmax1
*       call flush(6)
*       end if
*
*       if(rmin1.gt.rmmx)rmmx=rmax1
        iquot = int(rmmx/rmin1)
        ismax = min(50*iquot,500)
*
   70 s = ranf()
*
        if(s.lt.0.00001d0)s=0.00001d0
        if(s.gt.0.99999d0)s=0.99999d0
*
        s = -1.d0 + 2.d0*s
*
          rn1 = 0.5*(rmin1 + rmax1) +
     &                              0.25*(rmax1 - rmin1)*(3.0*s - s**3)
*
        if(rn1.gt.r(nj))goto 70
*       if(rn1.gt.rmmx)goto 70
*
        irn1 = imin1 - 1
        do 200 i = imin1,nj
        if(r(i).gt.rn1)goto 201
        irn1 = i
 200    continue
 201    continue
*
        dphdr = (phi(irn1+1)+phibin(irn1+1)-phi(irn1)-phibin(irn1))/
     *          (r(irn1+1)-r(irn1))
        u1 = phi(irn1) - phtid + phibin(irn1) + dphdr*(rn1-r(irn1))
*
      q0 = qmax*ranf()*dmax1(gmin1,gmax1)
*
      v1 = 2.0*(e1 - u1) - a1*a1/rn1/rn1
*
      if(v1.le.0.0) then
      ivel = ivel + 1
*
      if(ivel.gt.10000)then
      print*, ' ivel 10000 orbit'
      stop
      end if
*
      go to 70
      end if
      v1 = sqrt(v1)
      ivel = 0
*
      drds = 0.75*(rmax1 - rmin1)*(1.0 - s*s)
      q1 = drds/v1
*
      if(q1.gt.qmax*dmax1(gmin1,gmax1)) then
      write(6,*) 'o q1 > q  iseed,q1,gmin1,gmax1=',iseed,q1,gmin1,gmax1
        qmax = 2.0d0*qmax
        go to 70
      endif
*
      if(q0.gt.q1) go to 70
*
      isum = isum + 1
*
      xmsave = 0.d0
      pos(isum) = rn1
      vel(isum) = v1
      jn(isum) = isum
*
 450  continue
*
      if(isum.lt.ismax) go to 70
*
      call sort1(ismax,pos,jn)
*
      do 250 i = 1,nj
      if(r(i).lt.rmin1)irmin = i
 250  continue
*
      inj = irmin + 1
      if(irmin.eq.0)irmin=1
      ibsta(im1) = irmin
      jt = jn(1)
      period = 2.d0*(pos(1)-rmin1)/vel(jt)
      itold = 0
*
*     print*,' orbit: im1,rmin1,rmax1,imin1,imax1,irmin=',
*    &  im1,rmin1,rmax1,imin1,imax1,irmin
      do 300 it = 2,ismax-1
      it1  = it - 1
      jt = jn(it)
      jt1 = jn(it1)
      period = period + 2.d0*(pos(it)-pos(it1))/(vel(jt)+vel(jt1))
*
      if(r(inj).lt.pos(it))then
      xfrac = dble(it-itold)/dble(ismax)
      xmsave = xmsave + xmbin*xfrac
*
      if(lcheck.and.xmsave.lt.dexp(x(imr,inj)))then
      dmbin(im1,inj) = xmsave
*
      xmsave = 0.d0
*
      lcheck=.false.
*
      else
*
      if(.not.lcheck)then
      dmbin(im1,inj) = xmsave
*
      xmsave = 0.d0
      end if
*
      end if
*
      itold = it
      inj = inj + 1
      end if
*
 300  continue
*
      if(itold.gt.0)then
      xfrac = dble(ismax-itold)/dble(ismax)
      dmbin(im1,inj) = xmbin*xfrac
      else
      dmbin(im1,inj) = xmbin
      end if
      ilen(im1) = inj - irmin
*
      jt = jn(ismax)
      period = period + 2.d0*(rmax1-pos(ismax))/vel(jt)
*
      period = 2.d0*period
*
      if(idiff.gt.0)then
*      torb temp used in binsto/relax May 98
      torb(im1) = period
      else
*
      jt = jn(ismax/2)
      rn1 = pos(ismax/2)
*
      irn1 = imin1 - 1
        do 400 i = imin1,nj
        if(r(i).gt.rn1)goto 401
        irn1 = i
 400    continue
 401    continue
*
      dphdr = (phi(irn1+1)+phibin(irn1+1)-phi(irn1)-phibin(irn1))/
     *          (r(irn1+1)-r(irn1))
      u1 = phi(irn1) - phtid + phibin(irn1) + dphdr*(rn1-r(irn1))
*
      vt(im1) = 2.d0*(e1-u1) - vel(jt)**2
      vt(im1) = dsqrt(vt(im1))
      vr(im1) = vel(jt)
*      torb not used for binsto/relax Febr 99
      torb(im1) = 2.d0*pi*rn1/vt(im1)
*     print*,' orbit: nearly circular idiff,im1,p,torb=',
*    * idiff,im1,period,torb(im1)
      end if
*
      dmttot = 0.d0
      dmtold = 0.d0
      xmrbin(1) = 0.d0
*
      do 500 i=2,nj
      dmttot = dmttot + dmbin(im1,i)
      dmtold = dmtold + dmold(i)
      xmrbin(i) = xmrbin(i) + dmttot - dmtold
*         take care of numerically clean xmrbin inside
      if(xmrbin(i).lt.xmmin)xmrbin(i) = 0.d0
*         ...outside
      if(dabs(xmrbin(i)-ibin*xmbin).lt.xmmin)xmrbin(i) = ibin*xmbin
*         and inside the binary region
      if(dabs(xmrbin(i)-xmrbin(i-1)).lt.xmmin)xmrbin(i)=xmrbin(i-1)
 500  continue
*
*      TEST - determine torb as alocal dynamical time
*
*     if(rb(im1).lt.15.d0*rcore) then
*       torb(im1) = 10.d0*rb(im1)/sqrt(vr(im1)**2 + vt(im1)**2)
*     else
*       torb(im1) = 20.d0*rcore
*     endif
*
      return
      end
