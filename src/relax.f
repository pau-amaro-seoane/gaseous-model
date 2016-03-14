      subroutine relax(im1,im2,ipt1,ipt2,vr2,vt2,sm2,xnpart)
      implicit real*8(a-h,o-z),integer(i-k,m-n),logical(l)
*
*
*
*       determine relaxation process for stars between lmin and lmax
*       ------------------------------------------------------------
*       with time-step dt
*       -----------------
*
*
      include 'compar.f'
      PARAMETER(NTAB=32)
      INTEGER IDUM2,IV,IY
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
*
*     real*8 sm1,sm2,ww,den,a12,sin2b2,dt1,beta,dt,dt0,rx1,rx2,ek,ek0,
*    &       ep,ep1,ep2,voo
*
      real ranf
      external ranf
*
*     integer i,im1,lmin,lmax,nup,n,nc2,ibound,m,k
*
*      determine deflection angle and new velocities
*
*
        if(nmod.eq.0)goto 9999
*
         sm1 = body1(im1) + body2(im1)
*      take tbin determined in binsto as dt for single-binary relaxation
        dtrx = tbin(im1)
*      decide whether bin-single or bin-bin relaxation is done
*
         if(im2.gt.0)then
             sm2 = body1(im2) + body2(im2)
             vr2 = vr(im2)
             vt2 = vt(im2)
*       average of relaxation interval of both binaries
             dtrx = c12*(tbin(im1)+tbin(im2))
         end if
         if(dtrx.lt.1.d0/brem)dtrx=1.d0/brem
*
*       determine relative velocity of stars i and i+1
*
         if(ranf().lt.0.5) vr(im1) = -vr(im1)
*
         call relvel(im1,ww,vr2,vt2)
*
*       determine the deflection angle
*
         a12 = 2.d0*pi*grav*grav*(sm1 + sm2)**2/ww**3*xnpart*xctot
*
         ipihalf = 0
*        if(dt0.le.0.d0)dt0=dtrx
         dt0 = dtrx
*
         sin2b2 = a12*dt0
         defang = sin2b2
*
*       return in case of large deflection angle for test
*        if(sin2b2.gt.2.5d-3)then
*        ipt1=-1
*        return
*        end if
*
*       print*,' Binary ',im1,' defl. angle=',sin2b2
*       if beta greater than PI form a sequence of interaction
*
*        if(dt1.le.0.d0)dt1=dtrx
          dt1 = dtrx
          iddt1 = 0
   20    continue
*        if(sin2b2.gt.0.5) then
         if(sin2b2.gt.0.25) then
*        if(sin2b2.gt.2.74E-03) then
*
*          dt2 = 2.74E-03/a12
           dt2 = 0.25/a12
*          dt2 = 0.5/a12
*          dt2 = 1.0/a12
           dt2 = dt1 - dt2
*
*       find new velocities after interaction with deflection angle = PI/2
*
           beta = pi/3.d0
*          beta = 0.5*pi
*          beta = pi
           ipihalf = ipihalf + 1
           call newvel(im1,beta,vr2,vt2,vrn2,vtn2,sm2)
*       do not let binary with min rmin go down further
           if(im1.lt.0)return
*
           vr2 = vrn2
           vt2 = vtn2
*
           call relvel(im1,ww,vr2,vt2)
*
*       determine new deflection angle
*
         a12 = 2.d0*pi*grav*grav*(sm1 + sm2)**2/ww**3*xnpart*xctot
*       use approx. pi/8 as emergency deflection angle
           sin2b2 = a12*dt2
           dt1 = dt2
           iddt1 = iddt1 + 1
           if(iddt1.gt.10000) then
             sin2b2 = 0.25d0
             write(6,*) ' iddt1 > 10000 ','a12,ww=',a12,ww
             write(6,*) 'iddt1 grav,sm1,dt2=', dt2,grav,sm1
             write(6,*) 'iddt1 xnpart,xctot = ',xnpart,xctot
             call flush(6)
             go to 21
           else
             go to 20
           endif
         endif
*
 21      continue
*
         beta = 2.0*atan(sqrt(sin2b2/(1.0-sin2b2)))
*
*        print*,' im2,ib,n=',im2,im1,nameb(im1),
*    *    ' xnpart,t,dt=',xnpart,time,1.d0/brem,' ipih=',ipihalf,
*    *      ' beta=',beta,' defang, dt0, a12=',defang, dt0, a12
*       determine velocities after interaction
*
         call newvel(im1,beta,vr2,vt2,vrn2,vtn2,sm2)
*       do not let binary with min rmin go down further
           if(im1.lt.0)return
*
         vr2 = vrn2
         vt2 = vtn2
*
         if(im2.gt.0)then
             vr(im2) = vrn2
             vt(im2) = vtn2
         end if
*
*       determine new positions of two interacting stars
 9999 continue
*
      lnewp = .true.
*
      if(lnewp)then
         call newpos(im1,ipt1,defang)
         if(im2.gt.0)call newpos(im2,ipt2,defang)
      end if
*
      return
*
      end
      subroutine relvel(im1,ww,vr2,vt2)
      implicit real*8(a-h,o-z),integer(i-k,m-n),logical(l)
*
*
*       determine relative velocity of two interacting stars
*       ----------------------------------------------------
*
*
      include 'compar.f'
      PARAMETER(NTAB=32)
      INTEGER IDUM2,IV,IY
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
      common/relaxb/cosfi,sinfi
*
*
*     real*8 ww,vr1,vr2,vtx1,vtx2,vty1,vty2,sigr2,sigt2,
*    &       wx,wy,wz
*
      real*4 ranf
*
*     integer im1,irad
*
*     compute radial and tangential velocity for a star from the shell
*     by use of Schwarzschild Boltzmann distribution function
*
      vr1 = vr(im1)
*
      fi = 2.d0*pi*ranf()
      sinfi = sin(fi)
      cosfi = cos(fi)
*
      vtx1 = vt(im1)
      vtx2 = vt2*cosfi
      vty1 = 0.0
      vty2 = vt2*sinfi
*
      wx = vtx2 - vtx1
      wy = vty2 - vty1
      wz = vr2 - vr1
      ww = wx*wx + wy*wy + wz*wz
      ww = sqrt(ww)
*
*
      return
*
      end
      subroutine newvel(im1,beta,vr2,vt2,vrn2,vtn2,sm2)
      implicit real*8(a-h,o-z),integer(i-k,m-n),logical(l)
*
*
*       determine new velocities of two interacting stars after encounter
*       -----------------------------------------------------------------
*       M. HENON 1971, Astrophysics and Space Science Vol. 14, 151
*       ----------------------------------------------------------
*
*
      include 'compar.f'
      PARAMETER(NTAB=32)
      INTEGER IDUM2,IV,IY
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
      common/relaxb/cosfi,sinfi
*
*
*     real*8 vr1,vr2,vrn2,vtn2,vtx1,vtx2,vty1,vty2,beta,wx,wy,
*    &       wz,ww,wx1,wy1,wz1,wx2,wy2,wz2,wp,psi,sinpsi,cospsi,
*    &       sbeta,cbeta,wnx,wny,wnz,vnx1,vnx2,vny1,vny2,sm1,sm2,
*    &       sm12,sm112,sm212
*
      real*4 ranf
*
*     integer im1
*
      sm1 = body1(im1)+body2(im1)
*
      sm12 = sm1 + sm2
      sm112 = sm1/sm12
      sm212 = sm2/sm12
      vr1 = vr(im1)
      vtx1 = vt(im1)
      vtx2 = vt2*cosfi
      vty1 = 0.0
      vty2 = vt2*sinfi
*
      wx = vtx2 - vtx1
      wy = vty2 - vty1
      wz = vr2 - vr1
      wp = wx*wx + wy*wy
      ww = wp + wz*wz
      ww = sqrt(ww)
      wp = sqrt(wp)
*
      wx1 = wy*ww/wp
      wy1 = -wx*ww/wp
      wz1 = 0.0
      wx2 = -wx*wz/wp
      wy2 = -wy*wz/wp
      wz2 = wp
*
      psi = 2.d0*pi*ranf()
      sinpsi = sin(psi)
      cospsi = cos(psi)
      sbeta = sin(beta)
      cbeta = cos(beta)
*
      wnx = wx*cbeta + wx1*sbeta*cospsi + wx2*sbeta*sinpsi
      wny = wy*cbeta + wy1*sbeta*cospsi + wy2*sbeta*sinpsi
      wnz = wz*cbeta + wz1*sbeta*cospsi + wz2*sbeta*sinpsi
*
      vnx1 = vtx1 - sm212*(wnx - wx)
      vnx2 = vtx2 + sm112*(wnx - wx)
      vny1 = vty1 - sm212*(wny - wy)
      vny2 = vty2 + sm112*(wny - wy)
      vr(im1) = vr1 - sm212*(wnz - wz)
      vrn2 = vr2 + sm112*(wnz - wz)
      vt(im1) = sqrt(vnx1*vnx1 + vny1*vny1)
      vtn2 = sqrt(vnx2*vnx2 + vny2*vny2)
*
*       do not allow binary with min rmin to sink further
*     if(nameb(im1).eq.imcrit)then
*     deltae = vr(im1)**2 + vt(im1)**2 - vr1**2 - vtx1**2
*
*        if(deltae.lt.0.d0)then
*        vr(im1) = vr1
*        vt(im1) = vtx1
*        im1 = -im1
*        end if
*     end if
*
      return
*
      end
*
*
      subroutine newpos(im1,ipt1,defang)
      implicit real*8(a-h,o-z),integer(i-k,m-n),logical(l)
*
*
*       determine new positions of two interacting stars after encounter
*       ----------------------------------------------------------------
*       M. HENON 1971, Astrophysics and Space Science Vol. 14, 151
*       ----------------------------------------------------------
*
*
      include 'compar.f'
      PARAMETER(NTAB=32)
      INTEGER IDUM2,IV,IY
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
      dimension dmold(njo)
      common/ttx/dmtold,dmttot
*
*
*     real*8 e1,e2,a1,a2,q1,dz,b1,c1,gmin1,gmax1,gmin2,gmax2,
*    &       rmin1,rmin2,rmax1,rmax2,rmaxz,s,u1,u2,v1,drds,q0,
*    &       rr1,rr2
*
      real*4 ranf
*
*     integer k,im1,n,i,ipt1,ipt2,irm1,irm2,ipp,nup,
*    &        m,imi1,imi2,ipo,kmin1,kmax1,kmin2,kmax2,
*    &        n95,ipp1,ipp2,ichang,imk1,imk2,maxk
*
*
      xmmin = 1.d-13
      k = ish(im1)
      j = ico(im1)
      xmbin = body1(im1) + body2(im1)
*
      ipt1 = 0
      qmax = 1.9d0
*
*       determine total energy and angular momentum for interacting stars
*
      k1 = k - 1
      dphdr=(phi(k)+phibin(k)-phi(k1)-phibin(k1))/(r(k)-r(k1))
      u1 = phi(k1) - phtid + phibin(k1) + dphdr*(rb(im1)-r(k1))
 999  continue
      e1 = u1 + 0.5*(vr(im1)**2 + vt(im1)**2)
      a1 = rb(im1)*vt(im1)
*       enforce that initial binaries are bound even in tidal field
      if(nmod.eq.0.and.e1.gt.0)then
       write(6,*) "newpos  ",' IB, RB,E1,U1= ', im1,rb(im1),e1,u1 
       vr(im1) = vr(im1)*0.999D0
       vt(im1) = vt(im1)*0.999D0
       goto 999
      end if
*
      phief1 = c12*(a1/r(k1))**2 + phi(k1) - phtid + phibin(k1)
      phief2 = c12*(a1/r(k))**2 + phi(k) - phtid + phibin(k)
*
*       take negative name as a tag
      if(e1.gt.0.0) then
          ipt1 = 1
          go to 80
      endif
*
*       find rmin and rmax for interacting stars
*
      imax=0
      imin=0
      phimin=1.d30
*
      lmin=.false.
*
      do 100 i = 2,nj
      phieff = c12*(a1/r(i))**2 + phi(i) - phtid + phibin(i)
*
      if(phieff.lt.phimin)then
      phimin=phieff
      ipmin=i
      end if
*
      if(.not.lmin.and.phieff.le.e1)then
      lmin= .true.
      imin = i
      end if
      if(.not.lmin.and.phieff.gt.e1)goto 100
      if(lmin.and.phieff.gt.e1)goto 101
      lmax = .true.
      imax = i
 100  continue
 101  continue
*
*       return in case of circular orbit
*
      phimin = 1.d30
      if(imax.eq.0)then
      imax1 = ipmin + 1
      imin1 = ipmin - 1
      lmin = .false.
      lmax = .false.
*
      radxm = r(imin1)
      phixm = phi(imin1) - phtid + phibin(imin1)
      phixp = phi(imax1) - phtid + phibin(imax1)
      phiefm = c12*(a1/radxm)**2 + phixm
      dr = r(imax1) - r(imin1)
      dphdr = (phixp - phixm)/dr
      dr = dr/1.d2
*
      do 150 ix=1,100
      radx = r(imin1)+dr*dble(ix)
*
      phix = phi(imin1)- phtid +phibin(imin1)+dphdr*(radx-r(imin1))
      phieff = c12*(a1/radx)**2 + phix
*
      if(phieff.lt.phimin)then
      phimin=phieff
      phxmin = phix
      rmin1 = radx
      rmax1 = radx
      end if
*
      if(.not.lmin.and.phieff.le.e1)then
      lmin=.true.
      rmin1 = radx
      imin = ix
      drdph2 = (radx-radxm)/(phieff - phiefm)
      end if
*
      if(lmin.and.(.not.lmax).and.phieff.gt.e1)then
      rmax1 = radx
      imax = ix
      drdph1 = (radx-radxm)/(phieff - phiefm)
      lmax = .true.
      goto 140
      end if
*
      radxm = radx
      phiefm = phieff
 150  continue
*
 140  continue
*
*       first treat case of nearly circular orbit
      if(lmax)then
c
      iivelx = 0
 160  continue
      s = ranf()
*
      if(s.lt.0.1d0)s=0.1d0
      if(s.gt.0.9d0)s=0.9d0
*
      rn1 = rmin1 + s*(rmax1 - rmin1)
      phix = phixm + dphdr*(rn1 - r(imin1))
      vt(im1) = a1/rn1
      velx = 2.d0*(e1 - phix) - vt(im1)**2
      iivelx = iivelx + 1
      if(velx.le.0.d0) then
        if(iivelx.gt.10000) then
          write(6,*) ' iivelx > 10000 '
          call flush(6)
          velx = 0.d0
        else
          go to 160
        endif
      endif
      vr(im1) = dsqrt(velx)
      else
*
      rn1 = rmin1
*     print*,' check: im1,lmin,lmax,rmin1,rmax1,imax,imin,imax1,imin1,
*    *        ipmin,rn1,phxmin,phimin,e1,u1,a1,vr,vt,rb =',
*    *        im1,lmin,lmax,rmin1,
*    *        rmax1,imax,imin,imax1,imin1,ipmin,rn1,phxmin,phimin,
*    *        e1,u1,a1,vr(im1),vt(im1),rb(im1)
      if(lmin) then
        vt(im1) = dsqrt(2.d0*(e1-phxmin))
        vr(im1)= 0.0d0
      else
        rn1 = rb(im1)
        vt(im1) = a1/rb(im1)
        ipmin = k
        rmin1 = rn1
        rmax1 = rn1
        aaa = (2.d0*(e1-u1) - vt(im1)*vt(im1))
        if(aaa.gt.0.d0) then
        vr(im1) = sqrt(aaa)
        else
        print*, ' aaa <  0    aaa= ',aaa
        endif
      endif
*
      end if
*       save quantities for later use
      rb(im1) = rn1
      ish(im1) = ipmin
      rbmin(im1) = rmin1
      if(lmax)then
      rbmax(im1) = rmax1
      else
      rbmax(im1) = rmin1
      end if
*
*       restrict number of mass updates for primordial binaries
*
*       update binary mass distribution
        dmold(1) = 0.d0
        dmbin(im1,1) = 0.d0
*
      do 251 i = 2,nj
        dmold(i) = dmbin(im1,i)
        dmbin(im1,i) = 0.d0
 251    continue
*
      if(lmax)then
*
      dvol = rmax1**3 - rmin1**3
      if(rmax1.lt.r(ipmin))then
         ibsta(im1) = ipmin
         ilen(im1) = 1
         dmbin(im1,ipmin) = xmbin
      else 
         if(rmin1.gt.r(ipmin))then
             ibsta(im1) = ipmin + 1
             ilen(im1) = 1
             dmbin(im1,ipmin+1) = xmbin
         else
             ibsta(im1) = ipmin
             ilen(im1) = 2
             dmbin(im1,ipmin) = xmbin*(r(ipmin)**3 - rmin1**3)/dvol
             dmbin(im1,ipmin+1) = xmbin*(rmax1**3 - r(ipmin)**3)/dvol
         end if
      end if
*
      else
*
      ibsta(im1) = ipmin
      ilen(im1) = 1
      dmbin(im1,ipmin) = xmbin
*
      end if
*
      dmttot = 0.d0
      dmtold = 0.d0
      xmrbin(1) = 0.d0
*
      do 400 i=2,nj
      dmttot = dmttot + dmbin(im1,i)
      dmtold = dmtold + dmold(i)
      xmrbin(i) = xmrbin(i) + dmttot - dmtold
      if(xmrbin(i).lt.xmmin)xmrbin(i) = 0.d0
      if(dabs(xmrbin(i)-ibin*xmbin).lt.xmmin)xmrbin(i) = ibin*xmbin
      if(dabs(xmrbin(i)-xmrbin(i-1)).lt.xmmin)xmrbin(i)=xmrbin(i-1)
 400  continue
*
      if(lmax)then
      print*,' Binary ',im1,' N=',nameb(im1),
     *    ' nearly circular orbit at ',
     *    ' rmin1,R,rmax1=',rmin1,rn1,rmax1
      else
      print*,' Binary ',im1,' N=',nameb(im1),
     *    ' circular orbit at ',
     *    ' rmin1,R,rmax1=',rmin1,rn1,rmax1
      end if
*
      isum = -1
      goto 500
*
      else
*
*
*       Interpolate PHIEFF for rmax to get better value
*
*
c      if(time.gt.690.d0) then
c        print*, ' imod,time,im1,imax,imin = ',
c     *             imod,time,im1,imax,imin
c        call flush(6)
c      endif
      imax1 = imax
      imin1 = imin
      i1 = imax
      i2 = imax + 1
      phief1 = c12*(a1/r(i1))**2 + phi(i1) - phtid + phibin(i1)
      if(imax.eq.nj)then
      xtmass = dexp(x(imr,nj)) + xmrbin(nj)
      rmax1 = -c12*xtmass/e1*(1.d0 + sqrt(1.d0+2.d0*a1**2*e1/xtmass**2))
      drdph1=rmax1**2/(xtmass - a1**2/rmax1)
      print*,' rmax>rnj, rmax, phi, phibin, phieff, xtmass,=',rmax1,
     *   -dexp(x(imr,nj))/rmax1,-xmrbin(nj)/rmax1,c12*(a1/rmax1)**2,
     *    xtmass
      print*,' im1, imod, imax, imin, e1=',im1,imod,imax,imin,e1
      print*,' -m/e1, a1**2e1/m**2)=',-xtmass/e1,a1**2*e1/xtmass**2
*      set position to maximum for balance
      rb(im1) = r(nj)
      rbmax(im1) = rmax1
      rbmin(im1) = r(imin)
      ibsta(im1) = imin
      ilen(im1) = imax - imin + 1
      ish(im1) = nj
      u11 = -grav*(dexp(x(imr,nj))+xmrbin(nj))/r(nj)
      v1 = 2.0*(e1 - u11) - a1*a1/r(nj)/r(nj)
      vt(im1) = 0.d0
      vr(im1) = dsqrt(v1)
      ipt1=2
      print*,' im1,ibsta,ilen=',im1,ibsta(im1),ilen(im1),rmax1,rmin1
      goto 80
      else
      phief2 = c12*(a1/r(i2))**2 + phi(i2) - phtid + phibin(i2)
      drdph1=(r(i2)-r(i1))/(phief2-phief1)
      rmax1=r(i1)+drdph1*(e1-phief1)
      end if
*
      is1 = i1
      is2 = i2
*
      i1 = imin - 1
      i2 = imin
      if(i1.eq.1)then
      i1 = i1 + 1
      i2 = i2 + 1
      end if
      ri1=r(i1)
*
 777  continue
*
      phief1 = c12*(a1/ri1)**2 + phi(i1) - phtid + phibin(i1)
      phief2 = c12*(a1/r(i2))**2 + phi(i2) - phtid + phibin(i2)
      drdph2=(r(i2)-ri1)/(phief2-phief1)
      rmin1=ri1+drdph2*(e1-phief1)
*
      if(phief1.lt.e1)then
      ri1=ri1/1.25d0
      goto 777
      end if
*
      end if
*
      if(-c32*(rmax1-rmin1)*drdph2.le.0.d0)then
      print*,' gmin1 error im1,rmin1,rmax1,drdph2=',
     *    im1,rmin1,rmax1,drdph2
      print*,' i1,ri1,phi,phibin,phieff=',i1,ri1,
     *      phi(i1),phibin(i1),c12*(a1/r(i1))**2
      print*,' i2,r,phi,phibin,phieff=',i2,r(i2),phi(i2),phibin(i2),
     *      c12*(a1/r(i2))**2
      call flush(6)
      end if
*     gmin1 = 2.0*a1*a1/rmin1**3 + 2.0*b1/rmin1**2  (dQ/dr|rmin)
      gmin1 = sqrt(-c32*(rmax1-rmin1)*drdph2)
*
      if(c32*(rmax1-rmin1)*drdph1.le.0.d0)then
      print*,' gmax1 error im1,rmin1,rmax1,drdph1=',
     *    im1,rmin1,rmax1,drdph1
      print*,' is1,r(is1),phi,phibin,phieff=',is1,r(is1),
     *      phi(is1),phibin(is1),c12*(a1/r(is1))**2
      call flush(6)
      print*,' is2,r(is2),phi,phibin,phieff=',is2,r(is2),
     *      phi(is2),phibin(is2),c12*(a1/r(is2))**2
      call flush(6)
      end if
*     gmax1 = 2.0*a1*a1/rmax1**3 + 2.0*b1/rmax1**2  (dQ/dr|rmax)
      gmax1 = sqrt(c32*(rmax1-rmin1)*drdph1)
*
*       determination of the new position
*
        isum = 0
        ivel = 0
*
   70 s = ranf()
*
        s = -1.d0 + 2.d0*s
*
          rn1 = 0.5*(rmin1 + rmax1) + 
     &                              0.25*(rmax1 - rmin1)*(3.0*s - s**3)
*
        do 200 i = 2,nj
        if(r(i).gt.rn1)goto 201
        irn1 = i
 200    continue
 201    continue
* 
*
*       check whether there are enough stars inside rn1
*      if(dexp(x(imr,imax1))/xmind(j).ge.4.d0)then
*      if(dexp(x(imr,irn1))/xmind(j).lt.4.d0)then
*        TEST   May 98
       if(dexp(x(imr,imax1))/xmind(j).ge.1.d-5)then
       if(dexp(x(imr,irn1))/xmind(j).lt.1.d-5)then
       isum = isum + 1
       if(isum.gt.10000) then
       print*,' rmin1,rmax1=',rmin1,rmax1
       print*,' imin1,imax1=',imin1,imax1
       call flush(6)
       go to 71
*       stop 'isum > 10000 '
       end if
       goto 70
       end if
       end if
*
 71    continue
*
        dphdr = (phi(irn1+1)+phibin(irn1+1)-phi(irn1)-phibin(irn1))/
     *          (r(irn1+1)-r(irn1))
        u1 = phi(irn1) - phtid + phibin(irn1) + dphdr*(rn1-r(irn1))
*
      q0 = qmax*ranf()*dmax1(gmin1,gmax1)
*
      v1 = 2.0*(e1 - u1) - a1*a1/rn1/rn1
*
      if(v1.le.0.0)then
*
      ivel = ivel + 1
*
      if(ivel.gt.10000)then
      print*,'ivel 10000 in newpos'
      call flush(6)
      v1 = 0.d0
      goto 450
      end if
*
      go to 70
      end if
*
      v1 = sqrt(v1)
*
      drds = 0.75*(rmax1 - rmin1)*(1.0 - s*s)
      q1 = drds/v1
*
      if(q1.gt.qmax*dmax1(gmin1,gmax1)) then
        write(6,*) 'r q1 > q  iseed,q1,gmin1,gmax1=',iseed,q1,gmin1,gmax1
        call flush(6)
        qmax = 2.0d0*qmax
        go to 70
      endif
*  
      if(q0.gt.q1) go to 70
*
 450  continue
*
      vt(im1) = a1/rn1
      vr(im1) = v1
      rb(im1) = rn1
      ish(im1) = irn1 + 1
*       save rmin, rmax for further use
      rbmin(im1) = rmin1
      rbmax(im1) = rmax1
*
*    NO!
*     if(.not.ls(5).or.mod(imod,10).eq.0)
*
      call orbit(im1,imin1,imax1,rmax1,rmin1,gmin1,gmax1,e1,a1)
*
*
 500  continue
*
      if(iev(im1).ne.0)then
      ekt = xmind(1)*dexp(x(i20+1,1)-x(i00+1,1))
      write(60,600)ibint,ibin,nameb(im1),nmod,time,rmin1,rmax1,rn1,
     *   e1,a1,defang,imin1,imax1,isum,rcore,eb(im1),ekt,ecc(im1),
     *   body1(im1)+body2(im1)
      call flush(60)
      end if
 600  format(2i5,' bin',i5,1p,i8,7d12.4,3i5,5d12.4)
*
*      write(65,650)ibint,ibin,nameb(im1),nmod,time,rmin1,rmax1,rn1,
*    *    torb(im1),dtrx,tbin(im1),vr(im1),vt(im1),iev(im1),
*    *    imin,imax
*650  format(3i5,1p,i8,9d12.4,3i3)
*
 80   continue
*
      return
*
      end
*
*
*
      subroutine kick(im1,vb)
      implicit real*8(a-h,o-z),integer(i-k,m-n),logical(l)
*
*      compute components of the 'kick' velocity (gained during
*      --------------------------------------------------------
*      interaction between binary and field star) of binary or
*      ------------------------------------------------------- 
*      star and new velocity of interacting objects
*      --------------------------------------------
*
      include 'compar.f'
      PARAMETER(NTAB=32)
      INTEGER IDUM2,IV,IY
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
*
*     real*8 q1,q2,costeta,sinteta,cosffi,sinffi,vb,b1,c1,dv,vkk
*
      real*4 ranf
*
*     integer im1
*
      bodyb = body1(im1) + body2(im1)
*
*      used formulae from Stodolkiewicz Acta Astronomica 1986, 36, p19
*
      q1 =ranf() - 0.5
      costeta = -2.d0*q1
      sinteta = sqrt(1.d0 - costeta**2)
*
      q2 = 2.d0*pi*ranf()
      sinffi = sin(q2)
      cosffi = cos(q2)
      b1 = 2.d0*(vr(im1)*costeta + vt(im1)*sinteta*cosffi)
      c1 = -2.d0*vb/bodyb
      dv = 0.5d0*(-b1 + sqrt(b1*b1 - 4.d0*c1))
      vr(im1) = vr(im1) + dv*costeta
      b1 = vt(im1) + dv*sinteta*cosffi
      c1 = dv*sinteta*sinffi
      vt(im1) = sqrt(b1*b1 + c1*c1)
*
      return
      end
*
      subroutine kicks(vb,vr1,vt1,xm)
      implicit real*8(a-h,o-z),integer(i-k,m-n),logical(l)
*
*      compute components of the 'kick' velocity (gained during
*      --------------------------------------------------------
*      interaction between binary and field star) of binary or
*      -------------------------------------------------------
*      star and new velocity of interacting objects
*      --------------------------------------------
*
      include 'compar.f'
      PARAMETER(NTAB=32)
      INTEGER IDUM2,IV,IY
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
*
*     real*8 q1,q2,costeta,sinteta,cosffi,sinffi,vb,b1,c1,dv,vkk
*
      real*4 ranf
*
*     integer im1
*
*      used formulae from Stodolkiewicz Acta Astronomica 1986, 36, p19
*
      q1 =ranf() - 0.5
      costeta = -2.d0*q1
      sinteta = sqrt(1.d0 - costeta**2)
*
      q2 = 2.d0*pi*ranf()
      sinffi = sin(q2)
      cosffi = cos(q2)
      b1 = 2.d0*(vr1*costeta + vt1*sinteta*cosffi)
      c1 = -2.d0*vb/xm
      dv = 0.5d0*(-b1 + sqrt(b1*b1 - 4.d0*c1))
      vr1 = vr1 + dv*costeta
      b1 = vt1 + dv*sinteta*cosffi
      c1 = dv*sinteta*sinffi
      vt1 = sqrt(b1*b1 + c1*c1)
*
      return
      end
*
      SUBROUTINE PHIMAS(IPOT)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*
*
*         compute smooth potential for all particles
*         ------------------------------------------
*
*
      INCLUDE 'compar.f'
      REAL*8 RAUX(NBINO)
*       Flag IPOT for single star potential calculation in
*       case of binary formation
      IF(IPOT.EQ.1)THEN
*       Initialization in case of no binary
*
*          Start at outer boundary
         PHI(NJ)=-GRAV*DEXP(X(IMR,NJ))/R(NJ)-PHTID
*
         DO 340 I=1,NJ-1
 340     PHI(I)=0.D0
*
         DO 350 J=1,NCOMP
         FM(J,1)=0.D0
         DO 360 I=2,NJ
         DER3=R(I)**3-R(I-1)**3
         FM(J,I)=FM(J,I-1)+DEXP(X(I00+J,I))*DER3
 360     CONTINUE
*
         DER3=R(NJ)**3-R(NJ-1)**3
         RAV=(R(NJ)+R(NJ-1))/2.D0
         FP(J,NJ)=DEXP(X(I00+J,NJ))*DER3/RAV
         DO 370 I=NJ-1,2,-1
         DER3=R(I)**3-R(I-1)**3
         RAV=(R(I)+R(I-1))/2.D0
         FP(J,I)=FP(J,I+1)+DEXP(X(I00+J,I))*DER3/RAV
 370     CONTINUE
*
         PHI(1)=PHI(1)-PI43*GRAV*FP(J,2)-PHTID
         RFAC=R(3)/R(2)
         IF(LS(19))PHI(1)=PHI(1)-GRAV*XMHOLE/R(2)*RFAC
         DO 380 I=2,NJ-1
         RAV=(R(I)+R(I-1))/2.D0
         PHI(I)=PHI(I)-PI43*GRAV*(FM(J,I)/R(I)+FP(J,I+1))-PHTID
         IF(LS(19))PHI(I)=PHI(I)-GRAV*XMHOLE/R(I)
 380     CONTINUE
*
 350     CONTINUE
*
         END IF
*
      IF(IBIN.EQ.0)THEN
      DO 1 I = 1,NJ
      XMRBIN(I) = 0.D0
 1    PHIBIN(I) = 0.D0
      RETURN
      END IF
*
      IF(IBIN.EQ.1)THEN
         INAME(1) = 1
      ELSE
*
         DO 2 IB = 1,IBIN
         RAUX(IB) = RB(IB)
 2       INAME(IB) = IB
*
         CALL SORT1(IBIN,RAUX,INAME)
      END IF
*
*       Determination of smooth potential for gaseous model
*
      PHIBIN(NJ) = -GRAV*XMRBIN(NJ)/R(NJ)
*
      DO 70 I = NJ-1,2,-1
*
      IP = I+1
      PHIBIN(I) = PHIBIN(IP) - GRAV*XMRBIN(I)*(1.D0/R(I)-1.D0/R(IP))
 70   CONTINUE
*
      PHIBIN(1) = PHIBIN(2)
*     
      RETURN
      END
*
      SUBROUTINE SORT1(N,RA,RB)
*
*
*       Heapsort method (Press p. 231).
*       -------------------------------
*
      INTEGER  RB,RRB
      REAL*8 RA(N),RRA
      DIMENSION  RB(N)
*
*
      L = N/2+1
      IR=N
   10 CONTINUE
      IF(L.GT.1)THEN
	  L=L-1
	  RRA=RA(L)
	  RRB=RB(L)
      ELSE
          RRA=RA(IR)
	  RRB=RB(IR)
	  RA(IR)=RA(1)
	  RB(IR)=RB(1)
          IR=IR-1
	  IF(IR.EQ.1)THEN
	      RA(1)=RRA
	      RB(1)=RRB
	      RETURN
          END IF
      END IF
      I=L
      J=L+L
   20 IF(J.LE.IR)THEN
	  IF(J.LT.IR)THEN
	      IF(RA(J).LT.RA(J+1))J=J+1
          END IF
	  IF(RRA.LT.RA(J))THEN
	      RA(I)=RA(J)
	      RB(I)=RB(J)
	      I=J
	      J=J+J
          ELSE
	      J=IR+1
          END IF
          GO TO 20
      END IF
      RA(I)=RRA
      RB(I)=RRB
      GO TO 10
*
      END
