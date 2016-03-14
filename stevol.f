      subroutine stevol
      implicit real*8 (a-h,o-z),integer(i-k,m-n),logical(l)
      include 'compar.f'
      dimension vxmind(ncompo),vxm(ncompo),xm(ncompo)
*
*  rescale time from calculation units to years       
      t=time*cunt/year
*
*     if (t.lt.0.1) return
*
      do 900 j=1,ncomp
      xm(j)=xmind(j)*cunm/sunm
 900  vxm(j)=xm(j)
* input parameter maximum stellar lifetime tupin
      tup = tupin
*
      do 10 i=1,ncomp
      if (t.gt.tau(i).and.t.lt.tup) xm(i)=xmrem(i) +
     &   (xm(i)-xmrem(i))*(tup-t) / (tup-tau(i))
      if (t.ge.tau(i)) xm(i)=xmrem(i)
      print*,' i,t,tau=',i,t,tup,tau(i)
      print*,' xm,vxm=',i,xm(i)
      tup = tau(i)
   10 continue
*
      do 950 j=1,ncomp
      xmind(j)=xm(j)*sunm/cunm
 950  vxmind(j)=vxm(j)*sunm/cunm
*
*     if(i.eq.10)then
*     write(6,6)t,(xm(j),j=1,ncomp)
*   6 format(1x,'    t, xm=',1pd12.3,0p,10(1x,f5.2))
*     write(6,7)time,(xmind(j),j=1,ncomp)
*   7 format(1x,' t, xmind=',1pd12.3,1p,10(1x,d9.2))
*     write(6,8)time,(xmind(j)-vxmind(j),j=1,ncomp)
*   8 format(1x,' t, dmind=',1pd12.3,1p,10(1x,d9.2))
*     end if
*
      print*,' AFTER MASSES 2 '
      do 1000 j=1,ncomp
*
*  calculate dln m /dt ; brem = 1/deltat
*
      dlnmdt=(1.d0-xmind(j)/vxmind(j))*brem
*
*     if(iwarn.lt.1000.and.dabs(1.d0-xmind(j)/vxmind(j)).gt.corr)then
      iwarn=iwarn+1
      print*,' timestep warning t=',t,' years comp. ',j,' i=',i,' dt=',
     *  cunt/year/brem,' dm=',(xmind(j)-vxmind(j))*cunm/sunm,
     *  ' dlnmdt=',dlnmdt
      print*,' xmind,vxmind=',j,xmind(j),vxmind(j)
      print*,' xm,vxm=',j,xm(j),vxm(j)
*     end if
*
*  equation for dlnrho/dt; dlnpr/dt; dlnpt/dt
*
      g(i00+j)=g(i00+j)-dlnmdt
      g(i20+j)=g(i20+j)-dlnmdt
      g(i02+j)=g(i02+j)-dlnmdt
*
 1000 continue
*
      return
*
      end
