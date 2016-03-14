      SUBROUTINE PRINT(IPI)
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
      CHARACTER*5 AB
      CHARACTER*64 FSTR
      CHARACTER*64 MY_FSTR, MY_FSTR_HEAD
      CHARACTER*10 TCHAR,BLANK
         INCLUDE 'compar.f'
         INCLUDE 'equiv.f'
      PARAMETER(NTAB=32)
      INTEGER IDUM2,K,IV,IY
      REAL*8 DMSAVE(NBINO*50)
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
*
      CHARACTER*1 MC(NCOMPO,NJO)
      REAL XYZ(NJO/10,2),STAT(3)
      REAL  RES(NJO/10,2),RFL
      DIMENSION NFR(NJO/10),NFUR(NJO/10),M(NJO),TCHAR(8,NCOMPO)
C
          IF(LS(10))IZ00=1
          IF(IZ00.EQ.0)THEN
          IF(LS(3))THEN
          OPEN(1,FILE='fort.1',FORM='UNFORMATTED',
     *   STATUS='UNKNOWN',IOSTAT=IERR)
          IF(IERR.NE.0)PRINT*,' OPEN 1 IERR=',IERR
          END IF
          OPEN(2,FILE='fort.2',FORM='UNFORMATTED',
     *   STATUS='UNKNOWN',IOSTAT=IERR)
          IF(IERR.NE.0)PRINT*,' OPEN 2 IERR=',IERR
          OPEN(4,FILE='fort.4',FORM='UNFORMATTED',
     *   STATUS='UNKNOWN',IOSTAT=IERR)
          IF(IERR.NE.0)PRINT*,' OPEN 4 IERR=',IERR
          IZ00=1
          END IF
C
          IF(IPI.LE.0)GOTO 10
C
          IEI(1)=MODA
          IEI(2)=NMOD
          IEI(8)=IBINT
          IEI(9)=IBIN
          IEI(10)=IREC
*         IEI(11)=IBREC see below
          IEI(12)=ITER
          IEI(13)=NJ
          IEI(14)=NCOMP
          IEI(15)=IMCRIT
C
      AEI(2)=TIME
      AEI(4)=TTRX
      AEI(3)=BREM
C         Note that AEI(4) to AEI(8) are used in ZOMBI
      AEI(9)=XN3TOT
      AEI(10)=XMHOLE
      AEI(11)=DMTOT
      AEI(12)=XN3B
      AEI(13)=CS(17)
      AEI(14)=DPTOT
      AEI(15)=XMERR
      AEI(16)=XEERR
      AEI(17)=PHTID
      AEI(18)=RTID0
      AEI(19)=XMTOT0
      AEI(20)=XMBIN0
      AEI(21)=XMERR2
      AEI(22)=XEERR2
      AEI(23)=XMERR3
      AEI(24)=XEERR3
      AEI(25)=XMERR4
      AEI(26)=XEERR4
C
      IF(LS(3))THEN
      IF(IBIN.GT.0.AND.MOD(MODA,5).EQ.0)THEN
      IBREC=IBREC+1
      WRITE(1,ERR=998,IOSTAT=IWRERR)MODA,IBREC
      WRITE(1,ERR=998,IOSTAT=IWRERR)(BEI(K),K=1,NBEI)
      WRITE(1,ERR=998,IOSTAT=IWRERR)(IBEI(K),K=1,NIBEI)
      CALL FLUSH(1)
      REWIND 11
      WRITE(11,ERR=998,IOSTAT=IWRERR)MODA,IBREC
      WRITE(11,ERR=998,IOSTAT=IWRERR)(BEI(K),K=1,NBEI)
      WRITE(11,ERR=998,IOSTAT=IWRERR)(IBEI(K),K=1,NIBEI)
      CALL FLUSH(11)
      ILENT = 0
      DO 2999 IB=1,IBIN
 2999 ILENT = ILENT + ILEN(IB)
      IF(ILENT.GT.50*NBINO)THEN
      PRINT*,'ILENT too large'
      PRINT*,' ILEN=',(ILEN(IB),IB=1,IBIN)
      PRINT*,' ISTA=',(IBSTA(IB),IB=1,IBIN)
      CALL FLUSH(6)
      STOP
      END IF
*
      ILENT=0
      EKT=XMIND(1)*DEXP(X(I20+1,1)-X(I00+1,1))
      ITEN = 0
      I60 = 0
      IONE = 0
      WRITE(61,599)MODA,IBIN
 599  FORMAT(' -----------------MODA=',I8,' IBIN=',I8,'----------')
*
      DO 1000 IB=1,IBIN
      DO 1005 I=1,ILEN(IB)
      IBIX = IBSTA(IB)+I
      DMSAVE(ILENT+I)=DMBIN(IB,IBIX)
*     PRINT*,' write to DMSAVE(',ILENT+I,') DMBIN(',IB,',',IBIX,')',
*    *  ' value =',DMBIN(IB,IBIX),' bef=',DMBIN(IB,IBIX-1), 
*    *  ' aft=',DMBIN(IB,IBIX+1)
*     CALL FLUSH(6)
 1005 CONTINUE
*      Write Binary data on separate file
      RMINB=R(IBSTA(IB)+1)
      RMAXB=R(IBSTA(IB)+ILEN(IB))
      K=ISH(IB)
      K1 = K - 1
      DPHDR=(PHI(K)+PHIBIN(K)-PHI(K1)-PHIBIN(K1))/(R(K)-R(K1))
      U1 = PHI(K1) + PHIBIN(K1) + DPHDR*(RB(IB)-R(K1))
      E1 = U1 + 0.5D0*(VR(IB)**2 + VT(IB)**2)
      A1 = RB(IB)*VT(IB)
      DTIME = 1.D0/BREM
      IF(TORB(IB).GT.DTIME)THEN
      IONE=IONE+1
      IF(TORB(IB).GT.1.D1*DTIME)ITEN=ITEN+1
      IF(TORB(IB).GT.6.D1*DTIME)I60=I60+1
      END IF
*
      WRITE(61,600)IB,NAMEB(IB),TIME,RBMIN(IB),RBMAX(IB),RB(IB),
     * E1,A1,U1,RCORE,
     * EB(IB),EKT,ECC(IB),SEMIA(IB),BODY1(IB),BODY2(IB),SIZE1(IB),
     * SIZE2(IB),TORB(IB)*BREM
 600  FORMAT(2I8,1P,18D12.4)
 1000 ILENT = ILENT + ILEN(IB)
      CALL FLUSH(6)
      IF(ILENT.GT.50*NBINO)STOP 'ILENT too large'
      WRITE(1,ERR=998,IOSTAT=IWRERR)ILENT
      WRITE(1,ERR=998,IOSTAT=IWRERR)(DMSAVE(K),K=1,ILENT)
      WRITE(11,ERR=998,IOSTAT=IWRERR)ILENT
      WRITE(11,ERR=998,IOSTAT=IWRERR)(DMSAVE(K),K=1,ILENT)
      CLOSE(11)
      PRINT*,' Binaries saved IBIN=',IBIN,' IBREC=',IBREC
      PRINT*,' ILENT written=',ILENT,' IONE,I60=',IONE,I60
      END IF
*         For CRAY get random number sequence
      IFAIL = RANGET(IDUM2)
      WRITE(91,*)MODA,IDUM,IY,IDUM2,(IV(KK),KK=1,NTAB)
      CALL FLUSH(91)
          PRINT*,' Random Numbers Saved:       IDUM,IY,IDUM2=',
     *     IDUM,IY,IDUM2,' IV=',(IV(KK),KK=1,NTAB)
          END IF
*
          IEI(11) = IBREC
*
          WRITE(2,ERR=999,IOSTAT=IWRERR)(EI(K),K=1,NW/2)
          WRITE(2,ERR=999,IOSTAT=IWRERR)(EI(K),K=NW/2+1,NW)
*        Write last model on file for Restart
          REWIND 4
          WRITE(4,ERR=999,IOSTAT=IWRER4)(EI(K),K=1,NW/2)
          WRITE(4,ERR=999,IOSTAT=IWRER4)(EI(K),K=NW/2+1,NW)
C
      WRITE(6,1) MODA,NMOD,TIME,TTRX,
     *  TIME*CUNT/YEAR,1.D0/BREM,1.D0/BREM*CUNT/YEAR
 1    FORMAT(' PRINT: MODA=',I4,' NMOD=',I8,' TNB=',1P,D9.2,
     * ' TTRX=',D9.2,' T=',D9.2,' YRS.',/,
     * ' DTNB=',D9.2,' DT=',D9.2,' YRS.')
          IF(LS(19))WRITE(6,2)XMHOLE*CUNM/SUNM,DMTOT*CUNM/SUNM
 2    FORMAT(' PRINT WITH BLACK HOLE MASS=',1P,E9.2,
     *  ' DM=',E9.2,' SOL.MASSES ')
C
          IF(IPI.EQ.2)RETURN
C******************************************************************
C      END WRITE ON FILE/START WRITE OUTPUT
C******************************************************************
 10   CONTINUE
C*****************BEGIN PRINT LARGE MODEL OUTPUT*******************
       NS=NWRITE
       IF(LS(23))NS=2*NWRITE
C********Calculate total filling factor for each grid point*******
C
C********FP(I,J)**Mass of Comp. I inside R(J)******************
C********W(J)**Total Density at R(J)********************
      RFAC=R(3)/R(2)
      DO 8884 J=1,NCOMP
 8884 FP(J,1)=DLOG(PI43)+X(I00+J,1)+3.D0*DLOG(R(2)/RFAC)
*
      DO 8885 I=2,NJ
      W(I)=0.D0
*
      DO 8886 J=1,NCOMP
*
         IF(I.EQ.2)THEN
         DER=DLOG(RFAC)
         RAV=C12*DLOG(R(I)*R(I)/RFAC)
         RM=R(I)/RFAC
         ELSE
         DER=DLOG(R(I)/R(I-1))
         RAV=C12*DLOG(R(I)*R(I-1))
         RM=R(I-1)
         END IF
*
         IMASS=0
         AMAX=0.D0
 130     IMASS=IMASS+1
         IF(IMASS.EQ.1)THEN
         EFAC=DEXP(3.D0*DLOG(RM)-FP(J,I-1))
         VMASS=0.D0
         ELSE
         XMAV=C12*(FP(J,I)+FP(J,I-1))
         EFAC=DEXP(3.D0*RAV-XMAV)
         VMASS=FP(J,I)
         END IF
C
      FP(J,I)=FP(J,I-1)+PI4*DEXP(X(I00+J,I))*DER*EFAC
         XMCORR=FP(J,I)-VMASS
         IF(IMASS.GT.100)THEN
         PRINT*,' Failure of initial mass iteration at:'
         PRINT*,' J=',J,' I=',I
         XMCORR=0.D0
         END IF
         IF(DABS(XMCORR).GT.EPS)GOTO 130
*
         W(I)=W(I)+DEXP(X(I00+J,I))
 8886 CONTINUE
*
 8885 CONTINUE
*
      W(1)=0.D0
      DO 8888 J=1,NCOMP
 8888 W(1)=W(1)+DEXP(X(I00+J,1))
*
      RHOTOT=0.D0
      DO 8887 J=1,NCOMP
      RHOTOT=RHOTOT+DEXP(XOUT(J))
      HXP(I10+J)=0.D0
      HXP(I30+J)=0.D0
      HXP(I12+J)=0.D0
 8887 CONTINUE
*
      IRUN=0
      DO 8890 I=1,NJ
*       Prepare significant output for outer mass shells
      DO 8891 J=1,NCOMP
      FP(J,I)=DEXP(FP(J,I))
      MC(J,I)=' '
      IF(1.D0-FP(J,I).LT.1.D-2)THEN
      FP(J,I)=1.D0-FP(J,I)
      MC(J,I)='*'
      END IF
 8891 CONTINUE
*
      LWRITE=I.LT.10    .OR.
     *       MOD(I,NWRITE).EQ.0.AND.I.GE.10.AND.I.LT.NJ-NWRITE  .OR.
     *       I.GE.NJ-NJ/10
*
      IF(LWRITE)THEN
      IRUN=IRUN+1
      M(IRUN)=I
      END IF
*
 8890 CONTINUE
*
      MMAX=IRUN
C**********************************************************************
C
C      Output on 11 columns distributed for
C       NCOMP<=2 : 4 quantities per line  (900. formats)
C       NCOMP<=4 : 2 quantities per line  (1000. formats)
C       NCOMP>4 :  1 quantity per line (2000. formats)
C       NCOMP>9 :  splitting of lines occurs, not yet tested
C*****************************************************************
C Write Titel Strings in TCHAR-Variable
C
      DO 9888 K=1,NCOMP
      WRITE(TCHAR(1,K),121)K
 121  FORMAT(' MASS ',I2,2X)
      WRITE(TCHAR(2,K),122)K
 122  FORMAT(' RHO  ',I2,2X)
      WRITE(TCHAR(3,K),123)K
 123  FORMAT(' SIGR   ',I2)
      WRITE(TCHAR(4,K),124)K
 124  FORMAT(' SIGT   ',I2)
      WRITE(TCHAR(5,K),125)K
 125  FORMAT(' V-* ',I2,3X)
      WRITE(TCHAR(6,K),126)K
 126  FORMAT('LUMIN. ',I2,1X)
      WRITE(TCHAR(7,K),127)K
 127  FORMAT(' VR-ER ',I2,1X)
      WRITE(TCHAR(8,K),128)K
 128  FORMAT(' VR-ET ',I2,1X)
 9888 CONTINUE
      BLANK='          '
      FSTR='(1X,I4,1P,13(D11.4,A1))'
      write(MY_FSTR,3000) NCOMP+2
 3000 format('(1X,I4,1P,',I2,'(D11.4,A1))')
C*****************************************************************
C Compact Output if NCOMP.LE.4
C
      LCOMP=.FALSE.
      LCOMP2=.FALSE.
      IF(NCOMP.LE.4)LCOMP=.TRUE.
      IF(NCOMP.LE.1)LCOMP2=.TRUE.
C*****************************************************************
C Start Compact Output
C
      IF(LCOMP)THEN
C-----------------------------------------------------------------
C  Start Very Compact Output
C
      IF(LCOMP2)THEN
C
      NC=7-4*NCOMP
      NC=NC+7
      WRITE(6,140)' RADIUS   ',' TOTMASS  ',' XMRBIN   ',' VMRBIN   ',
     *  ' RHOBIN   ',(TCHAR(1,IK),IK=1,NCOMP),' TOTRHO   ',
     *  ((TCHAR(IJ,IK),IK=1,NCOMP),IJ=2,4),(BLANK,IK=1,NC)
 140  FORMAT(//,4X,20A10)
C-----------Masses, Densities, Radial and Tangential Vel.Dispersions---
      DO 9001 J=1,MMAX
 9001 WRITE(6,FSTR) M(J),R(M(J))/UNR,' ',DEXP(X(IMR,M(J)))/UNM,' ',
     *  XMRBIN(M(J))/UNM,' ',VMRBIN(M(J))/UNM,' ',RHOBIN(M(J))/UNRH,
     *  (MC(IJ,M(J)),FP(IJ,M(J))/UNM,IJ=1,NCOMP),' ',W(M(J))/UNRH,' ',
     *  (DEXP(X(I00+IJ,M(J)))/UNRH,' ',IJ=1,NCOMP),
     *  (DEXP((X(I20+IJ,M(J))-X(I00+IJ,M(J)))/2.D0)/UNU,' ',IJ=1,NCOMP),
     *  (DEXP((X(I02+IJ,M(J))-X(I00+IJ,M(J)))/2.D0)/UNU,' ',IJ=1,NCOMP)
*       Extrapolated outermost point
      RFAC=R(NJ)/R(NJ-1)
      WRITE(6,FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',DEXP(X(IMR,NJ))/UNM,
     *   ' ',0.D0,' ',0.D0,
     *  ('*',XMTOT(IJ)/UNM,IJ=1,NCOMP),' ',RHOTOT/UNRH,' ',
     *  (DEXP(XOUT(IJ))/UNRH,' ',IJ=1,NCOMP),
     *  (DEXP((PROUT(IJ)-XOUT(IJ))/2.D0)/UNU,' ',IJ=1,NCOMP),
     *  (DEXP((PTOUT(IJ)-XOUT(IJ))/2.D0)/UNU,' ',IJ=1,NCOMP)
C
      NC=11-2*NCOMP
      NC=NC+7
      WRITE(6,140)' RADIUS   ',
     *  ((TCHAR(IJ,IK),IK=1,NCOMP),IJ=5,8),(BLANK,IK=1,NC)
C---------Velocities,Luminosities,I30-and-I12-Fluxes-----------------
      DO 9002 J=1,MMAX
      IF(M(J)-1.GT.2)THEN
         SIG2=DEXP(X(I20+1,M(J))-X(I00+1,M(J)))*
     *               (1.D0+2.D0*DEXP(X(I02+1,M(J))-X(I20+1,M(J))))
         SIGD=DEXP(X(I20+1,M(J)-1)-X(I00+1,M(J)-1))*
     *               (1.D0+2.D0*DEXP(X(I02+1,M(J)-1)-X(I20+1,M(J)-1)))
         TRXX=(SIG2/3.D0)**C32/
     *    CS(60)/XMIND(1)/XCOUL(1)/DEXP(X(I00+1,M(J)))
      DSDR=(SIG2-SIGD)/(R(M(J))-R(M(J)-1))/3.D0
      VE  = -CS(21)/(4.D0*GRAV*DEXP(X(I00+1,M(J)))*TRXX)*DSDR
      TE = R(M(J))/VE
      ELSE
      VE=0.D0
      TE=0.D0
      END IF
 9002 WRITE(6,FSTR) M(J),R(M(J))/UNR,' ',
     * (X(I10+IJ,M(J))/UNU,' ',IJ=1,NCOMP),
     *(C12*(3.D0*DEXP(X(I20+IJ,M(J)))*X(I30+IJ,M(J))+
     * 2.D0*DEXP(X(I02+IJ,M(J)))*X(I12+IJ,M(J)))*PI4*R(M(J))**2/UNL,
     *          ' ',IJ=1,NCOMP),
     * (X(I30+IJ,M(J))/UNU,' ',IJ=1,NCOMP),
     * (X(I12+IJ,M(J))/UNU,' ',IJ=1,NCOMP),VE,' ',TE,' ',TRXX
*
      WRITE(6,FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',
     * (HXP(I10+IJ)/UNU,' ',IJ=1,NCOMP),
     *(C12*(3.D0*DEXP(PROUT(IJ))*HXP(I30+IJ)+
     * 2.D0*DEXP(PTOUT(IJ))*HXP(I12+IJ))*PI4*R(NJ+1)**2/UNL,
     *          ' ',IJ=1,NCOMP),
     * (HXP(I30+IJ)/UNU,' ',IJ=1,NCOMP),
     * (HXP(I12+IJ)/UNU,' ',IJ=1,NCOMP)
*
      ELSE
C
      NC=9-2*NCOMP
      NC=NC+7
      WRITE(6,140)' RADIUS   ',' TOTMASS  ',
     *  (TCHAR(1,IK),IK=1,NCOMP),' TOTRHO   ',
     *  (TCHAR(2,IK),IK=1,NCOMP),(BLANK,IK=1,NC)
C-----------Masses and Densities-------
      DO 10001 J=1,MMAX
10001 WRITE(6,FSTR) M(J),R(M(J))/UNR,' ',DEXP(X(IMR,M(J)))/UNM,
     *  (MC(IJ,M(J)),FP(IJ,M(J))/UNM,IJ=1,NCOMP),' ',W(M(J))/UNRH,' ',
     *  (DEXP(X(I00+IJ,M(J)))/UNRH,' ',IJ=1,NCOMP)
*       Extrapolated outermost point
      RFAC=R(NJ)/R(NJ-1)
      WRITE(6,FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',DEXP(X(IMR,NJ))/UNM,
     *  ('*',XMTOT(IJ)/UNM,IJ=1,NCOMP),' ',RHOTOT/UNRH,' ',
     *  (DEXP(XOUT(IJ))/UNRH,' ',IJ=1,NCOMP)
*
      NC=11-2*NCOMP
      NC=NC+7
      WRITE(6,140)' RADIUS   ',
     * ((TCHAR(IJ,IK),IK=1,NCOMP),IJ=3,4),(BLANK,IJ=1,NC)
C----------Radial/Tangential Velocity Disp.-------
      DO 10002 J=1,MMAX
10002 WRITE(6,FSTR) M(J),R(M(J))/UNR,' ',
     * (DEXP((X(I20+IJ,M(J))-X(I00+IJ,M(J)))/2.D0)/UNU,' ',IJ=1,NCOMP),
     * (DEXP((X(I02+IJ,M(J))-X(I00+IJ,M(J)))/2.D0)/UNU,' ',IJ=1,NCOMP)
*       Extrapolated outermost point
      WRITE(6,FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',
     *  (DEXP((PROUT(IJ)-XOUT(IJ))/2.D0)/UNU,' ',IJ=1,NCOMP),
     *  (DEXP((PTOUT(IJ)-XOUT(IJ))/2.D0)/UNU,' ',IJ=1,NCOMP)
*
      NC=11-2*NCOMP
      NC=NC+7
      WRITE(6,140)' RADIUS   ',
     * ((TCHAR(IJ,IK),IK=1,NCOMP),IJ=5,6),(BLANK,IJ=1,NC)
C-----------Velocities/Luminosities--------------
      DO 10003 J=1,MMAX
10003 WRITE(6,FSTR) M(J),R(M(J))/UNR,' ',
     * (X(I10+IJ,M(J))/UNU,' ',IJ=1,NCOMP),
     *(C12*(3.D0*DEXP(X(I20+IJ,M(J)))*X(I30+IJ,M(J))+
     * 2.D0*DEXP(X(I02+IJ,M(J)))*X(I12+IJ,M(J)))*PI4*R(M(J))**2/UNL,
     *          ' ',IJ=1,NCOMP)
*       Extrapolated outermost point
      WRITE(6,FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',
     * (HXP(I10+IJ)/UNU,' ',IJ=1,NCOMP),
     *(C12*(3.D0*DEXP(PROUT(IJ))*HXP(I30+IJ)+
     * 2.D0*DEXP(PTOUT(IJ))*HXP(I12+IJ))*PI4*R(NJ+1)**2/UNL,
     *          ' ',IJ=1,NCOMP)
C
      NC=11-2*NCOMP
      NC=NC+7
      WRITE(6,140)' RADIUS   ',
     * ((TCHAR(IJ,IK),IK=1,NCOMP),IJ=7,8),(BLANK,IJ=1,NC)
C-----------Energy Fluxes------------
      DO 10004 J=1,MMAX
10004 WRITE(6,FSTR) M(J),R(M(J))/UNR,' ',
     * (X(I30+IJ,M(J))/UNU,' ',IJ=1,NCOMP),
     * (X(I12+IJ,M(J))/UNU,' ',IJ=1,NCOMP)
*       Extrapolated outermost point
      RFAC=R(NJ)/R(NJ-1)
      WRITE(6,FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',
     * (HXP(I30+IJ)/UNU,' ',IJ=1,NCOMP),
     * (HXP(I12+IJ)/UNU,' ',IJ=1,NCOMP)
C
      END IF
C
C******************************************************************
C Start Long Output
C
      ELSE
C
      write(MY_FSTR_HEAD,3010) NCOMP+2
 3010 format('(//,4X,',I2,'(2X,A10))')

      WRITE(6,MY_FSTR_HEAD)' RADIUS   ',' TOTMASS  ',
     *     (TCHAR(1,IK),IK=1,NCOMP)
C--------------Masses----------------
      DO 20001 J=1,MMAX
20001 WRITE(6,MY_FSTR) M(J),R(M(J))/UNR,' ',DEXP(X(IMR,M(J)))/UNM,
     *  (MC(IJ,M(J)),FP(IJ,M(J))/UNM,IJ=1,NCOMP)
C
      WRITE(6,MY_FSTR_HEAD)' RADIUS   ',' TOTMASS  ',
     *     (TCHAR(2,IK),IK=1,NCOMP)
C-------------Densities--------------
      DO 20002 J=1,MMAX
20002 WRITE(6,MY_FSTR) M(J),R(M(J))/UNR,' ',W(M(J))/UNRH,' ',
     *  (DEXP(X(I00+IJ,M(J)))/UNRH,' ',IJ=1,NCOMP)
*       Extrapolated outermost point
      RFAC=R(NJ)/R(NJ-1)
      WRITE(6,MY_FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',RHOTOT/UNRH,' ',
     *  (DEXP(XOUT(IJ))/UNRH,' ',IJ=1,NCOMP)
C
      write(MY_FSTR_HEAD,3011) NCOMP+1
 3011 format('(//,4X,',I2,'(2X,A10))')

      WRITE(6,MY_FSTR_HEAD)' RADIUS   ',
     * (TCHAR(3,IK),IK=1,NCOMP)
C------------Radial Velocity Dispersions------------
      DO 20003 J=1,MMAX
20003 WRITE(6,MY_FSTR) M(J),R(M(J))/UNR,' ',
     * (DEXP((X(I20+IJ,M(J))-X(I00+IJ,M(J)))/2.D0)/UNU,' ',
     *                                       IJ=1,NCOMP)
*       Extrapolated outermost point
      WRITE(6,MY_FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',
     *  (DEXP((PROUT(IJ)-XOUT(IJ))/2.D0)/UNU,' ',IJ=1,NCOMP)
*
      WRITE(6,MY_FSTR_HEAD)' RADIUS   ',
     * (TCHAR(4,IK),IK=1,NCOMP)
C------------Tangential Velocity Dispersions--------
      DO 20004 J=1,MMAX
20004 WRITE(6,MY_FSTR) M(J),R(M(J))/UNR,' ',
     *(DEXP((X(I02+IJ,M(J))-X(I00+IJ,M(J)))/2.D0)/UNU,' ',
     *                                       IJ=1,NCOMP)
*       Extrapolated outermost point
      WRITE(6,MY_FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',
     *  (DEXP((PTOUT(IJ)-XOUT(IJ))/2.D0)/UNU,' ',IJ=1,NCOMP)
*
      WRITE(6,MY_FSTR_HEAD)' RADIUS   ',
     * (TCHAR(5,IK),IK=1,NCOMP)
C-------------Velocities-------------------------
      DO 20005 J=1,MMAX
20005 WRITE(6,MY_FSTR) M(J),R(M(J))/UNR,' ',
     * (X(I10+IJ,M(J))/UNU,' ',IJ=1,NCOMP)
*       Extrapolated outermost point
      WRITE(6,MY_FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',
     * (HXP(I10+IJ)/UNU,' ',IJ=1,NCOMP)
C
      WRITE(6,MY_FSTR_HEAD)' RADIUS   ',
     * (TCHAR(6,IK),IK=1,NCOMP)
C-------------Luminosities-------------------------
      DO 20006 J=1,MMAX
20006 WRITE(6,MY_FSTR) M(J),R(M(J))/UNR,' ',
     *(C12*(3.D0*DEXP(X(I20+IJ,M(J)))*X(I30+IJ,M(J))+
     * 2.D0*DEXP(X(I02+IJ,M(J)))*X(I12+IJ,M(J)))*PI4*R(M(J))**2/UNL,      
     *     ' ',IJ=1,NCOMP)
*       Extrapolated outermost point
      WRITE(6,MY_FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',
     *(C12*(3.D0*DEXP(PROUT(IJ))*HXP(I30+IJ)+
     * 2.D0*DEXP(PTOUT(IJ))*HXP(I12+IJ))*PI4*R(NJ+1)**2/UNL,
     *          ' ',IJ=1,NCOMP)
C
      WRITE(6,MY_FSTR_HEAD)' RADIUS   ',
     * (TCHAR(7,IK),IK=1,NCOMP)
C-------------Fluxes of Radial Energy----------------
      DO 20007 J=1,MMAX
20007 WRITE(6,MY_FSTR) M(J),R(M(J))/UNR,' ',
     * (X(I30+IJ,M(J))/UNU,' ',IJ=1,NCOMP)
*       Extrapolated outermost point
      WRITE(6,FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',
     * (HXP(I30+IJ)/UNU,' ',IJ=1,NCOMP)
C
      WRITE(6,MY_FSTR_HEAD)' RADIUS   ',
     * (TCHAR(8,IK),IK=1,NCOMP)
C-------------Fluxes of Tangential Energy--------------
      DO 20008 J=1,MMAX
20008 WRITE(6,MY_FSTR) M(J),R(M(J))/UNR,' ',
     * (X(I12+IJ,M(J))/UNU,' ',IJ=1,NCOMP)
*       Extrapolated outermost point
      WRITE(6,MY_FSTR) NJ+1,R(NJ)*RFAC/UNR,' ',
     * (HXP(I12+IJ)/UNU,' ',IJ=1,NCOMP)
C
      END IF
C********************************************************************
          DO 100 KKK=1,3*NCOMP
      IXX=I10+KKK-2*NCOMP
      IF(KKK.LE.2*NCOMP)IXX=I20+KKK-NCOMP
      IF(KKK.LE.NCOMP)IXX=I00+KKK
      NFIT=NJ/10
      NF=NJ
C
      DO 65 K=1,10
C
C REGRESSION AND STAND.DEV IS CALCULATED 10 TIMES
C
          IM=0
          N=NFIT
          NFUR(K)=NF
          NFR(K)=NF-N
C
C NFIT GOVERNS THE NUMBER OF POINTS FOR REGRESSION ANALYSIS
C
      DO 70 I=1,NFIT
      IR=NF-I
      RFL=REAL(R(IR))
      IF((RFL.LT.1.E-30).OR.(RFL.GT.1.E30))THEN
      XYZ(I,1)=0.E0
      ELSE
      XYZ(I,1)=ALOG(RFL)
      END IF
        IF(KKK.LE.2*NCOMP)THEN
          XYZ(I,2)=REAL(X(IXX,IR))
        ELSE
          RFL=REAL(DABS(X(IXX,IR)))
          IF((RFL.LT.1.E-30).OR.(RFL.GT.1.E30))THEN
          XYZ(I,2)=0.E0
          ELSE
          XYZ(I,2)=ALOG(RFL)
          END IF
        END IF
 70   CONTINUE
C
          CALL SLOPE(XYZ,NFIT,STAT)
C
          RES(K,1)=STAT(1)
          RES(K,2)=STAT(2)
C
          NF=NF-N
          IF(NF.EQ.NJ/10)NF=2+NFIT
C
 65       CONTINUE
C
       IF(KKK.EQ.1)WRITE(6,485)(NFR(IK),NFUR(IK),IK=1,10)
       AB='I10: '
       IF(KKK.LE.2*NCOMP)AB='I20: '
       IF(KKK.LE.NCOMP)AB='I00: '
       KMOD=MOD(KKK-1,NCOMP)+1
C
       WRITE(6,487) AB,KMOD,(RES(I,1),I=1,10)
       WRITE(6,488)(RES(I,2),I=1,10)
 100      CONTINUE
 485   FORMAT(14X,' FROM/TO I=:',10(I4,'/',I4,3X))
 487   FORMAT(1X,A5,' COMP.',I2,' REGRESSION:',1P,10(E9.2,1X))
 488   FORMAT(14X,' STAND.DEV.:',1P,10(E9.2,1X))
C
      RETURN
C
 998  CONTINUE
      PRINT*,' Error occurred during WRITE on UNIT 1 or 11'
      PRINT*,' IOSTAT2=',IWRERR,' ILENT,IBREC=',ILENT,IBREC
      PRINT*,' FORTRAN STOP '
      STOP
*
 999  CONTINUE
      PRINT*,' Error occurred during WRITE on UNIT 2 or 4'
      PRINT*,' IOSTAT2=',IWRERR,' IOSTAT4=',IWRER4
      PRINT*,' FORTRAN STOP '
      STOP
C
      END
