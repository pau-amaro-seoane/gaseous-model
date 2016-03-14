      SUBROUTINE READMO
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
         INCLUDE 'compar.f'
         INCLUDE 'equiv.f'
      PARAMETER(NTAB=32)
      INTEGER IDUM2,K,IV,IY
      REAL*8 DMSAVE(NBINO*50)
      COMMON/RANDO/IDUM,IY,IDUM2,IV(NTAB)
C
        IF(LS(3))OPEN(1,FILE='fort.1',FORM='UNFORMATTED',
     *   STATUS='UNKNOWN',IOSTAT=IERR)
          IF(IERR.NE.0)THEN
          PRINT*,' OPEN ERROR ON UNIT 1 IERR=',IERR
          STOP
          END IF
          OPEN(2,FILE='fort.2',FORM='UNFORMATTED',
     *   STATUS='UNKNOWN',IOSTAT=IERR)
          IF(IERR.NE.0)THEN
          PRINT*,' OPEN ERROR ON UNIT 2 IERR=',IERR
          STOP
          END IF
          OPEN(4,FILE='fort.4',FORM='UNFORMATTED',
     *   STATUS='UNKNOWN',IOSTAT=IERR)
          IF(IERR.NE.0)PRINT*,' OPEN 4 IERR=',IERR
C
 6662     CONTINUE
          READ(2,ERR=999,END=9998,IOSTAT=IREERR)(EI(K),K=1,NW/2)
          READ(2,ERR=999,END=9998,IOSTAT=IREERR)(EI(K),K=NW/2+1,NW)
C
          JMODA=IEI(1)
*
          IF(JMODA.EQ.MODA)THEN
          PRINT*,' Requested Model MODA=',MODA,' read '
          ELSE
          GOTO 6662
          END IF
*
          NMOD=IEI(2)
          IBINT=IEI(8)
          IBIN=IEI(9)
          IREC=IEI(10)
          IBREC=IEI(11)
          PRINT*,' IBREC=',IBREC
          ITER=IEI(12)
          NJ=IEI(13)
          NCOMP=IEI(14)
          IMCRIT=IEI(15)
C
          IF(LS(3))THEN
*
 6663     CONTINUE
          READ(91,*,ERR=9995,END=9995)JMODA,IDUM,IY,IDUM2,
     *           (IV(KK),KK=1,NTAB)
*
          IF(JMODA.EQ.MODA)THEN
          PRINT*,' Random number sequence for MODA=',MODA,' read'
          ELSE
          GOTO 6663
          END IF

*         For CRAY set random number sequence
          IFAIL = RANSET(IDUM2)
          PRINT*,' Random Numbers Initialized: IDUM,IY,IDUM2=',
     *     JMODA,IDUM,IY,IDUM2,' IV=',(IV(KK),KK=1,NTAB)
*
          IF(JMODA.NE.MODA)THEN
          PRINT*,' Error in reading random number sequence '
          STOP
          END IF
*
          END IF
*
          IF(LS(3).AND.IBIN.GT.0)THEN
*
 6664 CONTINUE
      READ(1,ERR=9994,END=9994,IOSTAT=IREERR)JMODA,JBREC
      READ(1,ERR=9994,END=9994,IOSTAT=IREERR)(BEI(K),K=1,NBEI)
      READ(1,ERR=9993,END=9993,IOSTAT=IREERR)(IBEI(K),K=1,NIBEI)
      READ(1,ERR=9992,END=9992,IOSTAT=IREERR)ILENT
      READ(1,ERR=9991,END=9991,IOSTAT=IREERR)(DMSAVE(K),K=1,ILENT)
*
      PRINT*,' Binary record ',JBREC,' read from fort.1 '
      PRINT*,' Read Binary number IBIN=',IBIN
      PRINT*,' belonging to MODA=',JMODA,' ILENT read =',ILENT

      IF(JMODA.EQ.MODA)THEN
      ELSE
      GOTO 6664
      END IF
*
      ILENT = 0
      DO 1000 IB=1,IBIN
      DO 1005 I=1,ILEN(IB)
      IBIX = IBSTA(IB)+I
      DMBIN(IB,IBIX)=DMSAVE(ILENT+I)
*     PRINT*,' read from DMSAVE(',ILENT+I,') DMBIN(',IB,',',IBIX,')',
*    *  ' value =',DMBIN(IB,IBIX),' bef=',DMBIN(IB,IBIX-1),
*    *  ' aft=',DMBIN(IB,IBIX+1)
*     CALL FLUSH(6)
 1005 CONTINUE
 1000 ILENT = ILENT + ILEN(IB)
*
          PRINT*,' Binaries read IBIN=',IBIN,' IBREC=',IBREC
*
      TIME=AEI(2)
      IC50=0
 2000 CONTINUE
      READ(50,*,ERR=2999,END=2998)XTIME
      IF(XTIME.GT.TIME)THEN
      GOTO 2998
      ELSE
      IC50=IC50+1
      GOTO 2000
      END IF
*
 2999 CONTINUE
      PRINT*,' WARNING'
      PRINT*,' No Escapers yet to be read from Unit 50' 
      PRINT*,' or read error from Unit 50'
*
 2998 CONTINUE
      PRINT*,IC50,' records read from Unit 50'
      BACKSPACE 50
      CALL FLUSH(6)
*
      IC60=0
 3000 CONTINUE
      READ(60,600,ERR=3999,END=3998)IDDM,IDDM2,IDDM3,IDDM4,XTIME
 600  FORMAT(2I5,4X,I5,1P,I8,D12.4)
      IF(XTIME.GT.TIME)THEN
      GOTO 3998
      ELSE
      IC60=IC60+1
      GOTO 3000
      END IF
*
 3999 CONTINUE
      PRINT*,' WARNING'
      PRINT*,' No Binaries yet to be read from Unit 60'
      PRINT*,' or read error from Unit 60'
*
 3998 CONTINUE
      PRINT*,IC60,' records read from Unit 60'
      BACKSPACE 60
      CALL FLUSH(6)
*
      END IF

          IF(ITIME.GE.0)THEN
          BREM=AEI(3)
          IF(ITIME.EQ.0)PRINT*,
     *       ' Timestep unchanged - read from data file'
          IF(ITIME.GT.0)THEN
          TIME=AEI(2)
          TTRX=AEI(4)
          PRINT*,' Timestep taken from input control data'
          XN3TOT=AEI(9)
          PRINT*,' Restart with XN3TOT=',XN3TOT
          PRINT*,' CLOCK SET BY FILE'
          END IF
          END IF
          IF(ITIME.LT.0)THEN
          TIME=0.D0
          TTRX=0.D0
          PRINT*,
     *      ' Timestep after restart determined in CHOOSE'
          END IF
C----------Note that AEI(4)-AEI(8) are used in ZOMBI
          XMHOLE=AEI(10)
          DMTOT=AEI(11)
          XN3B=AEI(12)
          CS(17)=AEI(13)
          DPTOT=AEI(14)
          XMERR=AEI(15)
          XEERR=AEI(16)
          PHTID=AEI(17)
          RTID0=AEI(18)
          XMTOT0=AEI(19)
          XMBIN0=AEI(20)
          XMERR2=AEI(21)
          XEERR2=AEI(22)
          XMERR3=AEI(23)
          XEERR3=AEI(24)
          XMERR4=AEI(25)
          XEERR4=AEI(26)
C
          WRITE(6,102)
 102      FORMAT(1X,'-------------START MESSAGE ----------------')
          WRITE(6,103) MODA,NMOD,TIME,TTRX,
     *      TIME*CUNT/YEAR,1.D0/BREM,1.D0/BREM*CUNT/YEAR
 103      FORMAT(' START: MODA=',I4,' NMOD=',I8,' TNB=',1P,D9.2,
     * ' TTRX=',D9.2,' T=',D9.2,' YRS.',/,
     * ' DTNB=',D9.2,' DT=',D9.2,' YRS.')
*      Adjust stellar masses if stellar evolution is present
       IF(LS(21))CALL STEVOL
*      Check whether switch on central black hole is due
*
       IF(LEHOLE.AND.TTRX.GT.TBINI)THEN
       LS(19)=.TRUE.
       LEHOLE=.FALSE.
       IF(LLOSS)LS(14)=.TRUE.
       LLOSS=.FALSE.
       IF(LELOSS)LS(18)=.TRUE.
       LELOSS=.FALSE.
       PRINT*,'*SUCCESSFUL RESTART AFTER BLACK HOLE SWITCH ON TIME'
       END IF
       IF(LEHOLE.AND.TTRX.LE.TBINI)THEN
       PRINT*,'*SUCCESSFUL RESTART BEFORE BLACK HOLE SWITCH ON TIME'
       END IF
*
          IF(LS(19))WRITE(6,104)XMHOLE*CUNM/SUNM,DMTOT*CUNM/SUNM
 104      FORMAT(' START WITH BLACK HOLE MASS=',1P,E9.2,
     *  ' DM=',E9.2,' SOL.MASSES ')
C
          WRITE(6,107)
 107      FORMAT(//,'****END START INFORMATIONS****',//)
C
C NECESSARY PARAMETERS (CS,IS, ETC.) STORED ON FILE
C HAVE TO BE RECOVERED HERE
C
C
          RETURN
C
 999      CONTINUE
          PRINT*,' Error occurred during READ from Unit 2 '
          PRINT*,' IOSTAT=',IREERR
          PRINT*,' MODA stands at MODA=',IMODA
          PRINT*,' FORTRAN STOP '
          STOP
C
 9998     PRINT*,' End of File in READ Unit 2'
          STOP
 9995     PRINT*,' Error 9995 during READ Unit 91'
          STOP
 9994     PRINT*,' Error 9994 during READ Unit 1'
          STOP
 9993     PRINT*,' Error 9993 during READ Unit 1'
          STOP
 9992     PRINT*,' Error 9992 during READ Unit 1'
          STOP
 9991     PRINT*,' ILENT=',ILENT,' Error 9991 during READ Unit 1'
          STOP
 6660     PRINT*,' Error 6660 during OPEN Unit 60'
          STOP
          END
