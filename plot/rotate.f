       IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
C
       IROT=1 
       PRINT*,' Rotation Axis IROT=',IROT,' (density) '
C
C       PRINT*,' Enter Rotation angle (degree):'
C       READ*,ANG
C
       PI=4.D0*DATAN(1.D0)
C
       DO 100 I=1,10
       REWIND 50
C
       ANG=DBLE(18*I)
       PRINT*,' I=',' Angle=',ANG
       XANG=(PI/180.D0)*ANG
C
       CANG=DCOS(XANG)
       SANG=DSIN(XANG)
       PRINT*,' ANG,XANG,COS,SIN=',ANG,XANG,CANG,SANG
C
       IREC=0
       SUM=0.D0
C
       PRINT*,' Data to be read from File Unit 50'
C
 1999  IREC=IREC+1
C
       READ(50,*,ERR=700,END=700)TIME,TTRX,RHO,SIG,XI
C
C       IF(MOD(IREC,500).EQ.0)PRINT*,IREC,' Records read '
C
       RHOP=RHO
       SIGP=CANG*SIG+SANG*XI
       XIP=-SANG*SIG+CANG*XI
       IF(IREC.EQ.1)VXIP=XIP
       SUM=SUM+(VXIP-XIP)
C
       VXIP=XIP
C
C       WRITE(51,5005)TIME,TTRX,RHOP,SIGP,XIP
C 5005 FORMAT(1P,5(D14.7,1X))
C
       GOTO 1999
C
 700   PRINT*,IREC-1,' Records transformed ...'
C
       PRINT*,' Angle=',ANG,' SUM=',SUM
C
 100   CONTINUE
C
       STOP
       END
