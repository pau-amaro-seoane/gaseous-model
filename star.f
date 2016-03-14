      PROGRAM STAR
C     MAIN PROGRAM "SPEDI" ("star" for hist. reasons)
C     
C     Numerical Simulation of Multi-Component Spherical Stellar Systems
C     by solution of moment equations
C     with an implicit Newton-Henyey-Raphson scheme
C     
C     R.Spurzem
C     Institute of Theoretical Physics and Astrophysics, 
C     Univ. Kiel
C     
C     Explanation of Parameters given in COMPAR
C-------------------------------------------------------------------
C     Various Common blocks are included see there
C     
C---------------------------------------------------------------------
C     Hints for execution on VAX/VMS-systems:
C     
C     INCLUDE-Statements need logical name  for the
C     directory containing the source codes.
C     
C     This directory should contain the following COM-Files:
C     1) FORALL.COM       compile all on Directory SC: 
C     2) LINKALL.COM      create Library and link;
C     module (.EXE) and library (.OLB) are
C     created on the logical directory
C     named SC:, its name must be user-predefined
C     3) STAR.COM     command file for running the program, 
C     expects parameters P1=Data File Number
C     P2=Module
C     P3=Additional Data File key
C     (P1: necessary; P2 default=STAR; P3 def. =empty)
C     
C     STAR called with P1 P2 P3 uses the input data files
C     "DKUG"+P1+P3+".DAT"
C     "MESS"+P1+P3+".DAT"
C     and the output files
C     "CWRIT"+P1+P3+".DAT"
C     "INFO"+P1+P3+".DAT"      
C     The input files are expected on directory SPLO 
C     The output files are created on directory RUNDIRn ,
C     where n is the parameter P1 above
C     Execution of the run happens also on directory RUNDIRn.
C     
C     Summary: The logical directory names SPLO, SC, RUNDIRn have
C     to be defined by user.
C     
C     4) STARS.COM   Same as STAR.COM, but for input  and
C     output Files like e.g. "DKUG"+P1+"S"+P3
C     or   "CWRIT"+P1"+"S"+P3
C     Here "S" means short; one has the opportunity
C     to execute on one RUNDIR - Directory
C     two runs which do not disturb the data
C     files for each other.
C     
C     5) STARn.COM for n=1,2,3... contains the call
C     of STAR with appropriate parameters for
C     use as a batch run (SUBMIT/... STARn.COM)
C     
C     6) STARnS.COM  Same as STARn.COM but calling STARS
C**********************************************************************
C     NOTES FOR INPUT/OUTPUT
C--------------------------------------------------------------
C     Input: DKUGn.DAT/DKUGnS.DAT is read as logical unit 15
C     MESSn.DAT/MESSnS.DAT is read as logical unit 14
C     
C     DKUG contains input data for the run;
C     MESS contains a heading message for the output
C     (Both are read and parttially printed in subr. READPM) 
C--------------------------------------------------------------
C     Output: CWRITn.DAT is written/read as logical unit 2
C     INFOn.DAT is written/read as logical unit 3
C     
C     CWRIT contains the full information of selected
C     models in form of 2 records of length NW REAL*8 (set below)
C     written sequentially, but unformatted 
C     INFO contains time-dependent quantities stored in 
C     vector VIT (see subroutine INFO) in form
C     of sequential records of length 100 REAL*8 
C     
C     CWRIT is written in subr. PRINT
C     INFO is written in subr. INFO
C     
C     In case of a RESTART run (start at some model already
C     written by another run)
C     the records of CWRIT and INFO are first read sequentially
C     to put the file pointer on a position after
C     which new information can be written;
C     this is done by subrs. READMO and INFO, respectively.
C     
C------------------------------------------------------------
C     
C     General Warning: Many working variables are stored in
C     common blocks which are always included
C     by the INCLUDE statement.
C     If you introduce new variables: make sure that a variable
C     of the same name is not already defined
C     in the COMMON-blocks!
C*******************************************************************
C     
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-K,M-N),LOGICAL(L)
*     
      INCLUDE 'compar.f'
*     
      COMMON/TIMER/TCOMP0
      COMMON/TEST/DEBIN
      CHARACTER*9 ZS
      REAL*4  TARRAY(2),ETIME
      COMMON/INFOHELP/ZS(IVDIM)
      DIMENSION NF(50),VST(50,IVDIM)
      DIMENSION DDPHI(NCOMPO)
C     Initialize Timer
      TCOMP0 = ETIME(TARRAY)
*     
      C13=1.D0/3.D0
      C23=2.D0*C13
      C53=5.D0*C13
      C56=C53/2.D0
      C12=5.D-1
      C32=3.D0*C12
      C52=2.5D0
      PI4=4.D0*PI
      RHOC=SUNM/PC**3
      CLFAC=CLICHT/2.997925D5
      CLFAC2=CLFAC*CLFAC
      OPEN(90,FILE='fort.90')
      OPEN(14,FILE='fort.14')
      OPEN(15,FILE='fort.15')
C---------------------------------------------------------------------
C     Initialize Counters: IMOD: counts models of this run
C     NMOD: counts models including previous runs
C     (is read from file by restart)
C     IRIT: governs interval of INFO-writing
C---------------------------------------------------------------------
      IMOD=0
      NMOD=0
      INFNUM=0
      IBREC=0
      IFRFIN=0
      IRIT=0
      TIME=0.D0
*     Initialize central trx counter
      TTRX=0.D0
      ITER=0
      ITERS=0
      IERR=0
C     ITER=0 Leaves Timestep unaffected by CHOOSE***************************
C     - - - - - - Read Initial Parameters - - - - - - - - - - - - - - - -
*     Locate Variables
      IMR=1
      I00=1
*     
      CALL READPM
*     
      I20=NCOMP+1
      I02=2*NCOMP+1
      I10=3*NCOMP+1
      I30=4*NCOMP+1
      I12=5*NCOMP+1
*     
*     Determine equation number etc
      NG = NCOMP*NEQC + NADD*NEQA + NEQX
      NPOS = NCOMP*NPOSC + NADD*NPOSA + NPOSX
*     
      CALL PRINTO
C***********************************************************************
C     Construct grid system and initial model or READ start model
      IF(LS(10))THEN
         CALL READMO
         MODA=MODA+1
      ELSE
         CALL RGRID
         CALL INIT
      END IF
C**************Call Black Hole Initialisation if desired***************
      IF(LS(19))THEN
         CALL AHOLE
C     Initialise loss-cone
         IF(LS(14))CALL LHOLE
         IF(LS(18))CALL LHOLEN
      END IF
C**************Print Every Model's Information*************************
C**************Compute PHI and primordial binary data if necessary*****
      CALL INFORM
      CALL PRINTS
      CALL INFO
C**************Print Out Large Info for First Model********************
      IIP=-1
      IF(MAXMOD.EQ.0)IIP=1
      CALL PRINT(IIP)
      IF(MAXMOD.EQ.0)STOP
C***********************************************************************
C     WORKING LOOP
C***********************************************************************
 60   CONTINUE
C     ********Choose Timestep, Call Henyey-Iteration************************
      CALL CHOOSE
C*********Stochastic Binary Energy Generation**************************
*     Test case only
*     IF(LREM.OR.IMOD.EQ.49)THEN
*     III=III+1
*     LREM=.FALSE.
*     IF(III.LT.10)LREM=.TRUE.
*     DO 7002 I=1,NJ
*     IF(III.EQ.1)THEN
*     DDPHI(I)=XBIN2B*PHI(I)
*     END IF
*     7002    PHIBIN(I)=PHIBIN(I)+DDPHI(I)
*     
*     XMRBIN(1)=0.D0
*     DO 7001 I=2,NJ
*     IM=I-1
*     DER=R(I)-R(IM)
*     XMRBIN(I)=R(I)**2/GRAV*(PHIBIN(I)-PHIBIN(I-1))/DER
*     7001    CONTINUE
*     
*     PRINT*,' Artificial Binary mass ',XMRBIN(NJ),' created'
*     
*     END IF
*     
*     
      IF(LS(3))CALL BINSTO
*     
C*********Late Switch on Binary Energy - reduce timestep****************
      IF(.NOT.LEBIN.AND.(LS(1).OR.LS(2)).AND.TIME.GT.TBINI)THEN
         LEBIN=.TRUE.
         BREM=1.D1*BREM
      END IF
C***********************************************************************
      CALL HENYEY
*     DESUM=DESUM+DEBIN
*     IF(MOD(NMOD,NCOL).EQ.0)
*     * PRINT*,'   star-DESUM=',DESUM,' DEBIN=',DEBIN
C****************************************************************
      IF(ITER.LE.ITMAX)ITERS=ITERS+ITER
      IMOD=IMOD+1
      NMOD=NMOD+1
*     
*     Check whether switch on central black hole is due
*     
      LGROW=.FALSE.
      IF(LEHOLE.AND.TTRX.GT.TBINI)THEN
         LGROW=.TRUE.
         LS(19)=.TRUE.
         LEHOLE=.FALSE.
         IF(LLOSS)LS(14)=.TRUE.
         LLOSS=.FALSE.
         IF(LELOSS)LS(18)=.TRUE.
         LELOSS=.FALSE.
         BREM=BREM*XBIN2B
         XMHFIN=XMHOLE
         XMHOLE=XMIND(1)
      END IF
*     Adiabatic growth of lately inserted central black hole
      IF(LGROW.AND.XMHOLE.LT.XMHFIN)THEN
         XMHOLE=XMHOLE*(1.D0+CORR)
         PRINT*,' Slow Growth of MHOLE=',XMHOLE
         BREM=BREM*FACT
         IF(XMHOLE.GE.XMHFIN)THEN
            XMHOLE=XMHFIN
            LGROW=.FALSE.
         END IF
      END IF  
C**************Call Black Hole Update if desired*****************
      IF(LS(19))THEN
         CALL AHOLE
C     loss cone
         IF(LS(14))CALL LHOLE
         IF(LS(18))CALL LHOLEN
      END IF
C****************************************************************
      CALL INFORM
      CALL PRINTS
C**** Control Information Printing********************************
      IRIT=IRIT+1
      IF(IRIT.EQ.NRIT.OR.NMOD.LT.20)THEN
         CALL INFO
         IRIT=0
      END IF
C------------Store VIT-INFO for final Printout--------------------
      INFINT=MAXMOD/50
      IF(MAXMOD.LT.50)INFINT=1
      IF((MOD(IMOD,INFINT).EQ.0).AND.(INFNUM.LT.50))THEN
         INFNUM=INFNUM+1
         DO 2500 J=1,IVDIM
            NF(INFNUM)=NMOD
            VST(INFNUM,J)=VIT(J)
 2500    CONTINUE
      END IF
C------------------------------------------------------------------
      IF(ITER.GE.ITMAX)GOTO 900
      IFALSE=0
      IF(IMOD.GT.MAXMOD)GOTO 4811
      IF((IMOD.LT.MODM).AND.(IMOD.LT.MODF))GOTO 60
      IPF=-1
      IF(IMOD.EQ.MODF)IPF=2
      IF((IMOD.EQ.MODF).AND.(IMOD.EQ.MODM))IPF=1
      CALL PRINT(IPF)
      IF(IPF.GT.0)MODA=MODA+1
      IF(IPF.LT.2)MODM=MODM+MODFAC
      IF  (IMOD.LT.MODF)GOTO 60
      MODF=MODF+MODFI
      GOTO 60
C***********************************************************************
 4811 WRITE(6,2000)
      CALL PRINT(-1)
      PRINT*,' This run has totally done ',ITERS,' iterations'
      PRINT*,' and totally written ',IREC,' INFO-records'
      GOTO 5010
C     
 900  WRITE(6,2002)
      WRITE(6,2003) TIME/YEAR,1.0D0/BREM/YEAR,ITER
      WRITE(6,2004) TTRX,DTTRX
      IF(IFALSE.GT.2)GOTO 960
      IF(LS(7))THEN
         WRITE(6,75)(JG(I),I=1,NG)
 75      FORMAT(1X,8(I3,1X))
         CALL PRINT(-1)
      END IF
      IERR=IERR+1
      IF(IERR.GT.ITMAX)GOTO 960
      TIME=TIME-1.D0/BREM
      TTRX=TTRX-DTTRX
      NMOD=NMOD-1
      BREM=BREM/FAA
      DTTRX=DTTRX*FAA
      DO 951 K=1,NG
         DO 950 J=1,NJ
            X(K,J)=VX(K,J)
            VX(K,J)=VVX(K,J)
 950     CONTINUE
 951  CONTINUE
      IF(IMOD.GT.1)IFALSE=IFALSE+1
      IF(IFALSE.LT.3)THEN
         ITER=0
         GOTO 60
      END IF
 960  CONTINUE
      CALL PRINT(-1)
      WRITE(6,2001)
 5010 CONTINUE
*     
      ICOUNT = 0
      DO 5050 IV = 1,IVDIM/13+1
*     
         ISTART = ICOUNT + 1
         ICOUNT = ICOUNT + 13
         IF(ICOUNT.GT.IVDIM)ICOUNT=IVDIM
*     
         INUM = ICOUNT - ISTART + 1
*     
         WRITE(6,5000)(ZS(K),K=ISTART,ICOUNT)
         WRITE(6,5001)ISTART,ICOUNT
         INUM = ICOUNT - ISTART + 1
         IF(INUM.EQ.13)THEN
            WRITE(6,5011)(NF(J),(VST(J,I),I=ISTART,ICOUNT),J=1,50)
         ELSE IF(INUM.LT.13)THEN
            WRITE(6,5011)(NF(J),(VST(J,I),I=ISTART,ICOUNT),
     *           (0.D0,I=1,13-INUM),J=1,50)
         END IF
*     
 5050 CONTINUE
*     
      STOP
C     
 5000 FORMAT(//,' HELP=',13A9)
 5001 FORMAT(' NMOD=; VIT(',I5,'-',I5,')= ',/)
 5011 FORMAT(50(1X,I5,1P,13D9.2,/))
 2000 FORMAT(1X,' Normal Stop in Main Program AGMULT')
 2001 FORMAT(1X,' Oops! Convergence not possible in Subroutine HENYEY')
 2002 FORMAT(/,1X,' TIMESTEP TOO LARGE, BOY!  DECREASED BY  *CS14')
 2003 FORMAT(1X,' TIME=' ,1P,E9.2,' DTIME=',1P,E9.2,' ITER=',I4)
 2004 FORMAT(1X,' TTRX=' ,1P,E9.2,' DTTRX=',1P,E9.2)
C     DEBUG SUBCHK
      END

      BLOCK DATA
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      INCLUDE 'compar.f'
C     
C     constants in cgs units *************************
C******mean molecular weight for fully H-ionized matter/solar composition
      DATA XMUE/6.D-1/
      DATA HMASS/1.6726D-24/
      DATA PI/3.141592654D0/
      DATA SUNR/6.96D10/
      DATA PI43/4.188790203D0/
C******xmk = hmass/Boltzmann constant
      DATA XMK/1.212162D-8/
      DATA HWIRK/6.626196D-27/
      DATA CLICHT/2.997925D10/
C******mean molecular weight for neutral matter/solar composition
      DATA XMUEI/1.505D0/
      DATA GRAV/6.6732D-8/
      DATA PC/3.0857D18/
      DATA SUNM/1.989D33/
      DATA SUNL/3.82D33/
C********year in seconds
      DATA YEAR/3.15576D7/
C*******Stefan-Boltzmann constant, beware of different definitions,
C*******here it is always:
C     Planckfunction B= SIGSBO *CLICHT *T**4/4PI
      DATA SIGSBO/7.565232619D-15/
C     
      END

