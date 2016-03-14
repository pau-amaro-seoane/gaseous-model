C Compar for 3-Point Henyey Program C*******************************************************************|
C NCOMP= number of Components                                       |
C NEQC= Number of Equations per Component                           |
C NPOSC= Number of Equations for positive Variables per Component   |
C*******************************************************************|
C NADD= Number of additional Components                             |
C NEQA= Number of Equations for additional component                |
C NPOSA= Number of additional Equations for positive Variable       |
C                                                                   |
C       Additional means here any component physically different    |
C       from a stellar component                                    |
C    NOTE: In this version NO additional component is implemented   |
C*******************************************************************|
C NEQX= Number of Component independent equations,                  |
C       here e.g. equation for total mass                           |
C NPOSX= Number of positive variables thereof                       !
C*******************************************************************|
C NG=Total Number of Eqs.=NCOMP*NEQ+NADD*NEQA+NEQX                  |
C NPOS= Total Number of Equations for positive Variables            |
C NPOS= NCOMP*NPOSC+NADD*NPOSA+NPOSX
C NJ= Number of Grid Points                                         |
C*******************************************************************|
C This is a parameter choice for first dynamical test calculations  |
C with 1 stellar dynamical component                                |
C*******************************************************************|
C                                                                   |
C  IMR, I00, I20 ... etc. are Index-Variables to refer easily       |
C     the physical quantities in the vectors X, VX; for example:    |
C                                                                   |
C       X(I30+J,I) refers to the 30-Moment of component J at grid I |
C      correspondingly moments 12, 20, 02, 00                       |
C      a little inconsistency: 10 denotes the velocity itself       |
C                              (10-moment would be rho*velocity)    |
C       X(IMR,I)  refers to the total mass at grid i                |
C*******************************************************************|
      INCLUDE 'params.f'
C       Save statement for static local variables
      SAVE
C*******************************************************************
C    Variables of Common-Block Work:
C    A,G,DG,JG,HY,C,D,E used by subr. HENYEY (see there)
C
C    HX,HXM,HXP contain lambda*xnew + (1-lambda)*xold for
C               semi-implicit scheme
C    X,VX,VVX variables at present/previous/pre-previous timestep
C    (VVX is useful for restauration after unsuccessful iteration)
C    Y loss cone filling factor
C    R radial independent coordinate
C    AEI, IEI working variables to be saved on file CWRIT
C    FP,FM used in Equations (e.g. I10eq) to store sign of velocities
C    F used in equations to switch off single equations
C    W working space; used in HENYEY to set minimal values for
C                     determination of relative corrections
C    VIT stores time-dependent information to be saved on file INFO
C********************************************************************
      COMMON/WORK/A(NGO*(2*NGO+1)),G(NGO),DG(NGO),JG(NGO),
     *   HY(NGO*(NGMAX+1)*NJO),
     *   C(NGMAX,NGO),D(NGMAX,NGO),E(NGMAX,NGO),
     *   HX(NGO),HXM(NGO),HXP(NGO),
     *   VVX(NGMAX,NJO),VX(NGMAX,NJO),X(NGMAX,NJO),Y(NCOMPO,NJO),
     *   Y2(NCOMPO,NJO),R(NJO),AEI(100+20*NCOMPO),IEI(20),
     *   TRX(NCOMPO,NCOMPO),TRXP(NCOMPO,NCOMPO),
     *   CR2B(NCOMPO,NCOMPO),ALP2B(NCOMPO,NCOMPO),
     *   SUMI(NCOMPO),SUMIP(NCOMPO),
     *   FP(NCOMPO,NJO),FM(NCOMPO,NJO),FX(NGO),W(NJO+NGO),PHI(NJO),
     *   VIT(IVDIM),TVEC(20)
C*******************************************************************
C    Common Block LOSS
C    containing all quantities for loss cone accretion of black hole
       COMMON/LOSS/ULC,VLC,VLCTH,VLCPH,
     *             DOMEGA,D0DOM,D1DOM,D2DOM,DTDOM,
     *             DOMEG1,D0DOM1,D1DOM1,D2DOM1,DTDOM1,
     *             EPSR,D0DER,D1DER,D2DER,DTDER,
     *             EPST,D0DET,D1DET,D2DET,DTDET,
     *             EPSR1,D0DER1,D1DER1,D2DER1,DTDER1,
     *             EPST1,D0DET1,D1DET1,D2DET1,DTDET1
C*******************************************************************
C    Common Block ITRE
C    IMOD count models of present run (i.e. calls of HENYEY)
C    NMOD count models including previous runs in case of restart
C    (Note: for unsucessful iteration IMOD is increased, NMOD not)
C    ITER returned by HENYEY with number of iterations or
C         with error condition (ITER>100)
C    MODA count models written on file CWRIT
C    ITMIN/ITMAX minimal/maximal number of iterations
C    BREM = 1/timestep
C    TIME = clock
C    (Note: BREM and TIME can be set by input parameters even
C     for restarts by use of the parameter ITIME)
C    EPS/FVEC = convergence limit for relative corrections
C                of positive /indefinite quantities, respectively.
C    FHEN = factor reducing the corrections applied in HENYEY
C           (see FACDX in HENYEY; necessary FHEN <= 1
C    FACT = factor by which timestep in CHOOSE is adjusted
C           (note: for extreme values of CS(17)=DMAX FACT is overriden)
C    CORR = target value for relative corrections of positive variables
C    FAA =  for unsuccessful iteration retry with DTIME=FAA*DTIME
C    IPR1, IPR2 ?
C    NW  :  length of records in CWRIT is NW/2
C    JDUM = dummy for IBM systems to keep everything in multiples
C           of double-word-lengths
C    STEU1-STEU3 parameters used in RGRID to construct radial grid
C    MODM, MODFAC : the first model is fully written on the output
C                   file after MODM models, thereafter in intervals of
C                   MODFAC models
C    MODF, MODFI :  as MODM, MODFAC for output on file CWRIT
C    MAXMOD : if IMOD=MAXMOD the run stops
C    IRIT,NRIT : after each NRIT models the VIT-vector is written
C                on file INFO; IRIT is the local counter in main prog.
C    NCOL: each NCOL models the medium large print out is done
C          by subroutine INFORM
C    NWRITE: for full model output: grid points 1-10, NJ, and therein
C                                   each NWRITE-th point printed
C    IREC: counter for the records on file INFO, saved on CWRIT
C********************************************************************
      COMMON/ITRE/
     * IMOD,NMOD,NJ,NCOMP,NG,NPOS,
     * J3DUM,ITER,MODA,ITMIN,ITMAX,BREM,TIME,EPS,FHEN,FACT,
     * CORR,FAA,FVEC,IPR1,IPR2,JDUM,J2DUM,STEU1,STEU2,STEU3,
     * MODM,MODF,MODFAC,MODFI,MAXMOD,IRIT,NRIT,NCOL,NWRITE,IREC
      COMMON/PHYS/PI,PI43,HWIRK,CLICHT,XMUEI,GRAV,PC,SUNM,SUNL,YEAR
     *    ,XMK,SUNR,XMUE,SIGSBO,HMASS,CONGAS
C*******************************************************************
C*****With Units Parameters some output can be switched to
C*****other units (see printo.f)
C
      COMMON/UNITS/UNR,UNM,UNRH,UNU,UNT,UNL
C*****With Cunits Parameters calculation units are determined -
C*****read in from Parameter data in CGS units
      COMMON/CUNITS/CUNR,CUNM,CUNRH,CUNU,CUNT,CUNL
C*******************************************************************
C General Working space REAL*8 (CS), Logical (LS), Integer (IS)
C Internally used are:
C     CS(13), CS(19), CS(20), CS(23) for move grid (LS(26), not supported!)
C                            (see e.g. GINIT)
C     CS(17), CS(18) max rel. corrections (HENYEY)
C     CS(21)  Lambda in stellardynamical heat flux equation
C     CS(22)  Factor for escape rate
C     CS(60) parameter to compute trelax easy (see READPM)
C     CS(n) for n>= 200 (see EQUIVALENCES below)
C
C     IS(n) for n<=13   (see EQUIVALENCES below)
C     LS(n) for n<=40   (input variables see READPM)
C**********************************************************************
      COMMON/PARAM/CS(300+55*NCOMPO),LS(100),IS(100)
C**********************************************************************
C    Vector IC governs diagnostic printout (LS(25)=TRUE)
C     each tripel (IC(I+K),K=1,3) for I=0,3,6,9,12,15
C     is interpreted as IMOD, ITER, grid number for GIDINI
C     to print out equation values and much more
C     IC is an input vector read by READPM
C**********************************************************************
      COMMON/CNTROL/IC(20)
C**********************************************************************
C    LEQ used to switch off equations, see PRINTO
C**********************************************************************
      COMMON/LEQC/LEQ(NGO)
C-----------------------------------------------------------------------
C    Note: all previous quantities used after HENYEY to simplify
C          expressions in equations, defined in GINIT
C-------------------------------------------------------------------
      DIMENSION AFAC(NCOMPO),ATFAC(NCOMPO),XBINY(NCOMPO),XMASSY(NCOMPO)
      DIMENSION XFBINY(NCOMPO),XSBINY(NCOMPO)
      DIMENSION UAV(NCOMPO),UAVP(NCOMPO),VRAV(NCOMPO),VTAV(NCOMPO)
      DIMENSION PAVR(NCOMPO),PAVT(NCOMPO),PAVRM(NCOMPO),PAVTM(NCOMPO)
      DIMENSION XAVR(NCOMPO),XAVRM(NCOMPO),XOUT(NCOMPO),PROUT(NCOMPO)
      DIMENSION PTOUT(NCOMPO),XINDHS(NCOMPO),XINDHB(NCOMPO)
C-------------------------------------------------------------------
C     Common block BINARY
C     contains binary data for stochastic 3-B binaries
      COMMON/BINARY/BODY1(NBINO),BODY2(NBINO),SIZE1(NBINO),SIZE2(NBINO),
     *ECC(NBINO),SEMIA(NBINO),RB(NBINO),TORB(NBINO),VR(NBINO),VT(NBINO),
     *EB(NBINO),RBMIN(NBINO),RBMAX(NBINO),XIN(NBINO),OMEG(NBINO),
     &XMRBIN(NJO),PHIBIN(NJO),
     *PHIOLD(NJO),EREM(NCOMPO),XHEAT(NCOMPO),XFEED(NCOMPO),
     *XKICK(NCOMPO),XMASS(NCOMPO),TBIN(NBINO),RHOBIN(NJO),VMRBIN(NJO)
      COMMON/BININT/INAME(NBINO),INAMEX(NBINO),
     *ISH(NBINO),ICO(NBINO),IEV(NBINO),NAMEB(NBINO),
     *IBSTA(NBINO),ILEN(NBINO)
      COMMON/BTMP/DMBIN(NBINO,NJO)
C-------------------------------------------------------------------
C    The following quantities are input parameters read from DKUG:
C     XMTOT(J)/XMIND(J) total mass/mass of individual for the component j
C     (note: READPM provides to give only one XMTOT(1) and
C            an index of mass spectrum in XMTOT(2)!)
C     RTOT(J)/RSC(J)/RIND(J) total/scaling/individual radius for comp. j
C     NOL : chooses initial density profile in DENS
C     EXP : chooses logarithmic density gradient (see DENS)
C     UINIT: initial velocity
C     TIND: internal temperature of component individuals
C     NOTE*******: The previous parameters can contain redundant
C                  information, for example:
C           TIND is meaningless for stellar components
C           EXP is meaningful only for one choice of NOL, and so on
C
C---------------------------------------------------------------------
C     XCOUL contains Coulomb-Logarithm (see READPM)
C---------------------------------------------------------------------
      DIMENSION XMTOT(NCOMPO),XMIND(NCOMPO),XMREM(NCOMPO),TAU(NCOMPO)
      DIMENSION DMESC(NCOMPO),DEESC(NCOMPO)
      DIMENSION DRHESC(NCOMPO),DERESC(NCOMPO),DETESC(NCOMPO)
      DIMENSION TMESC(NCOMPO),TEESC(NCOMPO),XMTOTI(NCOMPO)
      DIMENSION DMHOLE(5,NCOMPO),DELTM(5,NCOMPO)
      DIMENSION RTIDE(NCOMPO),XBIND(NCOMPO),RHOSTA(NCOMPO),RCRIT(NCOMPO)
      DIMENSION RTOT(NCOMPO),RSC(NCOMPO),XCOUL(NCOMPO)
      DIMENSION NOL(NCOMPO),EXP(NCOMPO)
      DIMENSION UINIT(NCOMPO),RIND(NCOMPO),TIND(NCOMPO)
      EQUIVALENCE (I,IS(1)),(IM,IS(2)),(IP,IS(3)),(L1,IS(4))
      EQUIVALENCE (ITIME,IS(5)),(L2,IS(6)),(LNJ,IS(7)),(ITVEC,IS(8))
      EQUIVALENCE (IUNIT,IS(9)),(ICTR,IS(10)),(NPR,IS(11))
      EQUIVALENCE (IPMAX,IS(12)),(ICOMP,IS(13)),(LEHOLE,IS(14))
      EQUIVALENCE (LLOSS,IS(15)),(LELOSS,IS(16)),(NRAND,IS(17))
      EQUIVALENCE (IBIN,IS(18)),(IEVENT,IS(19)),(IBINT,IS(20))
      EQUIVALENCE (IBREC,IS(21))
      EQUIVALENCE (IMR,IS(22)),(I00,IS(23)),(I10,IS(24)),(I20,IS(25))
      EQUIVALENCE (I02,IS(26)),(I30,IS(27)),(I12,IS(28))
      EQUIVALENCE (IMCRIT,IS(29)),(IMSAVE,IS(30)),(IFRFIN,IS(31))
C*********************************************************************
C*   CS(n) for n>200 is general working space
C*
C
      EQUIVALENCE (C23,CS(201)),(C56,CS(202)),(C13,CS(203))
      EQUIVALENCE (C53,CS(204)),(C12,CS(205)),(C32,CS(206))
      EQUIVALENCE (F2,CS(207)),(FNJ,CS(208)),(PI4,CS(209))
      EQUIVALENCE (CLFAC,CS(210)),(CLFAC2,CS(211)),(RHOC,CS(212))
      EQUIVALENCE (XCTOT,CS(213)),(XN3TOT,CS(214)),(ASPI,CS(215))
      EQUIVALENCE (AHEG,CS(216)),(BINFR,CS(217)),(PHTID,CS(218))
      EQUIVALENCE (W0,CS(219)),(RTIDAL,CS(220))
      EQUIVALENCE (XLR,CS(221)),(DER,CS(222)),(DERP,CS(223))
      EQUIVALENCE (RAV,CS(224)),(RAVI,CS(225)),(RAVP,CS(226))
      EQUIVALENCE (DERAV,CS(227)),(DYA,CS(228)),(C52,CS(229))
      EQUIVALENCE (CRX,CS(239)),(CBIN,CS(240)),(XTAN,CS(241))
C* Parameters for central black hole: XMHOLE,EFFI,RGRAV,XALPHA,XBETA
      EQUIVALENCE (XIMPL,CS(242)),(XBIN,CS(243)),(XMHOLE,CS(244))
      EQUIVALENCE (EFFI,CS(245)),(RGRAV,CS(246)),(DMTOT,CS(247))
      EQUIVALENCE (XALPHA,CS(248)),(XBETA,CS(249)),(ELL,CS(250))
      EQUIVALENCE (THETA,CS(251)),(DTTRX,CS(252)),(TTRX,CS(253))
      EQUIVALENCE (TBINI,CS(254)),(AEXP,CS(255)),(XLFAC,CS(256))
      EQUIVALENCE (RCORE,CS(257)),(XTEQ,CS(258)),(GAMMA,CS(259))
      EQUIVALENCE (XBIN2B,CS(260)),(FTRX,CS(261)),(XN3B,CS(262))
      EQUIVALENCE (DPTOT,CS(263)),(XMERR,CS(264)),(XEERR,CS(265))
      EQUIVALENCE (RTID0,CS(266)),(XMTOT0,CS(267)),(EBINI,CS(268))
      EQUIVALENCE (EMIN,CS(269)),(EMAX,CS(270)),(XMBIN0,CS(271))
      EQUIVALENCE (XMERR2,CS(272)),(XMERR3,CS(273)),(XMERR4,CS(274))
      EQUIVALENCE (XEERR2,CS(275)),(XEERR3,CS(276)),(XEERR4,CS(277))
      EQUIVALENCE (TUPIN,CS(278)),(RHOBAV,CS(279)),(XMBTOT,CS(280))
      EQUIVALENCE (XXESC,CS(281))
      EQUIVALENCE (UAV(1),CS(300)),(UAVP(1),CS(300+1*NCOMPO))
      EQUIVALENCE (VRAV(1),CS(300+2*NCOMPO)),(VTAV(1),CS(300+3*NCOMPO))
      EQUIVALENCE (PAVR(1),CS(300+4*NCOMPO)),(PAVRM(1),CS(300+5*NCOMPO))
      EQUIVALENCE (PAVT(1),CS(300+6*NCOMPO)),(PAVTM(1),CS(300+7*NCOMPO))
      EQUIVALENCE (XAVR(1),CS(300+8*NCOMPO)),(XAVRM(1),CS(300+9*NCOMPO))
      EQUIVALENCE (PROUT(1),CS(300+10*NCOMPO))
      EQUIVALENCE (PTOUT(1),CS(300+11*NCOMPO))
      EQUIVALENCE (XOUT(1),CS(300+12*NCOMPO))
      EQUIVALENCE (XMTOT(1),CS(300+13*NCOMPO))
      EQUIVALENCE (XMIND(1),CS(301+14*NCOMPO))
      EQUIVALENCE (XMREM(1),CS(301+15*NCOMPO))
      EQUIVALENCE (TAU(1),CS(301+16*NCOMPO))
      EQUIVALENCE (RTOT(1),CS(301+17*NCOMPO)),(RSC(1),CS(301+18*NCOMPO))
      EQUIVALENCE (NOL(1),CS(301+19*NCOMPO)),(EXP(1),CS(301+20*NCOMPO))
      EQUIVALENCE (UINIT(1),CS(301+21*NCOMPO))
      EQUIVALENCE (XCOUL(1),CS(301+22*NCOMPO))
      EQUIVALENCE (RIND(1),CS(301+23*NCOMPO))
      EQUIVALENCE (TIND(1),CS(301+24*NCOMPO)) 
      EQUIVALENCE (RTIDE(1),CS(301+25*NCOMPO))
      EQUIVALENCE (RHOSTA(1),CS(301+26*NCOMPO))
      EQUIVALENCE (XBIND(1),CS(301+27*NCOMPO))
      EQUIVALENCE (RCRIT(1),CS(301+28*NCOMPO))
      EQUIVALENCE (DMHOLE(1,1),CS(301+29*NCOMPO))
      EQUIVALENCE (DELTM(1,1),CS(301+34*NCOMPO))
      EQUIVALENCE (TMESC(1),CS(301+39*NCOMPO))
      EQUIVALENCE (TEESC(1),CS(301+40*NCOMPO))
      EQUIVALENCE (DMESC(1),CS(301+41*NCOMPO))
      EQUIVALENCE (DEESC(1),CS(301+42*NCOMPO))
      EQUIVALENCE (AFAC(1),CS(301+43*NCOMPO))
      EQUIVALENCE (ATFAC(1),CS(301+44*NCOMPO))
      EQUIVALENCE (XBINY(1),CS(301+45*NCOMPO))
      EQUIVALENCE (XFBINY(1),CS(301+46*NCOMPO))
      EQUIVALENCE (XSBINY(1),CS(301+47*NCOMPO))
      EQUIVALENCE (XMASSY(1),CS(301+48*NCOMPO))
      EQUIVALENCE (XINDHS(1),CS(301+49*NCOMPO))
      EQUIVALENCE (XINDHB(1),CS(301+50*NCOMPO))
      EQUIVALENCE (DRHESC(1),CS(301+51*NCOMPO))
      EQUIVALENCE (DERESC(1),CS(301+52*NCOMPO))
      EQUIVALENCE (DETESC(1),CS(301+53*NCOMPO))
      EQUIVALENCE (XMTOTI(1),CS(301+54*NCOMPO))
C* Values until CS(300+55*NCOMPO) occupied
