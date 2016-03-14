*     Parameter File   - change here for number of components       *
C*******************************************************************|
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
C NJMAX= Number of Grid Points                                      |
C NW= Dimension of Vector for Data storage on file                  |
C       (see main program)                                          |
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
      PARAMETER (NCMPMX=1,NEQC=6,NPOSC=3,NADD=0,NEQA=0,NPOSA=0)
      PARAMETER (NEQX=1,NPOSX=1)
      PARAMETER (NGMAX=NCMPMX*NEQC+NADD*NEQA+NEQX)
      PARAMETER (NPOSMX=NCMPMX*NPOSC+NADD*NPOSA+NPOSX,NJMAX=200)
C Parameters NCOMPO,NGO,NJO ensure efficient alignment in COMMON
      PARAMETER (NCOMPO=2*(NCMPMX/2+1),NGO=2*(NGMAX/2+1),
     *           NJO=2*(NJMAX/2+1))
C Parameter IVDIM is length of vector for time dependent data storage VIT
      PARAMETER (IVDIM=350*NCOMPO)
      PARAMETER (NBIN=30100)
      PARAMETER (NBINO=2*(NBIN/2+1))
C-----NW=Length of Storage Vector EI (see equiv.f and readmo.f/print.f)
      PARAMETER (NW=(2*NGO+2*NCOMPO+1)*NJO+20*NCOMPO+120)
C Parameters for length of binary save common data
      PARAMETER (NBEI=NBINO*16 + 5*NJO + 5*NCOMPO)
      PARAMETER (NIBEI=8*NBINO)
