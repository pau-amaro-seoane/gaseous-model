      =============================================================================
       Lame notes to use "spedi" (gas dynamical models of spherical star clusters)
      =============================================================================

    +-------------------------------------------------------------------------+
    | @run1 is a tcsh-script required to start the run                        |
    |                                                                         |
    | @run1r is the same but it can be used to restart ("r") a simulation     |
    | without deleting previous data. In practise I was using this script all |
    | the time                                                                |
    |                                                                         |
    | kug1.dat contains what you want to be the initial system                |
    |                                                                         |
    | To start a simulation: nohup @run1 spedi &                              |
    |                                                                         |
    | where "spedi" is the executable binary you get out when you compile the |
    | main code with "make"                                                   |
    +-------------------------------------------------------------------------+

 New style: Use launch_spedi.sh in utils instead of @run1
 ========

(0) General information: Refer to star.f and params.f

(1) see Run_SglMass_Expl for a directory ready to run single-mass model (use
    /@run1rstar_sglmass)
    see Run_MultiMass_Expl for a directory ready to run multi-component model (use
    star_Ncomp)

To run a multicomponent model, one has to compile spedi with
NCMPMX=number_of_comp in params.f

Use XMINF<0 in kug1.dat (IMF exponent) to generate automatically the mass
component with a given IMF exponent. The minimum mass corresponds to the first
component in kug1.dat. The max mass corresponds to the second component in the
file. The other components are not used but have to be present in the file
(dummies). Beware of the format of the file: strictly respect the number of
characters for each field. For instance "XIMF=-2.35000+000" is correct but not
"XIMF= -2.35000+000"...

(2) In Simula_BH_test/ you'll find an example to start a simulation of a star-accreting BH

(3) In plot/ you'll find scripts to process/ analyse data relative to
    ndcode/ (check plot/README)

(4) In kug1_examples/ you'll find "good" kug1.dat examples

(5) utils/ and smac/ are tools for plotting etc
   (smac contains supermongo macros and in utils you'll find scripts)
