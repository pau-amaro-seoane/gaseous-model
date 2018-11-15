#!/usr/bin/env zsh

# Define the run number (output files will have this appended)
export RUN=1

# Copy
cp kug$RUN.dat fort.15

# Backup, in case these files are around
if [ -f star$RUN.inf  ] ; then ln star$RUN.inf fort.3            ; fi
if [ -f star$RUN.cwr  ] ; then ln star$RUN.cwr fort.2            ; fi
if [ -f star$RUN.bin  ] ; then ln star$RUN.bin fort.1            ; fi
if [ -f star$RUN.rand ] ; then ln star$RUN.rand fort.91          ; fi
if [ -f fort.61       ] ; then mv fort.61 star$RUN.bindat.$$.old ; fi

# Some log at the head of the star file
echo "Run module used: " `date` ../spedi >> star$RUN.list.$$

# Launch it with the proviso that some files are not around.
# If they are, backup them and exit.

../spedi >> star$RUN.list.$$ ||\
            if [ -f fort.4  ] ; then mv fort.4 star$RUN.last       ; fi
            if [ -f fort.11 ] ; then mv fort.11 star$RUN.lastbin   ; fi
            if [ -f fort.61 ] ; then mv fort.61 star$RUN.bindat.$$ ; fi
            if [ -f fort.85 ] ; then mv fort.85 fort.85.$$         ; fi
            if [ -f fort.86 ] ; then mv fort.86 fort.86.$$         ; fi
            if [ -f fort.87 ] ; then mv fort.87 fort.87.$$         ; fi
            if [ -f fort.88 ] ; then mv fort.88 fort.88.$$         ; fi

# Backup new results
if [ -f fort.4  ] ; then mv fort.4 star$RUN.last       ; fi
if [ -f fort.11 ] ; then mv fort.11 star$RUN.lastbin   ; fi
if [ -f fort.61 ] ; then mv fort.61 star$RUN.bindat.$$ ; fi
if [ -f fort.85 ] ; then mv fort.85 fort.85.$$         ; fi
if [ -f fort.86 ] ; then mv fort.86 fort.86.$$         ; fi
if [ -f fort.87 ] ; then mv fort.87 fort.87.$$         ; fi
if [ -f fort.88 ] ; then mv fort.88 fort.88.$$         ; fi
