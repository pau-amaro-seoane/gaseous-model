#!/bin/ksh

tmp=/tmp/tmp_$$

gunzip --force --to-stdout | head -n 200 > $tmp

# Get gamma value (Coulomb log)

gamma=$(cat $tmp | grep GAMMA | head -n 1 | sed 's/.*GAMMA= *//; s/ *XBIN2.*//; s/D/e/')

# Get star number

Nstar=$(cat $tmp | grep XNTOTA | head -n 1 | sed 's/.*XNTOTA= *//; s/ *XIMF.*//; s/D/e/')

print 'Nstar= '$Nstar' gamma= '$gamma

rm $tmp

