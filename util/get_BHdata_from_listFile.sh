#!/bin/zsh

tmp=/tmp/tmp_$$

gunzip --force --to-stdout > $tmp

print '# Data from spedi simulation'

# Get gamma value (Coulomb log)

gamma=$(cat $tmp | grep GAMMA | head -n 1 | sed 's/.*GAMMA= *//; s/ *XBIN2.*//; s/D/e/')

# Get star number

Nstar=$(cat $tmp | grep XNTOTA | head -n 1 | sed 's/.*XNTOTA= *//; s/ *XIMF.*//; s/D/e/')

print '# gamma= '$gamma' Nstar= '$Nstar
print '# 1: t_nb 2: mbh_nb 3: dmbh_nb'

cat $tmp | gawk ' 
    /^ *T=.*DT=.*RHO=/ {
	t=$2; rho=$6
    }
    /^ *XMHOLE=.*DMTOT=.*DMHOLE\(1-5,1\)=/{
	mbh=$2
        dmbh=$8
	print t,mbh,dmbh
    }'

rm $tmp
