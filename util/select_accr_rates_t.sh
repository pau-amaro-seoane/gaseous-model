#!/bin/zsh 
 
tmp=/tmp/tmp_$$
 
cat > $tmp
 
print '# Data from spedi simulation'
 
# Get gamma value (Coulomb log)
 
gamma=$(cat $tmp | grep GAMMA | head -n 1 | sed 's/.*GAMMA= *//; s/ *XBIN2.*//; s/D/e/')
 
# Get star number
 
Nstar=$(cat $tmp | grep XNTOTA | head -n 1 | sed 's/.*XNTOTA= *//; s/ *XIMF.*//; s/D/e/')
 
print '# gamma='$gamma' Nstar='$Nstar
print '# 1: t_nb 2: m_acc_nb'
 
cat $tmp | gawk '
    /^ *T=.*DT=.*RHO=/ {
        t=$2; rho=$6
    }
    /^ *XMHOLE=.*DMTOT=.*DMHOLE\(1-5,1\)=/{
        m_acc=$8
        print t, m_acc
    }'
 
rm $tmp

