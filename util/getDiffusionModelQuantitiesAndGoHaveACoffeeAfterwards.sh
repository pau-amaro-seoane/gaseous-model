#!/bin/zsh

gawk '
    BEGIN{Radius[0]=-1;T=0}
    /^ *T= *.*DT=/ { T=$2 }
    /^ *RADIUS  *TOTMASS/ {
	#print "> Found radius information"
	getline
	i=0
	while (NF>0) {
	    Radius[i++]=$2
	    getline
	}
    }
    /^ *\*\*JC\*\*I\=K\(T\)\=/ {
	#print "> Found loss cone information", Radius[0]
	if (Radius[0]>=0) {
	    filename="DiffusionModelQuantities_" iFile++ ".asc"
	    print "> Creating " filename
	    print "# Time = " T " N-body units" > filename
	    print "# 0: Radius 1: K        2: MDOTLC  3: XER     4: XET     5: PF      6: TH_D    7: TH_LC   8: THLC_OM 9: UL_SR   10: VLC_ST 11: TIN    12: TOUT   13: DOMEGA" > filename
	    getline
	    i=0
	    while (NF==15) {
		print "  " Radius[i++],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15 > filename
		getline
	    }
	    Radius[0]=-1
	}
    }
'
