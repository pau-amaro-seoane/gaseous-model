#/usr/bin/ksh

for file in $(find . -name "star1.list.gz"); do
	realpath $file | tr '/' '\012' | grep ':' | sed 's/..*://; s/p/./' | tr '\012' ' ' # alpha n_comp delta m_max
	echo
	gunzip --to-stdout $file | grep 'T=..* DT=' | tail -n1 | sed 's/^ *T=//; s/^  *//; s/ .*//' # Tend
	gunzip --to-stdout $file | head -n 200 | grep '^ *I=' |\
	gawk '{
		if (Mmin==0) { Mmin=$8 }
		Mmax=$8
		Ntot=Ntot+$4; Mtot=Mtot+$6
	    } 
	    END{print Mmax, Mmin, Mtot/Ntot}' # Mmax/Mmin, M_avrg
done | gawk 'BEGIN{
	Npart=1e6; gamma_Coul=0.11
	print "# Tend in FP units (what else?)"
	print "# 1: alpha 2: n_comp 3: delta 4: mmax_over_min_nom  5: mmax_over_min 6: mmax_over_mavrg 7: Tend"}
	{ 
		alpha=$1; n_comp=$2; delta=$3; mmax_nom=$4;
		getline
		Tend=$1
		getline
		mmax=$1; mmin=$2; m_avrg=$3
		print sprintf("%5.2f",alpha)" "sprintf("%2.0f",n_comp)" "sprintf("%5.2f",delta)" "sprintf("%8.2f",mmax_nom)" "sprintf("%8.2f",mmax/mmin)" "sprintf("%8.2f",mmax/m_avrg)" "sprintf("%10.3e",Tend/Npart*log(gamma_Coul*Npart))
	}'
	
