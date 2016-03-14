#/usr/bin/ksh -x

cd ~/GAS_MODEL/spedi-heid03

Ncomp_List="$@"

[[ -z $Ncomp_List ]] && Ncomp_List="1 2 5 10 15 20 25 30 40 50"

[[ -f params.f.SVG ]] || cp params.f params.f.SVG

for Ncomp in $Ncomp_List; do
    cat params.f.SVG | sed 's/NCMPMX= *[1-9][0-9]* */NCMPMX='$Ncomp'/' > params.f
    make && mv star exec/star_$Ncomp""comp || {
	print -u2 '!! Error when compiling spedi for '$Ncomp' components !!'
	exit 1
    }
done
