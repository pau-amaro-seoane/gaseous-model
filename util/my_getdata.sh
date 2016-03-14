#!/usr/bin/zsh

echo ""
echo "-----my_getdata.sh Shellscript for extracting data from star1.encinf-----"
echo "              (20.2.2004)                                    "
echo "     <<< Must be used in simulation directory >>>"
echo ""

QuantityCodes=( "$@" ) # Put list of (up to 5) quantities to extract as options on the command line
                       # look for codes in fort.90

[[ -z $QuantityCodes ]] && {
    print -u2 "!!! You must provide list of (up to 5) quantities to extract as options on the command line !!!"
    print -u2 "    Look for codes in fort.90"
    exit 1
}

[[ -f star1.encinf ]] || {
    print -u2 "> File star1.encinf not found; running ndcode..."
    [[ -f ./my_ndcode.sh ]] && ./my_ndcode.sh || ~/GAS_MODEL/util/my_ndcode.sh
}

[[ -f getdata ]] || {
    print -u2 "!!! File getdata not found !!!"
    exit 1
}

# 150  FORMAT(' IW=',I1,' LEG=',L2,' LT5=',L2,' KUPA=',L2,' SCALE=',L2,' N1=',I4,' N2=',I4,' N3=',I4,' N4=',I4,' N5=',I4)
# IW=_ LEG=__ LT5=__ KUPA=__ SCALE=__ N1=____ N2=____ N3=____ N4=____ N5=____

typeset -R4 N1
typeset -R4 N2
typeset -R4 N3
typeset -R4 N4
typeset -R4 N5
data_file="data"

call_getdata=0
Nquant=$(echo $QuantityCodes | wc -w | sed 's/ //g')

N2=-2
N3=-2
N4=-2
N5=-2

prov=prov.$$
prov2=prov2.$$

i=1
j=1

for qq in $QuantityCodes; do
    eval N$i=$qq
    i=$(( i + 1 ))
    j=$(( j + 1 ))
    (( i > 5 )) && call_getdata=1
    (( j > Nquant )) && call_getdata=1
    data_file=$data_file"_"$qq

    if (( call_getdata == 1 )); then
	print -u2
        print -u2 "> Calling getdata for quantities N1=$N1 N2=$N2 N3=$N3 N4=$N4 N5=$N5"
	print -u2
	cat <<EOF > splo.plot999
1
1000
'XXX'
'XXX'
'XXX'
0.,2000000.,-4.,8.
 LALF= T. LT4= F. LDEC= T. LPR= F. LCL= T.
 NPLOT=  0 NTIME=  1 NMOMIN=   1 NMOMAX=9999 NINT=1 NCHR= 4 NASCII= 50
 IW=1 LEG=T. LT5=F. KUPA=F. SCALE=F. N1=$N1 N2=$N2 N3=$N3 N4=$N4 N5=$N5
 XLEG= 1.00-001 YLEG= 1.00-002 XT5= 0.00+000 YT5= 0.00+000 LSMO= F.
 UN1-5= 1.00+000 1.00+000 1.00+000 1.00+000 1.00+000

EOF

	print 99\\n999\\nstar1\\n | ./getdata | tee getdata.msg.$$
	print -u2 "> getdata exited"
	rm splo.plot999

	[[ -f plsp99.dat ]] || {
	    print -u2 "!!! getdata hasn't created the data file (plsp99.dat) !!!"
	    exit 1
	}

	[[ -z nlines ]] || nlines=$(
            grep '^ *From..*written' < getdata.msg.$$ |\
            sed 's/ *written *$//; s/^.* //' |\
            gawk '{sum=sum+$1}END{print sum}'
        )

	head -n $nlines plsp99.dat |\
	if [[ -f $prov ]]; then
	    paste $prov - 
	else
	    cat
	fi > $prov2
	mv $prov2 $prov
	rm plsp99.dat getdata.msg.$$

	N1=-2
	N2=-2
	N3=-2
	N4=-2
	N5=-2
        i=1
	call_getdata=0
    fi
done
    
data_file=$data_file".txt"
mv $prov $data_file
ln -sf $data_file my_getdata.out
rm plsp99.cmd