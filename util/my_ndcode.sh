#!/usr/bin/zsh

echo "                                                           "
echo "-----my_ndcode Shellscript for encoding unformatted data-----"
echo "              (20.2.2004)                                    "
echo "     200 Gridpoints assumed (?)                                "
echo "     <<< Must be used in simulation directory >>>"


[[ -f fort.2 ]] || {
	if [[ -f star1.cwr ]]; then
		ln -s star1.cwr fort.2
	else
		print -u2 "!!! file fort.2 or star1.cwr not found !!!"
		exit 1
	fi 
}
[[ -f fort.3 ]] || {
	if [[ -f star1.inf ]]; then
		ln -s star1.inf fort.3
	else
		print -u2 "!!! file fort.3 or star1.inf not found !!!"
		exit 1
	fi 
}

[[ -f ndcode.indat ]] || {
	cat <<EOF > ndcode.indat
 MODMIN=    1 MODINT=    1 MODMAX=99999 NCHR=  4 NASCII= 50 IPR=  0
 NMIN=  1 NINT=  2 NMAX=99999 NCHR=  4 NASCII= 50
EOF
}
cp ndcode.indat fort.15

[[ -f ndcode ]] || {
    print -u2 "!!! File ndcode not found !!!"
    exit 1
}

echo "...start encoding..."
./ndcode | tee ndcode.list
echo "...stop encoding..."

rm ndcode.indat fort.15
mv fort.10 star1.enccwr
mv fort.11 star1.encinf
