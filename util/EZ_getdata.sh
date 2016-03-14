#!/usr/bin/env zsh

inf2asc_path=~/bin/inf2asc 

echo ""
echo "-----EZ_getdata.sh Shellscript for extracting data from star1.inf/fort.3-----"
echo "                              (26.2.2004)                                    "
echo "               <<< Must be used in simulation directory >>>"
echo ""

QuantityCodes=( "$@" ) # Put list of quantities to extract as options on the command line
                       # look for codes in fort.90

[[ -z $QuantityCodes ]] && {
    print -u2 "!!! You must provide list of quantities to extract as options on the command line !!!"
    print -u2 "    Look for numeric codes in fort.90"
    exit 1
}

# Look for binary data file

[[ -f star1.inf ]] && input_file=star1.inf
[[ -z $input_file ]] && [[ -f fort.3 ]] && input_file=fort.3
[[ -z $input_file ]] && {
    print -u2 "!!! Binary \"inf\" file not found (star1.inf and/or fort.3) !!!"
    exit 1
}

# Determine number of components from kug1.dat

[[ -f kug1.dat ]] || {
    print -u2 "!!! kug1.dat not found !!!"
    exit 1
}

NCOMP=$(gawk '/^ *NCOMP=/ {print $2; exit}' < kug1.dat)

# Run the "translator" program

data_file="data"
for qq in $QuantityCodes; do
    data_file=$data_file"_"$qq
done
data_file=$data_file".txt"

if [[ -f $inf2asc_path ]]; then
    print -u2 "> Running $inf2asc_path $input_file $NCOMP $QuantityCodes"
    eval "$inf2asc_path $input_file $NCOMP $QuantityCodes > $data_file"
    ln -sf $data_file my_getdata.out
    print -u2 "> Data written into $data_file (with my_getdata.out linked to it)"
else
    print -u2 "!!! Executable $inf2asc_path not found !!!"
fi
