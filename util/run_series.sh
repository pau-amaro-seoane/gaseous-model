#/bin/ksh

DirList="$@"

[[ -z $DirList ]] && DirList=$(find . -name "kug1.dat" -print | sed 's/\/kug1\.dat//')

echo "> List of directories in which to run simulations : $DirList" 

base_dir=$(pwd)

for dir in $DirList; do
    if cd $dir; then
	if [[ -f star1.list ]] || [[ -f star1.list.gz ]]; then
	    echo "!!! star1.list(.gz) already exists in $dir !!!"
	else
	    echo "> Running spedi simulation in $(pwd)..."
	    nohup ./@run1r star &
	    wait
	    echo "< Simulation in $dir is seemingly over; Compressing files..."
            rm star1.list
	    gzip --force fort.* star1.list.* star1.last 
	    ln -s star1.list.*.gz star1.list.gz 
	fi
	cd $base_dir
    else
	echo "!!! cannot cd to $dir !!!"
    fi
done