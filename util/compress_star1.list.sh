#!/usr/bin/env zsh

base_dir=$1
[[ -z $base_dir ]] && base_dir=.

init_dir=$(pwd)
for directory in $(
	find $base_dir -name "star1.list*.gz" -or -name "star1.list.*" -print |\
	sed 's/\/[^\/]*$//' | sort -u
); do
	echo ">> " $directory
	cd $directory
	last_file=$(ls -1rt star1.list.* | grep 'star1.list.[1-9][0-9]*$' | tail -n1)
	[[ ! -z $last_file ]]  && mv -fv $last_file star1.list && gzip -v star1.list && {
		list=( $(ls -1rt star1.list.* | grep 'star1.list.[1-9][0-9]*$') )
		[[ -z $list ]] || rm -f $list
	}
	cd $init_dir
done
