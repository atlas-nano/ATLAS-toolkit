#!/bin/bash

if [ $# -lt 2 ]; then
	echo "usage: $0 bgf_file forcefield"
	exit 1
fi

if [ ! -e $1 ] || [ ! -r $1 ]; then
	echo "ERROR: Cannot access $1"
	exit 1
fi

if [ $(echo $2 | egrep -c '\s') -eq 0 ] && [ ! -e $2 ]; then
	if [ ! -e ~/ff/${2}.ff ]; then
		echo "ERROR: Cannot locate the forcefield $2"
		exit 1
	fi
fi
cell=(`grep CRYSTX $1 | tail -1 | awk '{print $2,$3,$4}'`)
if [ ${#cell[*]} -lt 3 ]; then
	echo "ERROR: Cannot determine cell info in $1"
	exit 1
fi
~/scripts/bgfmass.pl -b $1 -f "$2" | grep '^Addi' | awk -v c_str="${cell[*]}" 'BEGIN{split(c_str,cell," ")}{print $4/.6023/cell[1]/cell[2]/cell[3]}'
