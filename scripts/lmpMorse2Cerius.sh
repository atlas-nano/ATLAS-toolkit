#!/bin/bash

if [ $# -lt 3 ]; then
	echo "usage: $0 d alpha x0 [energy_units: 1=kcal/mol(default) 2=eV 3=hartree]"
	exit 1
fi

scale=1
if   [ $# -gt 3 ] && [ `echo $4 | egrep -c '^2$'` -eq 1 ]; then
	scale=23.06
elif [ $# -gt 3 ] && [ `echo $4 | egrep -c '^3$'` -eq 1 ]; then
	scale=627.509
fi

echo $1 $2 $3 | awk -v s=$scale '{print $2*$2*2*$1*s,$3,$1*s}'
