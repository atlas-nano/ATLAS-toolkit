#!/bin/bash

if [ $# -lt 3 ]; then
	echo "usage: $0 bgf_file ff_file target_dens [savename]"
	exit 1
fi

if [ ! -e $1 ] || [ ! -r $1 ]; then
	echo "ERROR: Cannot access $1"
	exit 1
fi
#if [ ! -e $2 ] || [ ! -r $2 ]; then
#	echo "ERROR: Cannot access forcefield file $2"
#	exit 1
#fi

c=$(echo $3 | awk '{if($1 ~ /^[0-9]+(\.[0-9]+){0,1}$/) print 1; else print 0}')
if [ $c -eq 0 ]; then
	echo "ERROR: Expected number for target_dens. Got '$3'"
	exit 1
fi
savename=$(basename $1)
savename=${savename%.*}
savename="${savename}.${3}dens.bgf"
if [ $# -gt 3 ]; then
	savename=$4
fi

#first remove any current box
/home/tpascal/scripts/removeBGFBox.pl -b $1 -s tmp.bgf > /dev/null
#add a box
/home/tpascal/scripts/addBoxToBGF.pl tmp.bgf tmp.bgf > /dev/null
#get mass
mass=$(/home/tpascal/scripts/bgfmass.pl -b tmp.bgf -f $2 | grep '^Adding' | awk '{print $4}')
#get correct cell
cell=($(grep CRYSTX tmp.bgf | awk -v d=$3 -v m=$mass '{cvol=$2*$3*$4; tvol=m/d/.6023; s=(tvol/cvol)**(1/3); print $2*s,$3*s,$4*s,$5,$6,$7}'))
#update bgf cell
/home/tpascal/scripts/updateBGFBox.pl -b tmp.bgf -s $savename -c "${cell[*]}" > /dev/null
rm -fr tmp.bgf
