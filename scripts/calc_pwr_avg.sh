#!/bin/bash

if [ $# -eq 0 ]; then
	echo "usage: $0 pwr_prefix [savename]"
	exit 1
fi

dirn=$(dirname $1)
flen=$(basename $1)

temp=($(find $dirn -name "${flen}.*.pwr"))
if [ ${#temp[*]} -eq 0 ]; then
	echo "ERROR: No valid files found while searching '${dirn}/${flen}.*.pwr'"
	exit 1
fi
if [ ${#temp[*]} -lt 2 ]; then
	echo "ERROR: Need at least 2 pwr files. Found only 1 while searching '${dirn}/${flen}.*.pwr'"
	exit 1
fi

plist=()
tlist=()

for i in $(seq 1 ${#temp[*]})
do
	pwr=${temp[i-1]}
	prefix=${pwr%.*}
	thermo=${prefix}.thermo
	if [ -e $thermo ]; then
		plist=(${plist[*]} $pwr)
		tlist=(${tlist[*]} $thermo)
	fi
done

if [ ${#plist[*]} -lt 2 ]; then
	echo "ERROR: Cannot find at least 2 pwr and matching thermo files while searching '${dirn}/${flen}.*.pwr'"
	exit 1
fi

sname=$(basename ${plist[0]})
sname=${sname%.*}
sname="${sname}.avg.dat"
if [ $# -gt 1 ]; then
	sname=$2
fi

echo "Averaging over ${#plist[*]} files and saving to $sname"

vals=($(awk '{if($1 ~ /^property/){nm[NF]++; nv[NF]=$NF}}END{n=asorti(nm,nm_sort); mx=nm_sort[n]; print mx,nv[mx]}' ${tlist[*]} | sed 's/ .*G\([0-9]\+\)\]/ \1/'))
ngrp_thermo=$(echo ${vals[1]} | awk '{print $1+0}')
ncolgrp_thermo=$(echo ${vals[0]} $ngrp_thermo | awk  '{print ($1-1)/$2}')
echo "ngrp_thermo: $ngrp_thermo ncolgrp_thermo: $ncolgrp_thermo"

vals=($(awk '{if($1 ~ /^freq/){nm[NF]++; nv[NF]=$NF}}END{n=asorti(nm,nm_sort); mx=nm_sort[n]; print mx,nv[mx]}' ${plist[*]} | sed 's/ .*G\([0-9]\+\)\]/ \1/'))
ngrp_pwr=$(echo ${vals[1]} | awk '{print $1+0}')
ncolgrp_pwr=$(echo ${vals[0]} $ngrp_pwr | awk  '{print ($1-1)/$2}')
echo "ngrp_pwr: $ngrp_pwr ncolgrp_pwr: $ncolgrp_pwr"

norm=($(awk -v ncg=$ncolgrp_thermo -v ngrp=$ngrp_thermo '{if($1 ~ /^nmol/ && NF==(1+ncg*ngrp)) { i=3; while(i<NF) { print $i; i+=ncg}}}' ${tlist[*]}))

echo "norm: ${norm[*]}"
gawk -v n_list="${norm[*]}" -v ncg=$ncolgrp_pwr -v ngrp=$ngrp_pwr '
BEGIN{
	split(n_list,norm," ");
	n=0
	cstart=-ngrp
}{
	t=NF
	if (FNR==3 && NF==(1+ngrp*ncg)) {
		n++
		cstart+=ngrp
	}
	if ($1 ~ /^[0-9]/) {
		c=cstart
		j=0
		for(i=2;i<=NF;i++) {
			if((j%ncg)==0)
				c++
			$i/=norm[c]
			val[$1][i]+=$i
			val2[$1][i]+=$i*$i
			j++
		}
	}
}
END{
	for (i in val) {
		printf "%-10.5f ",i
		for(j=2;j<=t;j++) {
			avg=val[i][j]/n
			avg2=val2[i][j]/n
			v=avg2-avg*avg
			std=0
			if(v>0)
				std=sqrt(avg2-avg*avg)
			printf "%10.5f %10.5f ",avg,std
		}
		printf "\n"
	}
}' ${plist[*]}  | sort -n > $sname
