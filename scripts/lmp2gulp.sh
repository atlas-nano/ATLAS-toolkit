#!/bin/bash

get_val_coeffs()
{
	#first get the valence types
	vstyle=""
	if [ $(grep -i "^${vt}_style" $lmp_ctl | sed 's/#.*//' | awk '{print NF-1}') -eq 1 ]; then
		vstyle=$(grep -i "^${vt}_style" $lmp_ctl | awk '{print $2}')
	fi
	#get coeffs from data file
	vname="${vt^}"
	vcoeffs=($(awk -v vnm=$vname -v vs=$vstyle '
BEGIN{
	valid = 0
}
{
	if($1 ~ vnm && $2 ~ /^Coeffs/)
		valid = 1
	else if (NF == 2 || (NF == 2 && $1 !~ vnm))
		valid = 0
	else if (valid == 1 && NF > 2) {
		printf "%d-->%s=",$1,vs;
		i=2;
		while(i<=NF) {
			if($i ~ /#/)
				i=NF+1;
			else {
				printf "%f-->",$i+0;
				i++;
			}
		}
		printf "\n";
	}
}' $lmp_data))
	#get coeffs from input/control file
	vcoeffs=(${vcoeffs[*]} $(grep -i "^${vt}_coeff" $lmp_ctl | awk -v vs=$vstyle '
{
	if(length(vs)>0) {
		vtype=vs; 
		i=3;
	} else { 
		vtype=$3; i=4; 
	} 
	printf "%d-->%s=",$2,vtype;
	while(i<=NF) {
		if($i ~ /#/)
			i=NF+1
		else {
			printf "%f-->",$i+0;
			i++;
		}
	}
	printf "\n"
}'))
	for i in $(seq 1 ${#vcoeffs[*]})
	do
		idxvals=($(echo ${vcoeffs[i-1]} | sed 's/\=/ /'))
		idx=($(echo ${idxvals[0]} | sed 's/-->/ /'))
		vtypes[${idx[0]}]=${idx[1]}
		vvals[${idx[0]}]=${idxvals[1]}
	done
	#now traverse the valence list and record unique types
	vlist=($(awk -v nfield=$n -v vnm=$vname -v a_str="${atoms[*]}"  -v vt_str="${vtypes[*]}" -v vv_str="${vvals[*]}" '
BEGIN {
	split(vt_str,vtypes," ")
	split(vv_str,vvars," ")
	split(a_str,atoms," ")
}
{
	if($1 ~ vnm && NF == 1)
		valid=1
	else if (NF == 1 || (NF == 1 && $1 !~ vnm))
		valid=0
	else if (valid == 1 && NF == (nfield+2))	{
		vt=vtypes[$2]
		vv=vvars[$2]
		ty_str=atoms[$3]
		found=0
		for(i=4;i<=NF;i++)
			ty_str = sprintf("%s-->%s",ty_str,atoms[$i])
		if(flist[ty_str]==1)
			found=1
		ty_str=atoms[$NF]
		for(i=NF-1;i>2;i--)
			ty_str = sprintf("%s-->%s",ty_str,atoms[$i])
		if(flist[ty_str]==1)
			found=1
		if(found==0) {
			flist[ty_str]=1
			printf "%s-->%s=%s\n",vt,ty_str,vv
		}
	}
}' $lmp_data | sed 's/-->$//' | sed 's/\=/-->/'))
}

print_val_coeffs()
{
	echo ${vals[*]} | awk -v n=$n -v s_str="${vscale[*]}" '
BEGIN {
	split(s_str,scale," ")
}
{
	for(i=1;i<=n;i++) {
		j=i+1
		printf "%s ",$j
	}
	j=1
	for(i=n+2;i<=NF;i++) 
		printf "%s ",$i*scale[j++]
	printf "0 10\n"
}'
}

if [ $# -lt 2 ]; then
	echo "usage: $0 lmp_control_file lmp_data_file > [gulp_file]"
	exit 1
fi

lmp_ctl=$1
lmp_data=$2
gin="$(echo $lmp_ctl | sed 's/^in.//').gin"
if [ $# -gt 2 ]; then
	gin=$3
fi

#first find masses line in data file and record the fftypes
ffids=()
ffnames=()
fftypes=()
aids=($(awk '
{
	if($1 ~ /^Masses/ && NF == 1)
		valid=1; 
	if(valid==1 && NF > 3) 
		printf "%s-->%s\n",$4,$1; 
	else if ($1 ~ /^Atoms/)
		valid=0
}' $lmp_data))
if [ ${#aids[*]} -eq 0 ]; then
	echo "ERROR: Cannot find Masses or fftype lines in $lmp_data"
	exit 1
fi
for i in $(seq 1 ${#aids[*]})
do
	vals=($(echo ${aids[i-1]} | sed 's/-->/ /'))
	fftypes[${vals[1]}]=${vals[0]}
	ffids=(${ffids[*]} ${vals[1]})
	ffnames=(${ffnames[*]} ${vals[0]})
done

#now find the Atoms line and record the fftype of each atom
atoms=()
aids=($(awk -v id_str="${ffids[*]}" -v nm_str="${ffnames[*]}" '
BEGIN{
	n=split(id_str,id," ")
	split(nm_str,nm," ")
	for(i=1;i<=n;i++) {
		fftype[id[i]]=nm[i]
	}
	valid=0	
}
{
	if($1 ~ /^Atoms/ && NF==1)
		valid=1
	else if (NF == 1)
		valid=0
	else if (valid == 1 && NF > 5)
		printf "%d-->%s\n",$1,fftype[$3]
	}' $lmp_data))
for i in $(seq 1 ${#aids[*]})
do
	vals=($(echo ${aids[i-1]} | sed 's/-->/ /'))
	atoms[${vals[0]}]=${vals[1]}
done

#BONDS
vt="bond"
n=2
get_val_coeffs
#now print bond info in gulp format
for i in $(seq 1 ${#vlist[*]})
do
	vals=($(echo ${vlist[i-1]} | sed 's/-->/ /g'))
	if [ "${vals[0]}" == "harmonic" ]; then
		echo "${vals[0]} bond kcal"
		vscale=(2 1)
		print_val_coeffs
	fi
done

#ANGLES
vt="angle"
n=3
get_val_coeffs
#now print angle info in gulp format
for i in $(seq 1 ${#vlist[*]})
do
	vals=($(echo ${vlist[i-1]} | sed 's/-->/ /g'))
	if [ "${vals[0]}" == "harmonic" ]; then
		echo "three bond regular kcal"
		vscale=(2 1)
		print_val_coeffs
	fi
done

#Dihedrals
vt="dihedral"
n=4
get_val_coeffs
#now print angle info in gulp format
for i in $(seq 1 ${#vlist[*]})
do
	vals=($(echo ${vlist[i-1]} | sed 's/-->/ /g'))
	if [ "${vals[0]}" == "harmonic" ]; then
		echo "torsion bond intra kcal"
		tmp=$(echo ${vals[7]} | awk '{print $1+0}')
		vals[7]=0
		if [ $(echo ${vals[6]} | awk '{print $1+0}') -eq -1 ]; then
			vals[7]=180
		fi	
		vals[6]=$tmp
		vscale=(1 1 1)
		print_val_coeffs
	fi
done
