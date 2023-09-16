#!/bin/bash

if [ $# -lt 2 ]; then
	echo "usage: $0 lammps_data_file lammps_trj_file [dq=0.005 (charge bin)] [save_prefix]"
	exit 1
fi

data_file=$1
trj_file=$2
dq=0.005
if [ $# -gt 2 ]; then
	dq=$3
fi
save_prefix=`basename $trj_file`
save_prefix=${save_prefix%.*}
if [ $# -gt 3 ]; then
	save_prefix=$4
fi

atoms_line=(`head -20 $trj_file | egrep '^ITEM: ATOMS'| head -1`)
i=0
q_col=0
for i in `seq 1 ${#atoms_line[*]}`
do
	j=$((i-1))
	if [ "${atoms_line[$j]}" == "q" ]; then
		q_col=$(($j-1))
	fi
done
if [ $q_col -eq 0 ]; then
	echo "ERROR: Cannot find field ' q ' in $trj_file"
	exit 1
fi
ml=$(($q_col-1))

type_str=(`awk 'BEGIN{start=0}{if($1 ~ /^Masses/){start=1;}else if(start&&$1 ~ /^[a-zA-Z]/){start=0;}else if(start&&NF>1){if($3 ~ /#/) print $NF; else print $1;}}' $data_file`)
cat <<DATA
trj_file:       $trj_file
data_file:      $data_file
dq:             $dq
prefix:         $save_prefix
q_col:          $q_col
min_line:       $ml
atom_types:     ${type_str[*]}
DATA
scripts_dir=`dirname $0`

gawk -v dq=$dq -v prefix=$save_prefix -v q=$q_col -v ml=$ml -v type_str="${type_str[*]}" -f ${scripts_dir}/get_q_distrib.awk $trj_file
