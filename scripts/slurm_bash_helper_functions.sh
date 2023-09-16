#!/bin/bash

get_rattle_str()
{
	unrattleStr=""
	rattleStr=""
	zero=0
	solv_shake=0
	solu_shake=0
	what_shake=""
	if [ -z "${solv_mol_constraints// }" ]; then
		solv_mol_constraints=0
	fi
	if [ -z "${solu_mol_constraints// }" ]; then
		solu_mol_constraints=0
	fi
	if [ ${hsolv} -eq 1 ] && [ ${solv_mol_constraints} -gt 0 ]; then
		solv_shake=1
	fi
	if [ ${hsolu} -eq 1 ] && [ ${solu_mol_constraints} -gt 0 ]; then
		solu_shake=1
	fi
	if [ ${solv_shake} -eq 1 ] && [ ${solu_shake} -eq 1 ]; then
		what_shake="all"
	elif [ ${solv_shake} -eq 1 ] && [ ${solu_shake} -eq 0 ]; then
		what_shake="solv"
		if [ ${hsolu} -eq 0 ]; then
			what_shake="all"
		fi
	elif [ ${solv_shake} -eq 0 ] && [ ${solu_shake} -eq 1 ]; then
		what_shake="solu"
		if [ ${hsolv} -eq 0 ]; then
			what_shake="all"
		fi
	fi

	if [ -n "${what_shake// }" ]; then	
	        unrattleStr="unfix rattleH"
		rattleAngle=""
		if [ ${solv_shake} -eq 1 ]; then
	      		rattleAngle=$(awk '{i=NF-3;j=NF-2;k=NF-1;if($NF ~ /^H/ && $k ~ /^O/ && $j ~ /^H/ && $i ~ /#/)print $1}' $dat_file | head -1)
		       	if [ -z "${rattleAngle// }" ]; then
		               	rattleAngle=$(awk '{i=NF-3;j=NF-2;k=NF-1;if($NF ~ /^H/ && $k ~ /^O/ && $j ~ /^H/ && $i ~ /#/)print $2}' $inp_file | head -1)
	       		fi
		fi
		if [ -n "${rattleAngle// }" ]; then
			rattleAngle="a ${rattleAngle}"
		fi	
	       	rattleStr="fix             rattleH ${what_shake} rattle 0.0001 20 500 m 1.008 ${rattleAngle}"
		zero=0
    	else
		zero=1
	fi
    	for i in ${inp_file}
    	do
	    echo "RATTLESTR UPDATING: $i"
	    sed -i "s/\"rattleStr_here\"/\"${rattleStr}\"/g" $i
	    sed -i "s/\${rattleStr}/${rattleStr}/g" $i
	    sed -i "s/\"unrattleStr_here\"/\"${unrattleStr}\"/g" $i
	    sed -i "s/\${unrattleStr}/${unrattleStr}/g" $i
	done
}

get_h_min_opt()
{
	#fix hydrogen vdw to prevent bad things during mininization for certain water models
	h_min_fix=""
	hh_min_fix=""
	h_min_unfix=""
	h_typeids=($(egrep -i '# (H|D|T)[0-9_A-Z]*\s*$' $dat_file | awk '{print $1}'))
	declare -A htypes
	for i in ${h_typeids[*]}
	do
		h_hyb=""
		h_eta=$(awk -v ht="^$i\$" 'BEGIN{valid=0}{if($1 ~ /^Pair/ && $2 ~ /^Coeff/)valid=1;else if(valid && $1 ~ ht) print $2; else if(valid && $2 ~ /^Coeff/) valid=0}' $dat_file | tail -1)
		h_sig=$(awk -v ht="^$i\$" 'BEGIN{valid=0}{if($1 ~ /^Pair/ && $2 ~ /^Coeff/)valid=1;else if(valid && $1 ~ ht) print $3; else if(valid && $2 ~ /^Coeff/) valid=0}' $dat_file | tail -1)
		if [ -z "${h_eta// }" ] || [ -z "${h_sig// }" ]; then
			#search control file
			h_eta=$(egrep "^\s*pair_coeff \s*${i} \s*${i} " ${inp_file} | awk '{if($4 ~ /[a-zA-Z]/ && $4 ~ /lj/)print $5; else if($4 !~ /[a-zA-Z]/)print $4}' | tail -1)
			h_sig=$(egrep "^\s*pair_coeff \s*${i} \s*${i} " ${inp_file} | awk '{if($4 ~ /[a-zA-Z]/ && $4 ~ /lj/)print $6; else if($4 !~ /[a-zA-Z]/) print $5}' | tail -1)
			h_hyb=$(egrep "^\s*pair_coeff \s*${i} \s*${i} " ${inp_file} | awk '{if($4 ~ /[a-zA-Z]/)print $4; }' | tail -1)
			echo "NOTE: Searched $inp_file for H VDW data and found ${h_hyb} type: $i eta: $h_eta sigma: $h_sig"
			if [ -z "${h_eta// }" ] || [ -z "${h_eta// }" ] && [ $(echo  $solv | egrep -ic '(mW|m3b|qeq)') -eq 0 ]; then
	           	echo "WARNING: Couldn't figure out the h_eta...Skipping"
			fi
		else
			echo "NOTE: Searched $dat_file and $inp_file for H VDW data and found type: $i eta: $h_eta sigma: $h_sig"
		fi
		h_eta_test=$(echo $h_eta | awk '{if($1 == 0) print 1; else print 0}')
		h_sig_test=$(echo $h_sig | awk '{if($1 == 0) print 1; else print 0}')

		echo "TESTING: h_typeid: ${i} h_eta_test: ${h_eta_test} solv: ${solv}"
		if [ -n "${i// }" ] && [ $h_eta_test -eq 1 ] && [ $(echo  $solv | egrep -ic '(mW|m3b|pqeq|rexpon)') -eq 0 ]; then
			echo "NOTE: fixing hydrogen atom type $i to prevent overlap during minimization..."
			echo
			if [ -z "${h_hyb// }" ]; then
				htypes["none"]=${htypes["none"]}"${i},"
			else
				htypes["${h_hyb}"]=${htypes["${h_hyb}"]}"${i},"
			fi
		fi
	done
	for i in ${!htypes[@]}
	do
		echo "hybrid: $i"
		hstr=""
		if [ $i != "none" ]; then
			hstr=$i
		fi
		istr=$(echo ${htypes[$i]}| sed 's/\,$//')
		indx=($(echo $istr | sed 's/\,/ /g'))
		for j in ${indx[*]}
		do
			h_min_fix="${h_min_fix} \npair_coeff $j * $hstr 0.1 3.000 "
			hh_min_fix="${hh_min_fix} \npair_coeff $j $j $hstr 0.0 3.000 "
			h_min_unfix="${h_min_unfix} \npair_coeff $j * $hstr 0.0 3.000 "
		done
	done
	echo "nhybd: $nhybd h_min_fix: $h_min_fix"
	for i in ${inp_file}
	do
	   	echo "H_MIN UPDATING: $i"
	   	sed -i "s;hh_min_fix_here;${hh_min_fix};g" $i
	   	sed -i "s;\${hh_min_fix};${hh_min_fix};g" $i
	   	sed -i "s;h_min_fix_here;${h_min_fix};g" $i
	   	sed -i "s;\${h_min_fix};${h_min_fix};g" $i
	   	sed -i "s;h_min_unfix_here;${h_min_unfix};g" $i
	   	sed -i "s;\${h_min_unfix};${h_min_unfix};g" $i
	done
}

get_hoh_min_opt()
{
    #make the oh bond and hoh angles stiff for minimization to avoid numerical errors
    oh_bond_str=""
    oh_bond_str_min=""
    hoh_angle_str=""
    hoh_angle_str_min=""
    #first get the bond info
    #search datafile
    oh_bond_str="$(awk '/# [H\|D\|T][_A-Z]\s*O[_A-Z]$/' ${dat_file} | tail -1)"
    if [ -n "${oh_bond_str// }" ]; then
        oh_bond_str="bond_coeff ${oh_bond_str}"
        if [ $(egrep -c 'bond_style\s*class2' ${inp_file}) -eq 1 ]; then
            oh_bond_str_min="bond_coeff $(echo $oh_bond_str | awk '{for(i=3;i<6;i++) $i*=100;print}')"
        else
            oh_bond_str_min="bond_coeff $(echo $oh_bond_str | awk '{$2=10000;print}')"
        fi
    else
    #search inputfile
        oh_bond_str="$(awk '/# [H\|D\|T][_A-Z]\s*O[_A-Z]$/' ${inp_file} | tail -1)"
        if [ -n "${oh_bond_str// }" ]; then
            oh_bond_str_min="$(echo $oh_bond_str | awk '{if($3 ~ /class/) {for(i=5;i<8;i++) $i*=100; } else if($3 ~ /[A-Za-z]/) $4=10000; else $3 = 10000; print}')"
        fi
    fi
    if [ -n "${oh_bond_str_min// }" ]; then
        #now get angle info
        #first search datafile
        hoh_angle_str="$(awk '/# [H\|D\|T][_A-Z]\s*O[_A-Z]\s*[H\|D\|T][_A-Z]$/' ${dat_file} | tail -1)"
        if [ -n "${hoh_angle_str// }" ]; then
            hoh_angle_str="angle_coeff ${hoh_angle_str}"
            hoh_angle_str_min="angle_coeff $(echo $hoh_angle_str | awk '{$2=1000; print}')"
        else
            #search inputfile
            hoh_angle_str="$(awk '/# [H\|D\|T][_A-Z]\s*O[_A-Z]\s*[H\|D\|T][_A-Z]$/' ${inp_file} | tail -1)"
            if [ -n "${hoh_angle_str// }" ]; then
                hoh_angle_str_min="$(echo $hoh_angle_str | awk '{if($3 ~ /[A-Za-z]/) $4=1000; else $3 = 1000; print}')"
            fi
        fi
    fi
    oh_bond_str="$(echo $oh_bond_str | sed 's/#.*//')"
    oh_bond_str_min="$(echo $oh_bond_str_min | sed 's/#.*//')"
    hoh_angle_str="$(echo $hoh_angle_str | sed 's/#.*//')"
    hoh_angle_str_min="$(echo $hoh_angle_str_min | sed 's/#.*//')"
}

get_2pt_mem()
{
	#get the memory on the node and set the 2PT memory
	mem=($(scontrol show job ${SLURM_JOB_ID} | grep mem | sed 's/.*mem=\([0-9\.]*\)\([a-zA-Z]*\),.*/\1 \2/'))
	if [ ${#mem[@]} -eq 0 ]; then
		mem=($(scontrol show job ${SLURM_JOB_ID} | grep MinMemoryNode | sed 's/.*MinMemoryNode=\([0-9\.]*\)\([a-zA-Z]*\) .*/\1 \2/'))
	fi
	fact=1000
	if [ ${#mem[@]} -eq 0 ]; then
		echo "WARNING: Cannot determine memory assigned to job. Setting to 1Gbyte"
		mem=(1024)
	elif [ ${mem[1]} == "M" ]; then
	    fact=1
	elif [ ${mem[1]} == "T" ]; then
	    fact=1000000
	fi
	export twoPT_mem=$(echo ${mem[0]} $fact | awk '{print $1*$2}')
}
