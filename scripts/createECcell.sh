#!/bin/bash

abort()                                                                                       
{                                                                                             
    echo >&2 '                                                                                
***************                                                                               
*** ABORTED ***                                                                               
***************                                                                               
'                                                                                             
    echo "An error occurred. Exiting..." >&2                                                  
    exit 1                                                                                    
} 
create_electrode()
{
	#assume that it's an element so look for bgf file
	is_valid=0
	ele_nm=$(echo $element | tr '[:upper:]' '[:lower:]')
	for i in ${scripts_dir}/dat/elements_bgfs/*.bgf
	do
		ele=$(basename $i | sed 's/\-.*//'| tr '[:upper:]' '[:lower:]')
		if [ $ele == $ele_nm ]; then
			${scripts_dir}/autoType.pl -i ${i} -s ${element}.xtal.bgf -f UFF > /dev/null
			echo "		Will use ${i}"
			is_valid=1
			break
		fi
	done
	if [ $is_valid -eq 0 ]; then
		echo "ERROR: Cannot locate any valid element structure file for $element"
		exit 1
	fi
	ele_bgf=${element}.xtal.bgf
	${scripts_dir}/replicate.pl -b $ele_bgf -u 1 -s $ele_bgf -d "1 1 4" > /dev/null
	eltrde_cell=($(grep -E CRYSTX $ele_bgf | awk '{print $2,$3,$4}'))
	etrde=($(echo ${eltrde_cell[*]} ${electrode_cell[*]} | awk '{print int($4/$1)+1,int($5/$2)+1,1}'))
	if [ ${etrde[0]} -gt 1 ] || [ ${etrde[1]} -gt 1 ]; then
		echo "		replicating $element xtal by ${etrde[*]}"
		${scripts_dir}/replicate.pl -b $ele_bgf -u 1 -s $ele_bgf -d "${etrde[*]}" > /dev/null
	fi
	${scripts_dir}/bondByDistance.pl -b $ele_bgf -s $ele_bgf -f UFF > /dev/null
	${scripts_dir}/decreaseDim.pl -b  $ele_bgf -s $ele_bgf -d z > /dev/null
}
scripts_dir=$(dirname $0)


if [ $# -lt 2 ]; then
	echo "usage: $0 electrode_element/electrode_structure_file solvent_file [save_name] [number_of_electrodes=1|2(default)] [center_in_cell=yes] [electrolyte thickness=60A(default)]"
	exit 1
fi

ele_bgf=$1

solv_bgf=$2
if [ ! -e $solv_bgf ]; then
	echo "ERROR: Cannot locate $solv_bgf"
	exit 1
fi

save_name=$(basename $ele_bgf)
save_name=${save_name%.*}
save_name="${save_name}.$(basename $solv_bgf)"
save_name=${save_name%.*}
save_name="${save_name}.ecCell.bgf"
if [ $# -gt 2 ]; then
	save_name=$3
fi

nelectrde=2
if [ $# -gt 3 ] && [ $(echo $4 | egrep -c '^1$') -eq 1 ]; then
	nelectrde=1
fi
echo "$nelectrde electrodes"
center=1
if [ $# -gt 4 ] && [ $(echo $5 | egrep -ic '0|no') -gt 0 ]; then
	center=0
fi
if [ ${center} -eq 1 ] && [ ${nelectrde} -eq 1 ]; then
	echo "Will center 1 electrode setup in unit cell"
elif [ ${nelectrde} -eq 1 ]; then 
	echo "Will not center 1 electrode setup in unit cell"
fi	
if [ $# -gt 5 ] && [ $(echo $6 | egrep -c '^[0-9]*$') -gt 0 ]; then
	dz=$6
else
	dz=60
fi
echo "solvent/electrolyte thickness: $dz A"

if [ $# -gt 5 ] && [ $(echo $6 | egrep -c '^[0-9]*$') -gt 0 ]; then
	nsolv=$6
else	
	nsolv=516
fi
#echo "${nsolv} solvent/electrolyte molecules" #turn off this option for now (not sure it makes sense...)

trap 'abort' 0                                                                                
set -e                                                                                        
STARTTIME=$(date +%s)
echo "Started at $(date)"
if [ $(egrep -c CRYSTX $solv_bgf) -eq 0 ]; then
	${scripts_dir}/addBoxToBGF.pl $solv_bgf $solv_bgf
fi
solv_cell=("$(grep CRYSTX $solv_bgf | awk '{print $2,$3,$4}')")
cp $solv_bgf ./__solv.bgf
solv_bgf="__solv.bgf"

eltrde_list=($1)
#let's do several cases in order to get the correct cell size
if [ ${#eltrde_list[*]} -gt 1 ] && [ -e "${eltrde_list[0]}" ] && [ -e "${eltrde_list[1]}" ]; then
#case 1: both electrode files specified, so get largest xy dimensions
	electrode_cell=("$(grep CRYSTX "${eltrde_list[*]}" | awk '
{
	x=$2;y=$4;z=$5; 
	i=8; 
	while(i<NF){
		j=i+1; k=i+2; 
		if(x<$i) x=$i; 
		if(y<$j) y=$j; 
		if(z<$k) z=$k; 
		i+=7
	} 
	print x,y,z
}')")
elif [ ${#eltrde_list[*]} -gt 1 ] && [ -e "${eltrde_list[1]}" ]; then
#case 2: electrode #2 is a bgf file
	electrode_cell=($(grep CRYSTX ${eltrde_list[1]} | awk '{print $2,$3,$4}'))
	echo "case 3: $(grep CRYSTX ${eltrde_list[1]} | awk '{print $2,$3,$4}')"
elif [ -e "${eltrde_list[0]}" ]; then
#case 3: electrode #1 is a bgf file
	echo "case 3: $(grep CRYSTX ${eltrde_list[0]} | awk '{print $2,$3,$4}')"
	electrode_cell=($(grep CRYSTX ${eltrde_list[0]} | awk '{print $2,$3,$4}'))
else
#case 4: get from solvent file
	echo "case 4: $(grep CRYSTX ${solv_bgf} | awk '{print $2,$3,$4}')"
	electrode_cell=($(grep CRYSTX ${solv_bgf} | awk '{print $2,$3,$4}'))		
fi
			
echo "electrode_cell: ${electrode_cell[*]} ec1: ${electrode_cell[0]} ec2: ${electrode_cell[1]} ec3: ${electrode_cell[2]}"
echo "solv_cell: ${solv_cell[*]}"
echo "Will create $save_name"

ele_bgf=${eltrde_list[0]}
echo "	electrode 1:"
if [ ! -e $ele_bgf ]; then
	element=$ele_bgf
	create_electrode
else 
	echo "		will use $ele_bgf"
fi
cp $ele_bgf __electrode1.bgf
ele_bgf="__electrode1.bgf"
${scripts_dir}/centerBGF.pl -b __electrode1.bgf -s __electrode1.bgf -c com_center > /dev/null	
#set electrode 1 zmin to 0
szmin=$(${scripts_dir}/getBounds.pl -b __electrode1.bgf | grep '^Z ' | awk '{$2*=-1; if($2>0) printf "+%s",$2; else print $2}')
${scripts_dir}/modifyAtomData.pl -s __electrode1.bgf -w __electrode1.bgf -a "index>0" -f "ZCOORD:${szmin}" > /dev/null

del=$(echo $dz 2.0 | awk '{print $1+$2}')
zmin=$(${scripts_dir}/getBounds.pl -b __electrode1.bgf | grep '^Z ' | awk '{$2*=-1; if($2>=0) printf "+%s",$2; else printf "%s",$2}')
${scripts_dir}/modifyAtomData.pl -s __electrode1.bgf -w __electrode1.bgf -a "index>0" -f "CHAIN:A ZCOORD:$zmin RESNAME:EL1" > /dev/null
off=$(${scripts_dir}/getBounds.pl -b __electrode1.bgf | grep '^Z ' | awk '{print $3+2.6}')
if [ $nelectrde -eq 2 ]; then
	echo "	electrode 2:"
	if [ ${#eltrde_list[@]} -gt 1 ]; then
		ele_bgf=${eltrde_list[1]}
		if [ ! -e $ele_bgf ]; then
			element=$ele_bgf
			create_electrode
		else
			echo "		Will use $ele_bgf"
		fi		
		cp $ele_bgf __electrode2.bgf
		ele_bgf="__electrode2.bgf"
	else
		echo "		Will use x-rotated $ele_bgf"
		${scripts_dir}/rotatemol.pl -b __electrode1.bgf -r "x:180" -s __electrode2.bgf -a "index>0" > /dev/null
	fi
	${scripts_dir}/centerBGF.pl -b __electrode2.bgf -s __electrode2.bgf -c com_center > /dev/null
	#set electrode 2 zmin to zero
	szmin=$(${scripts_dir}/getBounds.pl -b __electrode2.bgf | grep '^Z ' | awk '{$2*=-1; if($2>0) printf "+%s",$2; else print $2}')
	${scripts_dir}/modifyAtomData.pl -s __electrode2.bgf -w __electrode2.bgf -a "index>0" -f "ZCOORD:${szmin}" > /dev/null
	#offset electrode 2	
	eltrde2_zoffset=$(${scripts_dir}/getBounds.pl -b __electrode1.bgf | grep '^Z ' | awk -v dz=$del '{print $3+dz+5.2}')
	echo "		offsetting by ${eltrde2_zoffset}"
	${scripts_dir}/modifyAtomData.pl -s __electrode2.bgf -w __electrode2.bgf -a "index>0" -f "CHAIN:B ZCOORD:+${eltrde2_zoffset} RESNAME:EL2" > /dev/null 	
	${scripts_dir}/combineBGF.pl __electrode1.bgf __electrode2.bgf __electrode_final.bgf > /dev/null
	z1=$(${scripts_dir}/getBounds.pl -b __electrode1.bgf | grep '^Z' | awk '{print $3-$2}')
	z2=$(${scripts_dir}/getBounds.pl -b __electrode2.bgf | grep '^Z' | awk '{print $3-$2}')
	ztot=$(echo $z1 $z2 $dz | awk '{print $1+$2+$3}')
else
	ztot=$(${scripts_dir}/getBounds.pl -b __electrode1.bgf | grep '^Z ' | awk -v dz=$del '{print $3+dz}')
	mv __electrode1.bgf __electrode_final.bgf 
fi

echo "	solvent/electrolyte"
#now let's create the correct solvent cell
solv_cell_str=($(echo ${electrode_cell[*]} $dz ${solv_cell[*]} | awk '{print int($1/$5)+1,int($2/$6)+1,int($4/$7)+1}'))
#create larger solv cell if possible
if [ ${solv_cell_str[0]} -gt 1 ] || [ ${solv_cell_str[1]} -gt 1 ] || [ ${solv_cell_str[2]} -gt 1 ]; then
	nmol_solv=$(${scripts_dir}/countAtoms.pl -f $solv_bgf -m 1 | grep '^Found' | grep mole | awk '{print $2}')
	if [ $nmol_solv -gt 1 ]; then
		echo "	replicating solvent cell by ${solv_cell_str[*]}"
		${scripts_dir}/replicate.pl -b $solv_bgf -s $solv_bgf -d "${solv_cell_str[*]}" > /dev/null
	else
		echo "	1 molecule in solvent bgf, so creating amorphous structure will cell '"${electrode_cell[0]} ${electrode_cell[1]} $dz"'"
		nmol_solv=$(grep -E CRYSTX $solv_bgf | tail -1 | awk -v cell_str="${electrode_cell[0]} ${electrode_cell[1]} $dz" '
BEGIN{
	split(cell_str,cell," ")
}
{
	mol_vol=$2*$3*$4 ;
	cell_vol = cell[1]*cell[2]*cell[3];
	print int(cell_vol/mol_vol)+1
}');
		echo "		will add $nmol_solv molecules"
		#echo ${scripts_dir}/amorphousBuilder.pl -c "${electrode_cell[0]} ${electrode_cell[1]} $dz" -m $solv_bgf -n $nmol_solv	
		${scripts_dir}/amorphousBuilder.pl -c "${electrode_cell[0]} ${electrode_cell[1]} $dz" -m $solv_bgf -n $nmol_solv > /dev/null	
	fi
fi
echo "		trimming cell"
${scripts_dir}/centerBGF.pl -b __solv.bgf -s __solv.bgf -c com_center > /dev/null
${scripts_dir}/trimCell.pl -m 1 -b __solv.bgf -s __solv.bgf -c "${electrode_cell[0]} ${electrode_cell[1]} $dz" > /dev/null 
cSolv=$(${scripts_dir}/countAtoms.pl -f __solv.bgf -m 1  | grep '^Found ' | grep mole | awk '{print $2}')
if [ $cSolv -lt $nsolv ] && [ 1 -gt 2 ]; then
	echo "	adjusting number of solvent/electrolyte molecules to $nsolv"
	repz=$(echo $nsolv $cSolv | awk '{print int($1/$2)+1}')
	${scripts_dir}/replicate.pl -b __solv.bgf -s __solv.bgf -u 1 -d "1 1 $repz" > /dev/null
	cSolv=$(${scripts_dir}/countAtoms.pl -f __solv.bgf -m 1  | grep '^Found ' | grep mole | awk '{print $2}')
	zmax=$(echo $nsolv $cSolv $repz ${solv_cell[3]} | awk '{print $1*$4*$3/$2+3}')
	${scripts_dir}/trimCell.pl -m 1 -b __solv.bgf -s __solv.bgf -c "${electrode_cell[0]} ${electrode_cell[1]} $zmax" > /dev/null 
	cSolv=$(${scripts_dir}/countAtoms.pl -f __solv.bgf -m 1  | grep '^Found ' | grep mole | awk '{print $2}')
	nremove=$(echo $cSolv $nsolv | awk '{print $1-$2}')
	if [ $nremove -gt 1 ]; then
		echo "		removing $nremove molecules to achieve $nsolv"
		${scripts_dir}/removeAtoms.pl -a "index>0" -f __solv.bgf -s __solv.bgf -m 1 -r 1 -n $nremove > /dev/null
	fi
fi
zlen=$(echo ${solv_cell[3]} $dz | awk '{print $1+$2}')
#set the zmin of solvent to 0
szmin=$(${scripts_dir}/getBounds.pl -b __solv.bgf | grep '^Z ' | awk '{$2*=-1; if($2>0) printf "+%s",$2; else print $2}')
${scripts_dir}/modifyAtomData.pl -s __solv.bgf -w __solv.bgf -a "index>0" -f "CHAIN:X ZCOORD:${szmin}" > /dev/null
echo  "	offseting solvent by $off"
${scripts_dir}/modifyAtomData.pl -s __solv.bgf -w __solv.bgf -a "index>0" -f "ZCOORD:+${off}" > /dev/null

${scripts_dir}/combineBGF.pl __electrode_final.bgf __solv.bgf $save_name > /dev/null
if [ $nelectrde -eq 1 ]; then
	zlen=$(${scripts_dir}/getBounds.pl -b $save_name | grep '^Z ' | awk '{print $3-$2+2}')
else
	zlen=$(${scripts_dir}/getBounds.pl -b $save_name | grep '^Z ' | awk '{print $3-$2+52}')
fi

electrode_cell=(${electrode_cell[0]} ${electrode_cell[1]} $zlen)
${scripts_dir}/updateBGFBox.pl -b $save_name -s $save_name -c "${electrode_cell[*]} 90 90 90" > /dev/null
${scripts_dir}/splitAtomsByMol.pl -b $save_name -s $save_name -f resnum > /dev/null
echo "	centering atoms in cell"
if [ $nelectrde -eq 1 ] && [ ${center} -eq 1 ]; then
	${scripts_dir}/centerBGF.pl -b $save_name -s __test.bgf -c com_center -r "chain eq 'A'"  > /dev/null
	${scripts_dir}/imageAtomsIntoBox.pl -b __test.bgf -s $save_name -m 1 -o "chain ne 'A'" > /dev/null
else
	${scripts_dir}/centerBGF.pl -b $save_name -s __test.bgf -c com_center > /dev/null
fi
${scripts_dir}/centerBGF.pl -b $save_name -c com_center -d z -s $save_name > /dev/null
${scripts_dir}/splitAtomsByMol.pl -b $save_name -s $save_name -f resnum > /dev/null
${scripts_dir}/groupAtoms.pl -b $save_name -s $save_name -f chain -g resnum > /dev/null
rm -fr __electrode*.bgf __solv.bgf __test.bgf
echo "	final cell size: ${electrode_cell[*]}"
${scripts_dir}/bgfcharge.pl $save_name | awk '{print $1,$2,$3}'
echo "All tasks completed"

trap : 0
ENDTIME=$(date +%s)
echo "Ended at $(date)"
echo $STARTTIME $ENDTIME | awk '{diff=$2-$1; days=int(diff/86400); diff-= days*86400; hrs=int(diff/3600); diff-=hrs*3600; mins=int(diff/60); diff-=mins*60; printf "Took %d days %d hrs %d mins %d sec to complete this task..\n",days,hrs,mins,diff}'
