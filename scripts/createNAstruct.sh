#!/bin/bash
abort()                                                                                       
{                                                                                             
    echo >&2 '                                                                                
***************                                                                               
*** ABORTED ***                                                                               
***************                                                                               
'                                                                                             
    echo "An error occurred. Exiting..." >&2                                                  
    if [ -e ${curr_dir}/running/$slurm_file ]; then                                           
        mv ${curr_dir}/running/$slurm_file  ${curr_dir}/failed                                
        touch ${curr_dir}/failed/$slurm_file
    fi                                                                                        
    exit 1                                                                                    
}                                                                                             
scripts_dir=$(dirname $0)

trap 'abort' 0

set -e

if [ $# -lt 2 ]; then
	echo "usage: $0 'sequence' helix_type:[s=single strand, d = double strand] [savename=dna.bgf] [na_type=dna(default) or rna] [ds_na_form=b(default) or a/s/z)] [extra_lib(s)]"
	exit 1
fi

nseq=$1

if [ $(echo $2 | egrep -v -c -i '^[d|s]$') -gt 0 ]; then
	echo "ERROR: Expected s or d for helix_type. Got '$2'"
	exit 1
fi
htype=$(echo $2 | tr '[:upper:]' '[:lower:]')

ntype="d"
if [ $# -gt 3 ]; then
	if [ $(echo $4 | egrep -c '^r') -gt 0 ]; then
		ntype="r"
	fi
fi


nform="b"
if [ $# -gt 4 ]; then
	tmp=$(echo $5 | sed 's/^.*\(a\|s\|z\).*$/\1/i')
	if [ ${#tmp} -eq 1 ]; then
		nform=$(echo $tmp | tr '[:upper:]' '[:lower:]')
	fi
fi

ht="double"
if [ "$htype" == "s" ]; then
	ht="single"
fi
if [ -z "${AMBERHOME// }" ]; then
	export AMBERHOME=/home/tpascal/codes/amber18
fi
leaprc=${AMBERHOME}/dat/leap/cmd/leaprc.DNA.bsc1
if [ "$ntype" == "r" ]; then
	leaprc=${AMBERHOME}/dat/leap/cmd/leaprc.RNA.OL3
fi
cp $leaprc ./my_leaprc

savename="${htype}s-${nform}${ntype}na.bgf"
if [ $# -gt 2 ]; then
	savename=$3
fi

x_libs=""
if [ $# -gt 5 ]; then
	for i in $6
	do
		if [ ! -e $i ]; then
			continue
		fi
		suf="${i##*.}"
		if [ "$suf" == "lib" ]; then
			echo "loadoff $i" >> ./my_leaprc
			x_libs="$x_libs $i"
		else
			echo "loadAmberParams $i" >> ./my_leaprc
		fi
	done
fi

fd_helix="abdna"
if [ "$nform" == "a" ] && [ "$ntype" == "r" ]; then
	fd_helix="arna"
elif [ "$nform" == "s" ] && [ "$ntype" == "d" ]; then
	fd_helix="sbdna"
elif [ "$nform" == "a" ] && [ "$ntype" == "d" ]; then
	fd_helix="adna"
fi

cat <<DATA
OPTIONS
sqeuence:    $nseq
helix_type:  $ht strand
NA_type:     ${ntype}na
NA_form:     ${nform}${ntype}na
fd_helix:    $fd_helix
bgf_name:    $savename

AMBERHOME:   $AMBERHOME

DATA

module load cpu/0.15.4  gcc/10.2.0  mvapich2/2.3.6 amber/20.21 netcdf-c/4.7.4 &> /dev/null 2>&1

declare -A atm_list

if [ $(echo $nseq | egrep -c '\([a-zA-Z0-9\.]*=[A|T|C|G|U]\)') -gt 0 ]; then
	#non-natural residues, so store resname and corresponding atom_names
	nnaa_list=($(echo $nseq | sed 's/[^(]*(\([a-zA-Z0-9\._]*\)=[^)]*)[^(]*/\1 /g' | awk '
{
	for(i=1;i<=NF;i++) 
		val[$i]=1
}
END{
	for(i in val)
		print i
}'))
	for i in ${nnaa_list[*]}
	do
		#search all the libs provided and select all the atom names
		found=0
		for j in $x_libs
		do
			start_line=$(egrep -n '^!!index array str' $j | sed 's/:.*//' | head -1)
			stop_line=$(egrep -n '!entry.\w+.unit.atoms' $j | sed 's/:.*//' | awk '{if(FNR==1) print $1-1}')
			d=$(( $stop_line - $start_line ))
			entryid=$(head -${stop_line} $j | tail -${d} | egrep -n "\"$i\"" | sed 's/:.*//' | head -1)
			if [ -n "${entryid// }" ]; then
				start_line=$(egrep -n "^!entry.${i}.unit.atoms" $j | sed 's/:.*//' | head -1)
				stop_line=$(egrep -n "^!entry.${i}.unit.atomspertinfo" $j  | sed 's/:.*//' | awk '{if(FNR==1) print $1-1}')
				d=$(( $stop_line - $start_line ))
				#found the library match, so get atom names
				tmp_str=$(head -${stop_line} $j | tail -${d} | sed 's/^\s*\"\([^\"]*\)\".*/\1/' | awk '
{
	val[$1]=1
}
END{
	for (i in val) 
		printf "|%s",i
}')
				atm_list["$i"]="^(1$tmp_str)\$"
				found=1
				break
			fi
		done
		if [ $found -eq 0 ]; then
			echo "ERROR: Cannot locate any resname $i in libs $x_libs"
			exit 1
		fi
	done
fi
echo "Step 1: Running NAB builder"
#now we have to search the sequence string and extra any nonstandard residues
#expected format is xxx(lib_name=[A|T|G|C|U])xxx
tseq=$(echo $nseq | sed 's/(\([a-zA-Z0-9\._]*\)=\([A|C|G|T|U]\))/||\1||\2/g') #test sequence
tseq=$(echo $tseq | sed "s/5'-/||5p||/g")
tseq=$(echo $tseq | sed "s/-3'/||3p||/g")
res_ids=()
res_nms=()
res_atm=()
threep=()
fivep=()
while [ $(echo $tseq | grep -aob '||' | head -1 | wc -l) -gt 0 ]
do
	pos=$(echo $tseq | grep -aob '||' | head -1 | sed 's/:.*//' | awk '{print $1+1}')
	tstr=$(echo $tseq | sed 's/^[^|]*||\(5\|3\)p||.*/\1/')
	if [ "$tstr" == "5" ]; then
		fivep=(${fivep[*]} $pos)
		tseq=$(echo $tseq | sed 's/||5p||//')
	elif [ "$tstr" == "3" ]; then
		tmp=$(( $pos - 1 ))
		threep=(${threep[*]} $tmp)
		tseq=$(echo $tseq | sed 's/||3p||//')
	else
		b=$(echo $tseq | sed 's/[^|]*||\([a-zA-Z0-9\.]*\)||.*/\1/')
		tseq=$(echo $tseq | sed 's/||[a-zA-Z0-9\.]*||\([A|C|G|T|U]\)/\1/')
		atm_str=${atm_list["$b"]}
		if [ -z "${atm_str// }" ]; then
			echo "ERROR: Cannot find library file for resname $b"
			exit 1
		fi
		res_ids=(${res_ids[*]} $pos)
		res_nms=(${res_nms[*]} $b)
		res_atm=(${res_atm[*]} $atm_str)
	fi
done
nseq=$tseq

cat > nuc.nab <<DATA
molecule m;
m = fd_helix( "$fd_helix", "$nseq", "${ntype}na" );
putpdb( "nuc.pdb", m, "-wwpdb");
DATA
path=($AMBERHOME/bin:$path)
nab nuc.nab
./a.out > /dev/null
if [ ${#fivep[*]} -gt 0 ] || [ ${#threep[*]} -gt 0 ]; then
	if [ "$htype" == "d" ]; then
		t=${#nseq}
		for i in ${threep[*]}
		do
			nres=$(( $i + $t ))
			threep=(${threep[*]} $nres)
		done
		for i in ${fivep[*]}
		do
			nres=$(( $i + $t ))
			fivep=(${fivep[*]} $nres)
		done
	fi
	awk -v fp_str="${fivep[*]}" -v tp_str="${threep[*]}" '
BEGIN{
	n=split (fp_str,tmp1," ")
	for(i=1;i<=n;i++)
		prime[tmp1[i]]=5
	n=split (tp_str,tmp2," ")
	for(i=1;i<=n;i++)
		prime[tmp2[i]]=3
}
{
	if($1 ~ /^(HETATM|ATOM)/) {
		resnum=$5
		if(resnum in prime)
			$4=$4""prime[resnum]
	}
	print
}' nuc.pdb > tmp.pdb
	mv tmp.pdb nuc.pdb
fi

if [ "$htype" == "s" ]; then
	pos=$(grep -n 'TER' nuc.pdb | head -1 | sed 's/:.*//')
	sed -i "${pos},\$ d" nuc.pdb
fi

echo "Step 2: Fixing pdb for AMBER types"
if [ ${#res_ids[*]} -gt 0 ]; then
	awk -v id_str="${res_ids[*]}" -v nm_str="${res_nms[*]}" -v at_str="${res_atm[*]}" '
BEGIN{
	split(id_str,ids," ");
	split(at_str,atm," ");
	n=split(nm_str,nms," ");
	for(i=1;i<=n;i++) {
		nmap[ids[i]]=nms[i]
		amap[ids[i]]=atm[i]
	}
}
{
	if($1 ~ /^(HETATM|ATOM)/) {
		res_nm=$5
		if(res_nm in nmap) {
			atm_patrn=amap[res_nm]
			if($3 ~ atm_patrn) {
				$4=nmap[res_nm]
				print
			}
		} else {
			print
		}
	} else {
		print
	}	
}' nuc.pdb > __tmp.pdb
	mv __tmp.pdb test.pdb
else
	cp nuc.pdb test.pdb
fi
sed -i '/^.*H\s*$/d' test.pdb #strip all hydrogens
awk '
{
	if($1 ~ /^(HETATM|ATOM)/) {
		printf "%-6s %4d %-4s %-4s %4s %11.3f%8.3f%8.3f%6.2f%6.2f%12s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
	} else {
		print
	}
}' test.pdb > tmp.pdb
mv tmp.pdb test.pdb
		
echo "Step 3: Assigning AMBER atom types"
${scripts_dir}/pdb2amber.pl -p test.pdb -s test -l my_leaprc> /dev/null

echo "Step 4: Creating BGF file"
${scripts_dir}/amber2bgf.pl test.prmtop test.inpcrd ${savename} > /dev/null

rm -fr test.pdb test.inpcrd test.prmtop nuc.pdb nuc.nab nuc.c a.out tleap.log test.out leap.log tleap.out my_leaprc leaprc
trap : 0

echo
echo "All Done! Creating AMBER-type BGF file ${savename}"
