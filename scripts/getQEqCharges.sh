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

scripts_dir=$(dirname $0)
hdir=$(realpath ${scripts_dir}/../)
if [ $# -lt 2 ]; then
	echo "usage: $0 bgf_file 'ff_file(s)' [save_name] [total_charge=0] [qeq_data_file=${hdir}/ff/pqeq.dat] [lmp_executable=${hdir}/codes/bin/lmp_expanse]"
	exit 1
fi

bgf=$1
if [ ! -e $1 ] || [ ! -r $1 ]; then
	echo "ERROR: Cannot locate BGF file $1"
	exit 1
fi

ffs="$2"

if [ $# -lt 3 ]; then
	sname=$(basename $1)
	sname="${sname%.*}.qeq.bgf"
else
	sname=$3
fi

lmp=${hdir}/codes/bin/lmp_expanse.pqeq
if [ $# -gt 5 ]; then
	lmp=$6
fi
if [ ! -e $lmp ]; then
	echo "ERROR: Cannot locate lammps binary $lmp. Aborting"
	exit 1
fi

qeq_file=${hdir}/ff/pqeq.par
if [ $# -gt 4 ]; then
	qeq_file=$5
fi
if [ ! -e $qeq_file ]; then
	echo "ERROR: Cannot locate qeq data file: $qeq_file. Aborting"
	exit 1
fi

tot_q=0
if [ $# -gt 3 ]; then
	tot_q=$4
fi

cat <<DATA;
Will perform QEq charge calculation using the following:
bgf_file:	$bgf
ff_files:	$ffs
save_file:	$sname
lmp_exectuable:	$lmp
pqeq_data_file:	$qeq_file
total_charge:	$tot_q
DATA

trap 'abort' 0
set -e

#module load cpu slurm gcc openmpi cmake gsl intel-mkl &> /dev/null 2>&1
module load cpu/0.15.4  intel/19.1.1.217  gsl/2.5 &> /dev/null 2>&1
${scripts_dir}/createLammpsInput.pl -b $bgf -f "$ffs" -o "pqeq charge $tot_q" -s __test -q $qeq_file > /dev/null
${lmp} -in in.__test_singlepoint -screen none > /dev/null
${scripts_dir}/convertLammpsTrj.pl -b $bgf -l __test.lammps -t last -o bgf -s $sname > /dev/null
q=$(${scripts_dir}/bgfcharge.pl $sname | awk '{print $3}')
if [ $(echo $q $tot_q | awk '{if($1==$2)print 1; else print 0}') -eq 0 ]; then
	dq=$(echo $tot_q $q | awk '{dq=$1-$2; if(dq>0) printf "+%f",dq; else print dq}')
	${scripts_dir}/modifyAtomData.pl -s $sname -w $sname -a "index==1" -f "CHARGE:${dq}" > /dev/null
fi
rm -fr in.__test in.__test_singlepoint log.lammps log.cite __test.lammps __test.lammps.slurm data.__test
awk -v fname="__test.bgf" '
{
	if($1 ~ /^(ATOM|HETATM)/) {
		q=sprintf("%.3f",$NF); 
		split($0, a, " ", seps); 
		a[13]=sprintf("%.5f",q); 
		for (i=1;i<=NF;i++) 
			printf("%s%s", a[i], seps[i]) >> fname; 
		printf "\n" >> fname;
		tot+=$NF;
	}else 
		print >> fname;
}
END{
	print "sum of charges:",tot
}' $sname
mv __test.bgf $sname

trap : 0
echo "All done"
