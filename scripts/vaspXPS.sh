abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
	rm -fr ./running
    date
    exit 1
}

trap 'abort' 0

set -e

if [ $# -lt 2 ]; then
	echo "usage: $0 INCAR 'xps_atom(s)' (POSCAR_DIR)"
	exit 1
fi

INCAR=$1
xps_atom_list=($2)

if [ ! -e $INCAR ] || [ ! -r $INCAR ]; then
	echo "ERROR: Cannot access $INCAR"
	exit 1
fi

curr_dir=$SLURM_SUBMIT_DIR
func=`echo $INCAR | sed 's/INCAR.//'`

SCRATCH_DIR=$SCRATCH
if [ "${machine}" == "lrc" ]; then
	root_dir=/global/home/users/tpascal
	VSCRATCH=/clusterfs/vulcan/pscratch/$USER/
	ESCRATCH=/clusterfs/etna/pscratch/$USER/
	module unload intel
	module unload intel
	module load intel/2013_sp1.4.211 openmpi fftw mkl
	PARALLEL="mpirun -mca btl self,sm,openib"  
	vasp_cmd=${root_dir}/codes/bin/vasp
	vasp_pseudo_dir=${root_dir}/codes/vasp/Pseudo/paw-pbe
	if [ $( echo $SLURM_JOB_PARTITION | egrep -ic '(vulcan|nano)') -gt 0 ]; then
		SCRATCH_DIR=$VSCRATCH
	elif [ $( echo $SLURM_JOB_PARTITION | egrep -ic 'etna') -gt 0 ]; then
		SCRATCH_DIR=$ESCRATCH
	fi
elif [ "${machine}" == "nersc" ]; then
	root_dir=/global/homes/t/tpascal
	module unload intel
	module unload intel
	module load intel openmpi fftw vasp
	nprocs=`echo $SLURM_NNODES $SLURM_CPUS_ON_NODE | awk '{print $1*$2/2}'`
	PARALLEL="srun -n $nprocs -c2 --cpu_bind=cores"
	vasp_cmd=vasp_std
	vasp_pseudo_dir=${root_dir}/codes/vasp/Pseudo/paw-pbe
fi

scripts_dir=${root_dir}/scripts/
temp_dir=${SCRATCH_DIR}/vasp/xps/${prefix}/${func}/
results_dir=$curr_dir/results/${func}/
poscar_dir=$curr_dir/POSCARS
if [ $# -gt 2 ]; then
	poscar_dir=$3
fi

mkdir -p $results_dir

STARTTIME=$(date +%s)
cat <<DATA
VASP XPS calculation on $prefix using $func functional
started at `date`
root_dir:        $root_dir
curr_dir:        $curr_dir
temp_dir:        $temp_dir
scripts_dir:     $scripts_dir

results_dir:     $results_dir
poscar_dir:      $poscar_dir
vasp_pseudo_dir: $vasp_pseudo_dir

vasp_cmd:        $vasp_cmd
parallel_cmd:    $PARALLEL

xps_atom_list:   ${xps_atom_list[*]}

DATA

for i in ${poscar_dir}/*.POSCAR
do
	mol=`basename $i .POSCAR`
	echo "VASP GS of $mol"
	if [ -d ${temp_dir}/${mol} ] && [ -d ${temp_dir}/${mol}/GS ] && [ -e ${temp_dir}/${mol}/GS/running ]; then
		echo "	Appears to be currently running...Skipping to be safe"
	elif [ -e ${results_dir}/${mol}.GS.OUTCAR ]; then
		echo "	Appears to be already completed... Skipping"
	else
		mkdir -p ${temp_dir}/${mol}/GS
		cd ${temp_dir}/${mol}/GS
		echo "	curr_dir: $PWD"
		echo "1" > ./running
		cp $i ./POSCAR
	    ln -fs ~/codes/bin/vdw_kernel.bindat
		ln -fs $curr_dir/KPOINTS
	    cp $curr_dir/${INCAR} ./INCAR
		pstr=`head -1 POSCAR | awk -v dname="$vasp_pseudo_dir" '{for(i=1;i<=NF;i++) printf " %s/%s/POTCAR ",dname,$i}'`
		cat $pstr > POTCAR
		echo "		Started at "`date`
		$($PARALLEL $vasp_cmd > ${mol}.GS.OUT 2>&1) || true
		if [ `egrep -c '^   1 F=' ${mol}.GS.OUT ` -gt 0 ]; then
			cp ${mol}.GS.OUT ${results_dir}/
			cp OUTCAR ${results_dir}/${mol}.GS.OUTCAR
		else
			echo "		ERROR occurred..."
		fi
		rm -fr ./running
		echo "		Ended at "`date`
	fi
	for j in ${xps_atom_list[*]}
	do
		cd ${temp_dir}/${mol}
		xt=(`awk -v idx=$j -f ${scripts_dir}/vaspXPSposcar.awk $i`)
		xtype=${xt[0]}
		echo "	XPS of ${mol}.${xtype}${j}"
		if [ -d ${temp_dir}/${mol} ] && [ -d ${temp_dir}/${mol}/${xtype}${j}xps ] && [ -e ${temp_dir}/${mol}/${xtype}${j}xps/running ]; then
			echo "	Appears to be currently running...Skipping to be safe"
			continue
		elif [ -e ${results_dir}/${mol}.${xtype}${j}xps.OUTCAR ]; then
			echo "		Appears to be already completed... Skipping"
			continue
		fi
		mkdir -p ${temp_dir}/${mol}/${xtype}${j}xps
		cd ${temp_dir}/${mol}/${xtype}${j}xps
		xt=(`awk -v idx=$j -f ${scripts_dir}/vaspXPSposcar.awk $i`)
		echo "		curr_dir: $PWD"
		echo "1" > running
	    ln -fs ~/codes/bin/vdw_kernel.bindat
		ln -fs $curr_dir/KPOINTS
		cp $curr_dir/${INCAR} ./INCAR
		cat >> ./INCAR <<DATA
ICORELEVEL = 2
CLNT = ${xt[1]}
CLN = 1
CLL = 0
CLZ = 1
DATA
		pstr=`cat _xps_str | awk -v dname="$vasp_pseudo_dir" '{for(i=1;i<=NF;i++) printf " %s/%s/POTCAR ",dname,$i}'`
		cat $pstr > POTCAR
		cp $i _tmp
		sed -i '1d;6d' _tmp
		cat _xps_str _tmp > POSCAR
		sed -i '5r _xps_idx' POSCAR
		rm _tmp _xps_str _xps_idx
		echo "		Started at "`date`
		$($PARALLEL $vasp_cmd > ${mol}.${xtype}${j}xps.OUT 2>&1) || true
		if [ `egrep -c '^   1 F=' ${mol}.${xtype}${j}xps.OUT ` -gt 0 ]; then
			cp ${mol}.${xtype}${j}xps.OUT ${results_dir}/
			cp OUTCAR ${results_dir}/${mol}.${xtype}${j}xps.OUTCAR
		else
			echo "		ERROR occurred..."
		fi
		rm -fr ./running
		echo "		Ended at "`date`
	done
done

trap : 0

ENDTIME=$(date +%s)
echo "Ended at "`date`
echo $STARTTIME $ENDTIME | awk '{diff=$2-$1; days=int(diff/86400); diff-= days*86400; hrs=int(diff/3600); diff-=hrs*3600; mins=int(diff/60); diff-=mins*60; printf "Took %d days %d hrs %d mins %d sec to complete this task..\n",days,hrs,mins,diff}'
