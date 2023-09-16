abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    date
    exit 1
}

trap 'abort' 0

set -e

func=`basename $INCAR | sed 's/INCAR.//'`
SCRATCH_DIR=$SCRATCH
VSCRATCH=/clusterfs/vulcan/pscratch/$USER/
ESCRATCH=/clusterfs/etna/pscratch/$USER/
if [ $( echo $SLURM_JOB_PARTITION | egrep -ic '(vulcan|nano)') -gt 0 ]; then
	SCRATCH_DIR=$VSCRATCH
elif [ $( echo $SLURM_JOB_PARTITION | egrep -ic 'etna') -gt 0 ]; then
	SCRATCH_DIR=$ESCRATCH
fi

temp_dir=${SCRATCH_DIR}/vasp/sp/${prefix}/${func}/
results_dir=$curr_dir/results/${func}/
comp_dir=$curr_dir/completed/${func}/
poscar_dir=$curr_dir/POSCARS
vasp_pseudo_dir=~/codes/vasp/Pseudo/paw-pbe
vasp_pseudo_dir=/global/home/groups-sw/nano/software/sl-7.x86_64/vasp.5.4.4_vtst/Potentials/PBE/
vasp_exec=/global/home/groups-sw/nano/software/sl-7.x86_64/vasp.5.4.4_vtst/bin/vasp_std

bash /usr/Modules/init/bash
module unload intel
module unload intel
#module load intel/2013_sp1.4.211 openmpi fftw mkl
module load intel/2016.4.072 openmpi fftw mkl

mkdir -p $results_dir $comp_dir

STARTTIME=$(date +%s)
cat <<DATA
VASP singleponit calculation on $prefix using $func functional
started at `date`
curr_dir: $curr_dir
temp_dir: $temp_dir
results_dir: $results_dir
DATA

for i in ${poscar_dir}/*.POSCAR
do
	mol=`basename $i .POSCAR`
	echo "VASP calc of $mol"
	if [ -e ${results_dir}/${mol}.OSZICAR ]; then
		echo "	Already completed... Skipping"
		continue
	elif [ -e ${comp_dir}/${mol}.POSCAR ]; then
		echo "	Currently running... Skipping"
		continue
	fi
	rm -fr $temp_dir/${mol}
	mkdir -p ${temp_dir}/${mol}
	cd ${temp_dir}/${mol}
	ln -s $i ./POSCAR
	pstr="`head -1 POSCAR | awk -v dname="$vasp_pseudo_dir" '{for(i=1;i<=NF;i++) printf " %s/%s/POTCAR ",dname,$i}'`"
	cat $pstr > POTCAR
	ln -s ${i} ${comp_dir}/${mol}.POSCAR
        ln -s ~/codes/bin/vdw_kernel.bindat
        ln -s $curr_dir/KPOINTS
        ln -s $curr_dir/${INCAR} ./INCAR
        mpirun -mca btl self,sm,openib $vasp_exec  > ${mol}.OUT
	~/scripts/vasp/chgsum.pl AECCAR0 AECCAR2 > /dev/null
	~/codes/bin/bader -c bader -i chgcar CHGCAR -ref CHGCAR_sum > ${mol}.bader.out
	mv ACF.dat ${mol}.bader.atomic.charges.dat
	mv BCF.dat ${mol}.bader.volume.charges.dat
	cp *.charges.dat ${results_dir}
        #cp OUTCAR ${results_dir}/${mol}.OUTCAR
        cp OSZICAR ${results_dir}/${mol}.OSZICAR
done

trap : 0

ENDTIME=$(date +%s)
echo "Ended at "`date`
echo $STARTTIME $ENDTIME | awk '{diff=$2-$1; days=int(diff/86400); diff-= days*86400; hrs=int(diff/3600); diff-=hrs*3600; mins=int(diff/60); diff-=mins*60; printf "Took %d days %d hrs %d mins %d sec to complete this task..\n",days,hrs,mins,diff}'
