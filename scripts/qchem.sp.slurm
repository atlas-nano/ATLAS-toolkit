#!/bin/bash
#Job name:
#SBATCH --job-name=co2.sp.B3LYP_vdwDF2_aug-cc-pVTZ
#
#Partition:
#SBATCH --partition=vulcan
#
#Account
#SBATCH --account=nano
#
#QOS
#SBATCH --qos=normal
#
#Nodes
#SBATCH --nodes=1
#
#Processors:
#SBATCH --ntasks=8
#
#QUEUE
#SBATCH --constraint=vulcan
#
#Wall clock limit:
#SBATCH --time=10:0:0
#

module unload intel
module load qchem/5.1.2
temp_dir=/clusterfs/vulcan/pscratch/tpascal/sp/qchem/carbonates
temp_dir=/global/scratch/tpascal/sp/qchem/carbonates
mkdir -p $temp_dir
export QCSCRATCH=$temp_dir
export QCTMPDIR=$temp_dir
echo $SLURM_NODELIST > ${temp_dir}/${SLURM_JOB_ID}.machinefile
export QCMACHINEFILE=${temp_dir}/${SLURM_JOB_ID}.machinefile
export QCMPIRUN=mpirun
nprocs=8
NN=8

mkdir -p $SLURM_SUBMIT_DIR/results $SLURM_SUBMIT_DIR/running $SLURM_SUBMIT_DIR/completed $SLURM_SUBMIT_DIR/failed
cd $temp_dir
for i in $SLURM_SUBMIT_DIR/*.qcin
do
	if [ ! -e $i ]; then
		continue
	fi
	cp $i ./
	mv $i $SLURM_SUBMIT_DIR/running/
	prefix=$(basename $i .qcin)
	scontrol update jobid=$SLURM_JOB_ID name=${prefix}
	sed  -i '/^\s*CHELPG/a MEM_TOTAL      5000' ${prefix}.qcin
	qchem -np 8 -save ${prefix}.qcin ${prefix}.qcout ${prefix}
	if [ $(egrep -c 'Total energy in the final basis set' ${prefix}.qcout) -eq 0 ]; then
		mv $SLURM_SUBMIT_DIR/running/${prefix}.qcin $SLURM_SUBMIT_DIR/failed
	else	
		cp ${prefix}.qcout $SLURM_SUBMIT_DIR/results/
		mv $SLURM_SUBMIT_DIR/running/${prefix}.qcin $SLURM_SUBMIT_DIR/completed
	fi	
done
