#!/bin/bash
#Job name:
#SBATCH --job-name=fprefix_here.rtemp_hereK.cp2k
#
#Partition:
#SBATCH --partition=regular
#
#Account
#SBATCH --account=m1141
#
#Nodes
#SBATCH --nodes=nodes_here
#
#Processors:
#SBATCH --ntasks=nprocs_here
#
#Wall clock limit:
#SBATCH --time=36:0:0
#
#
fprefix=fprefix_herecell_p_here
nprocs=nprocs_here
rtemp=rtemp_here

cd $SLURM_SUBMIT_DIR

results_dir=$SLURM_SUBMIT_DIR/results
mkdir -p $SCRATCH/md/cp2k/${fprefix} ${results_dir}/

cd $SCRATCH/md/cp2k/${fprefix}
cp $SLURM_SUBMIT_DIR/${fprefix}.${rtemp}K.cp2k.in ./
cp $SLURM_SUBMIT_DIR/${fprefix}.cp2k.xyz ./
module load cp2k

srun -n $nprocs cp2k.popt -in ${fprefix}.${rtemp}K.cp2k.in > ${fprefix}.${rtemp}K.out 
#cp ${fprefix}.* $results_dir/
