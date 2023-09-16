#!/bin/bash
#Job name:
#SBATCH --job-name=fprefix_here.rtemp_hereK.cp2k
#
#Partition:
#SBATCH --partition=vulcan_c24
#
#Account
#SBATCH --account=vulcan
#
#QOS
#SBATCH --qos=normal
#
#Nodes
#SBATCH --nodes=5
#
#Processors:
#SBATCH --ntasks=24
#
#Wall clock limit:
#SBATCH --time=72:0:0
#

fprefix=fprefix_here
rtemp=rtemp_here
cp2k_exe=/global/home/users/tpascal/codes/bin/cp2k.popt

curr_dir=$SLURM_SUBMIT_DIR
results_dir=$curr_dir/results
temp_dir=$SCRATCH/md/cp2k/${fprefix}/${rtemp}K
cd $curr_dir

mkdir -p ${temp_dir} ${results_dir}/

cd ${temp_dir}
cp ${curr_dir}/${fprefix}.${rtemp}K.cp2k.in ./
cp ${curr_dir}/${fprefix}.cp2k.xyz ./

bash /usr/Modules/init/bash
module unload intel
module load intel/2015.6.233 openmpi fftw mkl

mpirun -mca btl self,sm,openib $cp2k_exe -in ${fprefix}.${rtemp}K.cp2k.in > ${fprefix}.${rtemp}K.screen.out
cp ${fprefix}.* $results_dir/
