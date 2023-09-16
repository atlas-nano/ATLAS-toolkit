#!/bin/bash
#Job name:
#SBATCH --job-name=fprefix_here
#
#Partition:
#SBATCH --partition=vulcan
#
#Account
#SBATCH --account=vulcan
#
#QOS
#SBATCH --qos=normal
#
#Nodes
#SBATCH --nodes=nodes_here
#
#Processors:
#SBATCH --ntasks=nprocs_here
#
#QUEUE
#SBATCH --constraint=vulcan
#
#Wall clock limit:
#SBATCH --time=72:0:0
#

module unload intel
module unload openmpi
module load intel/2015.0.090 fftw/3.3.4-intel mkl/2015.0.090 openmpi/1.8.3-intel lapack/3.5.0-intel
fprefix=fprefix_here

njob=nodes_here
ppp=ppn_here
SCRATCHDIR=$SCRATCH/
SCRIPT_HOME=/global/home/users/tpascal/scripts/
cp2k_cmd="mpirun -mca btl self,sm,openib /global/home/users/tpascal/codes/bin/cp2k.popt"
curr_dir=$SLURM_SUBMIT_DIR
results_dir=$curr_dir/results

mkdir -p $SCRATCHDIR/${fprefix}/sp
cd $SCRATCHDIR/${fprefix}/sp

echo "CP2K singlepoint calculation of ${fprefix}"
cp $curr_dir/${fprefix}.sp.cp2k.in ./
cp $curr_dir/${fprefix}.cp2k.xyz ./

$cp2k_cmd -in ${fprefix}.sp.cp2k.in > ${fprefix}.sp.cp2k.screen.out

