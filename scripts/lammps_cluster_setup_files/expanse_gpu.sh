echo "RUNNING on EXPANSE using $SLURM_GPUS GPUs"
export nprocs=$(( $SLURM_NTASKS_PER_NODE * $SLURM_NNODES ))
export ppn=$(( $SLURM_NTASKS_PER_NODE ))
export SCRATCH_DIR=$EXSCRATCH
export codes_dir=/home/tpascal/codes/
export scripts_dir=/home/tpascal/scripts/
export lmp_exec="${codes_dir}/bin/lmp_expanse"
export RSCRATCH=/expanse/lustre/projects/ddp381/$USER/results/
module purge
module load cpu/0.17.3b  gcc/10.2.0/npcyll4 python/3.8.12/7zdjza7 cmake/3.21.4/teqow32 intel-mpi/2019.10.317/kdx4qap fftw/3.3.10/sigimvh gsl/2.7/wtlsmyy intel-mkl/2020.4.304/ghfk3mu gnuplot/5.4.2/mfinpvw slurm
MPI_ROOT=/home/tpascal/codes/openmpi/openmpi-4.1.5_gpu
export PATH=${MPI_ROOT}/bin:$PATH
export LAMMPS_PARALLEL="mpirun -n $nprocs -mca btl vader,self"
export lmp_exec="${codes_dir}/bin/lmp_expanse_intel_mpi_gpu -sf gpu -pk gpu $SLURM_GPUS"
