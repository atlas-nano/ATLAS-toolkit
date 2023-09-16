echo "RUNNING on EXPANSE"
export nprocs=$(( $SLURM_NTASKS_PER_NODE * $SLURM_NNODES / 2 ))
export ppn=$(( $SLURM_NTASKS_PER_NODE ))
export SCRATCH_DIR=$EXSCRATCH
#export LAMMPS_PARALLEL="srun --mpi=pmi2 -n $nprocs --ntasks-per-node=$SLURM_NTASKS_PER_NODE"
export LAMMPS_PARALLEL="mpirun -n $nprocs -mca btl vader,self"
export codes_dir=/home/tpascal/codes/
export scripts_dir=/home/tpascal/scripts/
export lmp_exec="${codes_dir}/bin/lmp_expanse"
#export lmp_exec="${codes_dir}/bin/lmp_expanse_pqeq"
export RSCRATCH=/expanse/lustre/projects/csd626/$USER/results/
module purge
#module load cpu/0.15.4 intel/19.1.1.217 slurm mvapich2 cmake intel-mkl sdsc fftw gsl
#module load cpu sdsc slurm gcc/10.2.0 cmake  mvapich2/2.3.4 amdfftw/2.2 python gsl
#module load cpu/0.15.4 slurm gcc/10.2.0 openmpi cmake gsl intel-mkl amdfftw slurm sdsc gsl
#module load cpu/0.15.4 intel/19.1.1.217 cpu slurm openmpi cmake intel-mkl sdsc fftw
#module load cpu/0.17.3b  gcc/10.2.0/npcyll4 python/3.8.12/7zdjza7 cmake/3.21.4/teqow32 intel-mpi/2019.10.317/kdx4qap fftw/3.3.10/sigimvh gsl/2.7/wtlsmyy intel-mkl/2020.4.304/ghfk3mu gnuplot/5.4.2/mfinpvw slurm
module load cpu/0.17.3b  gcc/10.2.0/npcyll4 cmake/3.21.4/teqow32 openmpi/4.1.3/oq3qvsv fftw/3.3.10/bmvqsbr python/3.8.12/7zdjza7 intel-mkl/2020.4.304/ghfk3mu gsl/2.7/wtlsmyy gnuplot/5.4.2/mfinpvw sdsc slurm
#export LAMMPS_PARALLEL="mpirun -n $nprocs"
export lmp_exec="${codes_dir}/bin/lmp_expanse"
export MKL_ROOT=$INTEL_MKLHOME
export MKLROOT=$INTEL_MKLHOME
