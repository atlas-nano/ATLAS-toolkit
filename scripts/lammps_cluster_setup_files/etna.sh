echo "RUNNING on ETNA"
export SCRATCH_DIR=$ESCRATCH
export UCX_TLS=ud,sm,self
module purge
module load intel/2018.1.163 openmpi/3.0.1-intel mkl/2018.1.163 fftw/3.3.8-intel gnuplot
export nprocs=$(( $SLURM_NNODES * $SLURM_NPROCS ))
LAMMPS_PARALLEL="mpirun"
