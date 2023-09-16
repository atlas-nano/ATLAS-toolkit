export SCRATCH_DIR=$SHSCRATCH
export nprocs=$(( 32 * $SLURM_NNODES ))
echo "RUNNING on $nprocs processors on SHAHEEN"
export LAMMPS_PARALLEL="srun --hint=nomultithread --ntasks=$nprocs --ntasks-per-node=32 --ntasks-per-socket=16 --mem-bind=v,local --cpu-bind=cores"
export codes_dir=/project/k1493/x_pascalta/codes
export scripts_dir=/project/k1493/x_pascalta/scripts
export lmp_exec="${codes_dir}/bin/lammps_shaheen"
export RSCRATCH=$SHSCRATCH/results

module swap PrgEnv-cray PrgEnv-intel
module unload cray-libsci
module load gsl python/2.7.18-cdl
