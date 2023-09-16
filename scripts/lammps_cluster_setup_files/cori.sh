echo "RUNNING on CORI"
export SCRATCH_DIR=$CORISCRATCH
export nprocs=$(( $SLURM_NPROCS / 2 ))
export LAMMPS_PARALLEL="srun -c 1 --cpu_bind=cores -n $nprocs"
export codes_dir=/global/homes/t/tpascal/codes
export scripts_dir=/global/homes/t/tpascal/scripts
export lmp_exec="${codes_dir}/bin/lmp_cori"
export RSCRATCH=$CORISCRATCH
module load gnuplot
