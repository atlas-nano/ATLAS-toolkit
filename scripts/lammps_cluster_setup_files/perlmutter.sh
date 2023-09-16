export codes_dir=/global/homes/t/tpascal/codes
export scripts_dir=/global/homes/t/tpascal/scripts
module load cudatoolkit
module load craype-accel-nvidia80
export nprocs=$(( $SLURM_TASKS_PER_NODE * $SLURM_NNODES / 4 ))
echo "RUNNING on PERLMUTTER"
echo "using $SLURM_TASKS_PER_NODE tasks, $SLURM_NNODES nodes and $SLURM_CPUS_PER_TASK CPUs/tasks"
export LAMMPS_PARALLEL="srun  --cpu_bind=cores -C $SLURM_CPUS_PER_TASK -n $nprocs"
export LAMMPS_PARALLEL="srun  -n $nprocs"
export SCRATCH_DIR=$PMSCRATCH
export RSCRATCH=$PMSCRATCH
if [ -n "${SLURM_CPUS_PER_TASK// }" ]; then
	export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
	export OMP_NUM_THREADS=1
fi
if [ $(echo $SLURM_JOB_PARTITION | egrep -ic 'gpu') -gt 0 ]; then
	echo "GPU run activated!"
	export lmp_exec="${codes_dir}/bin/lmp_perlmutter -k on g 4 -sf kk -pk kokkos newton on neigh half binsize 2.8"
	export lmp_exec="${codes_dir}/bin/lmp_perlmutter -k on g 2 -sf kk"
else
	echo "RUNNING on CPUs only"
	export lmp_exec="${codes_dir}/bin/lmp_perlmutter -k on t $SLURM_CPUS_PER_TASK -sf kk"
fi
