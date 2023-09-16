#!/bin/bash

#scratch directories for various machines
VSCRATCH=/clusterfs/vulcan/pscratch/$USER/
VSCRATCH=/global/scratch/users/$USER/
ESCRATCH=/clusterfs/etna/pscratch/$USER/
CSCRATCH=/oasis/scratch/comet/$USER/temp_project
CORISCRATCH=/global/cscratch1/sd/$USER
SHSCRATCH=/scratch/$USER
EXSCRATCH=/expanse/lustre/scratch/$USER/temp_project
#EXSCRATCH=/scratch/$USER/job_${SLURM_JOB_ID}

#set defaults when using the LRC machines
export LAMMPS_PARALLEL="mpirun -mca btl vader,self"
export codes_dir=/global/home/users/tpascal/codes
export scripts_dir=/global/home/users/tpascal/scripts
export lmp_exec="${codes_dir}/bin/lmp_mftheory.march2018"
export lmp_exec="${codes_dir}/lammps/lammps-30Apr19/src/lmp_mftheory"
export SCRATCH_DIR=$SCRATCH
export RSCRATCH=$SCRATCH
export nprocs=$SLURM_NPROCS
export ppn=$SLURM_TASKS_PER_NODE

cluster=$SLURM_CLUSTER_NAME
if [ -z "${cluster// }" ] || [ $(echo ${cluster} | egrep -c perceus) -eq 1 ]; then
	cluster=$SLURM_JOB_PARTITION
fi
if [ -z "${cluster// }" ]; then
	echo "ERROR: Could not determine the name of the cluster... Searched \$SLURM_JOB_PARTITION and \$SLURM_CLUSTER_NAME environment variables"
	exit 1
fi
cluster=$(echo $cluster | tr '[[:upper:]]' '[[:lower:]]')
sdir=$(realpath $(dirname ${BASH_ARGV[0]}))

if [ ! -e ${sdir}/lammps_cluster_setup_files/${cluster}.sh ]; then
	echo "ERROR: Could not locate ${sdir}/lammps_cluster_setup_files/${cluster}.sh"
	exit 1
fi

. ${sdir}/lammps_cluster_setup_files/${cluster}.sh

#if [ $( echo $SLURM_JOB_PARTITION | egrep -ic '(vulcan|nano)') -gt 0 ]; then
#	echo "RUNNING on VULCAN"
#	export SCRATCH_DIR=$VSCRATCH
#   	module load intel/2018.1.163 openmpi/3.0.1-intel mkl/2018.1.163 fftw/3.3.8-intel gnuplot
#   	export nprocs=$(( $SLURM_NNODES * $SLURM_NPROCS ))
#elif [ $( echo $SLURM_JOB_PARTITION | egrep -ic 'etna') -gt 0 ]; then
#	echo "RUNNING on ETNA"
#   	export SCRATCH_DIR=$ESCRATCH
#   	module load intel/2018.1.163 openmpi/3.0.1-intel mkl/2018.1.163 fftw/3.3.8-intel gnuplot
#    LAMMPS_PARALLEL="mpirun"
#   	export nprocs=$(( $SLURM_NNODES * $SLURM_NPROCS ))
#elif [ $( echo $SLURM_CLUSTER_NAME | egrep -ic 'comet') -gt 0 ]; then
#	echo "RUNNING on COMET"
#   	export SCRATCH_DIR=$CSCRATCH
#   	export LAMMPS_PARALLEL="mpirun"
#   	export codes_dir=/home/tpascal/codes/
#   	export scripts_dir=/home/tpascal/scripts/
#   	export lmp_exec="${codes_dir}/bin/lmp_comet_intelmpi"
#   	export RSCRATCH=/oasis/projects/nsf/csd626/tpascal/results
#   	module purge
#   	module load gnu intel intelmpi mkl
#elif [ $( echo $SLURM_CLUSTER_NAME | egrep -ic 'expanse') -gt 0 ]; then
#   	echo "RUNNING on EXPANSE"
#   	nprocs=$(( $SLURM_NTASKS_PER_NODE * $SLURM_NNODES / 2 ))
#	ppn=$(( $SLURM_NTASKS_PER_NODE / 2 ))
#   	SCRATCH_DIR=$EXSCRATCH
#   	LAMMPS_PARALLEL="srun --mpi=pmi2 -n $nprocs"
#   	codes_dir=/home/tpascal/codes/
#   	scripts_dir=/home/tpascal/scripts/
#   	lmp_exec="${codes_dir}/bin/lmp_expanse"
#   	RSCRATCH=/expanse/lustre/projects/ddp381/$USER/results/
#   	#RSCRATCH=${SLURM_SUBMIT_DIR}/../results
#   	module purge
#	module load cpu/0.15.4 intel/19.1.1.217 slurm openmpi cmake intel-mkl sdsc fftw gcc gsl
#   	export MKL_ROOT=$INTEL_MKLHOME
#   	export MKLROOT=$INTEL_MKLHOME    
#elif [ $( echo $SLURM_CLUSTER_NAME | egrep -ic 'cori') -gt 0 ]; then
#	echo "RUNNING on CORI"
#   	export SCRATCH_DIR=$CORISCRATCH
#   	export nprocs=$(( $SLURM_NPROCS / 2 ))
#   	export LAMMPS_PARALLEL="srun -c 1 --cpu_bind=cores -n $nprocs"
#   	export codes_dir=/global/homes/t/tpascal/codes
#   	export scripts_dir=/global/homes/t/tpascal/scripts
#   	export lmp_exec="${codes_dir}/bin/lmp_cori"
#   	export RSCRATCH=$CORISCRATCH
#   	module load gnuplot
#elif [ $( echo $SLURM_CLUSTER_NAME | egrep -ic 'shaheen') -gt 0 ]; then
#   	export SCRATCH_DIR=$SHSCRATCH
#   	export nprocs=$(( 32 * $SLURM_NNODES ))
#   	echo "RUNNING on $nprocs processors on SHAHEEN"
#   	export LAMMPS_PARALLEL="srun --hint=nomultithread --ntasks=$nprocs --ntasks-per-node=32 --ntasks-per-socket=16 --mem-bind=v,local --cpu-bind=cores"
#   	export codes_dir=/project/k1493/x_pascalta/codes
#   	export scripts_dir=/project/k1493/x_pascalta/scripts
#   	export lmp_exec="${codes_dir}/bin/lmp_shaheen"
#   	export RSCRATCH=$SHSCRATCH/results
#
#   	module swap PrgEnv-cray PrgEnv-intel
#   	module unload cray-libsci
#	module load gsl python/2.7.18-cdl
#fi
