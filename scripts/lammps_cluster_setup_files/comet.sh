echo "RUNNING on COMET"
export SCRATCH_DIR=$CSCRATCH
export LAMMPS_PARALLEL="mpirun"
export codes_dir=/home/tpascal/codes/
export scripts_dir=/home/tpascal/scripts/
export lmp_exec="${codes_dir}/bin/lmp_comet_intelmpi"
export RSCRATCH=/oasis/projects/nsf/csd626/tpascal/results
module purge
module load gnu intel intelmpi mkl
