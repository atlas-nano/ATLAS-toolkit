#!/bin/bash

for i in *.qcin
do
	nm=${i%.*}
	mkdir -p results/${nm}/
	cat > ${nm}.qcscript <<DATA
#!/bin/bash
#Job name:
#SBATCH --job-name=${nm}
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
#SBATCH --nodes=1
#
#Processors:
#SBATCH --ntasks=8
#
#QUEUE
#SBATCH --constraint=vulcan
#
#Wall clock limit:
#SBATCH --time=10:0:0
#

module unload intel
module load qchem/5.1.2
nprocs=8
NN=8
prefix=${nm}
temp_dir=/clusterfs/vulcan/pscratch/tpascal/geopt/qchem/
cd \$SLURM_SUBMIT_DIR
mkdir -p \$temp_dir
export QCSCRATCH=\$temp_dir
echo \$SLURM_NODELIST > \${temp_dir}/\${SLURM_JOB_ID}.machinefile
export QCMACHINEFILE=\${temp_dir}/\${SLURM_JOB_ID}.machinefile
export QCMPIRUN=mpirun

cd \$temp_dir
cp \${SLURM_SUBMIT_DIR}/\${prefix}.qcin ./
qchem -np 8 -save \${prefix}.qcin \${prefix}.qcout \${prefix}
mv plot.mo \${prefix}.plot.mo 
rm -fr $QCMACHINEFILE
#cp -r results/plots/* \${SLURM_SUBMIT_DIR}/results/\${prefix}
#sed '1,/^======= MOLDEN-FORMATTED INPUT FILE FOLLOWS =======/d; /^======= END OF MOLDEN-FORMATTED INPUT FILE =======/,$ d' \${prefix}.qcout \${prefix}.molden
DATA
done
