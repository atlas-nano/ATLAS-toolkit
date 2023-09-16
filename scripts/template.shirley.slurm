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

fprefix=fprefix_here
cell_p=cell_p_here
c_a=c_a_here
c_b=c_b_here
c_c=c_c_here
c_alpha=`echo c_alpha_here | awk '{pi = atan2(1,1)*4; print cos($1*pi/180)}'`
c_beta=`echo c_beta_here | awk '{pi = atan2(1,1)*4; print cos($1*pi/180)}'`
c_gamma=`echo c_gamma_here | awk '{pi = atan2(1,1)*4; print cos($1*pi/180)}'`

njob=nodes_here
ppp=ppn_here
SCRATCHDIR=$SCRATCH/
SCRIPT_HOME=/global/home/users/tpascal/scripts/
cp2k_cmd="mpirun -mca btl self,sm,openib /global/home/users/tpascal/codes/bin/cp2k.popt"
curr_dir=$SLURM_SUBMIT_DIR
results_dir=$curr_dir/results

mkdir -p $SCRATCHDIR/${fprefix}/${cell_p}
cd $SCRATCHDIR/${fprefix}/${cell_p}

echo "Calculating XAS from ${fprefix}"
mkdir -p ./xas/xtal/ ${results_dir}/xas/xtal/
cd ./xas/xtal/
cp $curr_dir/${fprefix}.${cell_p}.Input_Block.in ./Input_Block.in
cp $curr_dir/${fprefix}.${cell_p}.xyz ./
cp ${SCRIPT_HOME}/arvid/*.sh ./

celldm_1=`echo $c_a | awk '{print $1*1.889725989}'`
celldm_2=`echo $c_a $c_b | awk '{print $2/$1}'`
celldm_3=`echo $c_a $c_c | awk '{print $2/$1}'`
sed -i 's/^A=.*$/CELLDM1='$celldm_1'/' ./Input_Block.in
sed -i 's/^B=.*$/CELLDM2='$celldm_2'/' ./Input_Block.in
sed -i 's/^C=.*$/CELLDM3='$celldm_3'/' ./Input_Block.in
sed -i 's/^COSBC=.*$/CELLDM4='$c_alpha'/' ./Input_Block.in
sed -i 's/^COSAC=.*$/CELLDM5='$c_beta'/' ./Input_Block.in
sed -i 's/^COSAB=.*$/CELLDM6='$c_gamma'/' ./Input_Block.in
sed -i 's/^MOLNAME=.*$/MOLNAME="'${fprefix}.${cell_p}'"/' ./Input_Block.in
sed -i 's/^NJOB=.*$/NJOB='$njob'/' ./Input_Block.in
sed -i 's/^PPP=.*$/PPP='$ppp'/' ./Input_Block.in

export SLURM_SUBMIT_DIR=$PWD
sh ./XAS-xyz.sh
sh ./XAS-xyz-ref.sh
sh ./XASAnalyse-xyz.sh

cp -r ./XAS/Spectrum-* $results_dir/xas/xtal/
