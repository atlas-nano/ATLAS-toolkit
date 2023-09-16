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

do_md=1
do_md_xas=1
do_xtal_xas=1
nsnaps=nsnaps_here
rtemp=rtemp_here
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

if [ $do_md -eq 1 ]; then
	echo "CP2K MD at ${rtemp} K"
	cd $SCRATCHDIR/${fprefix}/${cell_p}
	mkdir -p ./md/cp2k/${rtemp}K/ ${results_dir}/MD/
	cd ./md/cp2k/${rtemp}K/
	cp $curr_dir/${fprefix}.${cell_p}.${rtemp}K.cp2k.in ./
	cp $curr_dir/${fprefix}.${cell_p}.cp2k.xyz ./
	$cp2k_cmd -in ${fprefix}.${cell_p}.${rtemp}K.cp2k.in > ${fprefix}.${cell_p}.${rtemp}K.out 
	cp ${fprefix}.${cell_p}.* ${results_dir}/MD
fi

if [ $do_md_xas -eq 1 ]; then
	echo "Calculating Spectra from MD structures"
	cd $SCRATCHDIR/${fprefix}/${cell_p}
	mkdir -p ./xas/md/${rtemp}K/ ${results_dir}/xas/MD/${rtemp}K
	cd ./xas/md/${rtemp}K/
	cp $curr_dir/${fprefix}.${cell_p}.Input_Block.in ./Input_Block.in
	cp ${SCRIPT_HOME}/arvid/*.sh ./

	cell_a=`csh -f ${SCRIPT_HOME}/get_snaps.csh ${SCRATCHDIR}/${fprefix}/${cell_p}/md/cp2k/${rtemp}K/${fprefix}.${cell_p}.${rtemp}K $nsnaps ${fprefix}.${cell_p}.${rtemp}K | awk '{print $2}'`
	cell_b=`csh -f ${SCRIPT_HOME}/get_snaps.csh ${SCRATCHDIR}/${fprefix}/${cell_p}/md/cp2k/${rtemp}K/${fprefix}.${cell_p}.${rtemp}K $nsnaps ${fprefix}.${cell_p}.${rtemp}K | awk '{print $3}'`
	cell_c=`csh -f ${SCRIPT_HOME}/get_snaps.csh ${SCRATCHDIR}/${fprefix}/${cell_p}/md/cp2k/${rtemp}K/${fprefix}.${cell_p}.${rtemp}K $nsnaps ${fprefix}.${cell_p}.${rtemp}K | awk '{print $4}'`
	celldm_1=`echo $cell_a | awk '{print $1*1.889725989}'`
	celldm_2=`echo $cell_a $cell_b | awk '{print $2/$1}'`
	celldm_3=`echo $cell_a $cell_c | awk '{print $2/$1}'`
	sed -i 's/^A=.*$/CELLDM1='$celldm_1'/' ./Input_Block.in
	sed -i 's/^B=.*$/CELLDM2='$celldm_2'/' ./Input_Block.in
	sed -i 's/^C=.*$/CELLDM3='$celldm_3'/' ./Input_Block.in
	sed -i 's/^COSBC=.*$/CELLDM4='$c_alpha'/' ./Input_Block.in
	sed -i 's/^COSAC=.*$/CELLDM5='$c_beta'/' ./Input_Block.in
	sed -i 's/^COSAB=.*$/CELLDM6='$c_gamma'/' ./Input_Block.in
	sed -i 's/^MOLNAME=.*$/MOLNAME="'${fprefix}.${cell_p}.${rtemp}K'"/' ./Input_Block.in
	sed -i 's/^NJOB=.*/NJOB='$njob'/' ./Input_Block.in
	sed -i 's/^PPP=.*/PPP='$ppp'/' ./Input_Block.in

	export SLURM_SUBMIT_DIR=$PWD
	sh ./XAS-xyz.sh
	sh ./XAS-xyz-ref.sh
	sh ./XASAnalyse-xyz.sh

	cp -r ./XAS/Spectrum-* $results_dir/xas/MD/${rtemp}K/
	cd ../
fi

if [ $do_xtal_xas -eq 1 ]; then
	echo "Calculating spectra from Crystal structure"
	cd $SCRATCHDIR/${fprefix}/${cell_p}
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
	sed -i 's/^MOLNAME=.*$/MOLNAME="'${fprefix}.${cell_p}.crystal'"/' ./Input_Block.in
	sed -i 's/^NJOB=.*$/NJOB='$njob'/' ./Input_Block.in
	sed -i 's/^PPP=.*$/PPP='$ppp'/' ./Input_Block.in

	export SLURM_SUBMIT_DIR=$PWD
	sh ./XAS-xyz.sh
	sh ./XAS-xyz-ref.sh
	sh ./XASAnalyse-xyz.sh

	cp -r ./XAS/Spectrum-* $results_dir/xas/xtal/
fi
