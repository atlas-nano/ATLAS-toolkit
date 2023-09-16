#!/bin/bash

if [ $# -eq 0 ]; then
	echo "usage: $0 xyz_file(s)"
	exit 1
fi

list=($1)
for i in ${list[*]}
do
	mol=`basename $i .xyz`
	bgf=`echo $i | sed 's/.xyz/.bgf/'`
	if [ ! -e $bgf ]; then
		continue
	fi
	echo $mol
	cell=(`grep CRYSTX $bgf | awk '{print $2,$3,$4}'`)
	. ~tpascal/scripts/arvid/xyz2inp.sh $i angstrom > __test.sh
	. ./__test.sh
	cat > ${mol}.qe.charges.slurm <<DATA1;
#!/bin/bash -I
#SBATCH --job-name=${mol}.qe.charges
#SBATCH --partition=mako
#SBATCH --qos=mako_short
#SBATCH --nodes=8
#SBATCH --time=00:10:00

prefix=$mol
temp_dir=\$SCRATCH/charges/\${prefix}
curr_dir=\$SLURM_SUBMIT_DIR

#rm -fr \$temp_dir
mkdir -p \$temp_dir \${curr_dir}/results
cd \$temp_dir

cp /global/home/groups-sw/nano/software/sl-6.x86_64/shirley_QE4.3/pseudo/* ./

cat > \${temp_dir}/${mol}.scf.in <<DATA;
&control
    calculation='scf'
    prefix='\${prefix}'
    pseudo_dir='/global/home/groups-sw/nano/software/sl-6.x86_64/shirley_QE4.3/pseudo'
    outdir='./'
    wf_collect=.TRUE.
    disk_io='low'
    restart_mode='from_scratch'
    tprnfor=.true.
tstress=.true.
/ 
&system
    ibrav=14 , a=${cell[0]} , b=${cell[1]} , c=${cell[2]}
    nat=$NAT, ntyp=$NTYP, tot_charge=0.0
    nbnd=$NBND, occupations='smearing', degauss=0.0019
    ecutwfc=25, ecutrho=200
    
    
/
&electrons
    electron_maxstep=30
    conv_thr=1.0d-8
    mixing_beta=0.3
    
/ 
$ATOMIC_SPECIES
$ATOMIC_POSITIONS
DATA

cat > \${temp_dir}/${mol}.nscf.in <<DATA;
&control
    calculation='nscf'
    prefix='\${prefix}'
    pseudo_dir='/global/home/groups-sw/nano/software/sl-6.x86_64/shirley_QE4.3/pseudo'
    outdir='\${temp_dir}'
    wf_collect=.TRUE.
/ 
&system
    ibrav=14 , a=${cell[0]} , b=${cell[1]} , c=${cell[2]}
    nat=$NAT, ntyp=$NTYP, tot_charge=0.0
    nbnd=$NBND, occupations='smearing', degauss=0.0019
    ecutwfc=25, ecutrho=200
    
    
/
&electrons
    diagonalization='david',
    conv_thr=1.0d-8,
    
/ 
$ATOMIC_SPECIES
$ATOMIC_POSITIONS
DATA


cat > \${temp_dir}/${mol}.pdos.in <<DATA;
 &inputpp
    outdir='./'
    prefix='\${prefix}'
    io_choice='both'
    Emin=0.0, Emax=10.0, DeltaE=0.025
    ngauss=1, degauss=0.002
 /
DATA

cat > \${temp_dir}/${mol}.pp.in <<DATA;
&INPUTPP
	prefix='\${prefix}',
	outdir='./',
	filplot='charge_dens.dat',
	plot_num=0
/

&PLOT
	nfile=1,
	filepp(1)='./charge_dens.dat',
	iflag=3,
	output_format=6,
	fileout='\${prefix}.chg.cube',
/
DATA

	mpirun -mca btl self,sm,openib /global/home/groups-sw/nano/software/sl-6.x86_64/shirley_QE4.3/bin/pw.x < ${mol}.scf.in | tee ${mol}.scf.out
	mpirun -mca btl self,sm,openib /global/home/groups-sw/nano/software/sl-6.x86_64/shirley_QE4.3/bin/pw.x < ${mol}.nscf.in | tee ${mol}.nscf.out
	mpirun -mca btl self,sm,openib /global/home/groups-sw/nano/software/sl-6.x86_64/shirley_QE4.3/bin/projwfc.x < ${mol}.pdos.in | tee ${mol}.pdos.out
	mpirun -mca btl self,sm,openib /global/home/groups-sw/nano/software/sl-6.x86_64/shirley_QE4.3/bin/pp.x < ${mol}.pp.in | tee ${mol}.pp.out
	cp *.out \${curr_dir}/results
DATA1
done

