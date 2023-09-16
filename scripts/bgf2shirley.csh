#!/bin/tcsh
#!/bin/csh

if ($#argv == 0) then
	echo "usage: $0 bgf_file [num_procs] [ppn=8] [njobs]"
	exit(1)
endif

if !(-e $1) then
	echo "ERROR: Cannot locate $1"
	exit(1)
endif
set bgf_file = $1

set totatoms = `egrep "^CONECT " $bgf_file | tail -1 | awk '{print $2}'`
set nodes = 8
if ($nodes > $totatoms) set nodes = $totatoms
if ($#argv > 1) set nodes = $2
set ppn = 8
if ($#argv > 2) set ppn = $3
@ totprocs = $nodes * $ppn
@ njobs = $totatoms / $totprocs
if ($#argv > 3) set njobs = $4

if (`echo $totprocs $njobs | awk '{if ($1 % $2 > 0) print 1; else print 0;}'`) then
    echo "ERROR: Total number of processors: $totprocs not multiple of number of simultaneous jobs: $njobs"
    exit(1)
endif
set prefix = $bgf_file:r
set molname = `basename $prefix`
set xyzfile = ${molname}.xyz
babel -ibgf $bgf_file -oxyz $xyzfile >& /dev/null || goto error
set cell = (`grep CRYSTX $bgf_file | awk '{for (i=2;i<8;i++) print $i;}'`)
if ($#cell == 0) set cell = (20 20 20 90 90 90)

cat > ${prefix}.shirley.in <<DATA
# Version of the code
SHIRLEY_ROOT=/global/home/groups-sw/nano/software/sl-6.x86_64/shirley_QE4.3
# pseudopotentials
PSEUDO_DIR=/global/home/groups-sw/nano/software/sl-6.x86_64/shirley_QE4.3/pseudo

# Parallelization variables
# Number of atomic calculations that can run simultaneously
NJOB=$njobs
# procs per pool - used for parallelization of diag within Shirley
PPP=$ppn

PARA_PREFIX="mpirun -mca btl openib,tcp,self"
PARA_POSTFIX=""
PW_POSTFIX="-ntg \$PPP"

# where will output be dumped
TMP_DIR=./

# xc functional
PSEUDO_FUNCTIONAL=pbe
#pseudopotential postfix
PSEUDO_POT_ES_POST_FIX=pbe-van-dgp-1s1.UPF

# specific details for this calculation
MOLNAME="$molname"

#XAS info
XAS_ARG=3
CHAPPROX="XCH"

#Defines the variables from the 'system' namespace
IBRAV=14
A=$cell[1]
B=$cell[2]
C=$cell[3]
cosBC=$cell[4]
cosAC=$cell[5]
cosAB=$cell[6]
ECUT_WFC='45'
ECUT_RHO='360'
NBND_FAC=2
OCCUPATIONS='smearing'
SMEARING='fd'
DEGAUSS=0.0019

#Defines the variables from the 'electrons' namespace
DIAG='david'
DIAG_NSCF='david'
ELEC_CONV_THR='1.0d-8'
ELEC_MIXING_BETA='0.3'

# SCF k-points
K_POINTS="K_POINTS automatic
 1 1 1 0 0 0"

# coordinates
XYZFILES="$xyzfile"
XYZUNIT="angstrom"
XASELEMENTS="O"

#other options
FORCESTRESS="tprnfor=.true.
tstress=.true."
DATA

cat > ${prefix}.shirley.pbs <<DATA
#PBS -q vulcan_batch
#PBS -l nodes=${nodes}:ppn=${ppn}
#PBS -l walltime=48:00:00
#PBS -j oe

export prefix=${molname}
export results_dir=\$PBS_O_WORKDIR/results
export scratch_dir=\$SCRATCH/wat/md/cp2k/\${prefix}/spectra

mkdir -p \$scratch_dir \$results_dir

cd \$scratch_dir
cp \$PBS_O_WORKDIR/\${prefix}.shirley.in ./Input_Block.in
cp \$PBS_O_WORKDIR/${xyzfile} ./
cp \$HOME/scripts/arvid/*.sh ./

export PBS_O_WORKDIR=\$PWD
sh ./XAS-xyz.sh
DATA

exit:
echo "All tasks completed"
exit(0)

error:
echo "Error occurred"
exit(1)
