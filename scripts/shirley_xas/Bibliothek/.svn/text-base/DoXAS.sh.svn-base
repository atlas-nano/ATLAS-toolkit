#!/bin/bash

#Arg1: directory of the output of the XAS calculation
#Arg2: index of the excited atom
#Arg3: location of NodePool
TmpDir=`pwd`
OutDir=$1
XATOM=$2
XATOM_SYMBOL=$3
NODE_POOL=$4

cd $OutDir

. ./TMP_INPUT.in${PBS_JOBID}


#get hostfile information
XAS_HOSTFILE=`ls | grep "Node-${PBS_JOBID}" | sed -n '1p'`
PPN=`sed 's/,/ /g' "$XAS_HOSTFILE" | wc -w`
if [ -z $CRAY_SITE_LIST_DIR ]; then
  XAS_HOSTS=`cat "${XAS_HOSTFILE}"`
  TMP_PARA_PREFIX="$TMP_PARA_PREFIX -host ${XAS_HOSTS} -np $PPN"
else
  mppnppn=`qstat -f $PBS_JOBID | grep Resource_List.mppnppn | awk '{print $3}'`
  TMP_PARA_PREFIX="$TMP_PARA_PREFIX -n $PPN -N $mppnppn"
fi


function release {

if [ -z $BibDir ]
then
  ./ReleaseNode.sh "$NODE_POOL" "$OutDir"
else
  ${BibDir}/ReleaseNode.sh "$NODE_POOL" "$OutDir"
fi

}

#Check what actually needs to be done and define JOB
SCFOUT=$TMP_MOLNAME.scf.out
NSCFOUT=$TMP_MOLNAME.nscf.out
BASISOUT=$TMP_MOLNAME.basis.out
HAMOUT=$TMP_MOLNAME.ham.out
XASOUT=$TMP_MOLNAME.xas.out

JOB=''
if [ -f $SCFOUT ] && grep -q 'convergence has been achieved' $SCFOUT
then

  if [ -f $XASOUT ] && grep -q 'end shirley_xas' $XASOUT
  then
    echo
    echo " this atom seems to be completed already"
    echo " if this is an error then delete the output files in"
    echo " $OutDir"
    echo
    release
    exit 0
  fi

  if [ -f $HAMOUT ] && grep -q 'fft_scatter' $HAMOUT
  then
    JOB=1
  fi
  if [ -z $JOB ] && [ -f $BASISOUT ] && grep -q 'fft_scatter' $BASISOUT
  then
    JOB=2
  fi
  if [ -z $JOB ] && [ -f $NSCFOUT ] && grep -q 'fft_scatter' $NSCFOUT
  then
    JOB=3
  fi
  if [ -z $JOB ]
  then
    JOB=4
  fi

else
  JOB=5
fi
echo JOB = $JOB

#Check if XCH or FCH and change tot_charge
tot_charge=$TMP_tot_charge
if [[ -z $TMP_tot_charge ]]; then
  tot_charge=0.0
fi
if test "$TMP_CHAPPROX" == "FCH" ; then
tot_charge=1.0
fi
echo tot_charge = $tot_charge

#Get the number of atomic species
NTYP=`echo "$TMP_ATOMIC_POSITIONS"| awk '{print $1}'| sort | uniq -c |wc -l`
NTYP="$(( ${NTYP} - 1 ))"
#echo $NTYP



#Functions

function error_handler {
  exitstatus=$?
  echo =================================================================
  echo error_handler
  echo problem with an mpirun job
  echo exit status $exitstatus
  if [ $exitstatus -ne 0 ]
  then
    echo exiting
    echo =================================================================
    release
    exit 1
  fi
  echo =================================================================
}

function scf {
# ----------------------------------------------------------------------
if [[ -n $1 ]]; then
  TMP_startingpot=$1
else
  TMP_startingpot=atomic
fi

local PREFIX=$TMP_MOLNAME.scf
local INPUT=$PREFIX.in
local OUTPUT=$PREFIX.out

local LATTICE="ibrav=$TMP_IBRAV"
if [[ -n $TMP_A ]] ; then LATTICE="${LATTICE} , a=$TMP_A"; fi
if [[ -n $TMP_B ]] ; then LATTICE="${LATTICE} , b=$TMP_B"; fi
if [[ -n $TMP_C ]] ; then LATTICE="${LATTICE} , c=$TMP_C"; fi
if [[ -n $TMP_COSAB ]] ; then LATTICE="${LATTICE} , cosab=$TMP_COSAB"; fi
if [[ -n $TMP_COSAC ]] ; then LATTICE="${LATTICE} , cosac=$TMP_COSAC"; fi
if [[ -n $TMP_COSBC ]] ; then LATTICE="${LATTICE} , cosbc=$TMP_COSBC"; fi
if [[ -n $TMP_CELLDM1 ]] ; then LATTICE="${LATTICE} , celldm(1)=$TMP_CELLDM1"; fi
if [[ -n $TMP_CELLDM2 ]] ; then LATTICE="${LATTICE} , celldm(2)=$TMP_CELLDM2"; fi
if [[ -n $TMP_CELLDM3 ]] ; then LATTICE="${LATTICE} , celldm(3)=$TMP_CELLDM3"; fi
if [[ -n $TMP_CELLDM4 ]] ; then LATTICE="${LATTICE} , celldm(4)=$TMP_CELLDM4"; fi
if [[ -n $TMP_CELLDM5 ]] ; then LATTICE="${LATTICE} , celldm(5)=$TMP_CELLDM5"; fi
if [[ -n $TMP_CELLDM6 ]] ; then LATTICE="${LATTICE} , celldm(6)=$TMP_CELLDM6"; fi

cat > $INPUT <<EOFSCF
&control
    calculation='scf',
    prefix='$TMP_MOLNAME',
    pseudo_dir='$TMP_PSEUDO_DIR',
    outdir='$OutDir',
    wf_collect =$TMP_wf_collect,
    restart_mode='$TMP_restart_mode',
    $TMP_FORCESTRESS
/ 
&system
    $LATTICE
    nat=$TMP_NAT, ntyp=$NTYP, tot_charge=$tot_charge,
    nbnd=$TMP_NBND, occupations='smearing', smearing='$TMP_SMEARING', degauss=$TMP_DEGAUSS,
    ecutwfc =$TMP_ECUT_WFC, ecutrho=$TMP_ECUT_RHO,
    $TMP_SPIN
    $TMP_LDAU
/
&electrons
    diagonalization='$TMP_DIAG'
    conv_thr=$TMP_ELEC_CONV_THR,
    mixing_beta=$TMP_ELEC_MIXING_BETA,
    mixing_mode='$TMP_mixing_mode',
    electron_maxstep=$TMP_electron_maxstep,
    diago_thr_init=$TMP_diago_thr_init,
    startingpot='$TMP_startingpot',
/ 
$TMP_ATOMIC_SPECIES
$TMP_ATOMIC_POSITIONS
$TMP_K_POINTS
$TMP_CELL_PARAMETERS
EOFSCF
local PW="$TMP_PARA_PREFIX $TMP_BIN_DIR/pw.x $TMP_PARA_POSTFIX $TMP_PW_POSTFIX"
#$PW < $INPUT > $OUTPUT || error_handler
# To enable SCF restart
$PW < $INPUT > $OUTPUT


#save
if [ -f $TMP_MOLNAME.save ]
then
  cp $TMP_MOLNAME.save $PREFIX.save 
fi
if [ -f $TMP_MOLNAME.rho ]
then
  cp $TMP_MOLNAME.rho $PREFIX.rho 
fi

}

# NSCF calculation
# ----------------------------------------------------------------------
function nscf {
# ----------------------------------------------------------------------
local PREFIX=$TMP_MOLNAME.nscf
local INPUT=$PREFIX.in
local OUTPUT=$PREFIX.out

local NBND_SCF=`grep 'number of Kohn-Sham states=' $TMP_MOLNAME.scf.out | awk '{print $5}' | tail -1`
local NBND_NSCF=$(( NBND_SCF * TMP_NBND_FAC ))

local LATTICE="ibrav=$TMP_IBRAV"
if [[ -n $TMP_A ]] ; then LATTICE="${LATTICE} , a=$TMP_A"; fi
if [[ -n $TMP_B ]] ; then LATTICE="${LATTICE} , b=$TMP_B"; fi
if [[ -n $TMP_C ]] ; then LATTICE="${LATTICE} , c=$TMP_C"; fi
if [[ -n $TMP_COSAB ]] ; then LATTICE="${LATTICE} , cosab=$TMP_COSAB"; fi
if [[ -n $TMP_COSAC ]] ; then LATTICE="${LATTICE} , cosac=$TMP_COSAC"; fi
if [[ -n $TMP_COSBC ]] ; then LATTICE="${LATTICE} , cosbc=$TMP_COSBC"; fi
if [[ -n $TMP_CELLDM1 ]] ; then LATTICE="${LATTICE} , celldm(1)=$TMP_CELLDM1"; fi
if [[ -n $TMP_CELLDM2 ]] ; then LATTICE="${LATTICE} , celldm(2)=$TMP_CELLDM2"; fi
if [[ -n $TMP_CELLDM3 ]] ; then LATTICE="${LATTICE} , celldm(3)=$TMP_CELLDM3"; fi
if [[ -n $TMP_CELLDM4 ]] ; then LATTICE="${LATTICE} , celldm(4)=$TMP_CELLDM4"; fi
if [[ -n $TMP_CELLDM5 ]] ; then LATTICE="${LATTICE} , celldm(5)=$TMP_CELLDM5"; fi
if [[ -n $TMP_CELLDM6 ]] ; then LATTICE="${LATTICE} , celldm(6)=$TMP_CELLDM6"; fi

if [[ -z $TMP_DIAG_NSCF ]]; then TMP_DIAG_NSCF=cg; fi

cat > $INPUT <<EOFNSCF
&control
    calculation = 'nscf',
    prefix='$TMP_MOLNAME',
    pseudo_dir = '$TMP_PSEUDO_DIR',
    outdir='$OutDir',
    wf_collect =$TMP_wf_collect,
/ 
&system
    $LATTICE
    nat=$TMP_NAT, ntyp=$NTYP, tot_charge=$tot_charge,
    nbnd=$NBND_NSCF, occupations='smearing', smearing='fd', degauss=$TMP_DEGAUSS,
    ecutwfc=$TMP_ECUT_WFC, ecutrho=$TMP_ECUT_RHO,
    $TMP_SPIN
    $TMP_LDAU
/
&electrons
    diagonalization='$TMP_DIAG_NSCF',
    conv_thr=$TMP_ELEC_CONV_THR,
/ 
$TMP_ATOMIC_SPECIES
$TMP_ATOMIC_POSITIONS
  K_POINTS    automatic
 1 1 1 0 0 0
$TMP_CELL_PARAMETERS
EOFNSCF

local PW="$TMP_PARA_PREFIX $TMP_BIN_DIR/pw.x $TMP_PARA_POSTFIX $TMP_PW_POSTFIX"
$PW < $INPUT > $OUTPUT || error_handler

}


# SHIRLEY BASIS
# ----------------------------------------------------------------------
function basis {
# ----------------------------------------------------------------------
local PREFIX=$TMP_MOLNAME.basis
local INPUT=$PREFIX.in
local OUTPUT=$PREFIX.out

if [ -z $TMP_BASIS_TRACE_TOL ]
then
  TMP_BASIS_TRACE_TOL='1.d-8'
fi

cat > $INPUT <<EOFBASIS
&input
  prefix='$TMP_MOLNAME',
  outdir='$OutDir',
  trace_tol=$TMP_BASIS_TRACE_TOL,
/
EOFBASIS

local BASIS="$TMP_PARA_PREFIX $TMP_BIN_DIR/shirley_basis.x"

$BASIS < $INPUT > $OUTPUT || error_handler
}

# SHIRLEY HAMILTONIAN
# ----------------------------------------------------------------------
function ham {
# ----------------------------------------------------------------------
local PREFIX=$TMP_MOLNAME.ham
local INPUT=$PREFIX.in
local OUTPUT=$PREFIX.out

# Check for spin
local SPIN=`echo "$TMP_SPIN" | grep nspin | sed 's/nspin/nspin_ham/'`
cat > $INPUT <<EOFHAM
&input
  prefix='${TMP_MOLNAME}_opt',
  outdir='$OutDir',
  updatepp=.false.,
  ncpp=.false.,
  pseudo_dir='$TMP_PSEUDO_DIR',
  $SPIN
/
K_POINTS
  2 2 2 0 0 0
$UPDATEPP
EOFHAM
local HAM="$TMP_PARA_PREFIX $TMP_BIN_DIR/shirley_ham.x"
$HAM < $INPUT > $OUTPUT || error_handler

}

function xas {
# ----------------------------------------------------------------------
local nk=$1
local PREFIX=$TMP_MOLNAME.xas
local INPUT=$PREFIX.in
local OUTPUT=$PREFIX.out

local SUFFIX=`echo ${TMP_PSEUDO_POT_ES_POST_FIX} | sed 's/\.UPF.*$//'`

cat > $INPUT <<EOFXAS
&input
  prefix='${TMP_MOLNAME}_opt',
  outdir='$OutDir',
  outfile='$TMP_MOLNAME.xas.dump',
  readcorerep=.true.,
/
K_POINTS
  automatic
  $nk $nk $nk 0 0 0
COREREP
  $NTYP $TMP_NAT
  1
  $NTYP  $XATOM  '$TMP_PSEUDO_DIR/${XATOM_SYMBOL}.${SUFFIX}.pos'
EOFXAS

#Calcualting of the Delta
local XAS="$TMP_PARA_PREFIX $TMP_BIN_DIR/shirley_xas.x $PPP"
$XAS < $INPUT > $OUTPUT || error_handler

#local XASPARA="$TMP_PARA_PREFIX $TMP_BIN_DIR/xas_para.x $TMP_PARA_POSTFIX"
local XASPARA="$TMP_PARA_PREFIX $TMP_BIN_DIR/xas_para.x"
echo $XASPARA -30 40 1000 0.1 0 "$TMP_MOLNAME.xas.dump"
$XASPARA -30 40 1000 0.1 0 "$TMP_MOLNAME.xas.dump" || error_handler

rename dump "$nk" *dump*
#mv $TMP_MOLNAME.xas.dump.xas $TMP_MOLNAME.xas.$nk.xas

}


# Executable section

if [ $JOB -ge 5 ]; then
echo "$XATOM_SYMBOL$XATOM : Starting SCF..."
#Defaults
if [ -z $TMP_SMEARING ]
then
  TMP_SMEARING='fd'
fi
if [ -z $TMP_mixing_mode ]
then
  TMP_mixing_mode='plain'
fi
if [ -z $TMP_electron_maxstep ]
then
  TMP_electron_maxstep=30
fi
TMP_restart_mode='from_scratch'
if [ -z $TMP_DEGAUSS_STEP ]; then

  scf $TMP_startingpot
  # uh-oh, it might not have worked
  if grep -q 'convergence NOT achieved.* stopping' $TMP_MOLNAME.scf.out
  then
    count=0
    while [[ $count -lt 10 ]]; do
      count=$((count + 1))
      echo restarting unconverged SCF
      mv $TMP_MOLNAME.scf.out $TMP_MOLNAME.scf.out-$count
      scf file

      if grep -q 'convergence has been achieved' $TMP_MOLNAME.scf.out
      then
        break
      fi
    done

    if grep -q 'convergence NOT achieved.* stopping' $TMP_MOLNAME.scf.out
    then
      echo In directory `pwd`
      echo Seriously messed up SCF calculation
      echo Please investigate
      echo exiting
      exit
    fi
  
  fi
  
  echo "done"

else #TMP_DEGAUSS_STEP

  TMP_DEGAUSS_ORIG=$TMP_DEGAUSS
  TMP_DEGAUSS=`echo "$TMP_DEGAUSS $TMP_DEGAUSS_FAC $TMP_DEGAUSS_STEP" | awk '{print $1 * ($2 ^ ($3 - 1))}'`
  for g in `seq 1 $TMP_DEGAUSS_STEP`; do
    echo "SCF step $g of $TMP_DEGAUSS_STEP using degauss = $TMP_DEGAUSS"
    scf
    cp $TMP_MOLNAME.scf.out $TMP_MOLNAME.scf.out-$g
    TMP_DEGAUSS=`echo "$TMP_DEGAUSS $TMP_DEGAUSS_FAC" | awk '{print $1 / $2}'`
    TMP_restart_mode='restart'
  done
  TMP_DEGAUSS=$TMP_DEGAUSS_ORIG

fi #TMP_DEGAUSS_STEP

fi #JOB 5

if [ $JOB -ge 4 ]; then
echo "$XATOM_SYMBOL$XATOM : Starting NSCF..."
nscf
echo "done"
fi

if [ $JOB -ge 3 ]; then
echo "$XATOM_SYMBOL$XATOM : Starting BASIS..."
basis
echo "done"
fi

if [ $JOB -ge 2 ]; then
echo "$XATOM_SYMBOL$XATOM : Starting HAM..."
ham
echo "done"
fi

if [ $JOB -ge 1 ]; then
echo "$XATOM_SYMBOL$XATOM : Starting XAS..."
xas $TMP_XAS_ARG
echo "done"
fi

release


exit
