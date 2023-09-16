#!/bin/bash

#Arg1: directory of the output of the MD calculation
#Arg2: Indicator whether run should be parallel (1= parallel; 0= not parallel)

. ./TMP_INPUT.in${PBS_JOBID}

function F_MD {

local MD_OUTDIR=$1

local TmpDir=`pwd`
local NODE_POOL="${TmpDir}/Nodes"

local NAME=$TMP_MOLNAME
local PREFIX="$NAME"
local INPUT="$MD_OUTDIR/$PREFIX.in"
local OUTPUT="$MD_OUTDIR/$PREFIX.out"
local MyDir=`echo $0 | sed -e 's/MD.sh//g'`



#Check if run parallel
if [[ $2 -eq 1 ]] ; then

local MD_START_DIR=`pwd`
cd
cd "$MD_OUTDIR"
local MD_HOSTFILE=`ls Node*`
local PPN=`sed 's/,/ /g' "$MD_HOSTFILE" | wc -w`
cd
cd $MD_START_DIR


local MD_HOSTS=`cat "${MD_OUTDIR}/${MD_HOSTFILE}"` 
local PW_COMMAND="$TMP_PARA_PREFIX -host ${MD_HOSTS} -np $PPN $TMP_BIN_DIR/pw.x $TMP_PARA_POSTFIX"

else	#Not parallel
local PW_COMMAND="$TMP_PARA_PREFIX $TMP_BIN_DIR/pw.x $TMP_PARA_POSTFIX"
fi

#Check if XCH or FCH

if test "$TMP_CHAPPROX" == "XCH" ; then

local NELEC="$(( ${TMP_NELEC}+1 ))"

else
local NELEC="${TMP_NELEC}"
fi

#echo $NELEC

# Set ES values for the MD
local ES_INDICATOR=`echo "$TMP_ATOMIC_POSITIONS" | sed "1d" | awk '{print $1}'  |  grep "X"`

#echo $ES_INDICATOR

if test "$ES_INDICATOR" != "" ; then	#If ES
local MD_ES_BLOCK="nelec=$NELEC, nbnd=$TMP_NBND, occupations='smearing', degauss=$TMP_DEGAUSS"
else
local MD_ES_BLOCK="occupations='smearing', degauss=$TMP_DEGAUSS"
fi
#echo $MD_ES_BLOCK

#Get the number of atomic species
local NTYP=`echo "$TMP_ATOMIC_POSITIONS"| awk '{print $1}'| sort | uniq -c |wc -l`
NTYP="$(( ${NTYP} - 1 ))"
#echo $NTYP

cat > $INPUT << EOF
&control
        calculation='md'
        restart_mode='$TMP_RM',
        pseudo_dir='$TMP_PSEUDO_DIR/', 
        outdir='$MD_OUTDIR',
        dt=$TMP_DT
        nstep=$TMP_NSTEP
        wf_collect=$TMP_wf_collect 
        disk_io='$TMP_disk_io'
        tstress=.true.
        tprnfor=.true.
/
&system
        ibrav=$TMP_IBRAV, a=$TMP_A, b=$TMP_B, c=$TMP_C
        cosab=$TMP_COSAB, cosac=$TMP_COSAC, cosbc=$TMP_COSBC, 
        nat=$TMP_NAT, ntyp=$NTYP
        ecutwfc=$TMP_ECUT_WFC, ecutrho=$TMP_ECUT_RHO, nosym=.true.
        $MD_ES_BLOCK      
/
&electrons
        conv_thr=$TMP_ELEC_CONV_THR
        mixing_beta=$TMP_ELEC_MIXING_BETA
/
&ions
        pot_extrapolation='$TMP_POT_EXTRAPOLATION'
        wfc_extrapolation='$TMP_WFC_EXTRAPOLATION'
        ion_temperature='$TMP_ION_TEMPERATURE'
        tempw=$TMP_TEMPW
        tolp=$TMP_TOLP
        delta_t=$TMP_DELTA_T
        nraise=$TMP_NRAISE
/
$TMP_ATOMIC_SPECIES
$TMP_ATOMIC_POSITIONS
$TMP_K_POINTS
EOF
#echo $PW_COMMAND
$PW_COMMAND < $INPUT > $OUTPUT

if [[ $2 -eq 1 ]] ; then	# If MD parallel

${MyDir}ReleaseNode.sh "$NODE_POOL" "$MD_OUTDIR"

fi


}

F_MD $1 $2

exit
