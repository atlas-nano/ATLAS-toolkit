#!/usr/bin/env bash
# Title         : SCF.sh
# Description   : Runs a simple self-consistent field calculation
# Argument      : Directory of the output of the SCF calculation
# Argument      : Indicator whether run should be parallel (1 = parallel; 0 = not parallel)
# Argument      : Indicator whether to restart from file or not (1 = restart from file)
# Dependencies  : Requires that path to Shirley Bibliothek directory is defined in order to source Util.sh and Parallel.sh

source "$SHIRLEY_BIB/Util.sh"
source "$SHIRLEY_BIB/Parallel.sh"

DIR=$1
OPT_PARALLEL=$2
OPT_RESTART=$3

cd $DIR

source "$DIR/TMP_INPUT.in${PARA_JOB_ID}"

setup_parallel_prefix "$DIR"

LOCAL_PARA_PREFIX="$PARA_PREFIX"

run_scf ()
{
	local NAME=$TMP_MOLNAME
	local PREFIX="$NAME.scf"
	local INPUT="$DIR/$PREFIX.in"
	local OUTPUT="$DIR/$PREFIX.out"
	local PW_COMMAND="$LOCAL_PARA_PREFIX $TMP_BIN_DIR/pw.x $TMP_PARA_POSTFIX $TMP_PARA_PW_POSTFIX_SCF"

	if [ -n "$OPT_RESTART" ] && [ "$OPT_RESTART" -eq 1 ]; then
		STARTINGPOT="startingpot='file'"
	else
		STARTINGPOT=""
	fi

	# Check which core-hole approximation to use
	local tot_charge=$TMP_tot_charge
	if [ -z "$TMP_tot_charge" ]; then
		tot_charge=0.0
	fi
	if [ "$TMP_CHAPPROX" == "FCH" ]; then
		tot_charge=1.0
	fi
	if [ -z "$TMP_electron_maxstep" ]; then
		TMP_electron_maxstep=30
	fi
	if [ -z "$TMP_MIXING_MODE" ]; then
		TMP_MIXING_MODE='plain'
	fi

	# Get the number of atomic species
	local NTYP=$(echo "$TMP_ATOMIC_POSITIONS" | awk '{if(NR>1){print $1}}'| sort | uniq | wc -l)
	local LATTICE="ibrav=$TMP_IBRAV"
	if [ -n "$TMP_A" ];       then LATTICE="${LATTICE} , a=$TMP_A"; fi
	if [ -n "$TMP_B" ];       then LATTICE="${LATTICE} , b=$TMP_B"; fi
	if [ -n "$TMP_C" ];       then LATTICE="${LATTICE} , c=$TMP_C"; fi
	if [ -n "$TMP_COSAB" ];   then LATTICE="${LATTICE} , cosab=$TMP_COSAB"; fi
	if [ -n "$TMP_COSAC" ];   then LATTICE="${LATTICE} , cosac=$TMP_COSAC"; fi
	if [ -n "$TMP_COSBC" ];   then LATTICE="${LATTICE} , cosbc=$TMP_COSBC"; fi
	if [ -n "$TMP_CELLDM1" ]; then LATTICE="${LATTICE} , celldm(1)=$TMP_CELLDM1"; fi
	if [ -n "$TMP_CELLDM2" ]; then LATTICE="${LATTICE} , celldm(2)=$TMP_CELLDM2"; fi
	if [ -n "$TMP_CELLDM3" ]; then LATTICE="${LATTICE} , celldm(3)=$TMP_CELLDM3"; fi
	if [ -n "$TMP_CELLDM4" ]; then LATTICE="${LATTICE} , celldm(4)=$TMP_CELLDM4"; fi
	if [ -n "$TMP_CELLDM5" ]; then LATTICE="${LATTICE} , celldm(5)=$TMP_CELLDM5"; fi
	if [ -n "$TMP_CELLDM6" ]; then LATTICE="${LATTICE} , celldm(6)=$TMP_CELLDM6"; fi

	cat > $INPUT <<- EOFSCF
	&control
	    calculation='scf'
	    prefix='$NAME'
	    pseudo_dir='$TMP_PSEUDO_DIR'
	    outdir='$DIR'
	    wf_collect=$TMP_wf_collect
	    disk_io='$TMP_disk_io'
	    restart_mode='$TMP_RM'
	    $TMP_FORCESTRESS
	/
	&system
	    $LATTICE
	    nat=$TMP_NAT, ntyp=$NTYP, tot_charge=$tot_charge
	    nbnd=$TMP_NBND, occupations='smearing', degauss=$TMP_DEGAUSS
	    ecutwfc=$TMP_ECUT_WFC, ecutrho=$TMP_ECUT_RHO
	    $TMP_SPIN
	    $TMP_LDAU
	/
	&electrons
	    electron_maxstep=$TMP_electron_maxstep
	    conv_thr=$TMP_ELEC_CONV_THR
	    mixing_mode=$TMP_MIXING_MODE, mixing_beta=$TMP_ELEC_MIXING_BETA
	    $STARTINGPOT
	/
	$TMP_ATOMIC_SPECIES
	$TMP_ATOMIC_POSITIONS
	$TMP_K_POINTS
	$TMP_CELL_PARAMETERS
	EOFSCF

	print_time "Running pw.x as: $PW_COMMAND"
	$PW_COMMAND -inp "$INPUT" > "$OUTPUT"
}

run_scf "$DIR" "$OPT_PARALLEL"

count=0
while ! grep -q 'convergence has been achieved' "$DIR/$TMP_MOLNAME.scf.out"
do
	count=$((count + 1))
	print "Restarting SCF calculcation: restart $count"
	cp "$DIR/$TMP_MOLNAME.scf.out" "$DIR/$TMP_MOLNAME.scf.out-$count"
	run_scf "$DIR" "$OPT_PARALLEL" 1
	if [ "$count" -gt 2 ]; then
		break;
	fi
done

if [ -n "$OPT_PARALLEL" ] && [ "$OPT_PARALLEL" -eq 1 ]; then
	release_node "$DIR"
fi 

exit