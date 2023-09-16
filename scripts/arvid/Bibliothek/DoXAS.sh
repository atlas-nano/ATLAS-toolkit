#!/usr/bin/env bash
# Title         : DoXAS.sh
# Description   : Runs a Shirley XAS calculation
# Argument      : Directory of the output of the XAS calculation
# Argument      : Index of the excited atom
# Argument      : Element symbol of the excited atom
# Dependencies  : Requires that path to Shirley Bibliothek directory is defined in order to source Util.sh and Parallel.sh

source "$SHIRLEY_BIB/Util.sh"
source "$SHIRLEY_BIB/Parallel.sh"

DIR=$1
ATOM_INDEX=$2
ATOM_ELEMENT=$3

cd $DIR

source "$DIR/TMP_INPUT.in${PARA_JOB_ID}"

setup_parallel_prefix "$DIR"

LOCAL_PARA_PREFIX="$PARA_PREFIX"

# Get timestamp
DATE_START=$(date)
TIMESTAMP_START=$(date +"%s")

NATOMS=$TMP_NAT
LABEL_SHORT=$(printf "%s%d" $ATOM_ELEMENT $ATOM_INDEX)
LABEL_FULL=$(printf "%s%0${#NATOMS}d" $ATOM_ELEMENT $ATOM_INDEX)

print
print_line
print " Starting execution for $LABEL_SHORT at: $DATE_START"
print " Running on $PARA_NPROCS_PER_ATOM processors"
print_line
print_labeled $LABEL_FULL "Running directory $DIR"

# Check what actually needs to be done and define JOB
SCFOUT="$TMP_MOLNAME.scf.out"
NSCFOUT="$TMP_MOLNAME.nscf.out"
BASISOUT="$TMP_MOLNAME.basis.out"
HAMOUT="$TMP_MOLNAME.ham.out"
XASOUT="$TMP_MOLNAME.xas.out"

JOB=''
if [ -f "$SCFOUT" ] && grep -q 'convergence has been achieved' "$SCFOUT"; then

	if [ -f "$XASOUT" ] && grep -q 'end shirley_xas' "$XASOUT"; then
		print " this atom seems to be completed already"
		print " if this is an error then delete the output files in $DIR"
		release_node "$DIR"
		exit 0
	fi

	if [ -f "$HAMOUT" ] && grep -q 'fft_scatter' "$HAMOUT"; then
		JOB=1
	fi
	if [ -z $JOB ] && [ -f "$BASISOUT" ] && grep -q 'fft_scatter' "$BASISOUT"; then
		JOB=2
	fi
	if [ -z $JOB ] && [ -f "$NSCFOUT" ] && grep -q 'fft_scatter' "$NSCFOUT"; then
		JOB=3
	fi
	if [ -z $JOB ]; then
		JOB=4
	fi

else
	JOB=5
fi

# Check if XCH or FCH and change TOT_CHARGE
TOT_CHARGE=$TMP_tot_charge
if [ -z "$TMP_tot_charge" ]; then
	TOT_CHARGE=0.0
fi
if [ "$TMP_CHAPPROX" == "FCH" ]; then
	TOT_CHARGE=$(echo "$TOT_CHARGE + 1.0" | bc)
fi

# Get the number of atomic species
NTYP=$(echo "$TMP_ATOMIC_POSITIONS"| awk '{print $1}'| sort | uniq -c | wc -l)
NTYP="$(($NTYP - 1))"

run_scf ()
{
	if [ -n "$1" ]; then
		TMP_STARTINGPOT=$1
	else
		TMP_STARTINGPOT=atomic
	fi

	local PREFIX="$TMP_MOLNAME.scf"
	local INPUT="$PREFIX.in"
	local OUTPUT="$PREFIX.out"

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

	cat > "$INPUT" <<- EOFSCF
	&control
	    calculation='scf',
	    prefix='$TMP_MOLNAME',
	    pseudo_dir='$TMP_PSEUDO_DIR',
	    outdir='$DIR',
	    wf_collect =$TMP_wf_collect,
	    restart_mode='$TMP_restart_mode',
	    $TMP_FORCESTRESS
	/
	&system
	    $LATTICE
	    nat=$TMP_NAT, ntyp=$NTYP, tot_charge=$TOT_CHARGE,
	    nbnd=$TMP_NBND, occupations='smearing', smearing='$TMP_SMEARING', degauss=$TMP_DEGAUSS,
	    ecutwfc =$TMP_ECUT_WFC, ecutrho=$TMP_ECUT_RHO,
	    $TMP_SPIN
	    $TMP_LDAU
	/
	&electrons
	    diagonalization='$TMP_DIAG'
	    conv_thr=$TMP_ELEC_CONV_THR,
	    mixing_beta=$TMP_ELEC_MIXING_BETA,
	    mixing_mode=$TMP_MIXING_MODE,
	    electron_maxstep=$TMP_electron_maxstep,
	    diago_thr_init=$TMP_diago_thr_init,
	    startingpot='$TMP_STARTINGPOT',
	/
	$TMP_ATOMIC_SPECIES
	$TMP_ATOMIC_POSITIONS
	$TMP_K_POINTS
	$TMP_CELL_PARAMETERS
	EOFSCF
	local PW="$LOCAL_PARA_PREFIX $TMP_BIN_DIR/pw.x $TMP_PARA_POSTFIX $TMP_PARA_PW_POSTFIX_SCF"
	print_labeled $LABEL_FULL "Running pw.x as: $PW"
	$PW -inp "$INPUT" > "$OUTPUT"

	# Save
	if [ -f "$TMP_MOLNAME.save" ]; then
		cp "$TMP_MOLNAME.save" "$PREFIX.save"
	fi
	if [ -f "$TMP_MOLNAME.rho" ]; then
		cp "$TMP_MOLNAME.rho" "$PREFIX.rho"
	fi
}

run_nscf ()
{
	local PREFIX="$TMP_MOLNAME.nscf"
	local INPUT="$PREFIX.in"
	local OUTPUT="$PREFIX.out"

	local NBND_SCF=$(grep 'number of Kohn-Sham states=' "$TMP_MOLNAME.scf.out" | awk '{print $5}' | tail -1)
	local NBND_NSCF=$(printf "%.0f" $(echo "scale=2;$NBND_SCF * $TMP_NBND_FAC" | bc))
	if [ -n "$TMP_K_POINTS_NSCF" ]; then
		local K_POINTS="$TMP_K_POINTS_NSCF"
	else
		local K_POINTS="K_POINTS automatic
	1 1 1 0 0 0"
	fi

	local LATTICE="ibrav=$TMP_IBRAV"
	if [ -n "$TMP_A" ];         then LATTICE="${LATTICE} , a=$TMP_A"; fi
	if [ -n "$TMP_B" ];         then LATTICE="${LATTICE} , b=$TMP_B"; fi
	if [ -n "$TMP_C" ];         then LATTICE="${LATTICE} , c=$TMP_C"; fi
	if [ -n "$TMP_COSAB" ];     then LATTICE="${LATTICE} , cosab=$TMP_COSAB"; fi
	if [ -n "$TMP_COSAC" ];     then LATTICE="${LATTICE} , cosac=$TMP_COSAC"; fi
	if [ -n "$TMP_COSBC" ];     then LATTICE="${LATTICE} , cosbc=$TMP_COSBC"; fi
	if [ -n "$TMP_CELLDM1" ];   then LATTICE="${LATTICE} , celldm(1)=$TMP_CELLDM1"; fi
	if [ -n "$TMP_CELLDM2" ];   then LATTICE="${LATTICE} , celldm(2)=$TMP_CELLDM2"; fi
	if [ -n "$TMP_CELLDM3" ];   then LATTICE="${LATTICE} , celldm(3)=$TMP_CELLDM3"; fi
	if [ -n "$TMP_CELLDM4" ];   then LATTICE="${LATTICE} , celldm(4)=$TMP_CELLDM4"; fi
	if [ -n "$TMP_CELLDM5" ];   then LATTICE="${LATTICE} , celldm(5)=$TMP_CELLDM5"; fi
	if [ -n "$TMP_CELLDM6" ];   then LATTICE="${LATTICE} , celldm(6)=$TMP_CELLDM6"; fi
	if [ -z "$TMP_DIAG_NSCF" ]; then TMP_DIAG_NSCF="cg"; fi

	cat > "$INPUT" <<- EOFNSCF
	&control
	    calculation = 'nscf',
	    prefix='$TMP_MOLNAME',
	    pseudo_dir = '$TMP_PSEUDO_DIR',
	    outdir='$DIR',
	    wf_collect =$TMP_wf_collect,
	/
	&system
	    $LATTICE
	    nat=$TMP_NAT, ntyp=$NTYP, tot_charge=$TOT_CHARGE,
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
	$K_POINTS
	$TMP_CELL_PARAMETERS
	EOFNSCF

	local PW="$LOCAL_PARA_PREFIX $TMP_BIN_DIR/pw.x $TMP_PARA_POSTFIX $TMP_PARA_PW_POSTFIX_NSCF"
	print_labeled $LABEL_FULL "Running pw.x as: $PW"
	$PW -inp "$INPUT" > "$OUTPUT" || error_handler "$DIR"
}

run_basis ()
{
	local PREFIX="$TMP_MOLNAME.basis"
	local INPUT="$PREFIX.in"
	local OUTPUT="$PREFIX.out"

	if [ -z "$TMP_BASIS_TRACE_TOL" ]; then
		TMP_BASIS_TRACE_TOL='1.d-8'
	fi

	cat > "$INPUT" <<- EOFBASIS
	&input
	    prefix='$TMP_MOLNAME',
	    outdir='$DIR',
	    trace_tol=$TMP_BASIS_TRACE_TOL,
	/
	EOFBASIS

	local BASIS="$LOCAL_PARA_PREFIX $TMP_BIN_DIR/shirley_basis.x"
	print_labeled $LABEL_FULL "Running shirley_basis.x as: $BASIS"
	$BASIS < "$INPUT" > "$OUTPUT" || error_handler "$DIR"
}

run_ham ()
{
	local PREFIX="$TMP_MOLNAME.ham"
	local INPUT="$PREFIX.in"
	local OUTPUT="$PREFIX.out"

	# Check for spin
	local SPIN=$(echo "$TMP_SPIN" | grep nspin | sed -r 's/.*nspin\s*=\s*([0-9]).*/nspin_ham=\1/')
	cat > $INPUT <<- EOFHAM
	&input
	    prefix='${TMP_MOLNAME}_opt',
	    outdir='$DIR',
	    updatepp=.false.,
	    ncpp=.false.,
	    pseudo_dir='$TMP_PSEUDO_DIR',
	    $SPIN
	/
	K_POINTS
	    2 2 2 0 0 0
	    $UPDATEPP
	EOFHAM

	local HAM="$LOCAL_PARA_PREFIX $TMP_BIN_DIR/shirley_ham.x"
	print_labeled $LABEL_FULL "Running shirley_ham.x as: $HAM"
	$HAM < "$INPUT" > "$OUTPUT" || error_handler "$DIR"
}

run_xas ()
{
	local nk=$1
	local PREFIX="$TMP_MOLNAME.xas"
	local INPUT="$PREFIX.in"
	local OUTPUT="$PREFIX.out"

	local SUFFIX=$(echo ${TMP_PSEUDO_POT_ES_POST_FIX} | sed 's/\.UPF.*$//')

	cat > "$INPUT" <<- EOFXAS
	&input
	    prefix='${TMP_MOLNAME}_opt',
	    outdir='$DIR',
	    outfile='$TMP_MOLNAME.xas.dump',
	    readcorerep=.true.,
	/
	K_POINTS
	    automatic
	    $nk $nk $nk 0 0 0
	COREREP
	    $NTYP $TMP_NAT
	    1
	    $NTYP  $ATOM_INDEX  '$TMP_PSEUDO_DIR/${ATOM_ELEMENT}.${SUFFIX}.pos'
	EOFXAS

	local XAS="$LOCAL_PARA_PREFIX $TMP_BIN_DIR/shirley_xas.x $PARA_NPROCS_PER_POOL"
	print_labeled $LABEL_FULL "Running shirley_xas.x as: $XAS"
	$XAS < "$INPUT" > "$OUTPUT" || error_handler "$DIR"

	local XASPARA="$LOCAL_PARA_PREFIX $TMP_BIN_DIR/xas_para.x"
	print_labeled $LABEL_FULL "Running xas_para.x as: $XASPARA"
	$XASPARA -30 40 1000 0.1 0 "$TMP_MOLNAME.xas.dump" || error_handler "$DIR"

	rename dump "$nk" *dump*
}


if [ $JOB -ge 5 ]; then

	print_labeled $LABEL_FULL "Starting SCF..."

	# Defaults
	if [ -z "$TMP_SMEARING" ]; then
		TMP_SMEARING='fd'
	fi
	if [ -z "$TMP_mixing_mode" ]; then
		TMP_mixing_mode='plain'
	fi
	if [ -z "$TMP_electron_maxstep" ]; then
		TMP_electron_maxstep=30
	fi
	TMP_restart_mode='from_scratch'

	if [ -z "$TMP_DEGAUSS_STEP" ]; then

		run_scf $TMP_STARTINGPOT

		# SCF did not converge so retry
		if grep -q 'convergence NOT achieved.* stopping' "$TMP_MOLNAME.scf.out"; then
			count=0
			while [[ $count -lt 10 ]]; do
				count=$((count + 1))
				print_labeled $LABEL_FULL "restarting unconverged SCF"
				mv "$TMP_MOLNAME.scf.out" "$TMP_MOLNAME.scf.out-$count"
				run_scf file

				if grep -q 'convergence has been achieved' "$TMP_MOLNAME.scf.out"; then
					break
				fi
			done

			if grep -q 'convergence NOT achieved.* stopping' "$TMP_MOLNAME.scf.out"; then
				print_labeled $LABEL_FULL "Could not converge SCF - exiting"
				exit
			fi
		fi
		print_labeled $LABEL_FULL "SCF completed"

	else

		TMP_DEGAUSS_ORIG=$TMP_DEGAUSS
		TMP_DEGAUSS=$(echo "$TMP_DEGAUSS $TMP_DEGAUSS_FAC $TMP_DEGAUSS_STEP" | awk '{print $1 * ($2 ^ ($3 - 1))}')

		for g in $(seq 1 $TMP_DEGAUSS_STEP); do
			print_labeled $LABEL_FULL "SCF step $g of $TMP_DEGAUSS_STEP using degauss = $TMP_DEGAUSS"
			run_scf
			cp "$TMP_MOLNAME.scf.out" "$TMP_MOLNAME.scf.out-$g"
			TMP_DEGAUSS=$(echo "$TMP_DEGAUSS $TMP_DEGAUSS_FAC" | awk '{print $1 / $2}')
			TMP_restart_mode='restart'
		done

		TMP_DEGAUSS=$TMP_DEGAUSS_ORIG
	fi
fi

if [ $JOB -ge 4 ]; then
	print_labeled $LABEL_FULL "Starting NSCF..."
	run_nscf
	print_labeled $LABEL_FULL "NSCF completed"
fi

if [ $JOB -ge 3 ]; then
	print_labeled $LABEL_FULL "Starting BASIS..."
	run_basis
	print_labeled $LABEL_FULL "BASIS completed"
fi

if [ $JOB -ge 2 ]; then
	print_labeled $LABEL_FULL "Starting HAM..."
	run_ham
	print_labeled $LABEL_FULL "HAM completed"
fi

if [ $JOB -ge 1 ]; then
	print_labeled $LABEL_FULL "Starting XAS..."
	run_xas $TMP_XAS_ARG
	print_labeled $LABEL_FULL "XAS completed"
fi

# Get timestamp
DATE_STOP=$(date)
TIMESTAMP_STOP=$(date +"%s")

# Compute time elapsed
diff=$(($TIMESTAMP_STOP - $TIMESTAMP_START))
days=$(($diff / 86400))
rem=$(($diff % 86400))
hours=$(($rem / 3600))
rem=$(($rem % 3600))
mins=$(($rem / 60))
secs=$(($rem % 60))

print_line
print " Stopping execution for $LABEL_SHORT at: $DATE_STOP"
print " Wall time: ${days}d ${hours}h ${mins}m ${secs}s"
print_line
print

release_node "$DIR"

exit
