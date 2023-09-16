if [ -n "$SLURM_SUBMIT_DIR" ]; then
    PARA_JOB_DIR=$SLURM_SUBMIT_DIR
elif [ -n "$PBS_O_WORKDIR" ]; then
    PARA_JOB_DIR=$PBS_O_WORKDIR
else
    PARA_JOB_DIR=$(pwd)
fi

cd "$PARA_JOB_DIR"
source "$PARA_JOB_DIR/Input_Block.in"

DIR=$(pwd)
SHIRLEY_BIB="$SHIRLEY_ROOT/scripts/arvid/Bibliothek"

export SHIRLEY_BIB

source "$SHIRLEY_BIB/Util.sh"
source "$SHIRLEY_BIB/Parallel.sh"

print "____________________________________________________________________________________________________"
print "                ____ ___ ____ ___ ____    ___  ____ ____  _ ____ ____ ___ _ ____ _  _               "
print "                [__   |  |__|  |  |___    |__] |__/ |  |  | |___ |     |  | |  | |\ |               "
print "                ___]  |  |  |  |  |___    |    |  \ |__| _| |___ |___  |  | |__| | \|               "
print "____________________________________________________________________________________________________"
print

# Reset the number of simultaneous calculations
PARA_NJOB=1

setup_parallel_environment

# Chop up the total number of processors into chunks defined by PPN
$SHIRLEY_BIB/Chop.sh $PARA_NODEFILE $PARA_NPROCS_PER_ATOM

# Get timestamp
DATE_START=$(date)
TIMESTAMP_START=$(date +"%s")

print_line
print " Starting execution at: $DATE_START"
print " Running on $PARA_NPROCS_PER_ATOM processors"
print_line
print
print "Running directory $DIR"

ABORT=false
BIN_DIR="$SHIRLEY_ROOT/bin"
SCRIPT_DIR="$SHIRLEY_ROOT/scripts/arvid/Bibliothek"
STATES_DIR="$DIR/states"

# Check for directories
for dir in "$STATES_DIR" ; do
	[ -d "$dir" ] || mkdir -p "$dir"
done

# Check for binaries
for bin in "$BIN_DIR/pp.x" ; do
	[ -x "$bin" ] || { print "Required binary '$bin' not found or not executable"; ABORT=true; }
done

# Check for scripts
for script in "$SCRIPT_DIR/GetStates.sh" ; do
	[ -x "$script" ] || { print "Required script '$script' not found or not executable"; ABORT=true; }
done

print "Verifying existence of required directories, executables and resources..."
print

if [ "$ABORT" = true ]; then
	print "Aborting..."
	exit 1;
fi

$SHIRLEY_BIB/GetNode.sh "$STATES_DIR" "$PARA_NODEPOOL" || continue
ls "$STATES_DIR/Node-${PARA_JOB_ID}-*" > /dev/null

setup_parallel_prefix "$STATES_DIR"
LOCAL_PARA_PREFIX=$PARA_PREFIX
CMD_PP="$LOCAL_PARA_PREFIX $BIN_DIR/pp.x $PARA_PW_POSTFIX $PARA_POSTFIX"

print "Results directory:  $STATES_DIR"
print "Source directory:   $BIN_DIR"
print
print "Running pp.x as:    $CMD_PP"
print
print "Iterating over specified states"

structure=$(basename $(find "$DIR" -type f -name "*.xyz" | tail -1 | sed 's/.xyz//g'))
natoms=$(sed -n 1p "$DIR/${structure}.xyz" | tr -d '[[:space:]]')

# Check if states array has been defined otherwise invoke GetStates.sh
if [ ${#states[@]} -lt 1 ]; then
    states_string=$($SCRIPT_DIR/GetStates.sh -f)
    states=($states_string)
fi

for state in ${states[@]}; do

	atom_number=$(echo $state | awk 'BEGIN{FS="#"}{print $1}')
	band_number=$(echo $state | awk 'BEGIN{FS="#"}{print $2}')
	spin_number=$(echo $state | awk 'BEGIN{FS="#"}{print $3}')
	[ -z "$spin_number" ] && spin_number=1
	element_symbol=$(echo $atom_number | sed -n 's/\([A-Za-z]\+\).*/\1/p')
	element_number=$(echo $atom_number | sed -n 's/[A-Za-z]\+\([0-9]\+\)/\1/p')

	atom_short=$(printf "%s%d"  $element_symbol $element_number)
	atom_full=$(printf "%s%0${#natoms}d" $element_symbol $element_number)

	tmp_dir="$DIR/XAS/${structure}/${atom_short}"
	prefix="${MOLNAME}.${atom_full}-${CHAPPROX}"
	out_dir="$STATES_DIR/${atom_full}.${CHAPPROX}.${band_number}.${spin_number}"

	[ -d "$out_dir" ] || mkdir -p "$out_dir"

	state_in="$out_dir/${atom_full}.${CHAPPROX}.${band_number}.${spin_number}.in"
	state_out="$out_dir/${atom_full}.${CHAPPROX}.${band_number}.${spin_number}.out"
	[ -z "$state_type" ] && state_type='cube'
	case $state_type in 
		"xsf"		) format=5 ;; 
		"cube" | * 	) format=6; state_type='cube' ;;
	esac
	state_file="$out_dir/${atom_full}.${CHAPPROX}.${band_number}.${spin_number}.$state_type"

	if [ -f "$state_file" ]; then
		print "Output file for ${atom_short}#${band_number}#${spin_number} at '${state_file}' already exists"
		print "Assuming the calculation was already completed successfully previously: skipping..."
		continue
	fi

	cat > $state_in <<- EOF
	&inputpp
		prefix  = '$prefix'
		outdir  = '$tmp_dir/',
		plot_num= 7
		kpoint  = $spin_number
		kband   = $band_number
		lsign   = .true.
	/
	&plot
		iflag         = 3
		output_format = $format
		fileout       = '$state_file'
	/
	EOF
	print "Running post-processing for ${atom_short} and band number ${band_number} ${spin_number}..."
	$CMD_PP < $state_in > $state_out

done

rm tmp.pp

release_node "$STATES_DIR"

# Get timestamp
STOP=`date`
TIMESTAMP_STOP=$(date +"%s")

# Compute time elapsed
diff=$(($TIMESTAMP_STOP - $TIMESTAMP_START))
days=$(($diff / 86400))
rem=$(($diff % 86400))
hours=$(($rem / 3600))
rem=$(($rem % 3600))
mins=$(($rem / 60))
secs=$(($rem % 60))

print
print "#============================================================#"
print "# Stopping execution at: $STOP"
print "# Wall time            : ${days}d ${hours}h ${mins}m ${secs}s"
print "#============================================================#"
