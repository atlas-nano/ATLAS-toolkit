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
print "                        ____ _  _ _ ____ _    ____ _   _    _  _ ____ ____                          "
print "                        [__  |__| | |__/ |    |___  \_/      \/  |__| [__                           "
print "                        ___] |  | | |  \ |___ |___   |      _/\_ |  | ___]                          " 
print "____________________________________________________________________________________________________"
print

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

for xyzfile in $XYZFILES; do

    if [ -f "$xyzfile" ]; then
        print "Working on xyz-file $xyzfile"
    else
        print "Error: unable to find xz-file $xyzfile - skipping"
        continue
    fi

    # This defines the name of the directory where calculation results are stored
    CALC=$(echo "$xyzfile" | sed 's/\.xyz$//')

    if [ ! -d "$DIR/XAS/$CALC" ]; then
        mkdir -p "$DIR/XAS/$CALC"
    fi

    # Make a copy of Input_Block.in in the working directory
    input=$DIR/XAS/$CALC/Input_Block.in${PARA_JOB_ID}
    cp "$DIR/Input_Block.in" "$input"

    print "Converting $xyzfile to input"
    $SHIRLEY_BIB/../xyz2inp.sh "$xyzfile" $XYZUNIT $XASELEMENTS >> "$input"

    source "$input"

    for i in $(seq 1 $NAT) ; do

        if [ ${IND_EXCITATION[$i]} -eq 1 ]; then

            atom_dir="$DIR/XAS/$CALC/${ATOM_SYMBOL[$i]}${i}"

            if [ ! -d $atom_dir ]; then
                mkdir -p "$atom_dir"
            fi

            cp "$input" "$atom_dir/Input_Block.in"

            CATOM=$(seq -w $i $NAT | head -1)

            # Create file: TMP_INPUT.in
            atom_input="$atom_dir/TMP_INPUT.in${PARA_JOB_ID}"
            $SHIRLEY_BIB/ResetVariables.sh "$atom_dir"
            $SHIRLEY_BIB/VarPen.sh  "$atom_input" TMP_BIN_DIR="$SHIRLEY_ROOT/bin" TMP_MOLNAME=${MOLNAME}.${ATOM_SYMBOL[$i]}${CATOM}-${CHAPPROX} TMP_wf_collect=.TRUE. TMP_disk_io='high'
            $SHIRLEY_BIB/Labeler.sh $i "$atom_input"

            # Get nodes or skip atom if there was an error
            $SHIRLEY_BIB/GetNode.sh "$atom_dir" $PARA_NODEPOOL || continue

            # Refresh file system
            ls $atom_dir/Node-${PARA_JOB_ID}-* > /dev/null

            $SHIRLEY_BIB/DoXAS.sh "$atom_dir" $i ${ATOM_SYMBOL[$i]} &

            sleep 1

            # Refresh file system
            ls $PARA_NODEPOOL > /dev/null
        fi
    done
done

wait
exit