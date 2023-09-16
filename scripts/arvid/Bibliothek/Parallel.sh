#!/usr/bin/env bash
# Title         : Parallel.sh
# Description   : Defines functions that setup environment variables needed for parallel processing of Shirley
# Dependencies  : Requires that path to Shirley Bibliothek directory is defined in order to source Util.sh

source "$SHIRLEY_BIB/Util.sh"

setup_parallel_environment ()
{
    if [ -n "$SLURM_SUBMIT_DIR" ]; then
        PARA_SCHEDULER="slurm"
    elif [ -n "$PBS_O_WORKDIR" ]; then
        PARA_SCHEDULER="torque"
    fi

    if [ $PARA_SCHEDULER == "slurm" ]; then
        PARA_JOB_DIR=$SLURM_SUBMIT_DIR
        PARA_JOB_ID=$SLURM_JOB_ID
        PARA_NNODES=$SLURM_NNODES
        PARA_NPROCS=$SLURM_NPROCS
        PARA_NTASKS_PER_NODE=$SLURM_NTASKS_PER_NODE
        PARA_NODELIST=$(scontrol show hostnames $SLURM_JOB_NODELIST)
    elif [ $PARA_SCHEDULER == "torque" ]; then
        PARA_JOB_DIR=$PBS_O_WORKDIR
        PARA_JOB_ID=$PBS_JOBID
        PARA_NNODES=$PBS_NUM_NODES
        PARA_NPROCS=$PBS_NP
        PARA_NTASKS_PER_NODE=$PBS_NUM_PPN
        PARA_NODELIST=$(<$PBS_NODEFILE)
    fi

    PARA_NODEPOOL="$PARA_JOB_DIR/Nodes"
    PARA_NODEFILE="$PARA_JOB_DIR/Nodefile-$PARA_JOB_ID"
    PARA_NODELIST_COMMA=$(echo "$PARA_NODELIST" | tr '\n' ', ' | sed 's/,$//g')

    :> $PARA_NODEFILE
    for node in $PARA_NODELIST; do
        for task in $(seq 1 $PARA_NTASKS_PER_NODE); do
            echo $node >> $PARA_NODEFILE
        done
    done

    # Divide total allocated number of processors among atomic calculations
    PARA_NPROCS_PER_ATOM=$(echo $PARA_NPROCS / $PARA_NJOB | bc)
    PARA_NNODES_PER_ATOM=$(echo $PARA_NNODES / $PARA_NJOB | bc)

    if [ $((PARA_NNODES % PARA_NJOB)) -ne 0 ]; then
        print "Error: the total number of nodes $PARA_NNODES should be divisible by the number of simultaneous jobs PARA_NJOB = $PARA_NJOB"
        print "Error: the remainder is $((PARA_NNODES % PARA_NJOB))"
        exit
    fi

    if [ $((PARA_NPROCS_PER_ATOM % PARA_NPROCS_PER_POOL)) -ne 0 ]; then
        print "Error: the number of processors for each atom $PARA_NPROCS_PER_ATOM should be divisible by processors per pool $PARA_NPROCS_PER_POOL for shirley_xas.x"
        print "Error: the remainder is $((PARA_NPROCS_PER_ATOM % PARA_NPROCS_PER_POOL))"
        exit
    fi

    print_header " Parallel processing environment parameters"
    print " Job ID                        = $PARA_JOB_ID"
    print " Job hosts                     = $PARA_NODELIST_COMMA"
    print " Job submit directory          = $PARA_JOB_DIR"
    print_line '-'
    print " Number of nodes               = $PARA_NNODES"
    print " Number of processors          = $PARA_NPROCS"
    print " Number of tasks per node      = $PARA_NTASKS_PER_NODE"
    print " Number of simultaneous atoms  = $PARA_NJOB"
    print_line '-'
    print " Number of nodes per atom      = $PARA_NNODES_PER_ATOM"
    print " Number of processors per atom = $PARA_NPROCS_PER_ATOM"
    print " Number of processors per pool = $PARA_NPROCS_PER_POOL"
    print_line
    print

    export PARA_SCHEDULER
    export PARA_JOB_DIR
    export PARA_JOB_ID
    export PARA_NNODES
    export PARA_NPROCS
    export PARA_NTASKS_PER_NODE
    export PARA_NPROCS_PER_POOL
    export PARA_NPROCS_PER_ATOM
    export PARA_NNODES_PER_ATOM
    export PARA_NODEPOOL
    export PARA_NODEFILE
    export PARA_NODELIST
    export PARA_NODELIST_COMMA
    export PARA_COMMAND
    export PARA_COMMAND_FULLPATH
    export PARA_PREFIX_FORMAT
}

setup_parallel_prefix ()
{
    if [ -z "$PARA_COMMAND" ]; then
        print "Required variable PARA_COMMAND is not set"
        exit
    fi

    if [ -z "$PARA_COMMAND_FULLPATH" ]; then
        PARA_COMMAND_FULLPATH=$(which $PARA_COMMAND)
    fi

    if [ -z "$PARA_PREFIX_FORMAT" ]; then
        case "$PARA_COMMAND" in
            "srun")
                PARA_PREFIX_FORMAT="-n __PARA_NPROCS_PER_ATOM -N __PARA_NNODES_PER_ATOM" ;;
            "mpirun")
                PARA_PREFIX_FORMAT="-n __PARA_NPROCS_PER_ATOM -host __PARA_HOST_LIST" ;;
            *)
                print "Error: could not automatically determine the undefined but required PARA_PREFIX_FORMAT variable"
                exit
        esac
    fi

    # Optional argument to pass current working directory, look for node file
    if [ -n "$1" ]; then
        nodefiles=($(find "$1" -name "Node-*"))

        if [ "${#nodefiles[@]}" -eq 1 ]; then
            PARA_HOST_LIST=$(cat ${nodefiles[0]})
            print "Located nodefile ${nodefiles[0]}"
        else
            print "Warning: could not find nodefile in directory '$1'"
        fi
    fi

    PREFIX_FORMAT=$(echo "$PARA_PREFIX_FORMAT" | sed \
        -e "s/__PARA_NPROCS_PER_ATOM/$PARA_NPROCS_PER_ATOM/g" \
        -e "s/__PARA_NNODES_PER_ATOM/$PARA_NNODES_PER_ATOM/g" \
        -e "s/__PARA_HOST_LIST/$PARA_HOST_LIST/g")

    PARA_PREFIX="$PARA_COMMAND_FULLPATH $PREFIX_FORMAT"
}

release_node ()
{
    local dir=$1
    local node=$(find "$dir" -name "Node-${PARA_JOB_ID}-*")

    if [ -f "$node" ]; then
        mv "$node" "$PARA_NODEPOOL/"
    else
        print "Error: could not find a node file in '$dir' to release"
    fi
}

error_handler ()
{
    directory=$1
    exitstatus=$?

    print_line
    print "Error: problem with an mpirun job"
    print "Error: exit status $exitstatus"

    if [ $exitstatus -ne 0 ]; then
        print "exiting"
        print_line
        release_node "$1"
        exit 1
    fi

    print_line
}