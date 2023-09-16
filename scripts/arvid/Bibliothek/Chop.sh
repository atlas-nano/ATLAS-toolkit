#!/usr/bin/env bash

# This script decomposes the $PARA_NODEFILE into single nodes 
function chop
{
    local PARA_NODEFILE=$1
    local PARA_NPROCS_PER_ATOM=$2

    mkdir -p "$PARA_NODEPOOL"

    NPROCS=$(wc -l "$PARA_NODEFILE" | awk '{print $1}')
    i=0
    j=1
    while [ $i -lt $NPROCS ]; do
        ip1=$((i + 1 ))
        ipN=$((i + PARA_NPROCS_PER_ATOM))

        sed -n "${ip1},${ipN}p" "$PARA_NODEFILE" | tr '\n' ',' | sed 's/,$//' > "$PARA_NODEPOOL/Node-${PARA_JOB_ID}-${j}"
        j=$((j + 1))
        i=$ipN
    done
}

chop $1 $2

exit