#!/bin/bash

function WM {

local Checkdir=$1 #ARG1: directroy which should be checked by the watchman

local Startdir=`pwd` #starting dir

cd "$Checkdir"

while [ 1 ] ; do

Check=`ls | grep "Node-${PBS_JOBID}"`

if [[ -n $Check ]] ; then

break 

fi

sleep 1
done

cd "$Startdir"

}

WM $1

exit
