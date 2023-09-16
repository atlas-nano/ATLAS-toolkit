#!/bin/bash

function F_GetNode {
local Targetdir=$1 #ARG 1: directory where the next free node should be copied
local Sourcedir=$2 #ARG 2: directory where the nodes are stored

local MyDir=`echo $0 | sed -e 's/GetNode.sh//g'`	#Saves the directory where this script is saved (=Bibliothek)
#echo $MyDir

local Startdir=`pwd`

local dump=`ls $Targetdir | grep Node-*`
local SKIP=0
for file in $dump ; do
  local id=`basename $file | sed -e 's/Node-//' -e 's/\..*//'`
  local pbsid=`echo $PBS_JOBID | sed -e 's/\..*//'`
  if [[ $id != $pbsid ]] ; then
    if [[ -n $SLURM_JOBID ]] ; then
	local status=`squeue | grep $id | awk '{print $5}'`
    else
	local status=`qstat | grep $id | awk '{print $5}'`
    fi
    if [[ $status == 'R' ]] ; then
      echo
      echo "There is a running PBS_JOBID ($id) for this atom
            in directory $Targetdir"
      echo 
      SKIP=1
      break
    else
      rm -f $Targetdir/$file
    fi
  fi
done
if [[ $SKIP -ne 0 ]] ; then
  echo "skipping this atom to be safe"
  exit 1
fi


sleep 1

${MyDir}WatchMan.sh ${Sourcedir} #wait until the next node is free

cd $Sourcedir

local Node=`ls|grep "Node-${PBS_JOBID}" |sed -n '1p' `


mv $Node ${Targetdir}/${Node}

cd $Startdir
sleep 1 
}

F_GetNode $1 $2

exit

