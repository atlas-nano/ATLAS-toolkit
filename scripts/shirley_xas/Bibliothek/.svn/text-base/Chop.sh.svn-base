#!/bin/bash

#This script decomposes the $PBS_NODEFILE into single nodes 

function F_Chop {

local File=$1    #file to decompose (=PBS_NODEFILE)
local Nl=$2      #chunk size; For one node per file: Nl=PPN (processors per node)  


if [[ $# -eq 3 ]] ; then

local Dir=$3     # (optional) path where the directory for the created node files should be created (default=path of the calling script)

else

local Dir=`pwd`

fi


mkdir -p "${Dir}/Nodes"
NPROC=`wc -l $1 | awk '{print $1}'`

i=0
j=1
while [[ $i -lt $NPROC ]]
do
  ip1=$(( i + 1 ))
  ipN=$(( i + Nl ))
  
  sed -n "${ip1},${ipN}p" $File | tr '\n' ',' | sed 's/,$//' > ${Dir}/Nodes/Node-${PBS_JOBID}-${j}
  j=$(( j + 1 ))
  i=$ipN
done

}

F_Chop $1 $2 $3

exit

