#!/bin/bash

function F_ReleaseNode {

local Targetdir=$1 #ARG 1: Dir where all nodes are saved
local Sourcedir=$2 #ARG 2: dir where the node is currently saved

local Startdir=`pwd`
cd $Sourcedir

local Node=`ls|grep "Node-${PARA_JOB_ID}" |sed -n '1p' `
mv "${Sourcedir}/${Node}" "${Targetdir}/"

cd ${Startdir}
}

F_ReleaseNode $1 $2

exit

