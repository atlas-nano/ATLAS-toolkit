#!/bin/bash

#Function put File1 in File2 in a specific format
#Arg1: Specifies the kind of union (1=MD;2= )
#Arg2: Inital File, which contains the 'from_scratch' results
#Arg3: New File, which contains the 'restart' information and which will contain the combined information at the end of this executable


function F_Reverend {

local Husband=$2
local Wife=$3

if [[ $1 -eq 1 ]] ; then	#Union of MD files

#Check at which step the wife file starts

local WST=`less $Wife | grep Entering | sed -n '1p' | awk '{print $5}'`  

local HUB=`less $Husband | grep 'Entering' | wc -l`


if [[ $HUB -ge $WST ]] 	#Check whether the husband trajectory is long enough
then

local HFR=`grep -n "Entering Dynamics:    iteration = *${WST}$" $Husband| awk '{print $1}' | sed -e 's/://g'` #First row of the range, which should be deleted in the husband file

else
local HFR=`grep -n 'The maximum number of steps has been reached' $Husband | awk '{print $1}' | sed -e 's/://g'`

fi 
local HusbandCut=`less $Husband | sed -n "1,$(( ${HFR}-1 ))p"`








local WFR=`grep -n 'Entering' $Wife | sed -n "1p" | awk '{print $1}' | sed -e 's/://g' `
local WLR=`less $Wife | wc -l` #last row of the wife file

local WifeCut=`less $Wife | sed -n "${WFR},${WLR}p"`

cat > $Wife <<EOF
$HusbandCut

$WifeCut
EOF
rm $Husband

else

echo "Nothing"

fi 



}


F_Reverend $1 $2 $3 

exit
