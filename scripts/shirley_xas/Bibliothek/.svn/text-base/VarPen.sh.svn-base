#!/bin/bash

#This Function writes the argument in TMP_INPUT.in
#Arg: Operation (VARIABLE=VALUE)

function F_VarPen {

local i=''
local File=$1
shift

for i in "$@" ; do

local Variable=`echo "$i" | sed -n "1p" |sed -e 's/=/ /g' | awk '{print $1}'`

# delete old instances of Variable
if (grep -q "${Variable}=" $File) ; then

  local LN=`less $File |grep -n "${Variable}=" |  sed -e 's/:/ /g' | awk '{print $1}' `
  if (grep -q "${Variable}=\".*\" *$" $File) ; then
    OldLength=1
  elif (grep -q "${Variable}=\"" $File) ; then
    OldLength=`sed -n "/${Variable}=/,/\"$/p" $File | wc -l`
  else
    OldLength=1
  fi

  if [[ $OldLength -eq 1 ]] ; then
    sed -i "${LN}d" $File	#removes the old value
  else
    LM=$(( ${LN} + ${OldLength} - 1 ))
    sed -i "${LN},${LM}d" $File	#removes the old multiline value
  fi
fi

sed -i "/^$/d" $File	#removes all blank lines

local Length=`echo "$i" | wc -l`
if [[ $Length -gt 1 ]] ; then
  i=`echo "$i" | sed -e "1s/${Variable}=/${Variable}=\"/" -e "${Length}s/$/\"/"`
fi
cat >> $File <<EOF
$i
EOF

done


}

F_VarPen "$@"

exit
