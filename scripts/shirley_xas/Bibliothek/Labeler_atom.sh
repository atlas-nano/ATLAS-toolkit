#!/bin/bash

#ARG1: ID-Number of Atom to label

function F_Labeler {

local NAX=$1 #Arg1
local File=$2 

. $File

#label atomic positions

#echo "Labeler!";
#echo "$TMP_ATOMIC_POSITIONS" ;
#echo "$TMP_ATOMIC_POSITIONS" | sed -n "$(( $NAX + 1))p" ;
local LAtom=`echo "$TMP_ATOMIC_POSITIONS" | sed -n "$(( $NAX + 1))p" |awk '{print $1}'`
#echo $LAtom ;
local SEDSTR="$(( $NAX+1 ))s/${LAtom} /${LAtom}X /g"    #Label the excited atom
#echo $SEDSTR ;
local TMP_ATOMIC_POS=`echo "$TMP_ATOMIC_POSITIONS" |sed -e "$SEDSTR" `
#echo $TMP_ATOMIC_POS ;
#grep -n "TMP_ATOMIC_POSITIONS" $File;
local SL=`grep -n "TMP_ATOMIC_POSITIONS" $File | sed -e 's/:/ /g'| awk '{print $1}' `
#echo $SL ;
local Length=`echo "$TMP_ATOMIC_POSITIONS" | wc -l`
#echo $Length ;

sed -i "$SL,$(( $SL -1 + $Length ))d" $File
sed -i "/^$/d" $File	#removes all blank lines

cat >> $File <<EOF
TMP_ATOMIC_POSITIONS="$TMP_ATOMIC_POS"
EOF


#label atomic species

local AMass=`echo "$TMP_ATOMIC_SPECIES" | grep "^ *$LAtom " | awk '{print $2}'`
#local AMass=${ATOM_MASS[${NAX}]}
SL=`grep -n "TMP_ATOMIC_SPECIES" $File | sed -e 's/:/ /g'| awk '{print $1}' `


Length=`echo "$TMP_ATOMIC_SPECIES" | wc -l`
local NTYP=$(( $Length - 1 ))
local NTYP2=`echo "$TMP_ATOMIC_POS"| awk '{print $1}'| sort | uniq -c |wc -l`
NTYP2=$(( ${NTYP2} - 1  ))

sed -i "$SL,$(( $SL -1 + $Length ))d" $File
sed -i "/^$/d" $File	#removes all blank lines

if [[ $NTYP -eq $NTYP2 ]] ; then

local TMP=$TMP_ATOMIC_SPECIES

local TMPN=`echo -e "$TMP" | sed -e '1d' | grep -n "$LAtom " | sed -e 's/:/ /g' | awk '{print $1}'`

TMPN=$(( $TMPN + 1  ))

TMP=`echo "$TMP" | sed -e "${TMPN}d"`

cat >> $File <<EOF
TMP_ATOMIC_SPECIES="${TMP}
${LAtom}X ${AMass} ${LAtom}.${TMP_PSEUDO_POT_ES_POST_FIX}"

EOF

else

cat >> $File <<EOF
TMP_ATOMIC_SPECIES="$TMP_ATOMIC_SPECIES
${LAtom}X ${AMass} ${LAtom}.${TMP_PSEUDO_POT_ES_POST_FIX}"
EOF

fi

sed -i "/^$/d" $File	#removes all blank lines

}

F_Labeler $@

exit
