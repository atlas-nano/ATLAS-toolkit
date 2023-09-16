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

. $File

#echo $TMP_SPIN 
#Look for magnetized species
TMP_SPIN=`echo $TMP_SPIN | sed -e "s/,/ /g" | sed -e "s/start.*//"`;
NSPIN=`echo $TMP_SPIN | sed -e "s/nspin.*=//" | sed -e "s/^ *//"`;
#echo "myspin=$NSPIN" ;
if [[ $NSPIN -eq 2 ]] ; then
    i=0;
    for elem in `echo $TMP_MAG_SPEC` ; do
	i=$((i+1));
	MAGEL[$i]=$elem;
    done
    i=0;
    for val in `echo $TMP_MAG_VAL` ; do
	i=$((i+1));
	MAGV[$i]=$val;
    done
    NMAGE=$i;
    SPEC0=`echo $TMP_ATOMIC_SPECIES | sed -e "s/ATOMIC_SPECIES//"`;
#    echo "SPEC0=$SPEC0";
    ns=0;
    for element in ${SPEC0} ; do  
	if [[ "$element" =~ [a-zA-Z] ]] && [[ ! "$element" =~ [\.\-] ]] ; then
#	  echo $element; 
          ns=$((ns+1));
          SPEC[$ns]=`echo $element | awk '{print $1}'`;
	fi
    done	
#    OLD_IFS="$IFS";
#    IFS='UPF';
#    array=(${SPEC0}) ;
#    IFS="$OLD_IFS";
#    ns=0;
#    for j in "${array[@]}"
#    do
#	if [[ "$j" =~ [a-zA-Z] ]] ; then
#	ns=$((ns+1));
#	echo $ns;
#	echo $j ;
#	SPEC[$ns]=`echo $j | awk '{print $1}'`;
#	fi
#    done
#    echo $ns;
#    echo ${SPEC[*]};
    for ityp in `seq 1 $ns` ; do
	for i in `seq 1 $NMAGE` ; do
	   if [[ ${MAGEL[$i]} == ${SPEC[$ityp]} ]] || [[ ${MAGEL[$i]}X == ${SPEC[$ityp]} ]]  ; then
#	       echo "equal" , ${MAGEL[$i]}, ${SPEC[$ityp]};
	       TMP_SPIN="$TMP_SPIN , starting_magnetization($ityp)=${MAGV[$i]}";
	       break
	   fi
	done
    done

cat >> $File <<EOF
TMP_SPIN="$TMP_SPIN"
EOF
#echo SPIN=\"$TMP_SPIN\";
fi

#Set LDAU
if [[ -n $TMP_LDAU ]] ; then
TMP_LDAU=`echo $TMP_LDAU | sed -e "s/,/ /g" | sed -e "s/Hub.*//" | sed -e "s/^ *//"`;
    i=0;
    for elem in `echo $TMP_U_SPEC` ; do
	i=$((i+1));
	UEL[$i]=$elem;
    done
    i=0;
    for val in `echo $TMP_U_VAL` ; do
	i=$((i+1));
        UVAL[$i]=$val;
    done
    NUE=$i;

    SPEC0=`echo $TMP_ATOMIC_SPECIES | sed -e "s/ATOMIC_SPECIES//"`;
    ns=0;
    for element in ${SPEC0} ; do  
	if [[ "$element" =~ [a-zA-Z] ]] && [[ ! "$element" =~ [\.\-] ]] ; then
#	  echo $element; 
          ns=$((ns+1));
          SPEC[$ns]=`echo $element | awk '{print $1}'`;
	fi
    done	
#    OLD_IFS="$IFS";
#    IFS="UPF";
#    array=(${SPEC0}) ;
#    IFS="$OLD_IFS";
#    ns=0;
#    for j in "${array[@]}"
#    do
#	if [[ "$j" =~ [a-zA-Z] ]] ; then
#	ns=$((ns+1));
#	echo $ns;
#	echo $j ;
#	SPEC[$ns]=`echo $j | awk '{print $1}'`;
#	fi
#    done
#    echo $ns;
#    echo ${SPEC[*]};

    for ityp in `seq 1 $ns` ; do
	for i in `seq 1 $NUE` ; do
	   if [[ ${UEL[$i]} == ${SPEC[$ityp]} ]] || [[ ${UEL[$i]}X == ${SPEC[$ityp]} ]]  ; then
#	       echo "equal" , ${UEL[$i]}, ${SPEC[$ityp]};
	       TMP_LDAU="$TMP_LDAU , Hubbard_U($ityp)=${UVAL[$i]}";
	       break
	   fi
	done
    done
cat >> $File <<EOF
TMP_LDAU="$TMP_LDAU"
EOF
#echo LDAU=\"$TMP_LDAU\";
fi
sed -i "/^$/d" $File	#removes all blank lines

}

F_Labeler $@

exit
