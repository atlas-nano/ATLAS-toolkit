#!/bin/bash

#Defining the function
#Argument 1: step n for which the PWSCF.should be created 
#Argument 2: specifies the MD output file that should be read
#Argument 3: specifies the new file to write the out put (new pwscf.md)
#Argument 4: specifies the SCF file that should be read

function F_CreatePWSCF {

local n=$1
local MD_OUT_FILE=$2
local NEW_FILE=$3
local SCF_OUT_FILE=$4


#Create a new directory for the results of MD in excited states and a new pwscf.md file

:> $NEW_FILE

#read the original pwscf.md-file

# lattice parameter
local A=`grep 'lattice parameter' $MD_OUT_FILE | awk '{print $5}'`
#echo lattice parameter $A a.u.
local A=`echo $A*0.529177 | bc`
#echo lattice parameter $A Ang
# number of atoms
local NAT=`grep 'number of atoms/cell' $MD_OUT_FILE | awk '{print $5}'`
#echo $NAT atoms
# types
local NTYP=`grep 'number of atomic types' $MD_OUT_FILE | awk '{print $6}'`
#echo $NTYPE types
# dt
local DT=`grep 'Time step' $MD_OUT_FILE | awk '{print $4}'`

# type info
local masses=`grep -A$NTYP 'atomic species   valence    mass     pseudopotential' $MD_OUT_FILE | tail -$NTYP | awk '{print $3}'`
local species=`grep -A$NTYP 'atomic species   valence    mass     pseudopotential' $MD_OUT_FILE | tail -$NTYP | awk '{print $1}'`

local i=0
#echo $masses
local m=
for m in $masses
do
  i=$(( i + 1 ))
  mass_by_type[$i]=$m
  #echo $i $m ${mass_by_type[$i]}
done

i=0
local s=
for s in $species
do
  i=$(( i + 1 ))
  species_by_type[$i]=$s
  #echo $i $s ${species_by_type[$i]}
done

  

#total energy_________________________________________________
#local TOT_EN=`grep '!    total energy              =' $MD_OUT_FILE |sed -e 's/!    total energy              =  //' -e 's/Ry//' -e 's/ *$//'|awk "{if(NR==$n){print $NF}}"|tr -d ' '`

#echo $SCF_OUT_FILE
local TOT_EN=`grep '! *total energy' $SCF_OUT_FILE | awk '{print $5}'`
#echo $TOT_EN
echo "$TOT_EN" >> $NEW_FILE
echo "$n" >> $NEW_FILE

#atomic postions___________________________________________________________

# dump the positions of step n
# use grep options to help find the right step
#echo "# positions at step $n"

#grep -m $n -A $NAT 'ATOMIC_POSITIONS' $MD_OUT_FILE | tail -$(( NAT + 1 ))

local ATOM_POS=`grep -A $NAT 'ATOMIC_POSITIONS' $MD_OUT_FILE |sed -e '/^--/d' -e 's/--.*$//g'|sed -e '/^ATOMIC_POSITIONS/d' -e 's/ATOMIC_POSITIONS.*$//g'` #loscht die von grep automatisch eingefugten "--"
local ATOM_SPECIES=`echo "$ATOM_POS" |awk -v n=$n -v nat=$NAT '{if(NR>(n-2)*nat && NR<=(n-1)*nat){print }}' | awk '{print $1}'`
local ATOM_POS=`echo "$ATOM_POS" |awk -v n=$n -v nat=$NAT '{if(NR>(n-2)*nat && NR<=(n-1)*nat){print }}' | awk '{print $2, $3, $4}'`

#echo "$ATOM_POS"

i=0
local j=
#echo $ATOM_SPECIES
for s in $ATOM_SPECIES
do
  i=$(( i + 1 ))
  for j in `seq 1 $NTYP`
  do
    #echo $s $j ${species_by_type[$j]}
    if [ "$s" == "${species_by_type[$j]}" ]
    then
      species_by_atom[$i]=$j
      #echo $i ${species_by_atom[$i]}
      break
    fi
  done
done


for j in `seq 1 $NAT`; do

local TMP_POS=
local k=
for k in `seq 1 3`; do

TMP_POS=`echo "$ATOM_POS" |awk "{if(NR==$j) {print }}"`

if  [ $k -eq 1 ] ; then
TMP_POS=`echo "$TMP_POS" |awk '{print $1}'`
TMP_POS=`echo "scale=8; $TMP_POS/$A"|bc`
fi

if  [ $k -eq 2 ] ; then
TMP_POS=`echo "$TMP_POS" |awk '{print $2}'`
TMP_POS=`echo "scale=8; $TMP_POS/$A"|bc`
fi

if  [ $k -eq 3 ] ; then
TMP_POS=`echo "$TMP_POS" |awk '{print $3}'`
TMP_POS=`echo "scale=8; $TMP_POS/$A"|bc`
fi
echo "$TMP_POS" >> $NEW_FILE

done

done

echo "F" >> $NEW_FILE

#temperature___________________________________________________

local TEMP=`grep '^ *temperature *=' $MD_OUT_FILE | awk '{print $3}' | sed -n "${n}p"`
#echo temp $TEMP

echo "$TEMP" >> $NEW_FILE

#average temperaure____________________________________________
local TEMP_AV=0

for j in `seq 1 $n`; do
local tmp=`grep '^ *temperature *=' $MD_OUT_FILE | awk '{print $3}' | sed -n "${j}p"`
TEMP_AV=`echo "scale=6; $TEMP_AV+$tmp"| bc`
done
#echo temp_av $TEMP_AV
echo "$TEMP_AV" >> $NEW_FILE

#mass___________________________

local TOT_MASS=0
local MASS=
for j in `seq 1 $NAT`; do
#MASS=`echo "scale=6; ${ATOM_MASS[${j}]}*911.4442406"| bc`
#echo hello $j ${species_by_atom[$j]} ${mass_by_type[${species_by_atom[$j]}]}
MASS=`echo "scale=6; ${mass_by_type[${species_by_atom[$j]}]}*911.4442406"| bc`
TOT_MASS=`echo "scale=6; ${TOT_MASS}+${MASS}"| bc`

echo "$MASS" >> $NEW_FILE

done
echo "$TOT_MASS" >> $NEW_FILE

#elapsed_time____________________________________

local ETIME="0"`echo "scale=8; $DT*$n*0.00004837768652"|bc`
echo "$ETIME" >> $NEW_FILE

#Reference positions__________________________________
local REF_POS=`grep -A $NAT 'positions (a_0 units)' $MD_OUT_FILE |awk '{if(NR>1){print $7, $8, $9}}'`
echo "$REF_POS" >> $NEW_FILE


}

F_CreatePWSCF $1 $2 $3 $4
echo pwscf.md completed
exit
