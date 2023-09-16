#!/bin/bash
# System.sh - determines the system size
# number of atoms NAT
# number of types NTYP
# number of electrons NELEC
# and number of bands NBND

# Assumes access to Input_Block.in
. ./Input_Block.in

# Determine the number of atoms
NAT=`echo "$ATOMIC_POSITIONS" | wc -l`
NAT=$(( $NAT - 1 ))
# Determine the number of types
NTYP=`echo "$ATOMIC_SPECIES" | wc -l`
NTYP=$(( $NTYP - 1 ))
# Determine charges of ions
for ityp in `seq 1 $NTYP` ; do
  ityp1=$(( $ityp + 1 ))
  PP=`echo "$ATOMIC_SPECIES" | head -$ityp1 | tail -1 | awk '{print $3}'`
  if [ ! -f $PSEUDO_DIR/$PP ] ; then
    echo Pseudopotential not found: $PSEUDO_DIR/$PP
    exit 1
  fi
  Zion[$ityp]=`grep 'Z valence' $PSEUDO_DIR/$PP | awk '{print $1}'`
done
# Determine the number of electrons
NELEC=0.0
for i in `seq 1 $NAT` ; do
  NELEC=`echo "($NELEC + ${Zion[$ityp]})" | bc`
done
# Determine the number of bands
NBND=`echo "($NELEC*0.5*1.2)/1" | bc`

echo $NAT $NTYP $NELEC $NBND
