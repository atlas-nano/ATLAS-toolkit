#!/bin/bash

if [ $# -eq 0 ]; then
	echo "usage: $0 jaguar_input_file [qchem_input]"
	exit 1
fi

if [ ! -e $1 ] || [ ! -r $1 ]; then
	echo "ERROR: Cannot access $1"
	exit 1
fi

qchem_name=`basename $1`
qchem_name=${qchem_name%.*}
qchem_name="${qchem_name}.qcin"
if [ $# -gt 1 ]; then
	qchem_name=$2
fi

mcharge=0
multip=1

if [ `egrep -c '^molchg' $1` -gt 0 ]; then
	mcharge=`egrep '^molchg' $1 | sed 's/^.*=//'`
	echo "mcharge: $mcharge"
fi

cat > $qchem_name <<DATA
\$molecule
$mcharge $multip
`awk 'BEGIN{valid=0}{if($1 ~ /^\&zmat/)valid=1; else if($1 ~/^\&$/) valid=0; else if(valid==1 && NF==4)print}' $1`

\$rem
   METHOD             hf
   BASIS              6-311++G**
   IANLTY             200
\$end

\$plots
   Plot the HOMO and the LUMO on a line
  200   -15.0  15.0
  200   -15.0  15.0
  200  -15.0  15.0
   4   0   0   0
   1   2   3   4
\$end

DATA
