#!/bin/bash

if [ $# -eq 0 ]; then
	echo "usage: $0 lmp_input_file"
	exit 1
fi

if [ ! -e $1 ] || [ ! -r $1 ]; then
	echo "ERROR: Cannot access lammps input file $1"
	exit 1
fi

echo "harmonic bond kcal"
grep bond_coeff $1 | awk '{print $6,$7,$8,$3*2,$4,0,10}'
echo "three bond regular kcal"
grep angle_coeff $1 | awk '{print $6,$7,$8,$3*2,$4,0,10}' 
echo "torsion bond intra kcal"
grep dihedral_coeff $1 | awk '{a=0;if($4==-1)a=180;print $7,$8,$9,$10,$3,$5,a,0,10}'
cat <<DATA

omega  0.001 0.001 30000
omega_damping 0.1
broaden scale 1.0
omega_af rads 11.6 2.76 3615.0
temperature 300
DATA
