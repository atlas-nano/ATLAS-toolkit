#!/bin/bash

if [ $# -lt 2 ]; then
	echo "usage: xyz_coord xyz_vel [lmp_file] [velocity_scaling_factor] '[box_x] [box_y] [box_z]' [lmp_data_file]"
	exit 1
fi

if [ ! -e $1 ] || [ ! -r $1 ]; then
	echo "ERROR: Cannot open $1"
	exit 1
fi
xyz_coord=$1

if [ ! -e $2 ] || [ ! -r $2 ]; then
	echo "ERROR: Cannot open $2"
	exit 1
fi
xyz_vel=$2

lmp=$(basename $xyz_coord)
lmp="${lmp%.*}"
lmp="${lmp}.lammps"
if [ $# -gt 2 ]; then
	lmp=$3
fi

vel_scale=1
if [ $# -gt 3 ]; then
	vel_scale=$4
fi

has_box=0
box=()
if [ $# -gt 4 ]; then
	box=($5)
	if [ ${#box[*]} -eq 3 ]; then
		has_box=1
	fi
fi

has_lmp_data=0
lmp_data="blank"
if [ $# -gt 5 ] && [ -e $6 ] && [ -r $6 ]; then
	has_lmp_data=1
	lmp_data=$6
fi

cat <<DATA
XYZ_COORD               $xyz_coord
XYZ_VEL                 $xyz_vel
LAMMPS_FILE             $lmp
VELOCITY_SCALING_FACTOR $vel_scale
DATA

if [ $has_box -eq 1 ]; then
	echo "BOX_lengths             ${box[*]}"
fi

if [ $has_lmp_data -eq 1 ]; then
	echo "LAMMPS_DATA_FILE        $lmp_data"
fi

rm -fr $lmp
gawk  -v lmp_file=$lmp -v has_box=$has_box -v box_str="${box[*]}" -v hld=$has_lmp_data -v ld=$lmp_data -v vs=$vel_scale '
BEGIN{
	if(has_box)
		n=split(box_str,b," ");
	nfile=0	
	if (hld) { #read the lammps data file and save the fftype type map
		while(( getline line < ld) > 0) {
			n=split(line,lv," ")
			if(n == 4 && lv[1] ~ /^[0-9]/ && lv[2] ~ /^[0-9.]*/ && lv[3] ~ /#/ && lv[4] ~ /^[A-Za-z._]/) {
				etype[lv[4]]=lv[1]
				print "LMP_DATA: recorded typeid ",lv[1]," for fftype ",lv[4]
			}
		}
		close(ld)
	}
}
{
	if(FNR==1) {
		nfile++
		nsnap=0
		natom=0
		ntype=0
	}
	if($1 ~ /^[0-9]/ && NF==1) { #number of atoms line
		nsnap++ #snapshot counter
		natom=0 #reset atom counter
		#reset box
		if(has_box && nfile==1) {
			box[nsnap][1]=b[1]
			box[nsnap][2]=b[2]
			box[nsnap][3]=b[3]
		} else if(nfile==1) {
			bmin[1]=bmin[2]=bmin[3]=999999999
			bmax[1]=bmax[2]=bmax[3]=-999999999
		}
	} else if ($1 ~ /[A-Za-z]/ && NF == 4 && $2 ~ /[0-9]/ && $3 ~ /[0-9]/ && $4 ~ /[0-9]/) {
		natom++ #increase atom counter
		if(nfile==1) { #coords
			if (! (natom in atype)) { #store atom type info
				if(! ($1 in etype)) {
					ntype++
					etype[$1]=ntype
				}
				atype[natom]=etype[$1]
			}
			crd[nsnap][natom][1] = $2
			crd[nsnap][natom][2] = $3
			crd[nsnap][natom][3] = $4
			if (! has_box) {
				if($2<bmin[1]) bmin[1]=$2
				if($2>bmax[1]) bmax[1]=$2
				if($3<bmin[2]) bmin[2]=$3
				if($3>bmax[2]) bmax[2]=$3
				if($4<bmin[3]) bmin[3]=$4
				if($4>bmax[3]) bmax[3]=$4
				box[nsnap][1]=bmax[1]-bmin[1]	
				box[nsnap][2]=bmax[2]-bmin[2]	
				box[nsnap][3]=bmax[3]-bmin[3]	
			}
		} else { #vels
			vel[nsnap][natom][1] = $2*vs
			vel[nsnap][natom][2] = $3*vs
			vel[nsnap][natom][3] = $4*vs
		}
	}
}
END{
	n = asorti(etype, dest)
	for (i = 1; i <= n; i++) {
		print "label ",dest[i]," typeID ",etype[dest[i]]
	}
	for (i=1;i<=nsnap;i++) {
		print "ITEM: TIMESTEP" >> lmp_file
		print i >> lmp_file
		print "ITEM: NUMBER OF ATOMS" >> lmp_file
		print natom >> lmp_file
		print "ITEM: BOX BOUNDS pp pp pp" >> lmp_file
		print "0  ",box[i][1] >> lmp_file
		print "0  ",box[i][2] >> lmp_file
		print "0  ",box[i][3] >> lmp_file
		print "ITEM: ATOMS id type xu yu zu vx vy vz" >> lmp_file
		for (j=1;j<=natom;j++)
			print j,atype[j],crd[i][j][1],crd[i][j][2],crd[i][j][3],vel[i][j][1],vel[i][j][2],vel[i][j][3] >> lmp_file
	}
}' $xyz_coord $xyz_vel

			
